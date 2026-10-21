! The passive-scalar STATE and its constants (split out of scalar.f90 on
! 2026-10-21): scalar_type -- the per-scalar configuration and the device-
! mapped field arrays -- the [scalar] key constants, and the accessors that
! read nothing but the type. scalar_config.f90 fills it, scalar.f90 and
! scalar_conjugate.f90 advance it, scalar_stats.f90 reads it.
module scalar_types
    use, intrinsic :: iso_c_binding
    use :: init, only: VAR_U, VAR_P
    use :: boundary, only: boundary_type, NFACES, BC_DIRICHLET, BC_NEUMANN
    use :: turbulence, only: turb_type
    implicit none
    private

    ! Initial-condition profiles ([scalar.N] init_profile).
    integer(C_INT), parameter, public :: SC_INIT_UNIFORM = 0_C_INT
    integer(C_INT), parameter, public :: SC_INIT_LINEAR_Y = 1_C_INT
    ! Turbulent-Prandtl models ([scalar.N] prt_model): the per-scalar
    ! constant prt, or the Kays-Crawford correlation (prt_kays below).
    integer(C_INT), parameter, public :: SC_PRT_CONSTANT = 0_C_INT
    integer(C_INT), parameter, public :: SC_PRT_KAYS = 1_C_INT
    ! Volumetric source shape ([scalar.N] source_type).
    !   UNIFORM  s += source                       (the S0 behaviour)
    !   VELOCITY s += source * u_x                 KASAGI's fully-developed
    ! channel formulation (Kasagi, Tomita & Kuroda 1992). For a channel with a
    ! uniform wall flux the mean temperature rises linearly in x, so writing
    ! T = beta x + theta with theta homogeneous in x turns u.grad T into
    ! beta u + u.grad theta: the extra term is a source proportional to the
    ! INSTANTANEOUS streamwise velocity, not the mean. That is the only thing
    ! separating our bulk-heated channel from the reference DNS
    ! (tutorials/cht/channel), whose flux profile it reproduces exactly
    ! instead of approximating it by a uniform source.
    ! The homogeneous direction is [scalar.N] source_dir, DEFAULT x -- the
    ! channels' streamwise direction, so an ini that does not name one keeps
    ! the original arithmetic operand for operand. A pipe along z needs
    ! source_dir = z (the pipe campaign, docs/history/next_session_pipe_cht.md F1):
    ! stats_layout = plane averages over z and returns the (x,y) plane, i.e.
    ! the cross-section, which pins the axis to z from the other side.
    integer(C_INT), parameter, public :: SC_SRC_UNIFORM  = 0_C_INT
    integer(C_INT), parameter, public :: SC_SRC_VELOCITY = 1_C_INT
    ! Direction names, for the config echo. Index = the stored srcDir.
    character(len=1), parameter, public :: SRC_DIR_NAME(3) = ["x", "y", "z"]
    ! Immersed-body wall modes ([scalar.N] ibm_wall); used from S3.
    integer(C_INT), parameter, public :: SC_IBM_DIRICHLET = 0_C_INT
    integer(C_INT), parameter, public :: SC_IBM_ADIABATIC = 1_C_INT
    ! CONJUGATE heat transfer at the immersed interface (increment C1,
    ! docs/history/next_session_conjugate.md). The solid stops being a boundary
    ! condition and becomes a REAL unknown carrying its own conductivity
    ! ratio kappa_s = k_s/k_f and capacity ratio C_s = (rho c)_s/(rho c)_f
    ! (both 1 in the fluid by definition), and a face whose two cell centres
    ! straddle the interface takes the DISTANCE-WEIGHTED HARMONIC MEAN of the
    ! two materials' diffusivities. It is a THIRD branch, not a
    ! generalisation of the other two: S3's dirichlet shipped as an inline
    ! penalization statement and adiabatic as six face masks, both gated, and
    ! re-expressing them through this arithmetic could not be bit-exact
    ! (strategy doc Section 5).
    integer(C_INT), parameter, public :: SC_IBM_CONJUGATE = 2_C_INT

    ! Statistics layouts ([scalar] stats_layout, increment S4): wall-normal
    ! rows x-z averaged (the channel form) or rows of the global (x,y) plane
    ! z averaged (the boundary-layer form). See scalar_stats.f90.
    integer(C_INT), parameter, public :: SC_LAYOUT_PROFILE = 0_C_INT
    integer(C_INT), parameter, public :: SC_LAYOUT_PLANE = 1_C_INT

    integer, parameter, public :: SC_NAME_LEN = 32

    ! [scalar] convection: the factor f in conv = conv_div - f*s*(div u).
    integer(C_INT), parameter, public :: SC_CONV_DIV  = 0_C_INT   ! f = 0
    integer(C_INT), parameter, public :: SC_CONV_SKEW = 1_C_INT   ! f = 1/2
    integer(C_INT), parameter, public :: SC_CONV_ADV  = 2_C_INT   ! f = 1

    type, public :: scalar_type
        integer(C_INT) :: n = 0_C_INT
        ! [scalar] convection = divergence (default) | skew | advective.
        ! All three are `conv_div - f*s*(div u)|stencil` with f = 0, 1/2, 1,
        ! and the choice is a real trade, which is why it is a key:
        !   divergence (f=0)  exact global conservation (conserve.ini gates
        !                     it to round-off); a uniform scalar is NOT
        !                     preserved when div u is only projection-small.
        !   skew       (f=1/2) the MOMENTUM kernel's form: neutral in sum s^2
        !                     for any advecting field, but leaves 1/2 s div u
        !                     on a uniform field.
        !   advective  (f=1)  u.grad s: a uniform scalar is preserved exactly
        !                     for any advecting field; conservation given up.
        ! The scalar plan specified f=1 and CALLED it "skew", and the two
        ! operators shared `[flow] convection` with momentum until the 2026-09-25
        ! S3 lockdown. They are different operators, so the scalar's choice
        ! lives here under its own name.
        ! NOT MEASURED: every gate in validation/scalar runs the DEFAULT
        ! divergence form, so which form is best for a scalar is open -- see
        ! docs/history/next_session_port_finish.md.
        integer(C_INT) :: convMode = SC_CONV_DIV
        ! Was [scalar] convection written in the ini? A conjugate scalar
        ! resolves an UNSET key to advective (review finding F2, 2026-09-27):
        ! the cut-face convective mask drops a genuine fluid flux at fluid-
        ! centred cut faces, and only the advective form -- whose divergence
        ! is built from the SAME masked velocities -- keeps a uniform scalar
        ! uniform there (measured: drift exactly 0.0 vs 4.1e-5 in divergence
        ! form on validation/conjugate/wavy.ini, 1000 steps, flow spinning up).
        logical :: convSet = .false.
        ! Per-scalar configuration (all sized n; host + device).
        real(C_DOUBLE), allocatable :: pr(:), prt(:)
        integer(C_INT), allocatable :: prtModel(:)
        real(C_DOUBLE), allocatable :: source(:), initValue(:), ibmValue(:), inlet(:)
        integer(C_INT), allocatable :: ibmMode(:), initProfile(:), srcType(:)
        ! Which velocity component the VELOCITY source rides (1 = x, 2 = y,
        ! 3 = z); ignored by the UNIFORM source, and a source_dir on one is a
        ! config error rather than a silent no-op.
        integer(C_INT), allocatable :: srcDir(:)
        ! Conjugate body mode (increment C1): the solid's material ratios and
        ! its own initial value / volumetric source, per scalar. kappa and C
        ! are 1 in the FLUID by definition, so neither needs a field -- the
        ! sign of phi below selects the material pointwise.
        real(C_DOUBLE), allocatable :: solidK(:), solidC(:)
        real(C_DOUBLE), allocatable :: solidInit(:), solidSource(:), contactR(:)
        ! The solid's SECOND material band (the pipe jacket). solidDepth > 0
        ! splits the solid at the iso-surface phi = -solidDepth: a SHELL of
        ! that thickness wrapped around the body surface, carrying solidK /
        ! solidC / solidSource, and everything deeper carrying outerK / outerC
        ! and NO source. solidDepth = 0 (the default) is one solid material,
        ! i.e. today's scheme operand for operand. See band_face_diffusivity.
        real(C_DOUBLE), allocatable :: solidDepth(:), outerK(:), outerC(:)
        ! [scalar.N] tangential_correction (increment C2), 0/1 rather than a
        ! logical so it maps to the device exactly like ibmMode. DEFAULT OFF:
        ! the C1 baseline deliberately drops the tangential term of the exact
        ! cut-face flux, and whether adding it back is worth its cost is a
        ! MEASUREMENT (see validation/conjugate/README.md, C2).
        integer(C_INT), allocatable :: tangCorr(:)
        ! Host-only: which optional keys the ini actually set, so a solid_*
        ! key on a non-conjugate scalar (or an ibm_value on a conjugate one)
        ! is a hard config error rather than a silent no-op.
        logical, allocatable :: solidKeySet(:), solidInitSet(:), ibmValueSet(:)
        logical, allocatable :: srcDirSet(:)
        ! Signed distance to the immersed surface at the cell centres,
        ! GHOST-INCLUSIVE: phi = +dwall in the fluid, -dwall in the solid,
        ! the sign taken from the cell-centred IBM marker. This ONE field is
        ! the entire geometric input of the conjugate scheme (strategy doc
        ! Section 8: no new dataset, no moby_prepare change, no case-file
        ! format change -- dwall already exists at the VAR_P position for
        ! both geometry paths). A 1-cell dummy without a conjugate scalar,
        ! the nutNone idiom: mapped uniformly, every access mode-guarded.
        real(C_DOUBLE), allocatable :: phi(:,:,:,:)
        ! FLUID VOLUME FRACTION of each cell (increment C3), from the same phi
        ! and its central differences: the cell's volumetric capacity is then
        ! C_cell = f + (1 - f) C_s instead of the pointwise C1 value. It is a
        ! pure GEOMETRY field -- scalar-independent, built once at init beside
        ! phi, never touched again -- so it is stored rather than recomputed
        ! per substage, which also makes the number the transport kernel
        ! divides by, the number the time-step limiter uses and the number the
        ! checkers read (it rides the snapshot as "vfrac") THE SAME number.
        ! Same shape and same 1-cell-dummy idiom as phi.
        real(C_DOUBLE), allocatable :: vfrac(:,:,:,:)
        integer(C_INT) :: nConjugate = 0_C_INT
        ! The two transport kernels' scalar lists (1-based scalar indices),
        ! built once in validate_conjugate_config: plainList the dirichlet /
        ! adiabatic scalars (every scalar of a body-free run), conjList the
        ! conjugate ones. Sized max(1, count) -- the 1-cell dummy idiom, so
        ! the device maps stay uniform -- and a kernel launches only for a
        ! non-empty list.
        integer(C_INT), allocatable :: plainList(:), conjList(:)
        integer(C_INT) :: nPlain = 0_C_INT
        integer(C_INT) :: nTangential = 0_C_INT
        integer(C_INT) :: nBanded = 0_C_INT
        ! Per-face boundary rows (n, NFACES): BC_DIRICHLET / BC_NEUMANN with
        ! the face value / normal derivative. They live HERE, not in
        ! boundary_type, so bc keeps its fixed VAR_U:VAR_P shape and the
        ! restart bc_type/bc_value attributes keep their 24-entry layout.
        integer(C_INT), allocatable :: bcType(:,:)
        real(C_DOUBLE), allocatable :: bcValue(:,:)
        logical, allocatable :: bcTypeSet(:,:), bcValueSet(:,:)   ! host only
        ! q variable indices of the scalars, for the batched halo exchange.
        integer(C_INT), allocatable :: varList(:)
        ! Inverse centre-to-centre distances at the p position (the scalar
        ! diffusion's face metric). Mirrors turb%inv_d?(:,VAR_P,:), but the
        ! turbulence module returns early when the model is none, so a DNS
        ! run has no tables -- the scalar module owns its own.
        real(C_DOUBLE), allocatable :: invDx(:,:), invDy(:,:), invDz(:,:)  ! (0:nb+1,nBlocks)
        ! Inverse two-cell distances 1/(x(i+1) - x(i-1)) at the p position:
        ! the CENTRAL-difference metric the C2 tangential gradient needs (the
        ! face metric above is the one-cell arm). Same shape, same idiom.
        real(C_DOUBLE), allocatable :: cdx(:,:), cdy(:,:), cdz(:,:)        ! (0:nb+1,nBlocks)
        ! 1-cell dummy standing in for turb%nut when no turbulence model is
        ! active (that array is then not allocated at all): it gives the
        ! transport kernel something to map, and every access to it is
        ! guarded by useNut -- the turb_type "1-cell dummies, uniform device
        ! maps, all accesses model-guarded" idiom.
        real(C_DOUBLE), allocatable :: nutNone(:,:,:,:)
        ! Thermal wall function (S5a, [rans] wall_treatment = wall_function).
        ! Per scalar: Jayatilleke's P and the thermal sublayer thickness
        ! y+_T, both pure functions of (Pr, Pr_t) and computed once at init.
        real(C_DOUBLE), allocatable :: wfP(:), wfYpt(:)         ! (n)
        ! The wall-cell y+ field rans_wall_yplus fills (ghost-inclusive, zero
        ! away from wall cells). Full size ONLY under wall functions; a 1-cell
        ! dummy otherwise -- the nutNone idiom, so the device maps and the
        ! kernel map clauses stay uniform and every access is mode-guarded.
        real(C_DOUBLE), allocatable :: wfYplus(:,:,:,:)
        ! Host-only: dataset names and the config bookkeeping.
        character(len=SC_NAME_LEN), allocatable :: name(:)
        logical, allocatable :: sectionSeen(:)
        integer(C_INT) :: countKey = -1_C_INT    ! [scalar] count, -1 = unset
        ! Host-only statistics configuration ([scalar], increment S4). It is
        ! parsed here so ONE section configures everything about the scalars;
        ! the accumulator itself lives in scalar_stats.f90 (which uses this
        ! module, hence no state of its own can live in scalar_type). All
        ! intervals off by default: nothing in scalar_stats runs unless asked.
        integer(C_INT) :: statsLayout = SC_LAYOUT_PROFILE
        integer(C_INT) :: statsSample = -1_C_INT
        integer(C_INT) :: statsWrite = -1_C_INT
        integer(C_INT) :: heatInterval = -1_C_INT
        ! [scalar] indicator_interval (increment C2): how often the cut-face
        ! error indicator e_face = h_d s_t/(T_R - T_L) is reduced and printed.
        ! Off by default; it is a diagnostic of the BASELINE's dropped term,
        ! so it does not need (and does not imply) the correction.
        integer(C_INT) :: indicatorInterval = -1_C_INT
        character(len=256) :: statsFile = "scalar_stats.h5"
        character(len=256) :: heatFile = "scalar_heat.txt"
    end type scalar_type

    public :: scalar_conjugate_enabled, scalar_count, scalar_min_pr, scalar_min_prt, scalar_names
    public :: scalars_enabled

contains

    logical function scalars_enabled(sc)
        type(scalar_type), intent(in) :: sc

        scalars_enabled = sc%n > 0_C_INT
    end function scalars_enabled

    integer(C_INT) function scalar_count(sc) result(n)
        type(scalar_type), intent(in) :: sc

        n = sc%n
    end function scalar_count

    ! Host-only dataset-name table (io labels); empty when no scalars.
    function scalar_names(sc) result(names)
        type(scalar_type), intent(in) :: sc
        character(len=SC_NAME_LEN), allocatable :: names(:)

        if (allocated(sc%name)) then
            names = sc%name(1:int(sc%n))
        else
            allocate(names(0))
        end if
    end function scalar_names

    ! Smallest molecular Prandtl/Schmidt number (the binding Peclet limit).
    real(C_DOUBLE) function scalar_min_pr(sc) result(prmin)
        type(scalar_type), intent(in) :: sc

        prmin = 1.0d0
        if (sc%n > 0_C_INT) prmin = minval(sc%pr(1:int(sc%n)))
    end function scalar_min_pr

    logical function scalar_conjugate_enabled(sc)
        type(scalar_type), intent(in) :: sc

        scalar_conjugate_enabled = sc%nConjugate > 0_C_INT
    end function scalar_conjugate_enabled

    ! Smallest turbulent Prandtl number any scalar can reach (the binding
    ! eddy-diffusivity Peclet limit). Under prt_model = kays the correlation
    ! stays in [prt, 2 prt] -- prt is its Pe_t -> infinity asymptote -- so the
    ! configured prt is the floor for both models.
    real(C_DOUBLE) function scalar_min_prt(sc) result(prtmin)
        type(scalar_type), intent(in) :: sc

        prtmin = 1.0d0
        if (sc%n > 0_C_INT) prtmin = minval(sc%prt(1:int(sc%n)))
    end function scalar_min_prt

end module scalar_types
