!--------------------------!
!                          !
!     Passive scalars      !
!                          !
!--------------------------!
!
! Advected/diffused cell-centred scalars s_1..s_N (docs/next_session_scalar.md).
!
! DESIGN (increment S0/S1): a scalar is an EXTRA VARIABLE of blk%q, stored at
! the PRESSURE point, so
!   * the metric arrays never grow -- every scalar stencil reads the VAR_P
!     column of blk%x/d1x/...;
!   * the 2:1-interface halo exchange serves scalars for free: comm.f90 treats
!     any variable index > 3 as cell-centred (8-cell restrict average, injected
!     prolong, the (2 coarse + fine)/3 ghost blend), i.e. the pressure
!     treatment, and scalars ride the SAME MPI message as u,v,w,p;
!   * output/restart is one file, one collective write.
! `[scalar] count = 0` is bit-exact with a scalar-free build BY CONSTRUCTION:
! dns%nVar == NVAR reproduces every allocation shape, no kernel here is called,
! apply_bc is untouched and no extra exchange is issued.
!
! LANDMINE: the scalar index DIFFERS between the arrays -- q slot VAR_S0+is,
! qs/oldrhs slot SCR_S0+is (pressure has no scratch plane).
!
! Per-scalar state lives in scalar_type as flat ALLOCATABLE arrays mapped to
! the device once (enter_scalar_data), the turb_type/sst_type/bodyforce_type
! pattern. There is deliberately NO fixed bound on the scalar count.
module scalar
    use, intrinsic :: iso_c_binding
    use :: init, only: dns_type, VAR_U, VAR_V, VAR_W, VAR_P, &
        VAR_S0, SCR_S0, NVAR
    use :: blocks, only: block_set_type, FACE_CLOSED, FACE_COARSE, FACE_FINE
    use :: boundary, only: boundary_type, NFACES, boundary_face_id, &
        apply_scalar_bc_q, BC_DIRICHLET, BC_NEUMANN, &
        PATCH_GENERIC, PATCH_WALL, PATCH_INLET, PATCH_OUTLET
    use :: comm, only: comm_type, exchange_halos, comm_allreduce_sum, comm_allreduce_max
    use :: turbulence, only: turb_type, turbulence_is_enabled
    ! S5a: the log law of the wall is rans.f90's -- kappa and E are
    ! use-associated so the thermal wall function below and the momentum one
    ! (nut_wall_value / omega_wall_blend) can never drift apart.
    ! C1: the conjugate signed distance IS the RANS wall distance, so its two
    ! producers are use-associated rather than reimplemented (strategy doc
    ! Section 8: "no new computation, only a new trigger").
    use :: rans, only: WF_KAPPA, WF_E, fill_body_distance_analytic
    use :: ibmm, only: ibm_type, isInBody
    use :: io, only: read_dwall_blocks
    implicit none

    private

    ! Initial-condition profiles ([scalar.N] init_profile).
    integer(C_INT), parameter, public :: SC_INIT_UNIFORM = 0_C_INT
    integer(C_INT), parameter, public :: SC_INIT_LINEAR_Y = 1_C_INT
    ! Turbulent-Prandtl models ([scalar.N] prt_model): the per-scalar
    ! constant prt, or the Kays-Crawford correlation (prt_kays below).
    integer(C_INT), parameter, public :: SC_PRT_CONSTANT = 0_C_INT
    integer(C_INT), parameter, public :: SC_PRT_KAYS = 1_C_INT
    ! Immersed-body wall modes ([scalar.N] ibm_wall); used from S3.
    integer(C_INT), parameter, public :: SC_IBM_DIRICHLET = 0_C_INT
    integer(C_INT), parameter, public :: SC_IBM_ADIABATIC = 1_C_INT
    ! CONJUGATE heat transfer at the immersed interface (increment C1,
    ! docs/next_session_conjugate.md). The solid stops being a boundary
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

    ! Grazing-arm guard (strategy doc Section 8, LaTeX note Section 6.6). The
    ! level-set weight w is exact for a PLANE; for a curved interface the two
    ! cell centres measure to different nearest points and w carries a
    ! relative error O(curvature*h/a), with a = cos(angle between the arm and
    ! the interface normal). Arms nearly TANGENT to the surface therefore
    ! carry an unreliable w -- and little of the interface flux -- so below
    ! this direction cosine the face falls back to the plain harmonic mean.
    real(C_DOUBLE), parameter :: CONJ_MIN_COSINE = 5.0d-2
    ! Smallest magnitude a SOLID cell's signed distance is allowed to take,
    ! so that `phi < 0` is an EXACT material test even for a cell centre
    ! sitting exactly on the surface (see init_scalar_conjugate).
    real(C_DOUBLE), parameter :: CONJ_MIN_DISTANCE = 1.0d-300
    ! Statistics layouts ([scalar] stats_layout, increment S4): wall-normal
    ! rows x-z averaged (the channel form) or rows of the global (x,y) plane
    ! z averaged (the boundary-layer form). See scalar_stats.f90.
    integer(C_INT), parameter, public :: SC_LAYOUT_PROFILE = 0_C_INT
    integer(C_INT), parameter, public :: SC_LAYOUT_PLANE = 1_C_INT

    integer, parameter :: SC_NAME_LEN = 32

    type, public :: scalar_type
        integer(C_INT) :: n = 0_C_INT
        ! Per-scalar configuration (all sized n; host + device).
        real(C_DOUBLE), allocatable :: pr(:), prt(:)
        integer(C_INT), allocatable :: prtModel(:)
        real(C_DOUBLE), allocatable :: source(:), initValue(:), ibmValue(:), inlet(:)
        integer(C_INT), allocatable :: ibmMode(:), initProfile(:)
        ! Conjugate body mode (increment C1): the solid's material ratios and
        ! its own initial value / volumetric source, per scalar. kappa and C
        ! are 1 in the FLUID by definition, so neither needs a field -- the
        ! sign of phi below selects the material pointwise.
        real(C_DOUBLE), allocatable :: solidK(:), solidC(:)
        real(C_DOUBLE), allocatable :: solidInit(:), solidSource(:), contactR(:)
        ! Host-only: which optional keys the ini actually set, so a solid_*
        ! key on a non-conjugate scalar (or an ibm_value on a conjugate one)
        ! is a hard config error rather than a silent no-op.
        logical, allocatable :: solidKeySet(:), solidInitSet(:), ibmValueSet(:)
        ! Signed distance to the immersed surface at the cell centres,
        ! GHOST-INCLUSIVE: phi = +dwall in the fluid, -dwall in the solid,
        ! the sign taken from the cell-centred IBM marker. This ONE field is
        ! the entire geometric input of the conjugate scheme (strategy doc
        ! Section 8: no new dataset, no moby_prepare change, no case-file
        ! format change -- dwall already exists at the VAR_P position for
        ! both geometry paths). A 1-cell dummy without a conjugate scalar,
        ! the nutNone idiom: mapped uniformly, every access mode-guarded.
        real(C_DOUBLE), allocatable :: phi(:,:,:,:)
        integer(C_INT) :: nConjugate = 0_C_INT
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
        character(len=256) :: statsFile = "scalar_stats.h5"
        character(len=256) :: heatFile = "scalar_heat.txt"
    end type scalar_type

    public :: scalars_enabled, scalar_count
    public :: apply_scalar_config, validate_scalar_config, destroy_scalar
    public :: init_scalar, init_scalar_fields, enter_scalar_data, exit_scalar_data
    public :: scalar_sync, scalar_finish, scalar_transport
    public :: scalar_names, scalar_min_pr, scalar_min_prt
    public :: scalar_section_index, prt_kays, eddy_diffusivity
    ! Conjugate heat transfer (increment C1): the geometry trigger, the solid
    ! initial condition, the material-max diffusivity the Peclet limiter
    ! needs, and the one face coefficient the whole scheme consists of
    ! (public because scalar_stats.f90 must report the flux the transport
    ! kernel APPLIED, so it shares this expression rather than copying it).
    public :: scalar_conjugate_enabled, init_scalar_conjugate
    public :: init_scalar_solid_fields, scalar_conjugate_to_device
    public :: scalar_conjugate_peclet_rate, conjugate_face_diffusivity
    ! S5a thermal wall function: the face diffusivity scalar_stats.f90 needs
    ! (its wall flux must be the flux the transport kernel applied), and the
    ! two correlations, public for the host-side unit test.
    public :: wall_face_diffusivity, jayatilleke_p, thermal_yplus, wall_diffusivity

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

    ! Face diffusivity of the conjugate interface scheme -- the ONE new
    ! coefficient the whole increment consists of (strategy doc Section 1,
    ! LaTeX note Eq. (6.baseline)). phiL/phiR are the SIGNED distances at the
    ! two cell centres (positive in the fluid), invd the inverse
    ! centre-to-centre distance, dm = 1/(Re Pr) the FLUID molecular
    ! diffusivity and ks = kappa_s the solid-to-fluid conductivity ratio.
    !
    !   uncut  ->  kappa*dm: one material across the face, i.e. today's
    !              kernel line with kappa = 1 in the fluid;
    !   cut    ->  the DISTANCE-WEIGHTED HARMONIC MEAN dm/(w/k_L + (1-w)/k_R)
    !              on the level-set fraction w = phiL/(phiL - phiR),
    !              i.e. two resistances in series.
    !
    ! THE OBLIQUITY LEMMA is what makes w usable at any interface
    ! orientation: phi stores the PERPENDICULAR distance to the surface, not
    ! the distance along the arm, so both arms are shortened by the same
    ! direction cosine a = n.e_d and it CANCELS in the ratio. w is then
    ! exactly the fraction delta_L/h_d that the series resistance needs, with
    ! no normal ever computed and no per-arm cut record ever stored -- which
    ! is why this increment adds no dataset. Exact for a plane at any angle.
    !
    ! The cut test is the MATERIAL DISAGREEMENT of the two centres rather
    ! than the strategy doc's phiL*phiR < 0. They are the same test (the sign
    ! of phi IS the marker), but the marker form cannot be defeated by a cell
    ! centre lying exactly on the surface -- init_scalar_conjugate keeps a
    ! solid cell's phi strictly negative for precisely that reason.
    !
    ! A contact resistance rc (per unit area, same non-dimensionalisation)
    ! adds to the series: R = h_d(w/k_L + (1-w)/k_R)/dm + rc, hence the
    ! rc*dm*invd term. rc = 0 gives exactly the harmonic mean.
    real(C_DOUBLE) function conjugate_face_diffusivity(dm, phiL, phiR, ks, rc, invd) result(d)
        !$omp declare target
        real(C_DOUBLE), intent(in) :: dm, phiL, phiR, ks, rc, invd

        real(C_DOUBLE) :: w, a, kl, kr
        logical :: sl, sr

        sl = phiL < 0.0d0
        sr = phiR < 0.0d0
        kl = merge(ks, 1.0d0, sl)
        if (sl .eqv. sr) then
            d = dm*kl
            return
        end if
        kr = merge(ks, 1.0d0, sr)
        ! a = |phiL - phiR|/h_d is the direction cosine, a free by-product of
        ! the lemma; guard the grazing arms it identifies.
        a = abs(phiL - phiR)*invd
        if (a < CONJ_MIN_COSINE) then
            w = 0.5d0
        else
            w = phiL/(phiL - phiR)
        end if
        d = dm/(w/kl + (1.0d0 - w)/kr + rc*dm*invd)
    end function conjugate_face_diffusivity

    ! Smallest turbulent Prandtl number any scalar can reach (the binding
    ! eddy-diffusivity Peclet limit). Under prt_model = kays the correlation
    ! stays in [prt, 2 prt] -- prt is its Pe_t -> infinity asymptote -- so the
    ! configured prt is the floor for both models.
    real(C_DOUBLE) function scalar_min_prt(sc) result(prtmin)
        type(scalar_type), intent(in) :: sc

        prtmin = 1.0d0
        if (sc%n > 0_C_INT) prtmin = minval(sc%prt(1:int(sc%n)))
    end function scalar_min_prt

    ! Kays-Crawford turbulent Prandtl number (Kays 1994; Weigand, Ferguson &
    ! Crawford 1997), a pointwise function of the turbulent Peclet number
    ! Pe_t = (nu_t/nu) Pr:
    !
    !   1/Pr_t = 1/(2 Prt_inf) + C Pe_t/sqrt(Prt_inf)
    !            - (C Pe_t)^2 [1 - exp(-1/(C Pe_t sqrt(Prt_inf)))],   C = 0.3
    !
    ! Limits: Pe_t -> 0 gives Pr_t = 2 Prt_inf (the molecular limit, inert --
    ! it multiplies nu_t = 0); Pe_t -> infinity gives Pr_t -> Prt_inf, which
    ! is the per-scalar [scalar.N] prt (0.85 by default), so `constant` and
    ! `kays` share one configuration key.
    !
    ! NUMERICS: written with x = 1/(C Pe_t sqrt(Prt_inf)), the bracket is
    ! (C Pe_t)^2 [1 - exp(-x)] = a/b - 1/(2b^2) + ... and the leading a/b
    ! cancels the second term of the sum exactly -- catastrophically so once
    ! Pe_t is large. The small-x branch therefore evaluates the analytically
    ! equivalent series
    !
    !   1/Pr_t = (1/Prt_inf) [ 1 - x/3! + x^2/4! - x^3/5! + ... ]
    !
    ! which is also where the Pe_t -> infinity limit is exact by inspection.
    ! The two branches agree to round-off at the x = 1/2 crossover.
    real(C_DOUBLE) function prt_kays(pet, prtinf) result(prt)
        !$omp declare target
        real(C_DOUBLE), intent(in) :: pet      ! turbulent Peclet number
        real(C_DOUBLE), intent(in) :: prtinf   ! Pr_t at Pe_t -> infinity

        real(C_DOUBLE), parameter :: CKC = 0.3d0
        real(C_DOUBLE) :: a, b, x, u, sum, invprt
        integer :: m

        b = sqrt(prtinf)
        a = CKC*max(pet, 0.0d0)
        if (a*b < 1.0d-300) then
            prt = 2.0d0*prtinf              ! the molecular limit
            return
        end if
        x = 1.0d0/(a*b)

        if (x < 0.5d0) then
            sum = 1.0d0
            u = -x/6.0d0                    ! the m = 1 term, -x/3!
            do m = 1, 16                    ! machine precision for x < 1/2
                sum = sum + u
                u = -u*x/real(m + 3, C_DOUBLE)
            end do
            invprt = sum/prtinf
        else
            invprt = 0.5d0/prtinf + a/b - a*a*(1.0d0 - exp(-x))
        end if

        prt = 1.0d0/invprt
    end function prt_kays

    ! Face eddy diffusivity nu_t/Pr_t for one scalar. Called per face inside
    ! the transport kernel, so it is a declare-target pure function.
    real(C_DOUBLE) function eddy_diffusivity(nutFace, pr, prt, prtModel, re) result(d)
        !$omp declare target
        real(C_DOUBLE), intent(in) :: nutFace, pr, prt, re
        integer(C_INT), intent(in) :: prtModel

        real(C_DOUBLE) :: nt

        nt = max(nutFace, 0.0d0)
        if (prtModel == SC_PRT_KAYS) then
            ! Pe_t = (nu_t/nu) Pr = nu_t Re Pr.
            d = nt/prt_kays(nt*re*pr, prt)
        else
            d = nt/prt
        end if
    end function eddy_diffusivity

    !--------------------------------------------------------------------
    ! Thermal wall function (increment S5a)
    !--------------------------------------------------------------------
    !
    ! Under [rans] wall_treatment = wall_function the first cell sits in the
    ! log layer, so neither the scalar's own gradient nor the wall-cell nu_t
    ! the momentum wall function installs carries the right wall FLUX: the
    ! thermal sublayer has its own thickness (it scales with Pr, not with
    ! nu) and its own resistance. The classical closure is
    !
    !   theta+ = (theta_w - theta_P) u_tau* / q_w
    !          = Pr y+                              y+ <  y+_T   (conduction)
    !          = Pr_t [ln(E y+)/kappa + P(Pr/Pr_t)] y+ >= y+_T   (log layer)
    !
    ! with u_tau* = C_mu^(1/4) sqrt(k) (rans_wall_yplus supplies the same
    ! k-based y+ the momentum wall function uses) and Jayatilleke's P the
    ! extra sublayer resistance. It is delivered EXACTLY as T3 delivers the
    ! wall shear: as a wall-cell eddy diffusivity, so the ordinary
    ! face-flux discretisation reproduces q_w with no special-cased flux.
    ! With the wall-cell value copied into the no-slip ghost the wall face
    ! sees D = nu y+/theta+ and returns
    !   D (theta_P - theta_w)/y_P = u_tau* (theta_P - theta_w)/theta+,
    ! which IS the wall-function flux.
    !
    ! Pr_t here is the per-scalar CONSTANT [scalar.N] prt even under
    ! prt_model = kays: P and the log layer are defined with a constant
    ! Pr_t, and Kays-Crawford is a correlation for the RESOLVED interior
    ! (where it keeps running -- only wall cells take this branch).

    ! Jayatilleke's P (Jayatilleke 1969; the OpenFOAM Psmooth form): the
    ! conductive sublayer's extra resistance relative to the momentum one.
    ! P = 0 at Pr = Pr_t (the Reynolds analogy), positive for Pr > Pr_t.
    pure real(C_DOUBLE) function jayatilleke_p(prat) result(p)
        real(C_DOUBLE), intent(in) :: prat        ! Pr/Pr_t

        p = 9.24d0*(prat**0.75d0 - 1.0d0)*(1.0d0 + 0.28d0*exp(-0.007d0*prat))
    end function jayatilleke_p

    ! The thermal sublayer thickness y+_T: where the conductive branch
    ! Pr y+ meets the log branch. f(y) = Pr y - Pr_t[ln(E y)/kappa + P] has
    ! its minimum at y* = Pr_t/(kappa Pr) and, whenever a conductive
    ! sublayer exists at all, exactly one root beyond it -- the physical
    ! crossing (the root BELOW y* is spurious: there the log branch is the
    ! larger of the two). Bisection, host-side, once per scalar at init;
    ! OpenFOAM iterates the same fixed point.
    !
    ! Degenerate case f(y*) >= 0 (Pr well below Pr_t -- liquid metals): the
    ! log branch never rises above the conductive one, i.e. the thermal
    ! sublayer has no extent. y+_T = y* is then the switch, which keeps the
    ! log branch positive (it is increasing beyond y*) and reduces to
    ! conduction below it.
    real(C_DOUBLE) function thermal_yplus(pr, prt, p) result(ypt)
        real(C_DOUBLE), intent(in) :: pr, prt, p

        real(C_DOUBLE) :: ylo, yhi, ym, fm
        integer :: it

        ylo = prt/(WF_KAPPA*pr)                   ! the minimum of f
        if (pr*ylo - prt*(log(WF_E*ylo)/WF_KAPPA + p) >= 0.0d0) then
            ypt = ylo
            return
        end if
        yhi = ylo
        do it = 1, 200                            ! bracket the root above
            yhi = 2.0d0*yhi
            if (pr*yhi - prt*(log(WF_E*yhi)/WF_KAPPA + p) > 0.0d0) exit
        end do
        do it = 1, 200                            ! ~1e-60 of the bracket
            ym = 0.5d0*(ylo + yhi)
            fm = pr*ym - prt*(log(WF_E*ym)/WF_KAPPA + p)
            if (fm < 0.0d0) then
                ylo = ym
            else
                yhi = ym
            end if
        end do
        ypt = 0.5d0*(ylo + yhi)
    end function thermal_yplus

    ! Wall-cell eddy diffusivity of the thermal wall function: the total
    ! nu y+/theta+ minus the molecular nu/Pr, so the transport kernel adds
    ! it to its unchanged molecular term. EXACTLY zero on the conductive
    ! branch (theta+ = Pr y+ there), and continuous at y+_T by the
    ! definition of y+_T -- the structure of T3's nut_wall_value.
    real(C_DOUBLE) function wall_diffusivity(yplus, pr, prt, re, p, ypt) result(a)
        !$omp declare target
        real(C_DOUBLE), intent(in) :: yplus, pr, prt, re, p, ypt

        real(C_DOUBLE) :: tp

        if (yplus < ypt) then
            a = 0.0d0
        else
            tp = prt*(log(WF_E*yplus)/WF_KAPPA + p)
            a = max((yplus/tp - 1.0d0/pr)/re, 0.0d0)
        end if
    end function wall_diffusivity

    ! Cell eddy diffusivity for one scalar under wall functions: the thermal
    ! wall function at a wall cell (y+ > 0 marks them, rans_wall_yplus), the
    ! ordinary nu_t/Pr_t everywhere else.
    real(C_DOUBLE) function cell_diffusivity(nutc, yplus, pr, prt, prtModel, re, p, ypt) result(a)
        !$omp declare target
        real(C_DOUBLE), intent(in) :: nutc, yplus, pr, prt, re, p, ypt
        integer(C_INT), intent(in) :: prtModel

        if (yplus > 0.0d0) then
            a = wall_diffusivity(yplus, pr, prt, re, p, ypt)
        else
            a = eddy_diffusivity(nutc, pr, prt, prtModel, re)
        end if
    end function cell_diffusivity

    ! Face eddy diffusivity under wall functions: the 1/2 (L + R) average of
    ! the two CELL diffusivities. (Resolved mode averages nu_t first and
    ! divides by Pr_t afterwards; that is not available here, because a wall
    ! cell's diffusivity is not a function of its nu_t at all. The two agree
    ! identically for a constant Pr_t away from wall cells, and differ only
    ! in rounding -- but the wall-function branch is a separate code path,
    ! so resolved runs keep their arithmetic byte for byte.)
    real(C_DOUBLE) function wall_face_diffusivity(nutL, nutR, ypL, ypR, &
            pr, prt, prtModel, re, p, ypt) result(d)
        !$omp declare target
        real(C_DOUBLE), intent(in) :: nutL, nutR, ypL, ypR, pr, prt, re, p, ypt
        integer(C_INT), intent(in) :: prtModel

        d = 0.5d0*(cell_diffusivity(nutL, ypL, pr, prt, prtModel, re, p, ypt) &
                 + cell_diffusivity(nutR, ypR, pr, prt, prtModel, re, p, ypt))
    end function wall_face_diffusivity

    !--------------------------------------------------------------------
    ! Configuration
    !--------------------------------------------------------------------

    ! Index N of a [scalar.N] section, 0 for anything else (the
    ! grid_axis_index pattern).
    integer function scalar_section_index(section) result(idx)
        character(len=*), intent(in) :: section

        integer :: stat

        idx = 0
        if (len_trim(section) <= 7) return
        if (section(1:7) /= "scalar.") return
        read(section(8:len_trim(section)), *, iostat=stat) idx
        if (stat /= 0) idx = 0
    end function scalar_section_index

    ! One [scalar] (index 0) or [scalar.N] key. The arrays grow on demand, so
    ! there is no maximum scalar count anywhere.
    subroutine apply_scalar_config(index, key, value, sc, line_no, terminal)
        integer, intent(in) :: index
        character(len=*), intent(in) :: key, value
        type(scalar_type), intent(inout) :: sc
        integer, intent(in) :: line_no
        logical, intent(in) :: terminal

        integer :: is

        if (index == 0) then
            select case (trim(key))
            case ("count")
                is = int(read_int_value(value, key, line_no))
                if (is < 0) error stop "[scalar] count must be non-negative"
                sc%countKey = int(is, C_INT)
                call scalar_grow(sc, is)
            ! Statistics (increment S4, scalar_stats.f90). Off unless asked.
            case ("stats_sample_interval")
                sc%statsSample = read_int_value(value, key, line_no)
            case ("stats_write_interval")
                sc%statsWrite = read_int_value(value, key, line_no)
            case ("stats_file")
                sc%statsFile = trim(value)
            case ("stats_layout")
                select case (trim(value))
                case ("profile", "channel")
                    sc%statsLayout = SC_LAYOUT_PROFILE
                case ("plane", "boundary_layer", "xy")
                    sc%statsLayout = SC_LAYOUT_PLANE
                case default
                    if (terminal) print *, "error: [scalar] stats_layout must be profile or", &
                        " plane, input line", line_no
                    error stop "unknown [scalar] stats_layout"
                end select
            case ("heat_interval")
                sc%heatInterval = read_int_value(value, key, line_no)
            case ("heat_file")
                sc%heatFile = trim(value)
            case default
                if (terminal) print *, "warning: unknown [scalar] key on input line", &
                    line_no, ": ", trim(key)
            end select
            return
        end if

        if (index < 1) error stop "[scalar.N] section index must be >= 1"
        call scalar_grow(sc, index)
        sc%sectionSeen(index) = .true.
        is = index

        select case (trim(key))
        case ("name")
            sc%name(is) = trim(value)
        case ("pr", "prandtl", "schmidt", "sc")
            sc%pr(is) = read_real_value(value, key, line_no)
        case ("prt")
            sc%prt(is) = read_real_value(value, key, line_no)
        case ("prt_model")
            select case (trim(value))
            case ("constant")
                sc%prtModel(is) = SC_PRT_CONSTANT
            case ("kays", "kays_crawford", "kays-crawford")
                sc%prtModel(is) = SC_PRT_KAYS
            case default
                if (terminal) print *, "error: [scalar.N] prt_model must be constant or kays,", &
                    " input line", line_no
                error stop "unknown [scalar] prt_model"
            end select
        case ("initial", "init")
            sc%initValue(is) = read_real_value(value, key, line_no)
        case ("init_profile")
            select case (trim(value))
            case ("uniform")
                sc%initProfile(is) = SC_INIT_UNIFORM
            case ("linear_y", "linear-y")
                sc%initProfile(is) = SC_INIT_LINEAR_Y
            case default
                if (terminal) print *, "error: [scalar.N] init_profile must be uniform or", &
                    " linear_y, input line", line_no
                error stop "unknown [scalar] init_profile"
            end select
        case ("source")
            sc%source(is) = read_real_value(value, key, line_no)
        case ("ibm_wall")
            select case (trim(value))
            case ("dirichlet")
                sc%ibmMode(is) = SC_IBM_DIRICHLET
            case ("adiabatic", "neumann")
                sc%ibmMode(is) = SC_IBM_ADIABATIC
            case ("conjugate")
                sc%ibmMode(is) = SC_IBM_CONJUGATE
            case default
                if (terminal) print *, "error: [scalar.N] ibm_wall must be dirichlet,", &
                    " adiabatic or conjugate, input line", line_no
                error stop "unknown [scalar] ibm_wall"
            end select
        case ("ibm_value")
            sc%ibmValue(is) = read_real_value(value, key, line_no)
            sc%ibmValueSet(is) = .true.
        ! Conjugate solid properties (increment C1). All rejected outright on
        ! a non-conjugate scalar -- see validate_scalar_config.
        case ("solid_k", "solid_kappa")
            sc%solidK(is) = read_real_value(value, key, line_no)
            sc%solidKeySet(is) = .true.
        case ("solid_rhocp", "solid_capacity")
            sc%solidC(is) = read_real_value(value, key, line_no)
            sc%solidKeySet(is) = .true.
        case ("solid_init", "solid_initial")
            sc%solidInit(is) = read_real_value(value, key, line_no)
            sc%solidInitSet(is) = .true.
            sc%solidKeySet(is) = .true.
        case ("solid_source")
            sc%solidSource(is) = read_real_value(value, key, line_no)
            sc%solidKeySet(is) = .true.
        case ("contact_resistance")
            sc%contactR(is) = read_real_value(value, key, line_no)
            sc%solidKeySet(is) = .true.
        case ("inlet")
            sc%inlet(is) = read_real_value(value, key, line_no)
        case default
            call apply_scalar_face_key(is, key, value, sc, line_no, terminal)
        end select
    end subroutine apply_scalar_config

    ! Per-face rows: <dir>_<side>_type / _value, the [boundary] key style.
    subroutine apply_scalar_face_key(is, key, value, sc, line_no, terminal)
        integer, intent(in) :: is
        character(len=*), intent(in) :: key, value
        type(scalar_type), intent(inout) :: sc
        integer, intent(in) :: line_no
        logical, intent(in) :: terminal

        integer :: p1, p2, dir, side, face_id

        p1 = index(key, "_")
        p2 = 0
        if (p1 > 1) p2 = index(key(p1+1:), "_")
        if (p1 <= 1 .or. p2 <= 1) then
            if (terminal) print *, "warning: unknown [scalar.N] key on input line", &
                line_no, ": ", trim(key)
            return
        end if
        p2 = p1 + p2

        select case (key(:p1-1))
        case ("x"); dir = 1
        case ("y"); dir = 2
        case ("z"); dir = 3
        case default; dir = 0
        end select
        select case (key(p1+1:p2-1))
        case ("min"); side = 0
        case ("max"); side = 1
        case default; side = -1
        end select
        if (dir == 0 .or. side < 0) then
            if (terminal) print *, "warning: unknown [scalar.N] key on input line", &
                line_no, ": ", trim(key)
            return
        end if
        face_id = boundary_face_id(dir, side)

        select case (trim(key(p2+1:)))
        case ("type")
            select case (trim(value))
            case ("dirichlet", "0")
                sc%bcType(is,face_id) = BC_DIRICHLET
            case ("neumann", "1")
                sc%bcType(is,face_id) = BC_NEUMANN
            case default
                if (terminal) print *, "error: [scalar.N] boundary type must be dirichlet", &
                    " or neumann, input line", line_no
                error stop "unknown [scalar] boundary type"
            end select
            sc%bcTypeSet(is,face_id) = .true.
        case ("value")
            sc%bcValue(is,face_id) = read_real_value(value, key, line_no)
            sc%bcValueSet(is,face_id) = .true.
        case default
            if (terminal) print *, "warning: [scalar.N] boundary key must end in _type or", &
                " _value, input line", line_no
        end select
    end subroutine apply_scalar_face_key

    ! Grow the per-scalar arrays to at least n entries, preserving values
    ! (host-only ini-parse cost; the count is unknown until the file ends).
    subroutine scalar_grow(sc, n)
        type(scalar_type), intent(inout) :: sc
        integer, intent(in) :: n

        type(scalar_type) :: old
        integer :: nOld

        if (n <= int(sc%n)) return
        nOld = int(sc%n)
        if (nOld > 0) then
            old%pr = sc%pr; old%prt = sc%prt; old%prtModel = sc%prtModel
            old%source = sc%source; old%initValue = sc%initValue
            old%ibmValue = sc%ibmValue; old%inlet = sc%inlet
            old%ibmMode = sc%ibmMode; old%initProfile = sc%initProfile
            old%solidK = sc%solidK; old%solidC = sc%solidC
            old%solidInit = sc%solidInit; old%solidSource = sc%solidSource
            old%contactR = sc%contactR
            old%solidKeySet = sc%solidKeySet; old%solidInitSet = sc%solidInitSet
            old%ibmValueSet = sc%ibmValueSet
            old%bcType = sc%bcType; old%bcValue = sc%bcValue
            old%bcTypeSet = sc%bcTypeSet; old%bcValueSet = sc%bcValueSet
            old%name = sc%name; old%sectionSeen = sc%sectionSeen
            call deallocate_config(sc)
        end if

        allocate(sc%pr(n), sc%prt(n), sc%prtModel(n))
        allocate(sc%source(n), sc%initValue(n), sc%ibmValue(n), sc%inlet(n))
        allocate(sc%ibmMode(n), sc%initProfile(n))
        allocate(sc%solidK(n), sc%solidC(n), sc%solidInit(n), sc%solidSource(n))
        allocate(sc%contactR(n))
        allocate(sc%solidKeySet(n), sc%solidInitSet(n), sc%ibmValueSet(n))
        allocate(sc%bcType(n,NFACES), sc%bcValue(n,NFACES))
        allocate(sc%bcTypeSet(n,NFACES), sc%bcValueSet(n,NFACES))
        allocate(sc%name(n), sc%sectionSeen(n))

        ! Defaults. Neumann 0 on every face = adiabatic / zero gradient, which
        ! is also what the PATCH_WALL / PATCH_OUTLET / generic faces resolve to.
        sc%pr = 1.0d0
        sc%prt = 0.85d0
        sc%prtModel = SC_PRT_CONSTANT
        sc%source = 0.0d0
        sc%initValue = 0.0d0
        sc%ibmValue = 0.0d0
        sc%inlet = 0.0d0
        sc%ibmMode = SC_IBM_DIRICHLET
        sc%initProfile = SC_INIT_UNIFORM
        ! kappa_s = C_s = 1 is the inert conjugate solid (same material as the
        ! fluid), which is also the k_s = k_f case where the scheme's dropped
        ! tangential term vanishes identically -- a useful control.
        sc%solidK = 1.0d0
        sc%solidC = 1.0d0
        sc%solidInit = 0.0d0
        sc%solidSource = 0.0d0
        sc%contactR = 0.0d0
        sc%solidKeySet = .false.
        sc%solidInitSet = .false.
        sc%ibmValueSet = .false.
        sc%bcType = BC_NEUMANN
        sc%bcValue = 0.0d0
        sc%bcTypeSet = .false.
        sc%bcValueSet = .false.
        sc%name = ""
        sc%sectionSeen = .false.

        if (nOld > 0) then
            sc%pr(1:nOld) = old%pr; sc%prt(1:nOld) = old%prt
            sc%prtModel(1:nOld) = old%prtModel
            sc%source(1:nOld) = old%source; sc%initValue(1:nOld) = old%initValue
            sc%ibmValue(1:nOld) = old%ibmValue; sc%inlet(1:nOld) = old%inlet
            sc%ibmMode(1:nOld) = old%ibmMode; sc%initProfile(1:nOld) = old%initProfile
            sc%solidK(1:nOld) = old%solidK; sc%solidC(1:nOld) = old%solidC
            sc%solidInit(1:nOld) = old%solidInit
            sc%solidSource(1:nOld) = old%solidSource
            sc%contactR(1:nOld) = old%contactR
            sc%solidKeySet(1:nOld) = old%solidKeySet
            sc%solidInitSet(1:nOld) = old%solidInitSet
            sc%ibmValueSet(1:nOld) = old%ibmValueSet
            sc%bcType(1:nOld,:) = old%bcType; sc%bcValue(1:nOld,:) = old%bcValue
            sc%bcTypeSet(1:nOld,:) = old%bcTypeSet
            sc%bcValueSet(1:nOld,:) = old%bcValueSet
            sc%name(1:nOld) = old%name; sc%sectionSeen(1:nOld) = old%sectionSeen
            call deallocate_config(old)
        end if
        sc%n = int(n, C_INT)
    end subroutine scalar_grow

    ! Post-parse validation: contiguous sections, a count that agrees with
    ! the highest section index, positive Prandtl numbers and dataset names
    ! that do not collide with the reserved ones. Also derives dns%nVar.
    subroutine validate_scalar_config(sc, dns, terminal)
        type(scalar_type), intent(inout) :: sc
        type(dns_type), intent(inout) :: dns
        logical, intent(in) :: terminal

        integer :: is, js

        if (sc%countKey >= 0_C_INT .and. int(sc%countKey) /= int(sc%n)) then
            if (terminal) print *, "error: [scalar] count =", sc%countKey, &
                " disagrees with the highest [scalar.N] section", sc%n
            error stop "[scalar] count does not match the [scalar.N] sections"
        end if
        do is = 1, int(sc%n)
            if (.not. sc%sectionSeen(is)) then
                if (terminal) print *, "error: [scalar.", is, "] section missing;", &
                    " scalar sections must be contiguous from 1"
                error stop "missing [scalar.N] section"
            end if
            if (sc%pr(is) <= 0.0d0) error stop "[scalar.N] pr must be positive"
            if (sc%prt(is) <= 0.0d0) error stop "[scalar.N] prt must be positive"
            if (len_trim(sc%name(is)) == 0) write(sc%name(is), '(a,i0)') "s", is
            if (name_is_reserved(sc%name(is))) then
                if (terminal) print *, "error: [scalar.", is, "] name '", &
                    trim(sc%name(is)), "' is a reserved dataset name"
                error stop "[scalar.N] name collides with a reserved dataset"
            end if
            do js = 1, is - 1
                if (sc%name(js) == sc%name(is)) error stop "[scalar.N] duplicate name"
            end do
        end do

        call validate_conjugate_config(sc, dns, terminal)

        dns%nScalar = sc%n
        dns%nVar = NVAR + sc%n
    end subroutine validate_scalar_config

    ! Conjugate body mode (increment C1): the per-scalar solid properties and
    ! the case-level preconditions of docs/next_session_conjugate.md Section
    ! 12. Every one of these is a HARD ERROR rather than a warning: each
    ! silently produces a physically wrong answer that no gate in the passive
    ! set would flag.
    subroutine validate_conjugate_config(sc, dns, terminal)
        type(scalar_type), intent(inout) :: sc
        type(dns_type), intent(in) :: dns
        logical, intent(in) :: terminal

        integer :: is, nConj

        nConj = 0
        do is = 1, int(sc%n)
            if (sc%ibmMode(is) /= SC_IBM_CONJUGATE) then
                if (sc%solidKeySet(is)) then
                    if (terminal) print '(a,i0,a)', " error: [scalar.", is, &
                        "] solid_k / solid_rhocp / solid_init / solid_source /" // &
                        " contact_resistance need ibm_wall = conjugate"
                    error stop "[scalar.N] solid property without ibm_wall = conjugate"
                end if
                cycle
            end if
            nConj = nConj + 1
            ! The body temperature is an OUTCOME of the conjugate problem,
            ! not an input -- an ibm_value would be silently ignored.
            if (sc%ibmValueSet(is)) then
                if (terminal) print '(a,i0,a)', " error: [scalar.", is, &
                    "] ibm_value is meaningless with ibm_wall = conjugate:" // &
                    " the interface temperature is solved for, not imposed"
                error stop "[scalar.N] ibm_value with ibm_wall = conjugate"
            end if
            if (sc%solidK(is) <= 0.0d0) error stop "[scalar.N] solid_k must be positive"
            if (sc%solidC(is) <= 0.0d0) error stop "[scalar.N] solid_rhocp must be positive"
            if (sc%contactR(is) < 0.0d0) &
                error stop "[scalar.N] contact_resistance must be non-negative"
            ! Default solid initial value = the scalar's own `initial`, so a
            ! conjugate run with no solid_init starts uniform.
            if (.not. sc%solidInitSet(is)) sc%solidInit(is) = sc%initValue(is)
        end do
        sc%nConjugate = int(nConj, C_INT)
        if (nConj == 0) return

        if (.not. dns%ibm_enabled) then
            if (terminal) print *, "error: [scalar.N] ibm_wall = conjugate needs an", &
                " immersed body ([ibm] enabled = true)"
            error stop "[scalar] conjugate without an immersed body"
        end if
        ! The solid now carries a real temperature field, so a block buried
        ! inside the body is part of the SOLUTION domain. Removing it deletes
        ! the solid -- the same class of trap as the A3 penalization-force
        ! finding (validation/naca0012/README.md).
        if (dns%block_nb > 0_C_INT .and. dns%block_remove_solid &
                .and. .not. dns%block_refine_body) then
            if (terminal) print *, "error: ibm_wall = conjugate needs [blocks]", &
                " remove_solid = false -- buried blocks carry the solid temperature field"
            error stop "[scalar] conjugate with [blocks] remove_solid"
        end if
        if (dns%block_refine_body .and. .not. dns%block_keep_buried) then
            if (terminal) print *, "error: ibm_wall = conjugate needs [blocks]", &
                " keep_buried = true under refine_body -- buried leaves carry", &
                " the solid temperature field (and the case file must be", &
                " prepared with it too)"
            error stop "[scalar] conjugate without [blocks] keep_buried"
        end if
        ! Wall functions model an unresolved fluid sublayer; a conjugate
        ! interface needs the near-wall layer RESOLVED on both sides
        ! (strategy doc Section 6). Not validated -- so, an error.
        if (dns%rans_wall_treatment /= 0_C_INT) then
            if (terminal) print *, "error: ibm_wall = conjugate with [rans]", &
                " wall_treatment = wall_function is not supported"
            error stop "[scalar] conjugate with wall functions"
        end if
    end subroutine validate_conjugate_config

    logical function name_is_reserved(name)
        character(len=*), intent(in) :: name

        select case (trim(name))
        case ("un", "vn", "wn", "pn", "nut", "k", "omega", "gamma", "rethetat", &
              "fd", "x", "y", "z", "blocks")
            name_is_reserved = .true.
        case default
            name_is_reserved = .false.
        end select
    end function name_is_reserved

    real(C_DOUBLE) function read_real_value(value, key, line_no) result(x)
        character(len=*), intent(in) :: value, key
        integer, intent(in) :: line_no

        integer :: stat

        read(value, *, iostat=stat) x
        if (stat /= 0) then
            print *, "error: [scalar] ", trim(key), " needs a real value, input line", line_no
            error stop "could not parse a [scalar] value"
        end if
    end function read_real_value

    integer(C_INT) function read_int_value(value, key, line_no) result(n)
        character(len=*), intent(in) :: value, key
        integer, intent(in) :: line_no

        integer :: stat

        read(value, *, iostat=stat) n
        if (stat /= 0) then
            print *, "error: [scalar] ", trim(key), " needs an integer value, input line", line_no
            error stop "could not parse a [scalar] value"
        end if
    end function read_int_value

    !--------------------------------------------------------------------
    ! Initialisation
    !--------------------------------------------------------------------

    ! Resolve the patch-derived boundary rows and build the per-block metric
    ! tables. Runs after the config/restart metadata is final (the patch types
    ! and the block set must both exist).
    subroutine init_scalar(sc, blk, bc, wall_function, has_terminal)
        type(scalar_type), intent(inout) :: sc
        type(block_set_type), intent(in) :: blk
        type(boundary_type), intent(in) :: bc
        ! [rans] wall_treatment = wall_function: allocate the wall-cell y+
        ! field and report the thermal wall function's per-scalar constants.
        logical, intent(in) :: wall_function
        logical, intent(in) :: has_terminal

        integer :: is, dir, side, face_id, nx, ny, nz, i, b

        if (.not. scalars_enabled(sc)) return

        ! Patch-derived defaults, set-if-unset (the A0 face concept): a wall
        ! is adiabatic, an inlet holds the scalar's inlet value, an outlet
        ! and a generic patch are zero-gradient. BC_OUTFLOW is a
        ! face-staggered concept and does not apply to a cell-centred scalar.
        do dir = 1, 3
            do side = 0, 1
                face_id = boundary_face_id(dir, side)
                select case (bc%facePatchType(face_id))
                case (PATCH_INLET)
                    do is = 1, int(sc%n)
                        call resolve_scalar_row(sc, is, face_id, BC_DIRICHLET, &
                            sc%inlet(is), has_terminal, strict=.true.)
                    end do
                case (PATCH_OUTLET, PATCH_GENERIC)
                    do is = 1, int(sc%n)
                        call resolve_scalar_row(sc, is, face_id, BC_NEUMANN, 0.0d0, &
                            has_terminal, strict=.true.)
                    end do
                case (PATCH_WALL)
                    ! A wall is the ONE face kind where both scalar types are
                    ! physical: adiabatic (Neumann 0, the default) and
                    ! isothermal / fixed-concentration (Dirichlet) are equally
                    ! valid, so an explicit _type key is HONOURED here rather
                    ! than treated as a contradiction (Section 5 of
                    ! docs/next_session_scalar.md: "unless the ini gives a
                    ! type/value"). Inlet and outlet stay strict: an inlet must
                    ! impose a value, an outlet must be zero-gradient.
                    do is = 1, int(sc%n)
                        call resolve_scalar_row(sc, is, face_id, BC_NEUMANN, 0.0d0, &
                            has_terminal, strict=.false.)
                    end do
                end select
            end do
        end do

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))
        if (allocated(sc%invDx)) deallocate(sc%invDx, sc%invDy, sc%invDz)
        allocate(sc%invDx(0:nx+1,blk%nBlocks), sc%invDy(0:ny+1,blk%nBlocks), &
                 sc%invDz(0:nz+1,blk%nBlocks))
        do b = 1, int(blk%nBlocks)
            do i = 0, nx+1
                sc%invDx(i,b) = inv_delta(blk%x(i,VAR_P,b) - blk%x(i-1,VAR_P,b))
            end do
            do i = 0, ny+1
                sc%invDy(i,b) = inv_delta(blk%y(i,VAR_P,b) - blk%y(i-1,VAR_P,b))
            end do
            do i = 0, nz+1
                sc%invDz(i,b) = inv_delta(blk%z(i,VAR_P,b) - blk%z(i-1,VAR_P,b))
            end do
        end do

        if (allocated(sc%nutNone)) deallocate(sc%nutNone)
        allocate(sc%nutNone(0:0,0:0,0:0,1))
        sc%nutNone = 0.0d0

        ! Conjugate signed distance: full size (ghost-inclusive) only with a
        ! conjugate scalar, a 1-cell dummy otherwise -- the nutNone idiom, so
        ! the device maps stay uniform and every access is mode-guarded. It
        ! is allocated HERE (before enter_scalar_data maps it) but FILLED
        ! later by init_scalar_conjugate, which needs the IBM coefficients.
        if (allocated(sc%phi)) deallocate(sc%phi)
        if (scalar_conjugate_enabled(sc)) then
            allocate(sc%phi(0:nx+1,0:ny+1,0:nz+1,blk%nBlocks))
        else
            allocate(sc%phi(0:0,0:0,0:0,1))
        end if
        sc%phi = 0.0d0

        ! Thermal wall function (S5a): the per-scalar constants are pure
        ! functions of (Pr, Pr_t) and always cheap, so they are computed
        ! unconditionally; the y+ FIELD is full size only under wall
        ! functions (a 1-cell dummy otherwise, the nutNone idiom).
        if (allocated(sc%wfP)) deallocate(sc%wfP, sc%wfYpt)
        allocate(sc%wfP(int(sc%n)), sc%wfYpt(int(sc%n)))
        do is = 1, int(sc%n)
            sc%wfP(is) = jayatilleke_p(sc%pr(is)/sc%prt(is))
            sc%wfYpt(is) = thermal_yplus(sc%pr(is), sc%prt(is), sc%wfP(is))
        end do
        if (allocated(sc%wfYplus)) deallocate(sc%wfYplus)
        if (wall_function) then
            allocate(sc%wfYplus(0:nx+1,0:ny+1,0:nz+1,blk%nBlocks))
        else
            allocate(sc%wfYplus(0:0,0:0,0:0,1))
        end if
        sc%wfYplus = 0.0d0

        if (allocated(sc%varList)) deallocate(sc%varList)
        allocate(sc%varList(int(sc%n)))
        do is = 1, int(sc%n)
            sc%varList(is) = VAR_S0 + int(is, C_INT)
        end do

        if (has_terminal) then
            print *, "passive scalars:", sc%n
            do is = 1, int(sc%n)
                print '(a,i0,a,a,a,f8.4,a,es10.3)', "   s", is, " '", trim(sc%name(is)), &
                    "'  pr =", sc%pr(is), "  source =", sc%source(is)
                if (wall_function) print '(a,f8.4,a,f8.3)', &
                    "      thermal wall function: P =", sc%wfP(is), "  y+_T =", sc%wfYpt(is)
            end do
        end if
    end subroutine init_scalar

    ! Set-if-unset with the A0 contradiction rule: on a strict face an
    ! explicit _type key that disagrees with the declared patch is a hard
    ! config error; on a permissive one (a wall) it simply wins.
    subroutine resolve_scalar_row(sc, is, face_id, want, value, terminal, strict)
        type(scalar_type), intent(inout) :: sc
        integer, intent(in) :: is, face_id
        integer(C_INT), intent(in) :: want
        real(C_DOUBLE), intent(in) :: value
        logical, intent(in) :: terminal
        logical, intent(in) :: strict

        if (sc%bcTypeSet(is,face_id)) then
            if (strict .and. sc%bcType(is,face_id) /= want) then
                if (terminal) print '(a,i0,a,i0)', &
                    " error: explicit [scalar.N] _type key contradicts the declared " // &
                    "patch type on face ", face_id, ", scalar ", is
                error stop "[scalar] BC type contradicts the declared patch type"
            end if
        else
            sc%bcType(is,face_id) = want
        end if
        if (.not. sc%bcValueSet(is,face_id)) sc%bcValue(is,face_id) = value
    end subroutine resolve_scalar_row

    real(C_DOUBLE) function inv_delta(d) result(inv)
        real(C_DOUBLE), intent(in) :: d

        inv = 1.0d0/sign(max(abs(d), 1.0d-300), d)
    end function inv_delta

    ! Initial condition in the q scalar slots (host side, before the device
    ! map). Always called, including on restart: read_field then overwrites
    ! the scalars whose dataset the file carries, so an absent dataset (a
    ! newly added scalar, an older restart) keeps these values -- the RANS
    ! named-scalar restart precedent.
    subroutine init_scalar_fields(sc, dns, blk)
        type(scalar_type), intent(in) :: sc
        type(dns_type), intent(in) :: dns
        type(block_set_type), intent(inout) :: blk

        integer :: is, i, j, k, b, nx, ny, nz
        real(C_DOUBLE) :: s0, s1, y

        if (.not. scalars_enabled(sc)) return

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))

        do is = 1, int(sc%n)
            blk%q(:,:,:,VAR_S0+is,:) = sc%initValue(is)
            if (sc%initProfile(is) /= SC_INIT_LINEAR_Y) cycle
            ! linear_y: the straight line between the y faces' Dirichlet
            ! values (the conduction/heated-channel convenience IC).
            s0 = sc%bcValue(is, boundary_face_id(2,0))
            s1 = sc%bcValue(is, boundary_face_id(2,1))
            do b = 1, int(blk%nBlocks)
                do k = 0, nz+1
                    do j = 0, ny+1
                        y = blk%y(j,VAR_P,b)
                        do i = 0, nx+1
                            blk%q(i,j,k,VAR_S0+is,b) = s0 + (s1 - s0)*y/dns%leng(2)
                        end do
                    end do
                end do
            end do
        end do
    end subroutine init_scalar_fields

    !--------------------------------------------------------------------
    ! Conjugate interface geometry (increment C1)
    !--------------------------------------------------------------------

    ! Build the signed distance phi = +dwall (fluid) / -dwall (solid) at
    ! every cell centre INCLUDING THE GHOSTS, which is the entire geometric
    ! input of the conjugate scheme.
    !
    ! Nothing here is new computation: the magnitude is the RANS wall
    ! distance, use-associated from its two existing producers
    ! (fill_body_distance_analytic for the analytic indicator, the case
    ! file's dwall_blocks tiles for the file path), and the sign is the
    ! existing cell-centred IBM marker. That is the whole point of the
    ! baseline -- no new dataset, no moby_prepare change, no case-file format
    ! change (strategy doc Section 8 / Section 11).
    !
    ! THE MARKER: set_ibm_coeff writes SOLID/Re at every cell whose CENTRE is
    ! inside the body and a graded coefficient only at fluid-centred cells,
    ! so `|coef(VAR_P)| > SOLID_FACE_THRESHOLD` is exactly "this cell centre
    ! is in the solid" -- on both geometry paths, and ghost-inclusive on
    ! both, which is what lets a cut face on a block boundary see both signs.
    !
    ! HOST CODE writing a device-mapped array: the caller must push phi (and
    ! blk%q, if init_scalar_solid_fields ran) to the device afterwards. It
    ! must also have refreshed the HOST copy of ibm%coef first -- the
    ! analytic coefficients are computed on the device.
    subroutine init_scalar_conjugate(sc, dns, blk, bc, ibm, c)
        type(scalar_type), intent(inout) :: sc
        type(dns_type), intent(in) :: dns
        type(block_set_type), intent(in) :: blk
        type(boundary_type), intent(in) :: bc
        type(ibm_type), intent(in) :: ibm
        type(comm_type), intent(in) :: c

        real(C_DOUBLE), parameter :: SOLID_FACE_THRESHOLD = 1.0d20
        real(C_DOUBLE), allocatable :: dwall(:,:,:,:)
        integer :: i, j, k, b, nx, ny, nz
        logical :: found
        real(C_DOUBLE) :: counts(2)

        if (.not. scalar_conjugate_enabled(sc)) return

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))
        allocate(dwall(0:nx+1,0:ny+1,0:nz+1,blk%nBlocks))

        if (len_trim(dns%ibm_coeff_file) > 0) then
            call read_dwall_blocks(dwall, found, dns, blk, c%has_terminal)
            if (.not. found) then
                if (c%has_terminal) print *, "error: ibm_wall = conjugate needs the wall", &
                    " distance, but the case file has no dwall_blocks; re-run", &
                    " moby_prepare (it writes them by default): ", trim(dns%ibm_coeff_file)
                error stop "conjugate needs dwall_blocks"
            end if
        else
            call fill_body_distance_analytic(dwall, dns, blk, bc, ibm, &
                c%has_terminal, isInBody)
        end if

        ! A solid cell is given a STRICTLY negative phi, so that `phi < 0` is
        ! an exact material test even for a cell centre lying exactly on the
        ! surface (where dwall is 0). The kernel's cut test is that material
        ! disagreement, which is why the degenerate case cannot leak into it.
        counts = 0.0d0
        do b = 1, int(blk%nBlocks)
        do k = 0, nz+1
            do j = 0, ny+1
                do i = 0, nx+1
                    if (abs(ibm%coef(i,j,k,VAR_P,b)) > SOLID_FACE_THRESHOLD) then
                        sc%phi(i,j,k,b) = -max(dwall(i,j,k,b), CONJ_MIN_DISTANCE)
                        if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny &
                            .and. k >= 1 .and. k <= nz) counts(1) = counts(1) + 1.0d0
                    else
                        sc%phi(i,j,k,b) = dwall(i,j,k,b)
                    end if
                end do
            end do
        end do
        end do
        deallocate(dwall)

        counts(2) = real(count_cut_faces(sc, blk), C_DOUBLE)
        call comm_allreduce_sum(c, counts)
        if (c%has_terminal) print '(a,i0,a,i0,a)', &
            " conjugate interface: ", nint(counts(1)), " solid cells, ", &
            nint(counts(2)), " cut faces"

        call check_conjugate_refinement(sc, blk, c)
    end subroutine init_scalar_conjugate

    ! Interior cut faces of one rank (the two cell centres in different
    ! materials), for the init report. Counted on the low side of each cell
    ! so no face is counted twice.
    integer function count_cut_faces(sc, blk) result(n)
        type(scalar_type), intent(in) :: sc
        type(block_set_type), intent(in) :: blk

        integer :: i, j, k, b, nx, ny, nz
        logical :: s0

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))
        n = 0
        do b = 1, int(blk%nBlocks)
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    s0 = sc%phi(i,j,k,b) < 0.0d0
                    if ((sc%phi(i-1,j,k,b) < 0.0d0) .neqv. s0) n = n + 1
                    if ((sc%phi(i,j-1,k,b) < 0.0d0) .neqv. s0) n = n + 1
                    if ((sc%phi(i,j,k-1,b) < 0.0d0) .neqv. s0) n = n + 1
                end do
            end do
        end do
        end do
    end function count_cut_faces

    ! LANDMINE (strategy doc Section 12): refine_body's one-block
    ! 26-neighbour buffer is what keeps cut cells away from 2:1 block
    ! interfaces. The cut-face coefficient is a SAME-LEVEL two-point arm; on
    ! a coarse/fine block face the halo value is a restriction or a prolong
    ! of the other level, so neither phi nor the flux means what the scheme
    ! assumes -- and phi's own halo there is evaluated at THIS block's level's
    ! ghost coordinate, which is not a neighbouring cell centre at all. Away
    ! from the interface that is harmless (only the SIGN is read, and a point
    ! well inside one material has the right one), which is precisely what
    ! this check enforces. The buffer also keeps curvature*h small, which is
    ! what makes the level-set weight w trustworthy. CHECK it -- do not
    ! assume it.
    subroutine check_conjugate_refinement(sc, blk, c)
        type(scalar_type), intent(in) :: sc
        type(block_set_type), intent(in) :: blk
        type(comm_type), intent(in) :: c

        integer :: i, j, k, b, d, nx, ny, nz
        integer(C_INT) :: kind
        real(C_DOUBLE) :: bad(1)

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))
        bad = 0.0d0

        do b = 1, int(blk%nBlocks)
            do d = 1, 3
                kind = blk%physLow(d,b)
                if (kind == FACE_COARSE .or. kind == FACE_FINE) then
                    select case (d)
                    case (1)
                        do k = 1, nz
                            do j = 1, ny
                                if ((sc%phi(0,j,k,b) < 0.0d0) .neqv. &
                                    (sc%phi(1,j,k,b) < 0.0d0)) bad(1) = bad(1) + 1.0d0
                            end do
                        end do
                    case (2)
                        do k = 1, nz
                            do i = 1, nx
                                if ((sc%phi(i,0,k,b) < 0.0d0) .neqv. &
                                    (sc%phi(i,1,k,b) < 0.0d0)) bad(1) = bad(1) + 1.0d0
                            end do
                        end do
                    case default
                        do j = 1, ny
                            do i = 1, nx
                                if ((sc%phi(i,j,0,b) < 0.0d0) .neqv. &
                                    (sc%phi(i,j,1,b) < 0.0d0)) bad(1) = bad(1) + 1.0d0
                            end do
                        end do
                    end select
                end if
                kind = blk%physHigh(d,b)
                if (kind == FACE_COARSE .or. kind == FACE_FINE) then
                    select case (d)
                    case (1)
                        do k = 1, nz
                            do j = 1, ny
                                if ((sc%phi(nx,j,k,b) < 0.0d0) .neqv. &
                                    (sc%phi(nx+1,j,k,b) < 0.0d0)) bad(1) = bad(1) + 1.0d0
                            end do
                        end do
                    case (2)
                        do k = 1, nz
                            do i = 1, nx
                                if ((sc%phi(i,ny,k,b) < 0.0d0) .neqv. &
                                    (sc%phi(i,ny+1,k,b) < 0.0d0)) bad(1) = bad(1) + 1.0d0
                            end do
                        end do
                    case default
                        do j = 1, ny
                            do i = 1, nx
                                if ((sc%phi(i,j,nz,b) < 0.0d0) .neqv. &
                                    (sc%phi(i,j,nz+1,b) < 0.0d0)) bad(1) = bad(1) + 1.0d0
                            end do
                        end do
                    end select
                end if
            end do
        end do

        call comm_allreduce_sum(c, bad)
        if (bad(1) <= 0.0d0) return
        if (c%has_terminal) print '(a,i0,a)', &
            " error: ibm_wall = conjugate: ", nint(bad(1)), " interface cut faces sit on a" // &
            " 2:1 block face. The cut-face coefficient is a same-level arm and cannot" // &
            " read across a refinement interface. Use [blocks] refine_body = true (its" // &
            " one-block 26-neighbour buffer keeps the surface strictly inside the finest" // &
            " level), or refine the surface region explicitly."
        error stop "conjugate interface crosses a 2:1 block face"
    end subroutine check_conjugate_refinement

    ! The solid's own initial temperature ([scalar.N] solid_init, default =
    ! `initial`). COLD START ONLY: on a restart the solid field is part of
    ! the saved state and comes from the file like every other cell.
    ! Ghost-inclusive, so the halo layer is consistent before the first
    ! exchange. HOST code writing a device-mapped blk%q -- the caller pushes.
    subroutine init_scalar_solid_fields(sc, blk)
        type(scalar_type), intent(in) :: sc
        type(block_set_type), intent(inout) :: blk

        integer :: is, i, j, k, b, nx, ny, nz

        if (.not. scalar_conjugate_enabled(sc)) return

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))

        do is = 1, int(sc%n)
            if (sc%ibmMode(is) /= SC_IBM_CONJUGATE) cycle
            do b = 1, int(blk%nBlocks)
            do k = 0, nz+1
                do j = 0, ny+1
                    do i = 0, nx+1
                        if (sc%phi(i,j,k,b) < 0.0d0) &
                            blk%q(i,j,k,VAR_S0+is,b) = sc%solidInit(is)
                    end do
                end do
            end do
            end do
        end do
    end subroutine init_scalar_solid_fields

    ! The conjugate body's contribution to the explicit diffusive (Peclet)
    ! rate, max'ed into dns%peclet_rate by the caller once the interface
    ! exists.
    !
    ! TWO THINGS GATE 1 TAUGHT, both of which cost a NaN before they were
    ! understood, and neither of which a per-material bound can see:
    !
    ! 1. THE LIMIT IS NOT max OVER MATERIALS OF alpha = kappa/C. A cut face
    !    carries k_face up to max(kappa_L, kappa_R) -- that bound IS the
    !    strategy doc Section 7 argument that the interface is not stiff --
    !    but it feeds the cell on the OTHER side, whose capacity belongs to
    !    the other material. So a fluid cell against a kappa_s = 1000 solid
    !    sees 1000x the fluid rate even when alpha_s = alpha_f exactly.
    !    Building the rate from the ACTUAL face coefficients is both correct
    !    and far less conservative than max(kappa)/min(C): at w = 1/2 and
    !    kappa_s = 1000 the true penalty is 2x, not 1000x.
    !
    ! 2. THE CONVENTION MUST BE THE FULL DIAGONAL, AND A CUT CELL NEEDS THE
    !    GERSHGORIN FACTOR THE UNIFORM INTERIOR NEVER PAYS.
    !    precompute_peclet_rate reports alpha/h^2 per direction, i.e. HALF of
    !    one direction's diagonal; at a cut cell the two faces of a direction
    !    are not equal, so neither their max nor their mean is half the
    !    diagonal. Summing all six face terms and dividing by 6 reproduces
    !    alpha/h^2 exactly for a uniform isotropic cell, so pecletmax keeps
    !    its meaning -- that is the `diag/6` below.
    !    But the SPECTRAL RADIUS is bounded by 2*A_ii/C_i (Gershgorin, the
    !    operator has zero row sum), not by A_ii/C_i, and the existing
    !    convention is a factor ~1.9 short of that: uniform runs survive only
    !    because their extreme modes are never excited. An isolated cut cell
    !    is different -- its row is strongly ASYMMETRIC (a large k_face into a
    !    neighbour of the other material's capacity), so the worst mode IS
    !    local and IS attained, and this RK3's real-axis limit is 2.5.
    !    MEASURED on gate 1: (kappa_s, C_s) = (0.01, 0.01) at w = 0.95 blows
    !    up at pecletmax 0.3 and is stable at 0.2; (1000, 1000) at w = 0.80
    !    blows up at 0.4 and is stable at 0.2. Doubling the rate at cut cells
    !    -- diag/3 instead of diag/6 -- makes the nominal 0.4 behave as the
    !    measured-stable 0.2 in BOTH, which is exactly the factor Gershgorin
    !    predicts. Only interface cells pay it; the bulk of a conjugate run
    !    keeps today's step.
    real(C_DOUBLE) function scalar_conjugate_peclet_rate(sc, blk, dns, c) result(rate)
        type(scalar_type), intent(in) :: sc
        type(block_set_type), intent(in) :: blk
        type(dns_type), intent(in) :: dns
        type(comm_type), intent(in) :: c

        integer :: i, j, k, b, is, nx, ny, nz
        real(C_DOUBLE) :: dm, ks, rc, cc, phc, diag, share, r(1)
        logical :: solc, cut

        rate = 0.0d0
        if (.not. scalar_conjugate_enabled(sc)) return

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))
        r = 0.0d0

        do is = 1, int(sc%n)
            if (sc%ibmMode(is) /= SC_IBM_CONJUGATE) cycle
            dm = 1.0d0/(dns%re*sc%pr(is))
            ks = sc%solidK(is)
            rc = sc%contactR(is)
            do b = 1, int(blk%nBlocks)
            do k = 1, nz
                do j = 1, ny
                    do i = 1, nx
                        phc = sc%phi(i,j,k,b)
                        solc = phc < 0.0d0
                        cc = merge(sc%solidC(is), 1.0d0, solc)
                        cut = ((sc%phi(i-1,j,k,b) < 0.0d0) .neqv. solc) &
                         .or. ((sc%phi(i+1,j,k,b) < 0.0d0) .neqv. solc) &
                         .or. ((sc%phi(i,j-1,k,b) < 0.0d0) .neqv. solc) &
                         .or. ((sc%phi(i,j+1,k,b) < 0.0d0) .neqv. solc) &
                         .or. ((sc%phi(i,j,k-1,b) < 0.0d0) .neqv. solc) &
                         .or. ((sc%phi(i,j,k+1,b) < 0.0d0) .neqv. solc)
                        share = merge(3.0d0, 6.0d0, cut)
                        diag = (conjugate_face_diffusivity(dm, sc%phi(i-1,j,k,b), phc, &
                                    ks, rc, sc%invDx(i,b))*sc%invDx(i,b) &
                              + conjugate_face_diffusivity(dm, phc, sc%phi(i+1,j,k,b), &
                                    ks, rc, sc%invDx(i+1,b))*sc%invDx(i+1,b)) &
                                *blk%d1x(i,VAR_P,b) &
                             + (conjugate_face_diffusivity(dm, sc%phi(i,j-1,k,b), phc, &
                                    ks, rc, sc%invDy(j,b))*sc%invDy(j,b) &
                              + conjugate_face_diffusivity(dm, phc, sc%phi(i,j+1,k,b), &
                                    ks, rc, sc%invDy(j+1,b))*sc%invDy(j+1,b)) &
                                *blk%d1y(j,VAR_P,b) &
                             + (conjugate_face_diffusivity(dm, sc%phi(i,j,k-1,b), phc, &
                                    ks, rc, sc%invDz(k,b))*sc%invDz(k,b) &
                              + conjugate_face_diffusivity(dm, phc, sc%phi(i,j,k+1,b), &
                                    ks, rc, sc%invDz(k+1,b))*sc%invDz(k+1,b)) &
                                *blk%d1z(k,VAR_P,b)
                        r(1) = max(r(1), diag/(share*cc))
                    end do
                end do
            end do
            end do
        end do

        call comm_allreduce_max(c, r)
        rate = r(1)
    end function scalar_conjugate_peclet_rate

    ! Push the host-filled signed distance to the device (a no-op on the CPU
    ! build and without a conjugate scalar).
    subroutine scalar_conjugate_to_device(sc)
        type(scalar_type), intent(in) :: sc

        if (.not. scalar_conjugate_enabled(sc)) return

#ifdef USE_OPENMP_OFFLOAD
        !$omp target update to(sc%phi)
#endif
    end subroutine scalar_conjugate_to_device

    subroutine enter_scalar_data(sc)
        type(scalar_type), intent(inout) :: sc

        if (.not. scalars_enabled(sc)) return

#ifdef USE_OPENMP_OFFLOAD
        !$omp target enter data map(to: sc)
        !$omp target enter data map(to: sc%pr, sc%prt, sc%prtModel, sc%source, &
        !$omp& sc%initValue, sc%ibmValue, sc%inlet, sc%ibmMode, sc%initProfile, &
        !$omp& sc%bcType, sc%bcValue, sc%invDx, sc%invDy, sc%invDz, sc%nutNone, &
        !$omp& sc%wfP, sc%wfYpt, sc%wfYplus, sc%phi, sc%solidK, sc%solidC, &
        !$omp& sc%solidSource, sc%contactR)
#endif
    end subroutine enter_scalar_data

    subroutine exit_scalar_data(sc)
        type(scalar_type), intent(inout) :: sc

        if (.not. scalars_enabled(sc)) return

#ifdef USE_OPENMP_OFFLOAD
        !$omp target exit data map(delete: sc%pr, sc%prt, sc%prtModel, sc%source, &
        !$omp& sc%initValue, sc%ibmValue, sc%inlet, sc%ibmMode, sc%initProfile, &
        !$omp& sc%bcType, sc%bcValue, sc%invDx, sc%invDy, sc%invDz, sc%nutNone, &
        !$omp& sc%wfP, sc%wfYpt, sc%wfYplus, sc%phi, sc%solidK, sc%solidC, &
        !$omp& sc%solidSource, sc%contactR)
        !$omp target exit data map(delete: sc)
#endif
    end subroutine exit_scalar_data

    subroutine deallocate_config(sc)
        type(scalar_type), intent(inout) :: sc

        if (allocated(sc%pr)) deallocate(sc%pr)
        if (allocated(sc%prt)) deallocate(sc%prt)
        if (allocated(sc%prtModel)) deallocate(sc%prtModel)
        if (allocated(sc%source)) deallocate(sc%source)
        if (allocated(sc%initValue)) deallocate(sc%initValue)
        if (allocated(sc%ibmValue)) deallocate(sc%ibmValue)
        if (allocated(sc%inlet)) deallocate(sc%inlet)
        if (allocated(sc%ibmMode)) deallocate(sc%ibmMode)
        if (allocated(sc%initProfile)) deallocate(sc%initProfile)
        if (allocated(sc%solidK)) deallocate(sc%solidK)
        if (allocated(sc%solidC)) deallocate(sc%solidC)
        if (allocated(sc%solidInit)) deallocate(sc%solidInit)
        if (allocated(sc%solidSource)) deallocate(sc%solidSource)
        if (allocated(sc%contactR)) deallocate(sc%contactR)
        if (allocated(sc%solidKeySet)) deallocate(sc%solidKeySet)
        if (allocated(sc%solidInitSet)) deallocate(sc%solidInitSet)
        if (allocated(sc%ibmValueSet)) deallocate(sc%ibmValueSet)
        if (allocated(sc%bcType)) deallocate(sc%bcType)
        if (allocated(sc%bcValue)) deallocate(sc%bcValue)
        if (allocated(sc%bcTypeSet)) deallocate(sc%bcTypeSet)
        if (allocated(sc%bcValueSet)) deallocate(sc%bcValueSet)
        if (allocated(sc%name)) deallocate(sc%name)
        if (allocated(sc%sectionSeen)) deallocate(sc%sectionSeen)
    end subroutine deallocate_config

    subroutine destroy_scalar(sc)
        type(scalar_type), intent(inout) :: sc

        call deallocate_config(sc)
        if (allocated(sc%varList)) deallocate(sc%varList)
        if (allocated(sc%nutNone)) deallocate(sc%nutNone)
        if (allocated(sc%wfP)) deallocate(sc%wfP)
        if (allocated(sc%wfYpt)) deallocate(sc%wfYpt)
        if (allocated(sc%wfYplus)) deallocate(sc%wfYplus)
        if (allocated(sc%invDx)) deallocate(sc%invDx)
        if (allocated(sc%invDy)) deallocate(sc%invDy)
        if (allocated(sc%invDz)) deallocate(sc%invDz)
        if (allocated(sc%phi)) deallocate(sc%phi)
        sc%n = 0_C_INT
        sc%nConjugate = 0_C_INT
    end subroutine destroy_scalar

    !--------------------------------------------------------------------
    ! Transport
    !--------------------------------------------------------------------

    ! One RK3 substage of every scalar, in ONE fused kernel (collapse(4) over
    ! b,k,j,i with an inner scalar loop, so the GPU sees a single launch
    ! regardless of the scalar count). Reads the START-of-substage q and
    ! writes only the qs/oldrhs scalar slots.
    !
    ! Convection: 2nd-order central in DIVERGENCE form on the p cell's own six
    ! face velocities -- no velocity interpolation is needed on the staggered
    ! mesh, and the flux form telescopes, so sum(s dV) changes only by the
    ! boundary flux.
    !
    ! [flow] convection = skew subtracts s*(div u)|stencil built from the SAME
    ! face velocities (docs/next_session_scalar.md Section 2). NOTE what that
    ! buys, since it differs from the momentum kernel's skew term: the FULL
    ! subtraction is the ADVECTIVE form u.grad s, which preserves a UNIFORM
    ! scalar exactly for ANY advecting field (div u never enters), at the cost
    ! of the exact global conservation the divergence form has. Subtracting
    ! HALF would instead be the skew-symmetric form (the momentum kernel's
    ! choice, energy-neutral in sum s^2 for any advecting field) but leaves
    ! 1/2 s div u on a uniform field. The projection makes div u small either
    ! way; the trade is uniform-preservation vs. sum s^2 neutrality.
    !
    ! There is deliberately NO upwind option (module header of rans.f90 and
    ! docs/next_session_scalar.md Section 2): central/skew only until the
    ! shared TVD increment lands.
    !
    ! Diffusion: face-flux form on the module's own inverse centre-to-centre
    ! distances, with the face diffusivity
    !
    !   D_face = 1/(Re Pr) + 1/2 (nu_t,L + nu_t,R)/Pr_t(face)
    !
    ! reading turb%nut -- the ONE blended eddy viscosity that LES (WALE /
    ! Smagorinsky), RANS (SST) and IDDES all write, so a single code path
    ! serves all three models (increment S2). Pr_t is the per-scalar constant
    ! or the Kays-Crawford correlation, selected by [scalar.N] prt_model.
    !
    ! Immersed body (increment S3), per scalar via [scalar.N] ibm_wall:
    !   * dirichlet -- volume penalization toward [scalar.N] ibm_value with
    !     the implicit factor mu_s = 1/(1 + dt_gamma coef_p/Pr) formed inline
    !     (ibm.f90's implicit form, so no dt restriction; no mu array per
    !     scalar). The stored coefficient carries the 1/Re scaling, so
    !     coef_p/Pr is exactly the scalar's own 1/(Re Pr) -- ONE cell-centred
    !     coefficient array serves every scalar, each with its own Pr. Inside
    !     the body coef_p = SOLID/Re, so mu_s -> 0 and the cell holds
    !     ibm_value. Diffusive fluxes are NOT masked: the solid cell holds
    !     the wall value and delivers the flux. First-order/staircase, like
    !     the velocity penalization.
    !   * adiabatic -- no penalization; instead BOTH the convective and the
    !     diffusive flux are masked on each of the six p-cell faces whose
    !     staggered velocity coefficient is solid (the solw/sole/... test
    !     rans.f90 uses for k). The mask is symmetric across a face, so the
    !     flux form still telescopes and the fluid conserves sum(s dV)
    !     exactly.
    !   * conjugate (increment C1, docs/next_session_conjugate.md) -- the
    !     solid is a REAL unknown. No penalization either; instead
    !       - the face diffusivity becomes conjugate_face_diffusivity, which
    !         is kappa*dm within one material and the distance-weighted
    !         harmonic mean across a cut face. Every UNCUT FLUID face keeps
    !         today's expression verbatim, which is why the feature is
    !         untouched away from the body;
    !       - CONVECTION is hard-masked on every face whose staggered
    !         velocity node is solid AND on every cut face (u.n = 0 on the
    !         interface exactly, so this is consistent to second order).
    !         Relying on the penalised velocity being "small" is not
    !         acceptable here: the solid now carries a genuine temperature
    !         field that must not be advected. The skew term's divergence is
    !         then built from the SAME masked face velocities, so a uniform
    !         scalar is still preserved exactly;
    !       - the flux divergence is divided by the local capacity C
    !         (1 in the fluid, C_s = (rho c)_s/(rho c)_f in the solid), and
    !         the solid takes its own volumetric source;
    !       - nu_t enters NEITHER the solid NOR a cut face.
    !     The interface term is written as a FACE FLUX, never a cell source,
    !     which is what keeps sum(C s dV) conserved to round-off whatever the
    !     scheme's accuracy -- and what lets the C2 tangential correction and
    !     the C4 wedge model drop in without touching anything around them.
    ! With no immersed body (useIbm off) none of the branches is taken and
    ! the arithmetic is the S2 kernel's, byte for byte.
    !
    ! The nut ghost cells follow the momentum SGS convention: the halo
    ! exchange fills them across block and rank boundaries, and they stay at
    ! zero on a physical face (a resolved wall has nu_t -> 0 there anyway).
    !
    ! Under [rans] wall_treatment = wall_function (increment S5a) the face
    ! diffusivity at WALL CELLS comes from the thermal wall function instead
    ! (Kader/Jayatilleke, above) -- the momentum wall function's wall-cell
    ! nu_t is deliberately NOT reused there, because it carries the momentum
    ! sublayer's resistance, not the thermal one. Every other cell keeps
    ! nu_t/Pr_t, and with wall functions off the branch is not entered.
    subroutine scalar_transport(sc, blk, dns, turb, coef, dt_alpha, dt_beta, dt_gamma)
        type(scalar_type), intent(in) :: sc
        type(block_set_type), intent(inout) :: blk
        type(dns_type), intent(in) :: dns
        type(turb_type), intent(in) :: turb
        ! ibm%coef: with scalars on it always carries the cell-centred VAR_P
        ! column (init_ibm), so there is no dummy-array case here -- the
        ! body is switched off by dns%ibm_enabled alone.
        real(C_DOUBLE), intent(in) :: coef(0:,0:,0:,1:,1:)
        real(C_DOUBLE), intent(in) :: dt_alpha, dt_beta, dt_gamma

        if (.not. scalars_enabled(sc)) return

        ! With no turbulence model turb%nut does not exist; hand the kernel
        ! the 1-cell dummy and switch the eddy term off, which makes the
        ! arithmetic byte-identical to the molecular-only (S1) kernel.
        if (turbulence_is_enabled(turb) .and. allocated(turb%nut)) then
            call scalar_transport_kernel(sc, blk, dns, dt_alpha, dt_beta, dt_gamma, &
                turb%nut, .true., coef)
        else
            call scalar_transport_kernel(sc, blk, dns, dt_alpha, dt_beta, dt_gamma, &
                sc%nutNone, .false., coef)
        end if
    end subroutine scalar_transport

    subroutine scalar_transport_kernel(sc, blk, dns, dt_alpha, dt_beta, dt_gamma, &
            nut, useNut, coef)
        type(scalar_type), intent(in) :: sc
        type(block_set_type), intent(inout) :: blk
        type(dns_type), intent(in) :: dns
        real(C_DOUBLE), intent(in) :: dt_alpha, dt_beta, dt_gamma
        real(C_DOUBLE), intent(in) :: nut(0:,0:,0:,1:)
        logical, intent(in) :: useNut
        real(C_DOUBLE), intent(in) :: coef(0:,0:,0:,1:,1:)

        ! The solid-coefficient test the SGS and RANS kernels use (les.f90,
        ! rans.f90's SOLID_FACE_THRESHOLD): SOLID/Re is 1e30/Re.
        real(C_DOUBLE), parameter :: SOLID_FACE_THRESHOLD = 1.0d20

        integer :: i, j, k, b, is, nx, ny, nz, nBlocks, nScal, var, scr
        real(C_DOUBLE) :: ire, re, uw, ue, vs, vn, wb, wt, divu, divuse
        real(C_DOUBLE) :: s0, conv, diff, rhs, dm, fw, fe, ss, mus, ipr
        real(C_DOUBLE) :: ntw, nte, nts, ntn, ntb, ntt
        real(C_DOUBLE) :: dxw, dxe, dys, dyn, dzb, dzt
        ! Conjugate interface (C1): the six neighbour signed distances and
        ! this cell's, plus the per-scalar solid properties.
        real(C_DOUBLE) :: phc, phw, phe, phs, phn, phb, pht, ks, rc
        logical :: skew, useIbm, adiab, wallfn, anyConj, conjug, solc
        logical :: cutw, cute, cuts, cutn, cutb, cutt
        logical :: clw, cle, cls, cln, clb, clt
        logical :: solw, sole, sols, soln, solb, solt
        logical :: mw, me, ms, mn, mb, mt
        logical :: cmw, cme, cms, cmn, cmb, cmt

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))
        nBlocks = int(blk%nBlocks)
        nScal = int(sc%n)
        re = dns%re
        ire = 1.0d0/re
        skew = logical(dns%conv_skew)
        useIbm = logical(dns%ibm_enabled)
        ! S5a: the thermal wall function replaces the face eddy diffusivity
        ! AT WALL CELLS. Off (every resolved-wall and every non-RANS run) the
        ! branch is not entered and the S2/S3 arithmetic is reproduced byte
        ! for byte -- sc%wfYplus is then a 1-cell dummy, mapped but unread.
        wallfn = dns%rans_wall_treatment == 1_C_INT
        ! Conjugate branch (C1): off, sc%phi is a 1-cell dummy and nothing
        ! below reads it -- the S3 arithmetic is reproduced byte for byte.
        anyConj = sc%nConjugate > 0_C_INT

        !$omp target teams distribute parallel do collapse(4) &
        !$omp& map(to: ire, re, dt_alpha, dt_beta, dt_gamma, skew, useNut, useIbm, &
        !$omp& nScal, nx, ny, nz, anyConj, &
        !$omp& blk%q, blk%d1x, blk%d1y, blk%d1z, blk%physLow, blk%physHigh, nut, coef, &
        !$omp& sc%pr, sc%prt, sc%prtModel, sc%source, sc%invDx, sc%invDy, sc%invDz, &
        !$omp& sc%ibmMode, sc%ibmValue, sc%wfP, sc%wfYpt, sc%wfYplus, wallfn, &
        !$omp& sc%phi, sc%solidK, sc%solidC, sc%solidSource, sc%contactR) &
        !$omp& map(tofrom: blk%qs, blk%oldrhs) &
        !$omp& private(i,j,k,b,is,var,scr,uw,ue,vs,vn,wb,wt,divu,divuse, &
        !$omp& s0,conv,diff,rhs,dm,fw,fe,ss,mus,ipr,adiab,conjug,solc,ks,rc, &
        !$omp& phc,phw,phe,phs,phn,phb,pht, &
        !$omp& cutw,cute,cuts,cutn,cutb,cutt, &
        !$omp& clw,cle,cls,cln,clb,clt,solw,sole,sols,soln,solb,solt, &
        !$omp& mw,me,ms,mn,mb,mt,cmw,cme,cms,cmn,cmb,cmt, &
        !$omp& ntw,nte,nts,ntn,ntb,ntt,dxw,dxe,dys,dyn,dzb,dzt)
        do b = 1, nBlocks
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    ! The p cell's six face velocities: exactly the six
                    ! staggered components, no interpolation.
                    uw = blk%q(i,  j,  k,  VAR_U,b)
                    ue = blk%q(i+1,j,  k,  VAR_U,b)
                    vs = blk%q(i,  j,  k,  VAR_V,b)
                    vn = blk%q(i,  j+1,k,  VAR_V,b)
                    wb = blk%q(i,  j,  k,  VAR_W,b)
                    wt = blk%q(i,  j,  k+1,VAR_W,b)
                    divu = (ue - uw)*blk%d1x(i,VAR_P,b) &
                         + (vn - vs)*blk%d1y(j,VAR_P,b) &
                         + (wt - wb)*blk%d1z(k,VAR_P,b)

                    ! A FACE_CLOSED face borders a block removed inside the
                    ! immersed body: its halo is held at zero, so the
                    ! DIFFUSIVE flux must be masked (the convective one
                    ! already vanishes -- the face velocity is pinned to 0).
                    ! FACE_PHYS faces are NOT masked: their ghosts carry the
                    ! boundary condition.
                    clw = i == 1  .and. blk%physLow(1,b)  == FACE_CLOSED
                    cle = i == nx .and. blk%physHigh(1,b) == FACE_CLOSED
                    cls = j == 1  .and. blk%physLow(2,b)  == FACE_CLOSED
                    cln = j == ny .and. blk%physHigh(2,b) == FACE_CLOSED
                    clb = k == 1  .and. blk%physLow(3,b)  == FACE_CLOSED
                    clt = k == nz .and. blk%physHigh(3,b) == FACE_CLOSED

                    ! Solid staggered faces of this p cell (adiabatic body
                    ! mode only; scalar-independent, so hoisted out of the
                    ! scalar loop). rans.f90's test for k, verbatim.
                    solw = .false.; sole = .false.; sols = .false.
                    soln = .false.; solb = .false.; solt = .false.
                    if (useIbm) then
                        solw = abs(coef(i,  j,  k,  VAR_U,b)) > SOLID_FACE_THRESHOLD
                        sole = abs(coef(i+1,j,  k,  VAR_U,b)) > SOLID_FACE_THRESHOLD
                        sols = abs(coef(i,  j,  k,  VAR_V,b)) > SOLID_FACE_THRESHOLD
                        soln = abs(coef(i,  j+1,k,  VAR_V,b)) > SOLID_FACE_THRESHOLD
                        solb = abs(coef(i,  j,  k,  VAR_W,b)) > SOLID_FACE_THRESHOLD
                        solt = abs(coef(i,  j,  k+1,VAR_W,b)) > SOLID_FACE_THRESHOLD
                    end if

                    ! Conjugate interface geometry (C1): the six neighbour
                    ! signed distances and the cut flags. Pure geometry --
                    ! scalar-independent, so hoisted out of the scalar loop
                    ! (the level-set weight itself needs no normal at all,
                    ! the obliquity lemma; see conjugate_face_diffusivity).
                    phc = 0.0d0; phw = 0.0d0; phe = 0.0d0
                    phs = 0.0d0; phn = 0.0d0; phb = 0.0d0; pht = 0.0d0
                    solc = .false.
                    cutw = .false.; cute = .false.; cuts = .false.
                    cutn = .false.; cutb = .false.; cutt = .false.
                    if (anyConj) then
                        phc = sc%phi(i,  j,  k,  b)
                        phw = sc%phi(i-1,j,  k,  b)
                        phe = sc%phi(i+1,j,  k,  b)
                        phs = sc%phi(i,  j-1,k,  b)
                        phn = sc%phi(i,  j+1,k,  b)
                        phb = sc%phi(i,  j,  k-1,b)
                        pht = sc%phi(i,  j,  k+1,b)
                        solc = phc < 0.0d0
                        cutw = (phw < 0.0d0) .neqv. solc
                        cute = (phe < 0.0d0) .neqv. solc
                        cuts = (phs < 0.0d0) .neqv. solc
                        cutn = (phn < 0.0d0) .neqv. solc
                        cutb = (phb < 0.0d0) .neqv. solc
                        cutt = (pht < 0.0d0) .neqv. solc
                    end if

                    ! Face eddy viscosities: the same 1/2 (L + R) average the
                    ! momentum SGS correction uses for its face-normal terms.
                    ! Scalar-independent, so they are hoisted out of the
                    ! scalar loop; Pr_t is not (it may be per-scalar).
                    ntw = 0.0d0; nte = 0.0d0; nts = 0.0d0
                    ntn = 0.0d0; ntb = 0.0d0; ntt = 0.0d0
                    if (useNut) then
                        ntw = 0.5d0*(nut(i-1,j,k,b) + nut(i,j,k,b))
                        nte = 0.5d0*(nut(i,j,k,b) + nut(i+1,j,k,b))
                        nts = 0.5d0*(nut(i,j-1,k,b) + nut(i,j,k,b))
                        ntn = 0.5d0*(nut(i,j,k,b) + nut(i,j+1,k,b))
                        ntb = 0.5d0*(nut(i,j,k-1,b) + nut(i,j,k,b))
                        ntt = 0.5d0*(nut(i,j,k,b) + nut(i,j,k+1,b))
                    end if

                    do is = 1, nScal
                        var = VAR_S0 + is
                        scr = SCR_S0 + is
                        s0 = blk%q(i,j,k,var,b)

                        ! An adiabatic body seals this cell's solid faces:
                        ! convective AND diffusive flux masked. Every other
                        ! configuration leaves adiab .false., so the six
                        ! masks below collapse to the FACE_CLOSED flags and
                        ! the S2 expressions are reproduced exactly.
                        adiab = useIbm .and. sc%ibmMode(is) == SC_IBM_ADIABATIC
                        conjug = useIbm .and. sc%ibmMode(is) == SC_IBM_CONJUGATE
                        mw = clw .or. (adiab .and. solw)
                        me = cle .or. (adiab .and. sole)
                        ms = cls .or. (adiab .and. sols)
                        mn = cln .or. (adiab .and. soln)
                        mb = clb .or. (adiab .and. solb)
                        mt = clt .or. (adiab .and. solt)
                        ! The CONVECTIVE masks. They coincide with the
                        ! diffusive ones in every mode but conjugate, where
                        ! the solid conducts (diffusion is NOT masked) while
                        ! nothing at all is advected across a solid or cut
                        ! face.
                        cmw = mw .or. (conjug .and. (solw .or. cutw))
                        cme = me .or. (conjug .and. (sole .or. cute))
                        cms = ms .or. (conjug .and. (sols .or. cuts))
                        cmn = mn .or. (conjug .and. (soln .or. cutn))
                        cmb = mb .or. (conjug .and. (solb .or. cutb))
                        cmt = mt .or. (conjug .and. (solt .or. cutt))

                        if (adiab .or. conjug) then
                            conv = (merge(0.0d0, ue*0.5d0*(s0 + blk%q(i+1,j,k,var,b)), cme) &
                                  - merge(0.0d0, uw*0.5d0*(blk%q(i-1,j,k,var,b) + s0), cmw)) &
                                    *blk%d1x(i,VAR_P,b) &
                                 + (merge(0.0d0, vn*0.5d0*(s0 + blk%q(i,j+1,k,var,b)), cmn) &
                                  - merge(0.0d0, vs*0.5d0*(blk%q(i,j-1,k,var,b) + s0), cms)) &
                                    *blk%d1y(j,VAR_P,b) &
                                 + (merge(0.0d0, wt*0.5d0*(s0 + blk%q(i,j,k+1,var,b)), cmt) &
                                  - merge(0.0d0, wb*0.5d0*(blk%q(i,j,k-1,var,b) + s0), cmb)) &
                                    *blk%d1z(k,VAR_P,b)
                        else
                            conv = (ue*0.5d0*(s0 + blk%q(i+1,j,k,var,b)) &
                                  - uw*0.5d0*(blk%q(i-1,j,k,var,b) + s0))*blk%d1x(i,VAR_P,b) &
                                 + (vn*0.5d0*(s0 + blk%q(i,j+1,k,var,b)) &
                                  - vs*0.5d0*(blk%q(i,j-1,k,var,b) + s0))*blk%d1y(j,VAR_P,b) &
                                 + (wt*0.5d0*(s0 + blk%q(i,j,k+1,var,b)) &
                                  - wb*0.5d0*(blk%q(i,j,k-1,var,b) + s0))*blk%d1z(k,VAR_P,b)
                        end if
                        ! The skew term must subtract the divergence of the
                        ! SAME face velocities the convective flux used, or a
                        ! sealed solid cell would be "advected" by the
                        ! penalised velocity's residual divergence. With the
                        ! masked form a uniform scalar stays uniform exactly,
                        ! in the solid as in the fluid. (The adiabatic mode
                        ! keeps its shipped, gated arithmetic: divuse = divu.)
                        divuse = divu
                        if (skew .and. conjug) &
                            divuse = (merge(0.0d0, ue, cme) - merge(0.0d0, uw, cmw)) &
                                        *blk%d1x(i,VAR_P,b) &
                                   + (merge(0.0d0, vn, cmn) - merge(0.0d0, vs, cms)) &
                                        *blk%d1y(j,VAR_P,b) &
                                   + (merge(0.0d0, wt, cmt) - merge(0.0d0, wb, cmb)) &
                                        *blk%d1z(k,VAR_P,b)
                        if (skew) conv = conv - s0*divuse

                        ! Molecular diffusivity + the eddy part on each face.
                        ! With useNut off every face keeps dm exactly, so the
                        ! S1 arithmetic is reproduced bit-for-bit.
                        dm = ire/sc%pr(is)
                        dxw = dm; dxe = dm; dys = dm
                        dyn = dm; dzb = dm; dzt = dm
                        if (conjug) then
                            ! THE conjugate scheme: one face coefficient.
                            ! Within one material this is kappa*dm (kappa = 1
                            ! in the fluid, so a fluid-fluid face is today's
                            ! line unchanged); across a cut face it is the
                            ! distance-weighted harmonic mean built on the
                            ! level-set fraction of the arm.
                            ks = sc%solidK(is)
                            rc = sc%contactR(is)
                            dxw = conjugate_face_diffusivity(dm, phw, phc, ks, rc, sc%invDx(i,b))
                            dxe = conjugate_face_diffusivity(dm, phc, phe, ks, rc, sc%invDx(i+1,b))
                            dys = conjugate_face_diffusivity(dm, phs, phc, ks, rc, sc%invDy(j,b))
                            dyn = conjugate_face_diffusivity(dm, phc, phn, ks, rc, sc%invDy(j+1,b))
                            dzb = conjugate_face_diffusivity(dm, phb, phc, ks, rc, sc%invDz(k,b))
                            dzt = conjugate_face_diffusivity(dm, phc, pht, ks, rc, sc%invDz(k+1,b))
                            ! nu_t enters NEITHER the solid NOR a cut face
                            ! (strategy doc Section 12): at DNS resolution
                            ! nu_t -> 0 at the wall anyway, and ibm_aware
                            ! already zeroes it in solid cells, but the
                            ! interface coefficient is molecular by
                            ! construction and must not be polluted.
                            if (useNut .and. .not. solc) then
                                if (.not. cutw) dxw = dxw + eddy_diffusivity(ntw, sc%pr(is), &
                                    sc%prt(is), sc%prtModel(is), re)
                                if (.not. cute) dxe = dxe + eddy_diffusivity(nte, sc%pr(is), &
                                    sc%prt(is), sc%prtModel(is), re)
                                if (.not. cuts) dys = dys + eddy_diffusivity(nts, sc%pr(is), &
                                    sc%prt(is), sc%prtModel(is), re)
                                if (.not. cutn) dyn = dyn + eddy_diffusivity(ntn, sc%pr(is), &
                                    sc%prt(is), sc%prtModel(is), re)
                                if (.not. cutb) dzb = dzb + eddy_diffusivity(ntb, sc%pr(is), &
                                    sc%prt(is), sc%prtModel(is), re)
                                if (.not. cutt) dzt = dzt + eddy_diffusivity(ntt, sc%pr(is), &
                                    sc%prt(is), sc%prtModel(is), re)
                            end if
                        else if (useNut .and. wallfn) then
                            ! Thermal wall function: the face value is built
                            ! from the two CELL diffusivities, because a wall
                            ! cell's is not a function of its nu_t at all.
                            dxw = dm + wall_face_diffusivity(nut(i-1,j,k,b), nut(i,j,k,b), &
                                sc%wfYplus(i-1,j,k,b), sc%wfYplus(i,j,k,b), &
                                sc%pr(is), sc%prt(is), sc%prtModel(is), re, &
                                sc%wfP(is), sc%wfYpt(is))
                            dxe = dm + wall_face_diffusivity(nut(i,j,k,b), nut(i+1,j,k,b), &
                                sc%wfYplus(i,j,k,b), sc%wfYplus(i+1,j,k,b), &
                                sc%pr(is), sc%prt(is), sc%prtModel(is), re, &
                                sc%wfP(is), sc%wfYpt(is))
                            dys = dm + wall_face_diffusivity(nut(i,j-1,k,b), nut(i,j,k,b), &
                                sc%wfYplus(i,j-1,k,b), sc%wfYplus(i,j,k,b), &
                                sc%pr(is), sc%prt(is), sc%prtModel(is), re, &
                                sc%wfP(is), sc%wfYpt(is))
                            dyn = dm + wall_face_diffusivity(nut(i,j,k,b), nut(i,j+1,k,b), &
                                sc%wfYplus(i,j,k,b), sc%wfYplus(i,j+1,k,b), &
                                sc%pr(is), sc%prt(is), sc%prtModel(is), re, &
                                sc%wfP(is), sc%wfYpt(is))
                            dzb = dm + wall_face_diffusivity(nut(i,j,k-1,b), nut(i,j,k,b), &
                                sc%wfYplus(i,j,k-1,b), sc%wfYplus(i,j,k,b), &
                                sc%pr(is), sc%prt(is), sc%prtModel(is), re, &
                                sc%wfP(is), sc%wfYpt(is))
                            dzt = dm + wall_face_diffusivity(nut(i,j,k,b), nut(i,j,k+1,b), &
                                sc%wfYplus(i,j,k,b), sc%wfYplus(i,j,k+1,b), &
                                sc%pr(is), sc%prt(is), sc%prtModel(is), re, &
                                sc%wfP(is), sc%wfYpt(is))
                        else if (useNut) then
                            dxw = dm + eddy_diffusivity(ntw, sc%pr(is), sc%prt(is), &
                                sc%prtModel(is), re)
                            dxe = dm + eddy_diffusivity(nte, sc%pr(is), sc%prt(is), &
                                sc%prtModel(is), re)
                            dys = dm + eddy_diffusivity(nts, sc%pr(is), sc%prt(is), &
                                sc%prtModel(is), re)
                            dyn = dm + eddy_diffusivity(ntn, sc%pr(is), sc%prt(is), &
                                sc%prtModel(is), re)
                            dzb = dm + eddy_diffusivity(ntb, sc%pr(is), sc%prt(is), &
                                sc%prtModel(is), re)
                            dzt = dm + eddy_diffusivity(ntt, sc%pr(is), sc%prt(is), &
                                sc%prtModel(is), re)
                        end if

                        fw = merge(0.0d0, dxw*(s0 - blk%q(i-1,j,k,var,b)) &
                            *sc%invDx(i,b), mw)
                        fe = merge(0.0d0, dxe*(blk%q(i+1,j,k,var,b) - s0) &
                            *sc%invDx(i+1,b), me)
                        diff = (fe - fw)*blk%d1x(i,VAR_P,b)
                        fw = merge(0.0d0, dys*(s0 - blk%q(i,j-1,k,var,b)) &
                            *sc%invDy(j,b), ms)
                        fe = merge(0.0d0, dyn*(blk%q(i,j+1,k,var,b) - s0) &
                            *sc%invDy(j+1,b), mn)
                        diff = diff + (fe - fw)*blk%d1y(j,VAR_P,b)
                        fw = merge(0.0d0, dzb*(s0 - blk%q(i,j,k-1,var,b)) &
                            *sc%invDz(k,b), mb)
                        fe = merge(0.0d0, dzt*(blk%q(i,j,k+1,var,b) - s0) &
                            *sc%invDz(k+1,b), mt)
                        diff = diff + (fe - fw)*blk%d1z(k,VAR_P,b)

                        if (conjug) then
                            ! Divide the flux divergence by the LOCAL
                            ! volumetric capacity, and let the solid carry
                            ! its own volumetric source (Joule heating, ...).
                            ! C is pointwise here; the fluid-fraction-weighted
                            ! capacity of a cut cell is increment C3 (it
                            ! matters for transients, not at steady state).
                            rhs = (-conv + diff)/merge(sc%solidC(is), 1.0d0, solc) &
                                + merge(sc%solidSource(is), sc%source(is), solc)
                        else
                            rhs = -conv + diff + sc%source(is)
                        end if
                        ss = s0 + dt_alpha*rhs + dt_beta*blk%oldrhs(i,j,k,scr,b)
                        ! Implicit volume penalization toward the body value
                        ! (dirichlet mode). oldrhs keeps the UNpenalized rhs,
                        ! exactly as the momentum predictor does with mu.
                        ! The conjugate mode does NOT penalise at all: its
                        ! solid cells are unknowns, not boundary values.
                        if (useIbm .and. sc%ibmMode(is) == SC_IBM_DIRICHLET) then
                            ipr = 1.0d0/sc%pr(is)
                            mus = 1.0d0/(1.0d0 + dt_gamma*coef(i,j,k,VAR_P,b)*ipr)
                            ss = ss*mus + (1.0d0 - mus)*sc%ibmValue(is)
                        end if
                        blk%qs(i,j,k,scr,b) = ss
                        blk%oldrhs(i,j,k,scr,b) = rhs
                    end do
                end do
            end do
        end do
        end do
        !$omp end target teams distribute parallel do
    end subroutine scalar_transport_kernel

    ! Substage tail: qs -> q for the scalar slots, then the ghost/halo sync.
    ! Runs AFTER the pressure projection returns -- the scalar never enters
    ! the projection loop.
    subroutine scalar_finish(sc, blk, bc, c)
        type(scalar_type), intent(in) :: sc
        type(block_set_type), intent(inout) :: blk
        type(boundary_type), intent(in) :: bc
        type(comm_type), intent(inout) :: c

        integer :: i, j, k, b, is, nx, ny, nz, nBlocks, nScal

        if (.not. scalars_enabled(sc)) return

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))
        nBlocks = int(blk%nBlocks)
        nScal = int(sc%n)

        !$omp target teams distribute parallel do collapse(4) &
        !$omp& map(to: nScal, blk%qs) map(tofrom: blk%q) &
        !$omp& private(i,j,k,b,is)
        do b = 1, nBlocks
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    do is = 1, nScal
                        blk%q(i,j,k,VAR_S0+is,b) = blk%qs(i,j,k,SCR_S0+is,b)
                    end do
                end do
            end do
        end do
        end do
        !$omp end target teams distribute parallel do

        call scalar_sync(sc, blk, bc, c)
    end subroutine scalar_finish

    ! Physical ghosts + one batched halo exchange over all scalar variables
    ! (they ride the same entries, and the same message, as u,v,w,p).
    subroutine scalar_sync(sc, blk, bc, c)
        type(scalar_type), intent(in) :: sc
        type(block_set_type), intent(inout) :: blk
        type(boundary_type), intent(in) :: bc
        type(comm_type), intent(inout) :: c

        if (.not. scalars_enabled(sc)) return

        call apply_scalar_bc_q(blk, bc, int(sc%n), sc%bcType, sc%bcValue)
        call exchange_halos(c, blk, sc%varList)
    end subroutine scalar_sync

end module scalar
