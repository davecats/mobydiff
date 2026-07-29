module init
    use, intrinsic :: iso_c_binding
    implicit none

    integer(C_INT), parameter :: VAR_U = 1_C_INT
    integer(C_INT), parameter :: VAR_V = 2_C_INT
    integer(C_INT), parameter :: VAR_W = 3_C_INT
    integer(C_INT), parameter :: VAR_P = 4_C_INT
    integer(C_INT), parameter :: NVEL = 3_C_INT
    integer(C_INT), parameter :: NVAR = 4_C_INT
    integer(C_INT), parameter :: GRID_UNIFORM = 1_C_INT
    integer(C_INT), parameter :: GRID_COSINE  = 2_C_INT
    integer(C_INT), parameter :: GRID_TANH    = 3_C_INT
    integer(C_INT), parameter :: GRID_NATURAL = 4_C_INT
    integer(C_INT), parameter :: CFL_COURANT = 1_C_INT
    integer(C_INT), parameter :: CFL_PECLET  = 2_C_INT
    integer(C_INT), parameter :: NCFL = 2_C_INT

    ! Records which values the ini set explicitly, so config can win over the
    ! restart-file metadata for those fields (see read_restart_metadata). Lives
    ! here in the base module so both config (writer) and io (reader) can use it
    ! without an io<->config module cycle.
    type :: config_seen_type
        logical :: size(1:3) = .false.
        logical :: length(1:3) = .false.
        logical :: forcing(1:3) = .false.
        logical :: re = .false.
        logical :: dt = .false.
        logical :: nsteps = .false.
        logical :: t_final = .false.
        logical :: cflmax = .false.
        logical :: pecletmax = .false.
        logical :: dtmax = .false.
        logical :: pressure_niter = .false.
        logical :: pressure_sor = .false.
        logical :: turbulence_model = .false.
    end type config_seen_type

    ! Runtime/domain state shared by the solver modules.
    type :: dns_type
        integer(C_INT) :: globalSize(1:3) = 0_C_INT
        ! Per-direction local first index, last index, and count.
        integer(C_INT) :: localSize(1:3,0:2) = 0_C_INT
        integer(C_INT) :: step_current = 0_C_INT
        integer(C_INT) :: nsteps = 0_C_INT
        real(C_DOUBLE) :: leng(1:3) = 0.0d0
        real(C_DOUBLE) :: re = 0.0d0
        real(C_DOUBLE) :: dt = 0.0d0
        real(C_DOUBLE) :: t_final = 0.0d0
        real(C_DOUBLE) :: t_current = 0.0d0
        real(C_DOUBLE) :: cfl(1:NCFL) = 0.0d0
        real(C_DOUBLE) :: cflmax = 0.0d0
        real(C_DOUBLE) :: pecletmax = 0.0d0
        real(C_DOUBLE) :: peclet_rate = 0.0d0
        real(C_DOUBLE) :: dtmax = 0.0d0
        real(C_DOUBLE) :: forcing(1:3) = 0.0d0
        ! [flow] initial_u/v/w: uniform initial velocity (generic case).
        real(C_DOUBLE) :: initial_velocity(1:3) = 0.0d0
        real(C_DOUBLE) :: initial_noise = 0.0d0
        ! [flow] initial: analytic initial condition for the generic case --
        ! "uniform" (default), "beltrami" (3D ABC flow), "tgv" (2D Taylor-Green)
        ! or "tgv3d" (fully-3D product field, a manufactured momentum-operator
        ! test: every component varies in every direction, so the wall-normal
        ! velocity varies in the normal direction at all three interface
        ! orientations -- the property Beltrami lacks). tgv3d is NOT an exact
        ! NS solution; it is a manufactured momentum-operator test field.
        ! uniform/beltrami/tgv are exact incompressible-NS solutions used as gates.
        character(len=32) :: initial = "uniform"
        ! [blocks] nb: cubic block edge in cells; 0 = one block per rank box.
        integer(C_INT) :: block_nb = 0_C_INT
        ! [blocks] remove_solid: drop blocks buried inside the immersed body.
        logical(C_BOOL) :: block_remove_solid = .true.
        ! [blocks] refine = x0 x1 y0 y1 z0 z1 [level]: refine blocks
        ! intersecting this physical box to `level` (optional 7th value;
        ! default = refine_levels, i.e. the finest). Repeatable, up to 16
        ! boxes — nested per-level boxes build graded far-field decompositions
        ! (concentric wake-skewed rings, tutorials/naca B11).
        real(C_DOUBLE) :: block_refine_box(6,16) = 0.0d0
        integer(C_INT) :: block_refine_box_level(16) = -1_C_INT
        integer(C_INT) :: block_refine_nboxes = 0_C_INT
        ! [blocks] refine_levels: rounds of box refinement (max level).
        integer(C_INT) :: block_refine_levels = 1_C_INT
        ! [blocks] refine_dims = xyz (default octree) | xz (quadtree: blocks
        ! refine in x and z only; y keeps the single global node line at all
        ! levels). Stored as a per-direction mask (1 = the direction halves
        ! per level, 0 = fixed) so level scalings stay per-direction
        ! expressions; docs/next_session_refine2d.md.
        integer(C_INT) :: block_refine_mask(1:3) = 1_C_INT
        ! [blocks] refine_body: refine blocks whose dilated region meets the
        ! immersed surface to the finest level (+1 block buffer), and remove
        ! buried blocks at every level (analytic IBM).
        logical(C_BOOL) :: block_refine_body = .false.
        ! Keep buried refine_body blocks (mobygeom --keep-buried):
        ! LOAD-BEARING for penalization forces -- a removed solid core
        ! absorbs the body's pressure loading through its closed faces
        ! outside the coef bookkeeping (validation/naca0012/README.md).
        logical(C_BOOL) :: block_keep_buried = .false.
        logical(C_BOOL) :: ibm_enabled = .true.
        character(len=256) :: ibm_coeff_file = ""
        ! STL geometry (moby_prepare input only; the solver rejects it
        ! without a coeff_file). stl_file is repeatable -- one binary STL
        ! path per occurrence, so paths may contain spaces. The optional
        ! transform is v*scale + translate (mobygeom's convention).
        character(len=256) :: ibm_stl_file(8) = ""
        integer(C_INT) :: ibm_stl_count = 0_C_INT
        real(C_DOUBLE) :: ibm_stl_scale = 1.0d0
        real(C_DOUBLE) :: ibm_stl_translate(3) = 0.0d0
        ! [ibm] band_filter: optional 3-point low-pass on the predicted
        ! velocity in a thin near-body band (band_width cells, strength
        ! band_theta; theta = 1 annihilates the 2-cell mode per direction).
        ! Damps the cell-Reynolds parasite fan seeded by the staircase ring
        ! (validation/naca0012 README "LE fan root cause"/R1: the ring SEED
        ! controls the fan; refinement is the reference answer, the filter
        ! the production option). A separate correction pass — OFF means the
        ! pass is never called: bit-exact and zero cost by construction.
        logical(C_BOOL) :: ibm_band_filter = .false.
        integer(C_INT) :: ibm_band_width = 3_C_INT
        real(C_DOUBLE) :: ibm_band_theta = 0.5d0
        ! [force] optional spatially-varying volumetric body force
        ! (bodyforce.f90), added to the momentum predictor on top of the
        ! constant [flow] forcing_*. Disabled -> bit-exact with no force.
        !   type:    profile | file | custom
        !   profile: named analytic form (constant | sine)
        !   amp:     per-component amplitude
        !   wavenumber / dir: sine profile wavenumber and variation direction
        !   file:    HDF5 field (un/vn/wn = fx/fy/fz) for type = file
        logical(C_BOOL) :: force_enabled = .false.
        character(len=16) :: force_type = "profile"
        character(len=16) :: force_profile = "constant"
        real(C_DOUBLE) :: force_amp(1:3) = 0.0d0
        real(C_DOUBLE) :: force_wavenumber(1:3) = 0.0d0
        integer(C_INT) :: force_dir = 1_C_INT
        character(len=256) :: force_file = ""
        ! [force] type = trip: the Schlatter & Orlu (2012) random wall-normal
        ! trip forcing that triggers laminar->turbulent transition (see
        ! bodyforce.f90). x0/lx/ly place the streamwise+wall-normal Gaussian
        ! envelope; amp is the force amplitude; ts the temporal scale of the
        ! smooth-step random walk; nmodes the spanwise Fourier modes; seed the
        ! deterministic RNG seed.
        real(C_DOUBLE) :: trip_x0 = 0.0d0, trip_lx = 4.0d0, trip_ly = 1.0d0
        real(C_DOUBLE) :: trip_amp = 0.0d0, trip_ts = 4.0d0
        integer(C_INT) :: trip_nmodes = 16_C_INT, trip_seed = 1_C_INT
        ! [rans] section (rans.f90). In the T1 increment the section's mere
        ! presence builds the SST geometry state (wall distance + IBM wall
        ! cells) at init so it can be validated before any RANS transport
        ! exists; dump_geometry writes the diagnostic HDF5 file.
        logical(C_BOOL) :: rans_configured = .false.
        logical(C_BOOL) :: rans_dump_geometry = .false.
        ! Absolute tolerance of the geometry-agnostic analytic wall distance
        ! (walldist.f90 polish); gates sweep it to demonstrate convergence.
        real(C_DOUBLE) :: rans_dwall_tol = 1.0d-10
        ! T2 transport sub-model: 0 = none, 1 = k-omega SST. Wall treatment:
        ! 0 = resolved (y+_1 <~ 1), 1 = wall_function (phase T3). tu is the
        ! initial turbulence intensity in PERCENT of the local velocity;
        ! nut_ratio = nut/nu sets the initial omega = k/(nut_ratio*nu).
        integer(C_INT) :: rans_model = 0_C_INT
        integer(C_INT) :: rans_wall_treatment = 0_C_INT
        logical(C_BOOL) :: rans_transition = .false.
        real(C_DOUBLE) :: rans_tu = 5.0d0
        real(C_DOUBLE) :: rans_nut_ratio = 10.0d0
        ! [rans] ambient_sustain: add the Rumsey (2007) freestream-sustaining
        ! sources S_k = beta* k_inf omega_inf, S_omega = beta omega_inf^2 to
        ! the SST transport, making (k_inf, omega_inf) an exact fixed point:
        ! the ambient turbulence no longer free-decays on the way to the
        ! body, so a CONTROLLED tu arrives at ANY small nut_ratio (large
        ! nut_ratio ambients break the explicit eddy-diffusion dt bound in
        ! fine cells -- the C10 blow-up). k_inf/omega_inf from [rans] tu /
        ! nut_ratio and |initial_velocity| (must be nonzero when enabled).
        logical(C_BOOL) :: rans_ambient_sustain = .false.
        ! [rans] kpin_box / ktrip_box: OpenFOAM-matching forced-transition
        ! devices (their fvOptions scalarFixedValueConstraint +
        ! scalarSemiImplicitSource). Each kpin_box (x0 x1 y0 y1 z0 z1)
        ! pins k = 0 in its cells every substage (forced-laminar zone);
        ! each ktrip_box (x0 x1 y0 y1 z0 z1 rate) adds the constant
        ! volumetric source rate [k/time] to FLUID cells (trip strip).
        ! No boxes configured => the transport arithmetic is untouched.
        integer(C_INT) :: rans_n_kpin = 0_C_INT
        integer(C_INT) :: rans_n_ktrip = 0_C_INT
        real(C_DOUBLE) :: rans_kpin_box(6, 8) = 0.0d0
        real(C_DOUBLE) :: rans_ktrip_box(7, 8) = 0.0d0
        ! kpin_dwall = x0 x1 dmax: pin k = 0 in cells with x in [x0, x1]
        ! AND wall distance < dmax — a surface-following laminar band.
        ! Motivated by the C11 v1 blow-up (2026-07-22): a box pin spanning
        ! the freestream zeroes nut ON the 2:1 interface cascade and the
        ! undamped coarse-owns interface mode grows without bound; the
        ! dwall band keeps every interface in ambient-damped fluid.
        integer(C_INT) :: rans_n_kpin_dwall = 0_C_INT
        real(C_DOUBLE) :: rans_kpin_dwall(3, 8) = 0.0d0
        ! [rans] boostconv: steady-state residual-recombination
        ! accelerator (docs/next_session_boostconv.md); RANS-only.
        logical(C_BOOL) :: rans_boostconv = .false.
        integer(C_INT) :: rans_boostconv_interval = 25_C_INT
        integer(C_INT) :: rans_boostconv_capacity = 10_C_INT
        real(C_DOUBLE) :: rans_boostconv_tau = 1.0d-3
        real(C_DOUBLE) :: rans_boostconv_alpha = 0.02d0
        integer(C_INT) :: rans_boostconv_start = 0_C_INT
        character(len=256) :: field_prefix = ""
        integer :: field_interval = 0
        character(len=256) :: restart_file = ""
    end type dns_type

    ! Grid generation parameters and the global node lines. The staggered
    ! coordinates and finite-difference metrics live per block in
    ! block_set_type, sliced from these lines (slice_grid_direction).
    type :: grid_type
        integer(C_INT) :: distribution(1:3) = GRID_UNIFORM
        ! Build the n-point line by midpoint subdivision of the (n/2)-point
        ! line generated with the same parameters - bitwise identical to
        ! one refinement level of the coarser line (blocks.f90 does the
        ! same subdivision), so a uniformly fine reference run can share
        ! its grid exactly with a refined run's fine level.
        logical(C_BOOL) :: subdivided(1:3) = .false.
        real(C_DOUBLE) :: stretch(1:3) = 0.0d0
        real(C_DOUBLE) :: natural_dyw_plus(1:3) = 0.05d0
        logical(C_BOOL) :: natural_one_sided(1:3) = .false.
        real(C_DOUBLE), allocatable :: xNode(:), yNode(:), zNode(:)
    end type grid_type

contains

subroutine splash(has_terminal)
  logical, intent(in), optional :: has_terminal
  logical :: terminal

  terminal = .true.
  if (present(has_terminal)) terminal = has_terminal
  if (.not. terminal) return

  write(*,'(A)') "      __. - ~ ~ ~ - ."
  write(*,'(A)') "_   ,//           __  ' ,"
  write(*,'(A)') " \\  ||       __--  --    ,   ~~~~"
  write(*,'(A)') " , \\|\____---    o   \    ~~~    ~~~~"
  write(*,'(A)') ",   \ _            __/   ~~ ,  ~~~               mobyDiff"          
  write(*,'(A)') ",       \---/ / __--   ~~   ,~~                  commit: 7aa1c7b"
  write(*,'(A)') " ,          \/       ~~   ~~"
  write(*,'(A)') "  ,         ~~~ ~~~     ~~,"
  write(*,'(A)') "    ,    ~~~           , '"
  write(*,'(A)') "      ' - , _ _ _ ,  '"
  write(*,'(A)') ""
  write(*,'(A)') ""
end subroutine splash

subroutine init_grid(g, dns, periodic)
    type(grid_type), intent(inout) :: g
    type(dns_type), intent(inout)  :: dns
    logical(C_BOOL), intent(in)    :: periodic(1:3)

    call destroy_grid(g)

    allocate(g%xNode(0:int(dns%globalSize(1))))
    allocate(g%yNode(0:int(dns%globalSize(2))))
    allocate(g%zNode(0:int(dns%globalSize(3))))

    call build_node_line(g%xNode, dns%globalSize(1), dns%leng(1), &
        g%distribution(1), g%stretch(1), g%natural_one_sided(1), g%natural_dyw_plus(1), &
        g%subdivided(1))
    call build_node_line(g%yNode, dns%globalSize(2), dns%leng(2), &
        g%distribution(2), g%stretch(2), g%natural_one_sided(2), g%natural_dyw_plus(2), &
        g%subdivided(2))
    call build_node_line(g%zNode, dns%globalSize(3), dns%leng(3), &
        g%distribution(3), g%stretch(3), g%natural_one_sided(3), g%natural_dyw_plus(3), &
        g%subdivided(3))
end subroutine init_grid

! moby_prepare runs without the MPI Cartesian decomposition; give dns the
! whole-grid local size it would otherwise get from comm_init.
subroutine set_serial_local_size(dns)
    type(dns_type), intent(inout) :: dns
    integer :: dir

    do dir = 1, 3
        dns%localSize(dir,0) = 1_C_INT
        dns%localSize(dir,1) = dns%globalSize(dir)
        dns%localSize(dir,2) = dns%globalSize(dir)
    end do
end subroutine set_serial_local_size

subroutine destroy_grid(g)
    type(grid_type), intent(inout) :: g

    if (allocated(g%xNode)) deallocate(g%xNode)
    if (allocated(g%yNode)) deallocate(g%yNode)
    if (allocated(g%zNode)) deallocate(g%zNode)
end subroutine destroy_grid

recursive subroutine build_node_line(node, nGlobal, length, distribution, stretch, &
        natural_one_sided, natural_dyw_plus, subdivided)
    real(C_DOUBLE), intent(inout) :: node(0:)
    integer(C_INT), intent(in) :: nGlobal, distribution
    real(C_DOUBLE), intent(in) :: length, stretch, natural_dyw_plus
    logical(C_BOOL), intent(in) :: natural_one_sided
    logical(C_BOOL), intent(in), optional :: subdivided

    integer :: i, n
    real(C_DOUBLE) :: s
    real(C_DOUBLE), allocatable :: coarse(:)

    n = int(nGlobal)
    if (present(subdivided)) then
        if (subdivided) then
            ! Midpoint subdivision of the half-resolution line, exactly as
            ! blocks.f90 builds refinement-level lines.
            if (mod(n, 2) /= 0) error stop "subdivided grid needs an even point count"
            allocate(coarse(0:n/2))
            call build_node_line(coarse, int(n/2, C_INT), length, distribution, stretch, &
                natural_one_sided, natural_dyw_plus)
            do i = 0, n/2 - 1
                node(2*i) = coarse(i)
                node(2*i+1) = 0.5d0*(coarse(i) + coarse(i+1))
            end do
            node(n) = coarse(n/2)
            return
        end if
    end if
    do i = 0, n
        s = real(i, C_DOUBLE) / real(n, C_DOUBLE)
        node(i) = distribution_coordinate(s, length, distribution, stretch, &
            natural_one_sided, natural_dyw_plus, n)
    end do
    node(0) = 0.0d0
    node(n) = length
end subroutine build_node_line

! Sample the local, variable-staggered coordinates and the second-order
! finite-difference metrics for a window of nLocal cells starting at the
! 1-based global cell index `first` of a given global node line.
!
! This is shared by the rank-local grid setup (init_grid_direction above) and
! by the per-block metric setup in the blocks module, so the discrete
! operators are defined in exactly one place.
subroutine slice_grid_direction(node, coord, d1, lapM, lap0, lapP, nGlobal, first, nLocal, &
        length, periodic, dir)
    real(C_DOUBLE), intent(in) :: node(0:)
    real(C_DOUBLE), intent(inout) :: coord(-1:,:)
    real(C_DOUBLE), intent(inout) :: d1(0:,:)
    real(C_DOUBLE), intent(inout) :: lapM(0:,:), lap0(0:,:), lapP(0:,:)
    integer(C_INT), intent(in) :: nGlobal, first
    integer, intent(in) :: nLocal, dir
    real(C_DOUBLE), intent(in) :: length
    logical(C_BOOL), intent(in) :: periodic

    integer :: i, var, n, loCoord, hiCoord
    real(C_DOUBLE) :: hm, hp

    n = int(nGlobal)
    loCoord = lbound(coord,1)
    hiCoord = ubound(coord,1)

    do var = VAR_U, VAR_P
        ! Coordinates depend on both direction and variable because of staggering.
        do i = loCoord, hiCoord
            if (is_face_staggered(dir, var)) then
                coord(i,var) = face_at(node, n, length, int(first) + i - 2, periodic)
            else
                coord(i,var) = cell_center_at(node, n, length, int(first) + i - 1, periodic)
            end if
        end do

        d1(:,var) = 0.0d0
        lapM(:,var) = 0.0d0
        lap0(:,var) = 0.0d0
        lapP(:,var) = 0.0d0

        ! First-derivative inverse spacings connect the opposite staggering.
        do i = 0, nLocal+1
            if (is_face_staggered(dir, var)) then
                d1(i,var) = 1.0d0 / (cell_center_at(node, n, length, int(first) + i - 1, periodic) &
                                   - cell_center_at(node, n, length, int(first) + i - 2, periodic))
            else
                d1(i,var) = 1.0d0 / (face_at(node, n, length, int(first) + i - 1, periodic) &
                                   - face_at(node, n, length, int(first) + i - 2, periodic))
            end if
        end do

        ! Three-point second-derivative stencil on nonuniform spacing.
        do i = 1, nLocal
            hm = coord(i,var) - coord(i-1,var)
            hp = coord(i+1,var) - coord(i,var)
            lapM(i,var) = 2.0d0 / (hm * (hm + hp))
            lapP(i,var) = 2.0d0 / (hp * (hm + hp))
            lap0(i,var) = -(lapM(i,var) + lapP(i,var))
        end do
    end do
end subroutine slice_grid_direction

logical function is_face_staggered(dir, var)
    integer, intent(in) :: dir, var

    is_face_staggered = (dir == 1 .and. var == VAR_U) .or. &
                        (dir == 2 .and. var == VAR_V) .or. &
                        (dir == 3 .and. var == VAR_W)
end function is_face_staggered

real(C_DOUBLE) function distribution_coordinate(s, length, distribution, stretch, natural_one_sided, &
        natural_dyw_plus, n) result(x)
    real(C_DOUBLE), intent(in) :: s, length, stretch, natural_dyw_plus
    integer(C_INT), intent(in) :: distribution
    logical(C_BOOL), intent(in) :: natural_one_sided
    integer, intent(in) :: n
    real(C_DOUBLE), parameter :: pi = 3.1415926535897932384626433832795d0
    real(C_DOUBLE) :: a

    select case (distribution)
    case (GRID_COSINE)
        x = 0.5d0 * length * (1.0d0 - cos(pi*s))
    case (GRID_TANH)
        a = max(stretch, 1.0d-12)
        x = 0.5d0 * length * (1.0d0 + tanh(a*(2.0d0*s - 1.0d0))/tanh(a))
    case (GRID_NATURAL)
        x = natural_channel_coordinate(s, length, stretch, natural_one_sided, natural_dyw_plus, n)
    case default
        x = length*s
    end select
end function distribution_coordinate

real(C_DOUBLE) function natural_channel_coordinate(s, length, blend_index, one_sided, dy_wall_plus, n) result(x)
    real(C_DOUBLE), intent(in) :: s, length, blend_index, dy_wall_plus
    logical(C_BOOL), intent(in) :: one_sided
    integer, intent(in) :: n

    real(C_DOUBLE) :: half_length, j, jmax, yplus, yplus_max

    if (.not. one_sided) then
        half_length = 0.5d0*length
        j = min(s, 1.0d0 - s) * real(n, C_DOUBLE)
        jmax = 0.5d0 * real(n, C_DOUBLE)
        yplus = natural_wall_coordinate(j, blend_index, dy_wall_plus)
        yplus_max = natural_wall_coordinate(jmax, blend_index, dy_wall_plus)
        if (yplus_max <= 0.0d0) then
            x = length*s
        else if (s <= 0.5d0) then
            x = half_length * yplus/yplus_max
        else
            x = length - half_length * yplus/yplus_max
        end if
    else
        j = s * real(n, C_DOUBLE)
        jmax = real(n, C_DOUBLE)
        yplus = natural_wall_coordinate(j, blend_index, dy_wall_plus)
        yplus_max = natural_wall_coordinate(jmax, blend_index, dy_wall_plus)
        if (yplus_max <= 0.0d0) then
            x = length*s
        else
            x = length * yplus/yplus_max
        end if
    end if
end function natural_channel_coordinate

real(C_DOUBLE) function natural_wall_coordinate(j, blend_index, dy_wall_plus) result(yplus)
    real(C_DOUBLE), intent(in) :: j, blend_index, dy_wall_plus

    real(C_DOUBLE), parameter :: alpha = 1.25d0
    real(C_DOUBLE), parameter :: c_eta = 0.8d0
    real(C_DOUBLE) :: jb, blend, outer, dy_wall

    if (j <= 0.0d0) then
        yplus = 0.0d0
        return
    end if

    jb = merge(blend_index, 40.0d0, blend_index > 0.0d0)
    dy_wall = merge(dy_wall_plus, 0.05d0, dy_wall_plus > 0.0d0)
    blend = (j/jb)**2
    outer = (0.75d0*alpha*c_eta*j)**(4.0d0/3.0d0)
    yplus = (dy_wall*j + outer*blend)/(1.0d0 + blend)
end function natural_wall_coordinate

real(C_DOUBLE) function face_at(node, n, length, idx, periodic) result(x)
    real(C_DOUBLE), intent(in) :: node(0:)
    integer, intent(in) :: n, idx
    real(C_DOUBLE), intent(in) :: length
    logical(C_BOOL), intent(in) :: periodic

    ! Extend the node line into halos by wrapping or mirroring.
    if (idx < 0) then
        if (periodic) then
            x = node(idx + n) - length
        else
            x = 2.0d0*node(0) - node(-idx)
        end if
    else if (idx > n) then
        if (periodic) then
            x = node(idx - n) + length
        else
            x = 2.0d0*node(n) - node(2*n - idx)
        end if
    else
        x = node(idx)
    end if
end function face_at

real(C_DOUBLE) function cell_center_at(node, n, length, idx, periodic) result(x)
    real(C_DOUBLE), intent(in) :: node(0:)
    integer, intent(in) :: n, idx
    real(C_DOUBLE), intent(in) :: length
    logical(C_BOOL), intent(in) :: periodic

    x = 0.5d0 * (face_at(node, n, length, idx - 1, periodic) + &
                 face_at(node, n, length, idx, periodic))
end function cell_center_at

end module init
