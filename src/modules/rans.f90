!--------------------------!
!                          !
!      RANS (k-omega       !
!      SST) module         !
!                          !
!--------------------------!
!
! k-omega SST (Menter 2003 constants), resolved-wall mode (IDDES phase T2,
! docs/next_session_iddes.md), the T4 gamma-Re_thetat transition variant
! ([rans] transition = true; Langtry & Menter 2009 = OpenFOAM kOmegaSSTLM,
! resolved walls only): two extra transported scalars, the intermittency
! gamma and the transition-onset Reynolds number Re_thetat~, riding the
! same substage machinery as k/omega; the empirical correlations are
! transcribed VERBATIM from OpenFOAM kOmegaSSTLM.C (unit-tested by
! src/test_transition.f90); coupling into SST is P_k -> gamma_eff P_k and
! the k destruction scaled by min(max(gamma_eff, 0.1), 1), plus the LM F1
! modification max(F1, F3). gamma_eff = gamma — the separation-induced
! gamma_sep branch (Langtry-Menter Eq. 18) is a SEPARATE, not yet
! implemented increment. Every transition path is branch-gated on the
! flag, so transition = false stays bit-exact vs T3. Also the T3
! wall-function mode
! ([rans] wall_treatment = wall_function; Weber thesis Eqs. 4.39-4.42):
! the constrained-cell omega becomes the stepwise viscous/log blend on the
! k-based y+ (omega_wall_blend), the wall-cell nut is overwritten with the
! log-law value (nut_wall_value, copied into the domain-wall ghosts so the
! wall FACE carries it), and the log-branch wall-cell k production comes
! from the tangential velocity relative to the local wall normal
! (normalized grad dwall). Every wall-function path is branch-gated on the
! mode, so resolved mode stays bit-exact vs T2. The module owns the SST
! GEOMETRY state (phase T1: wall distance + IBM wall cells, below) and the
! TRANSPORT state (k, omega + their low-storage RK3 rhs history). The
! model is a producer of the one cell-centred turb%nut; everything
! downstream of nut (halo exchange, momentum correction, dt limit, io) is
! the untouched consumer chain.
!
! Physics note: the eddy-viscosity correction applies the DEVIATORIC
! Boussinesq stress; the -(2/3) k delta_ij part is absorbed into pressure,
! so in RANS runs the output p is a modified pressure.
!
! Discretization (per substage, rans_substage):
!   - omega is PINNED before the transport kernel reads neighbours (the
!     explicit-RK analogue of OpenFOAM's matrix constraint): IBM wall cells
!     and solid cells (Weber thesis 4.39-4.42, viscous limb — resolved
!     mode), and the first cell rows adjacent to no-slip physical domain
!     faces (Menter's omega wall condition).
!   - convection is FIRST-ORDER UPWIND: the doc's TVD van Leer limiter
!     needs a second upwind cell, which the single halo layer does not
!     carry — falling back near block edges would make results depend on
!     the block decomposition. T4 STEP-0 DECISION (2026-07-09): kept for
!     the transition scalars too. A velocity inlet CAN be composed from
!     the existing faces (Dirichlet velocity + Neumann pressure, zero-
!     gradient outlet), but no such case is validated. Since T5 STEP 0 an
!     inlet face CAN at least be classified correctly ([boundary]
!     <dir>_<side>_patch = patch overrides the tangential-Dirichlet
!     inference in domain_face_is_wall, boundary.f90), but the transported
!     scalars still have no Dirichlet-inlet ghost values (the applicator's
!     SCALAR_BC_VALUE mode is the hook). Hence the
!     T4 gates are channels, whose gamma fronts are wall-normal with ~0
!     mean cross-front velocity: the upwind numerical diffusion
!     |v| dy/2 across the front is orders of magnitude below the
!     physical diffusivity nu + nut/sigma_f (quantified on the T4 gate
!     fields by validation/rans_sst/t4_front_check.py), and the T4 gates
!     discriminate laminar vs turbulent correctly. Revisit (a second
!     scalar-only halo layer + TVD; block-edge fallbacks stay forbidden)
!     together with the flat-plate increment: inlet-vs-wall face
!     classification + scalar inlet values + outflow validation.
!   - diffusion is a face-averaged effective diffusivity (nu + sigma*nut,
!     arithmetic mean, nut lagged by one assembly) times the central
!     gradient; fluxes through solid staggered faces are masked (the
!     per-cell analogue of FACE_CLOSED), which makes k zero-gradient into
!     the IBM wall with no wall function.
!   - the destruction terms are POINT-IMPLICIT (Patankar): sinks linear in
!     the own variable are integrated by division, never subtracted
!     explicitly (near walls omega ~ 6 nu/(beta1 y^2) is huge); floors
!     k >= 0, omega >= OMEGA_MIN after every update.
!
! Phase T1 geometry state (unchanged):
!
! Wall distance (Weber thesis section 4.3, translated to the staggered
! penalization solver):
!   dwall = min(distance to the immersed surface, distance to the domain
!               wall faces), at CELL CENTRES, ghost layer included (every
!               value is evaluated pointwise from geometry, so no halo
!               exchange is ever needed).
!   yeff  = max(dwall, half the smallest local grid spacing). Without the
!           floor, cell centres grazing the staircase surface get dwall -> 0
!           and the omega wall condition ~ y^-2 explodes cell-to-cell.
! Sources for the immersed-surface part:
!   file-based IBM  the per-leaf dwall_blocks tiles written by
!                   `mobygeom.py block-table` (evaluated at each leaf's
!                   level, exactly like coef_blocks, cross-checked against
!                   the solver's leaf table at read);
!   analytic IBM    the geometry-agnostic walldist machinery (walldist.f90,
!                   IDDES phase T1b): surface point cloud from the isInBody
!                   indicator alone + kd-tree nearest + local polish to
!                   [rans] dwall_tol, queried at each block's own
!                   (level-sliced) cell centres.
! The domain-wall part comes from the global node-line ends of every
! non-periodic direction whose face is no-slip (Dirichlet on both
! tangential velocity components).
!
! IBM wall cells: fluid cells with at least one of their six staggered
! faces solid (the ibm_aware coefficient test les.f90 uses) — the
! cell-centred stand-in for a wall patch where T2 pins omega. Cells with
! all six faces solid are marked solid (benign values only).

module rans
    use, intrinsic :: iso_c_binding
    use :: init, only: dns_type, grid_type, VAR_U, VAR_V, VAR_W, VAR_P
    use :: blocks, only: block_set_type, FACE_PHYS
    use :: boundary, only: boundary_type, boundary_face_id, NFACES, &
        domain_face_is_wall, apply_scalar_bc, PATCH_UNSET, &
        SCALAR_BC_NONE, SCALAR_BC_COPY, SCALAR_BC_MIRROR
    use :: ibmm, only: ibm_type, isInBody
    use :: walldist, only: walldist_type, build_walldist, destroy_walldist, &
        walldist_distance
    use :: turbulence, only: turb_type
    use :: les_model, only: velocity_gradient_tensor
    use :: io, only: read_dwall_blocks, write_rans_geometry_file, read_scalar_field
    use :: comm, only: comm_type, comm_allreduce_sum, exchange_scalar_halos
    implicit none

    private
    public :: sst_type
    public :: init_rans_geometry, destroy_rans_geometry
    public :: enter_rans_data, exit_rans_data
    public :: write_rans_geometry
    public :: init_rans_transport, rans_prepare, rans_substage
    public :: rans_set_constrained_cells
    ! T4 transition correlations, public for the host-side unit tests
    ! (src/test_transition.f90).
    public :: lm_rethetac, lm_flength, lm_fonset, lm_fturb, lm_rethetat0, &
        lm_fthetat

    ! SST closure constants (Menter 2003). Set 1 = inner/near-wall,
    ! set 2 = outer; blended by F1.
    real(C_DOUBLE), parameter :: SST_A1 = 0.31d0
    real(C_DOUBLE), parameter :: SST_BETA_STAR = 0.09d0
    real(C_DOUBLE), parameter :: SST_SIGK1 = 0.85d0, SST_SIGK2 = 1.0d0
    real(C_DOUBLE), parameter :: SST_SIGW1 = 0.5d0, SST_SIGW2 = 0.856d0
    real(C_DOUBLE), parameter :: SST_BETA1 = 0.075d0, SST_BETA2 = 0.0828d0
    real(C_DOUBLE), parameter :: SST_ALPHA1 = 5.0d0/9.0d0, SST_ALPHA2 = 0.44d0
    ! Floors: k may touch zero; omega must stay positive (it divides).
    real(C_DOUBLE), parameter :: OMEGA_MIN = 1.0d-8

    ! Wall-function constants (T3, Weber thesis Eqs. 4.39-4.42 = the
    ! OpenFOAM omega/nutk wall functions): kappa, E, C_mu^(1/4), and the
    ! viscous/log switch y+_lam = the fixed point of y+ = ln(E y+)/kappa
    ! (OpenFOAM iterates the same fixed point; kappa 0.41, E 9.8 -> 11.53).
    ! y+ is built from the wall-cell k: y+ = C_mu^(1/4) sqrt(k) y_eff / nu.
    real(C_DOUBLE), parameter :: WF_KAPPA = 0.41d0
    real(C_DOUBLE), parameter :: WF_E = 9.8d0
    real(C_DOUBLE), parameter :: WF_CMU25 = 0.5477225575051661d0   ! 0.09**0.25
    real(C_DOUBLE), parameter :: WF_YPLUS_LAM = 11.530107402304532d0

    ! gamma-Re_thetat transition constants (T4; Langtry & Menter 2009, the
    ! OpenFOAM kOmegaSSTLM defaults).
    real(C_DOUBLE), parameter :: LM_CA1 = 2.0d0, LM_CA2 = 0.06d0
    real(C_DOUBLE), parameter :: LM_CE1 = 1.0d0, LM_CE2 = 50.0d0
    real(C_DOUBLE), parameter :: LM_CTHETAT = 0.03d0
    ! Diffusivity structure differs per equation (kOmegaSSTLM):
    ! gamma uses nu + nut/sigma_f; Re_thetat~ uses sigma_thetat*(nu + nut).
    real(C_DOUBLE), parameter :: LM_SIGMAF = 1.0d0
    real(C_DOUBLE), parameter :: LM_SIGMATHETAT = 2.0d0
    ! The Re_thetat,eq fixed point in the pressure-gradient parameter
    ! lambda: OpenFOAM iterates to lambdaErr with no hard cap (it only
    ! warns past maxLambdaIter = 10); a device kernel needs a bound, and
    ! the map contracts in a handful of iterations.
    real(C_DOUBLE), parameter :: LM_LAMBDA_ERR = 1.0d-6
    integer, parameter :: LM_MAX_LAMBDA_ITER = 32
    ! Us floor (OpenFOAM deltaU = small): keeps Tu, dU/ds and the
    ! relaxation time finite at stagnant cells; every source vanishes in
    ! that limit.
    real(C_DOUBLE), parameter :: LM_DELTAU = 1.0d-15

    ! Wall-cell marker values (per-cell byte).
    integer(C_SIGNED_CHAR), parameter, public :: WALL_CELL_FLUID = 0_C_SIGNED_CHAR
    integer(C_SIGNED_CHAR), parameter, public :: WALL_CELL_WALL  = 1_C_SIGNED_CHAR
    integer(C_SIGNED_CHAR), parameter, public :: WALL_CELL_SOLID = 2_C_SIGNED_CHAR

    ! "No wall anywhere" distance (fully periodic, no IBM): large but finite,
    ! which is also the correct SST limit (the wall-blend functions go to
    ! their far-field branch).
    real(C_DOUBLE), parameter :: NO_WALL_DISTANCE = 1.0d30
    ! Same staggered-coefficient threshold the LES ibm_aware test uses:
    ! solid faces carry coef = 1e30/re.
    real(C_DOUBLE), parameter :: SOLID_FACE_THRESHOLD = 1.0d20

    type :: sst_type
        logical(C_BOOL) :: geometry_built = .false.
        ! Cell-centred wall distance, ghost layer included (0:nb+1).
        real(C_DOUBLE), allocatable :: dwall(:,:,:,:)   ! raw distance
        real(C_DOUBLE), allocatable :: yeff(:,:,:,:)    ! floored (use this in the model)
        ! Interior-only per-cell marker (WALL_CELL_*).
        integer(C_SIGNED_CHAR), allocatable :: wallcell(:,:,:,:)  ! (1:nb,...,nBlocks)

        ! T2 transport state (allocated only for [turbulence] model = rans).
        logical(C_BOOL) :: transport_built = .false.
        real(C_DOUBLE), allocatable :: k(:,:,:,:), omg(:,:,:,:)      ! (0:nb+1,...)
        real(C_DOUBLE), allocatable :: ks(:,:,:,:), omgs(:,:,:,:)    ! substage scratch
        real(C_DOUBLE), allocatable :: koldrhs(:,:,:,:), omgoldrhs(:,:,:,:)
        ! Interior byte: 1 = first cell adjacent to a no-slip physical DOMAIN
        ! face (omega pinned there like Menter's wall condition; distinct
        ! from the IBM wallcell marker, which additionally zeroes P_k/nut).
        integer(C_SIGNED_CHAR), allocatable :: domwall(:,:,:,:)
        ! Unit wall-normal (normalized grad dwall) at the constrained
        ! (IBM-wall + domain-wall) cells, zero elsewhere; the T3
        ! wall-function log-branch production projects it out of the
        ! cell-centre velocity. Filled only for wall_treatment =
        ! wall_function (allocated always so the device maps are uniform).
        real(C_DOUBLE), allocatable :: wnorm(:,:,:,:,:)   ! (3,1:nb,...,nBlocks)
        ! Per-face no-slip flags for the cell-centred scalar BCs (index =
        ! boundary_face_id).
        integer(C_INT) :: facewall(NFACES) = 0_C_INT

        ! T4 transition state (gamma, Re_thetat~): full arrays only when
        ! [rans] transition = true, 1-cell dummies otherwise so the device
        ! maps and kernel map clauses stay uniform (the wnorm idiom) —
        ! every access is guarded by the transition flag.
        real(C_DOUBLE), allocatable :: gam(:,:,:,:), ret(:,:,:,:)       ! (0:nb+1,...)
        real(C_DOUBLE), allocatable :: gams(:,:,:,:), rets(:,:,:,:)     ! substage scratch
        real(C_DOUBLE), allocatable :: gamoldrhs(:,:,:,:), retoldrhs(:,:,:,:)
    end type sst_type

contains

    subroutine init_rans_geometry(sst, dns, g, blk, bc, ibm, c)
        type(sst_type), intent(inout) :: sst
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        type(block_set_type), intent(in) :: blk
        type(boundary_type), intent(in) :: bc
        type(ibm_type), intent(in) :: ibm
        type(comm_type), intent(in) :: c

        integer :: nx, ny, nz

        call destroy_rans_geometry(sst)

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))

        allocate(sst%dwall(0:nx+1,0:ny+1,0:nz+1,blk%nBlocks))
        allocate(sst%yeff(0:nx+1,0:ny+1,0:nz+1,blk%nBlocks))
        allocate(sst%wallcell(nx,ny,nz,blk%nBlocks))

        sst%dwall = NO_WALL_DISTANCE
        if (dns%ibm_enabled) then
            if (len_trim(dns%ibm_coeff_file) > 0) then
                call read_body_distance_file(sst, dns, blk, c%has_terminal)
            else
                call fill_body_distance_analytic(sst, dns, blk, bc, ibm, c%has_terminal)
            end if
        end if
        call min_in_domain_wall_distance(sst, dns, g, blk, bc)
        call apply_distance_floor(sst, blk)
        call classify_wall_cells(sst, dns, blk, ibm)

        sst%geometry_built = .true.
        call report_geometry(sst, blk, c)
        call report_patch_types(bc, c)
    end subroutine init_rans_geometry

    ! Print the resolved patch type of every non-periodic domain face, and
    ! whether it was declared ([boundary] <dir>_<side>_patch) or inferred
    ! from the tangential velocity BCs (T5 STEP 0 visibility).
    subroutine report_patch_types(bc, c)
        type(boundary_type), intent(in) :: bc
        type(comm_type), intent(in) :: c

        character(len=5), parameter :: face_names(NFACES) = &
            [character(len=5) :: "x_min", "x_max", "y_min", "y_max", "z_min", "z_max"]
        integer :: dir, side, face_id

        if (.not. c%has_terminal) return
        do dir = 1, 3
            if (bc%isPeriodic(dir)) cycle
            do side = 0, 1
                face_id = boundary_face_id(dir, side)
                print '(5A)', " RANS domain face ", trim(face_names(face_id)), ": ", &
                    trim(merge("wall ", "patch", domain_face_is_wall(bc, dir, side))), &
                    trim(merge(" (declared)", " (inferred)", &
                        bc%facePatchType(face_id) /= PATCH_UNSET))
            end do
        end do
    end subroutine report_patch_types

    subroutine destroy_rans_geometry(sst)
        type(sst_type), intent(inout) :: sst

        if (allocated(sst%dwall)) deallocate(sst%dwall)
        if (allocated(sst%yeff)) deallocate(sst%yeff)
        if (allocated(sst%wallcell)) deallocate(sst%wallcell)
        if (allocated(sst%k)) deallocate(sst%k, sst%omg, sst%ks, sst%omgs, &
            sst%koldrhs, sst%omgoldrhs, sst%domwall, sst%wnorm, &
            sst%gam, sst%ret, sst%gams, sst%rets, sst%gamoldrhs, sst%retoldrhs)
        sst%geometry_built = .false.
        sst%transport_built = .false.
    end subroutine destroy_rans_geometry

    subroutine enter_rans_data(sst)
        type(sst_type), intent(inout) :: sst

        if (.not. allocated(sst%dwall)) return

        !$omp target enter data map(to: sst)
        !$omp target enter data map(to: sst%dwall, sst%yeff, sst%wallcell)
        if (allocated(sst%k)) then
            !$omp target enter data map(to: sst%k, sst%omg, sst%ks, sst%omgs)
            !$omp target enter data map(to: sst%koldrhs, sst%omgoldrhs, sst%domwall, sst%wnorm)
            !$omp target enter data map(to: sst%gam, sst%ret, sst%gams, sst%rets)
            !$omp target enter data map(to: sst%gamoldrhs, sst%retoldrhs)
        end if
    end subroutine enter_rans_data

    subroutine exit_rans_data(sst)
        type(sst_type), intent(inout) :: sst

        if (.not. allocated(sst%dwall)) return

        if (allocated(sst%k)) then
            !$omp target exit data map(delete: sst%gamoldrhs, sst%retoldrhs)
            !$omp target exit data map(delete: sst%gam, sst%ret, sst%gams, sst%rets)
            !$omp target exit data map(delete: sst%koldrhs, sst%omgoldrhs, sst%domwall, sst%wnorm)
            !$omp target exit data map(delete: sst%k, sst%omg, sst%ks, sst%omgs)
        end if
        !$omp target exit data map(delete: sst%dwall, sst%yeff, sst%wallcell)
        !$omp target exit data map(delete: sst)
    end subroutine exit_rans_data

    ! Immersed-surface distance from the coefficient file: the per-leaf
    ! dwall_blocks tiles (mobygeom block-table). A coefficient file without
    ! them (legacy stl-ibm-coeff layout) cannot serve a RANS run.
    subroutine read_body_distance_file(sst, dns, blk, has_terminal)
        type(sst_type), intent(inout) :: sst
        type(dns_type), intent(in) :: dns
        type(block_set_type), intent(in) :: blk
        logical, intent(in) :: has_terminal

        logical :: found

        call read_dwall_blocks(sst%dwall, found, dns, blk, has_terminal)
        if (.not. found) then
            if (has_terminal) print *, "error: [rans] needs the wall distance, but the", &
                " coefficient file has no dwall_blocks; regenerate it with mobygeom block-table: ", &
                trim(dns%ibm_coeff_file)
            error stop
        end if
    end subroutine read_body_distance_file

    ! Immersed-surface distance for the analytic IBM: geometry-agnostic,
    ! computed from the isInBody indicator alone (walldist.f90, phase T1b) —
    ! isInBody is the one user-editable analytic-geometry hook, so dwall
    ! stays correct for ANY body defined there. The surface cloud samples
    ! the finest refinement level; queries run at each block's own
    ! (level-sliced) cell centres, so refined leaves are exact at their
    ! level for free.
    subroutine fill_body_distance_analytic(sst, dns, blk, bc, ibm, has_terminal)
        type(sst_type), intent(inout) :: sst
        type(dns_type), intent(in) :: dns
        type(block_set_type), intent(in) :: blk
        type(boundary_type), intent(in) :: bc
        type(ibm_type), intent(in) :: ibm
        logical, intent(in) :: has_terminal

        type(walldist_type) :: w
        integer :: i, j, k, b, nx, ny, nz, nf(3)
        real(C_DOUBLE) :: xA(1:3)

        nf = int(dns%globalSize)*2**(int(blk%nLevels) - 1)
        call build_walldist(w, isInBody, ibm, dns, &
            blk%lineX(0:nf(1), int(blk%nLevels)), &
            blk%lineY(0:nf(2), int(blk%nLevels)), &
            blk%lineZ(0:nf(3), int(blk%nLevels)), &
            nf, dns%leng, logical(bc%isPeriodic), has_terminal)
        ! No surface crossing (body outside the domain / below grid
        ! resolution): keep the no-wall distance, build_walldist warned.
        if (w%nCloud == 0) return

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))

        !$omp parallel do collapse(2) private(i,j,k,b,xA)
        do b = 1, int(blk%nBlocks)
        do k = 0, nz+1
            do j = 0, ny+1
                do i = 0, nx+1
                    xA(1) = blk%x(i, VAR_P, b)
                    xA(2) = blk%y(j, VAR_P, b)
                    xA(3) = blk%z(k, VAR_P, b)
                    sst%dwall(i,j,k,b) = walldist_distance(w, isInBody, ibm, dns, &
                        xA, dns%rans_dwall_tol)
                end do
            end do
        end do
        end do
        !$omp end parallel do

        call destroy_walldist(w)
    end subroutine fill_body_distance_analytic

    ! min in the distance to every no-slip domain face: non-periodic
    ! direction, Dirichlet BC on both tangential velocity components. The
    ! wall planes are the global node-line ends; ghost-cell centres outside
    ! the domain get the mirrored (positive) distance.
    subroutine min_in_domain_wall_distance(sst, dns, g, blk, bc)
        type(sst_type), intent(inout) :: sst
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        type(block_set_type), intent(in) :: blk
        type(boundary_type), intent(in) :: bc

        integer :: i, j, k, b, dir, side, nx, ny, nz
        real(C_DOUBLE) :: wall_coord, d
        logical :: is_wall(1:3, 0:1)

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))

        do dir = 1, 3
            do side = 0, 1
                is_wall(dir, side) = domain_face_is_wall(bc, dir, side)
            end do
        end do
        if (.not. any(is_wall)) return

        do dir = 1, 3
            do side = 0, 1
                if (.not. is_wall(dir, side)) cycle
                wall_coord = domain_wall_coordinate(g, dns, dir, side)
                do b = 1, int(blk%nBlocks)
                do k = 0, nz+1
                    do j = 0, ny+1
                        do i = 0, nx+1
                            select case (dir)
                            case (1)
                                d = abs(blk%x(i, VAR_P, b) - wall_coord)
                            case (2)
                                d = abs(blk%y(j, VAR_P, b) - wall_coord)
                            case default
                                d = abs(blk%z(k, VAR_P, b) - wall_coord)
                            end select
                            sst%dwall(i,j,k,b) = min(sst%dwall(i,j,k,b), d)
                        end do
                    end do
                end do
                end do
            end do
        end do
    end subroutine min_in_domain_wall_distance

    real(C_DOUBLE) function domain_wall_coordinate(g, dns, dir, side) result(x)
        type(grid_type), intent(in) :: g
        type(dns_type), intent(in) :: dns
        integer, intent(in) :: dir, side

        select case (dir)
        case (1)
            x = g%xNode(merge(int(dns%globalSize(1)), 0, side == 1))
        case (2)
            x = g%yNode(merge(int(dns%globalSize(2)), 0, side == 1))
        case default
            x = g%zNode(merge(int(dns%globalSize(3)), 0, side == 1))
        end select
    end function domain_wall_coordinate

    ! yeff = max(dwall, half the smallest local spacing), per cell from the
    ! block metric tables (blk%d1? hold the inverse cell widths).
    subroutine apply_distance_floor(sst, blk)
        type(sst_type), intent(inout) :: sst
        type(block_set_type), intent(in) :: blk

        integer :: i, j, k, b, nx, ny, nz
        real(C_DOUBLE) :: floor_dist

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))

        do b = 1, int(blk%nBlocks)
        do k = 0, nz+1
            do j = 0, ny+1
                do i = 0, nx+1
                    floor_dist = 0.5d0*min(1.0d0/blk%d1x(i, VAR_P, b), &
                                           1.0d0/blk%d1y(j, VAR_P, b), &
                                           1.0d0/blk%d1z(k, VAR_P, b))
                    sst%yeff(i,j,k,b) = max(sst%dwall(i,j,k,b), floor_dist)
                end do
            end do
        end do
        end do
    end subroutine apply_distance_floor

    ! Classify interior cells from the staggered IBM coefficients: solid
    ! faces carry coef = 1e30/re (the ibm_aware test). 0 of 6 solid = fluid,
    ! all 6 = solid, anything between = IBM wall cell. Reads the HOST copy
    ! of ibm%coef — the caller refreshes it from the device first when the
    ! analytic coefficients were computed there.
    subroutine classify_wall_cells(sst, dns, blk, ibm)
        type(sst_type), intent(inout) :: sst
        type(dns_type), intent(in) :: dns
        type(block_set_type), intent(in) :: blk
        type(ibm_type), intent(in) :: ibm

        integer :: i, j, k, b, nx, ny, nz, nSolid

        sst%wallcell = WALL_CELL_FLUID
        if (.not. dns%ibm_enabled) return

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))

        do b = 1, int(blk%nBlocks)
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    nSolid = count([abs(ibm%coef(i,  j,  k,  VAR_U, b)) > SOLID_FACE_THRESHOLD, &
                                    abs(ibm%coef(i+1,j,  k,  VAR_U, b)) > SOLID_FACE_THRESHOLD, &
                                    abs(ibm%coef(i,  j,  k,  VAR_V, b)) > SOLID_FACE_THRESHOLD, &
                                    abs(ibm%coef(i,  j+1,k,  VAR_V, b)) > SOLID_FACE_THRESHOLD, &
                                    abs(ibm%coef(i,  j,  k,  VAR_W, b)) > SOLID_FACE_THRESHOLD, &
                                    abs(ibm%coef(i,  j,  k+1,VAR_W, b)) > SOLID_FACE_THRESHOLD])
                    if (nSolid == 0) then
                        sst%wallcell(i,j,k,b) = WALL_CELL_FLUID
                    else if (nSolid == 6) then
                        sst%wallcell(i,j,k,b) = WALL_CELL_SOLID
                    else
                        sst%wallcell(i,j,k,b) = WALL_CELL_WALL
                    end if
                end do
            end do
        end do
        end do
    end subroutine classify_wall_cells

    subroutine report_geometry(sst, blk, c)
        type(sst_type), intent(in) :: sst
        type(block_set_type), intent(in) :: blk
        type(comm_type), intent(in) :: c

        real(C_DOUBLE) :: counts(2)
        integer :: nx, ny, nz

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))

        counts(1) = real(count(sst%wallcell == WALL_CELL_WALL), C_DOUBLE)
        counts(2) = real(count(sst%wallcell == WALL_CELL_SOLID), C_DOUBLE)
        call comm_allreduce_sum(c, counts)
        if (c%has_terminal) then
            print '(A,I0,A,I0,A)', " RANS geometry: ", nint(counts(1)), &
                " IBM wall cells, ", nint(counts(2)), " solid cells"
        end if
    end subroutine report_geometry

    ! Diagnostic dump ([rans] dump_geometry): one self-contained HDF5 file
    ! with the leaf table, interior dwall/yeff/wallcell and the per-block
    ! cell-centre coordinates (so gate scripts need no grid reconstruction).
    subroutine write_rans_geometry(sst, blk, dns, c)
        type(sst_type), intent(in) :: sst
        type(block_set_type), intent(in) :: blk
        type(dns_type), intent(in) :: dns
        type(comm_type), intent(in) :: c

        integer(C_INT), allocatable :: wallcell_int(:,:,:,:)
        real(C_DOUBLE), allocatable :: xc(:,:), yc(:,:), zc(:,:)
        character(len=280) :: file_name
        integer :: b, nx, ny, nz

        if (.not. sst%geometry_built) return

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))

        allocate(wallcell_int(nx,ny,nz,blk%nBlocks))
        wallcell_int = int(sst%wallcell, C_INT)
        allocate(xc(nx,blk%nBlocks), yc(ny,blk%nBlocks), zc(nz,blk%nBlocks))
        do b = 1, int(blk%nBlocks)
            xc(:,b) = blk%x(1:nx, VAR_P, b)
            yc(:,b) = blk%y(1:ny, VAR_P, b)
            zc(:,b) = blk%z(1:nz, VAR_P, b)
        end do

        if (len_trim(dns%field_prefix) > 0) then
            file_name = trim(dns%field_prefix)//"_ransgeom.h5"
        else
            file_name = "ransgeom.h5"
        end if
        if (c%has_terminal) print *, "writing RANS geometry: ", trim(file_name)
        call write_rans_geometry_file(file_name, blk, sst%dwall, sst%yeff, &
            wallcell_int, xc, yc, zc, c%has_terminal)

        deallocate(wallcell_int, xc, yc, zc)
    end subroutine write_rans_geometry

    !========================
    ! T2: k-omega SST transport (resolved wall mode)
    !========================

    ! Allocate the transport state and set the initial condition: per cell,
    ! k = 1.5 (tu/100 |u|)^2 and omega = k/(nut_ratio nu) (floored), then the
    ! constrained cells (solid k = 0, pinned omega). Requires the geometry
    ! state; runs on the HOST before enter_rans_data maps everything.
    subroutine init_rans_transport(sst, dns, blk, bc, ibm, has_terminal)
        type(sst_type), intent(inout) :: sst
        type(dns_type), intent(in) :: dns
        type(block_set_type), intent(in) :: blk
        type(boundary_type), intent(in) :: bc
        type(ibm_type), intent(in) :: ibm
        logical, intent(in) :: has_terminal

        integer :: i, j, k, b, nx, ny, nz, dir, side
        real(C_DOUBLE) :: umag2, uc, vc, wc, tu_frac, nu, kin
        logical :: found_k, found_omg, found_gam, found_ret

        if (.not. sst%geometry_built) error stop "RANS transport needs the geometry state"

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))

        allocate(sst%k(0:nx+1,0:ny+1,0:nz+1,blk%nBlocks))
        allocate(sst%omg(0:nx+1,0:ny+1,0:nz+1,blk%nBlocks))
        allocate(sst%ks(0:nx+1,0:ny+1,0:nz+1,blk%nBlocks))
        allocate(sst%omgs(0:nx+1,0:ny+1,0:nz+1,blk%nBlocks))
        allocate(sst%koldrhs(0:nx+1,0:ny+1,0:nz+1,blk%nBlocks))
        allocate(sst%omgoldrhs(0:nx+1,0:ny+1,0:nz+1,blk%nBlocks))
        allocate(sst%domwall(nx,ny,nz,blk%nBlocks))
        allocate(sst%wnorm(3,nx,ny,nz,blk%nBlocks))
        sst%wnorm = 0.0d0

        ! T4 transition scalars: full arrays only when the model is on,
        ! 1-cell dummies otherwise (uniform device maps; all accesses are
        ! flag-guarded).
        if (dns%rans_transition) then
            allocate(sst%gam(0:nx+1,0:ny+1,0:nz+1,blk%nBlocks))
            allocate(sst%ret(0:nx+1,0:ny+1,0:nz+1,blk%nBlocks))
            allocate(sst%gams(0:nx+1,0:ny+1,0:nz+1,blk%nBlocks))
            allocate(sst%rets(0:nx+1,0:ny+1,0:nz+1,blk%nBlocks))
            allocate(sst%gamoldrhs(0:nx+1,0:ny+1,0:nz+1,blk%nBlocks))
            allocate(sst%retoldrhs(0:nx+1,0:ny+1,0:nz+1,blk%nBlocks))
        else
            allocate(sst%gam(0:0,0:0,0:0,1), sst%ret(0:0,0:0,0:0,1), &
                sst%gams(0:0,0:0,0:0,1), sst%rets(0:0,0:0,0:0,1), &
                sst%gamoldrhs(0:0,0:0,0:0,1), sst%retoldrhs(0:0,0:0,0:0,1))
        end if
        ! IC: gamma = 1 (fully intermittent freestream, per Langtry-Menter);
        ! Re_thetat~ = the equilibrium correlation at the configured
        ! freestream Tu with zero pressure gradient (lambda = 0 makes it
        ! closed-form, independent of nu/U).
        sst%gam = 1.0d0
        sst%ret = lm_rethetat0(dns%rans_tu, 0.0d0, 1.0d0, 1.0d0)
        sst%gamoldrhs = 0.0d0
        sst%retoldrhs = 0.0d0

        ! First interior cell rows adjacent to a no-slip physical domain
        ! face: omega is pinned there (Menter wall condition).
        sst%domwall = 0_C_SIGNED_CHAR
        do dir = 1, 3
            do side = 0, 1
                sst%facewall(boundary_face_id(dir, side)) = &
                    merge(1_C_INT, 0_C_INT, domain_face_is_wall(bc, dir, side))
            end do
        end do
        do b = 1, int(blk%nBlocks)
            if (blk%physLow(1,b) == FACE_PHYS .and. sst%facewall(boundary_face_id(1,0)) == 1) &
                sst%domwall(1,:,:,b) = 1_C_SIGNED_CHAR
            if (blk%physHigh(1,b) == FACE_PHYS .and. sst%facewall(boundary_face_id(1,1)) == 1) &
                sst%domwall(nx,:,:,b) = 1_C_SIGNED_CHAR
            if (blk%physLow(2,b) == FACE_PHYS .and. sst%facewall(boundary_face_id(2,0)) == 1) &
                sst%domwall(:,1,:,b) = 1_C_SIGNED_CHAR
            if (blk%physHigh(2,b) == FACE_PHYS .and. sst%facewall(boundary_face_id(2,1)) == 1) &
                sst%domwall(:,ny,:,b) = 1_C_SIGNED_CHAR
            if (blk%physLow(3,b) == FACE_PHYS .and. sst%facewall(boundary_face_id(3,0)) == 1) &
                sst%domwall(:,:,1,b) = 1_C_SIGNED_CHAR
            if (blk%physHigh(3,b) == FACE_PHYS .and. sst%facewall(boundary_face_id(3,1)) == 1) &
                sst%domwall(:,:,nz,b) = 1_C_SIGNED_CHAR
        end do

        ! T3 wall functions need the wall-normal direction in every
        ! constrained cell (log-branch production projects it out of the
        ! velocity). Resolved mode never reads wnorm; keep it zero there.
        if (dns%rans_wall_treatment == 1_C_INT) &
            call compute_wall_normals(sst, dns, blk, ibm)

        tu_frac = dns%rans_tu/100.0d0
        nu = 1.0d0/dns%re
        sst%k = 0.0d0
        sst%omg = OMEGA_MIN
        do b = 1, int(blk%nBlocks)
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    uc = 0.5d0*(blk%q(i,j,k,VAR_U,b) + blk%q(i+1,j,k,VAR_U,b))
                    vc = 0.5d0*(blk%q(i,j,k,VAR_V,b) + blk%q(i,j+1,k,VAR_V,b))
                    wc = 0.5d0*(blk%q(i,j,k,VAR_W,b) + blk%q(i,j,k+1,VAR_W,b))
                    umag2 = uc*uc + vc*vc + wc*wc
                    kin = max(1.5d0*tu_frac*tu_frac*umag2, 1.0d-10)
                    sst%k(i,j,k,b) = kin
                    ! Wall-consistent omega: the freestream value blended
                    ! with the viscous wall profile 6 nu/(beta1 y^2) at every
                    ! cell's own wall distance. A flat freestream omega
                    ! against the pinned wall rows is a shock that the
                    ! explicit transport cannot survive.
                    sst%omg(i,j,k,b) = max(kin/(dns%rans_nut_ratio*nu), &
                        omega_wall_value(nu, sst%yeff(i,j,k,b)), OMEGA_MIN)
                end do
            end do
        end do
        end do
        sst%koldrhs = 0.0d0
        sst%omgoldrhs = 0.0d0

        ! Restarted runs carry k/omega in the restart file (T2 snapshots
        ! append them); absent datasets (older restarts) keep the tu/
        ! nut_ratio initialization above, with a warning.
        if (len_trim(dns%restart_file) > 0) then
            call read_scalar_field(blk, "k", dns%restart_file, sst%k, found_k, has_terminal)
            call read_scalar_field(blk, "omega", dns%restart_file, sst%omg, found_omg, has_terminal)
            if (.not. (found_k .and. found_omg)) then
                if (has_terminal) print *, "warning: restart file has no k/omega datasets;", &
                    " initializing from [rans] tu / nut_ratio"
            end if
            if (dns%rans_transition) then
                call read_scalar_field(blk, "gamma", dns%restart_file, sst%gam, &
                    found_gam, has_terminal)
                call read_scalar_field(blk, "rethetat", dns%restart_file, sst%ret, &
                    found_ret, has_terminal)
                if (.not. (found_gam .and. found_ret)) then
                    ! Partial presence reinitializes BOTH: the pair is one
                    ! coupled state.
                    if (has_terminal) print *, "warning: restart file has no gamma/rethetat", &
                        " datasets; initializing transition scalars from [rans] tu"
                    sst%gam = 1.0d0
                    sst%ret = lm_rethetat0(dns%rans_tu, 0.0d0, 1.0d0, 1.0d0)
                end if
            end if
        end if

        call set_constrained_cells_host(sst, dns, blk)
        sst%ks = sst%k
        sst%omgs = sst%omg
        sst%gams = sst%gam
        sst%rets = sst%ret
        sst%transport_built = .true.
    end subroutine init_rans_transport

    ! Menter viscous-limb wall omega at the effective distance.
    pure real(C_DOUBLE) function omega_wall_value(nu, y) result(w)
!$omp declare target
        real(C_DOUBLE), intent(in) :: nu, y

        w = 6.0d0*nu/(SST_BETA1*y*y)
    end function omega_wall_value

    ! T3 wall-function wall omega: stepwise viscous/log switch on the
    ! k-based y+ (Weber Eqs. 4.39-4.42). At k -> 0 (initial condition,
    ! laminarizing flow) the switch always takes the viscous limb, so the
    ! blend degrades continuously to the resolved-mode value.
    pure real(C_DOUBLE) function omega_wall_blend(nu, y, kv) result(w)
!$omp declare target
        real(C_DOUBLE), intent(in) :: nu, y, kv

        real(C_DOUBLE) :: sqrtk

        sqrtk = sqrt(max(kv, 0.0d0))
        if (WF_CMU25*sqrtk*y/nu < WF_YPLUS_LAM) then
            w = omega_wall_value(nu, y)
        else
            w = sqrtk/(WF_CMU25*WF_KAPPA*y)
        end if
    end function omega_wall_blend

    ! T3 wall-function wall-cell eddy viscosity (Weber Eq. 4.41 / OpenFOAM
    ! nutkWallFunction): zero on the viscous branch, nu (y+ kappa /
    ! ln(E y+) - 1) on the log branch (continuous at y+_lam by the fixed
    ! point). This is what carries the log-law wall shear into the
    ! momentum equation on a coarse near-wall grid.
    pure real(C_DOUBLE) function nut_wall_value(nu, yplus) result(nutw)
!$omp declare target
        real(C_DOUBLE), intent(in) :: nu, yplus

        if (yplus < WF_YPLUS_LAM) then
            nutw = 0.0d0
        else
            nutw = nu*(yplus*WF_KAPPA/log(WF_E*yplus) - 1.0d0)
        end if
    end function nut_wall_value

    !========================
    ! T4: gamma-Re_thetat transition correlations (Langtry & Menter 2009)
    !========================
    ! Transcribed VERBATIM from OpenFOAM kOmegaSSTLM.C — the piecewise fits
    ! are empirical, do NOT re-derive or "simplify" the coefficients. Pure
    ! device functions, unit-tested host-side by src/test_transition.f90
    ! against an independent transcription BEFORE they run in a kernel.

    ! Critical (onset) Reynolds number Re_thetac(Re_thetat~).
    pure real(C_DOUBLE) function lm_rethetac(ret) result(rc)
!$omp declare target
        real(C_DOUBLE), intent(in) :: ret

        if (ret <= 1870.0d0) then
            rc = ret - 396.035d-2 + 120.656d-4*ret - 868.230d-6*ret**2 &
               + 696.506d-9*ret**3 - 174.105d-12*ret**4
        else
            rc = ret - 593.11d0 - 0.482d0*(ret - 1870.0d0)
        end if
    end function lm_rethetac

    ! Transition length function F_length(Re_thetat~), blended to 40 inside
    ! the viscous sublayer on R_omega = y^2 omega / (500 nu) (the exponent
    ! groups it as y^2 omega/(200 nu), OpenFOAM's form).
    pure real(C_DOUBLE) function lm_flength(ret, y, wv, nu) result(fl)
!$omp declare target
        real(C_DOUBLE), intent(in) :: ret, y, wv, nu

        real(C_DOUBLE) :: fsub

        if (ret < 400.0d0) then
            fl = 398.189d-1 - 119.270d-4*ret - 132.567d-6*ret**2
        else if (ret < 596.0d0) then
            fl = 263.404d0 - 123.939d-2*ret + 194.548d-5*ret**2 &
               - 101.695d-8*ret**3
        else if (ret < 1200.0d0) then
            fl = 0.5d0 - 3.0d-4*(ret - 596.0d0)
        else
            fl = 0.3188d0
        end if
        fsub = exp(-(y*y*wv/(200.0d0*nu))**2)
        fl = fl*(1.0d0 - fsub) + 40.0d0*fsub
    end function lm_flength

    ! Onset trigger F_onset(Re_v, Re_thetac, R_T): the strain-rate Reynolds
    ! number against the critical value, shut off at high turbulent
    ! Reynolds number R_T = k/(nu omega).
    pure real(C_DOUBLE) function lm_fonset(rev, rethetac, rt) result(fon)
!$omp declare target
        real(C_DOUBLE), intent(in) :: rev, rethetac, rt

        real(C_DOUBLE) :: fon1, fon2, fon3

        fon1 = rev/(2.193d0*rethetac)
        fon2 = min(max(fon1, fon1**4), 2.0d0)
        fon3 = max(1.0d0 - (rt/2.5d0)**3, 0.0d0)
        fon = max(fon2 - fon3, 0.0d0)
    end function lm_fonset

    ! F_turb: disables the gamma destruction in fully turbulent regions.
    pure real(C_DOUBLE) function lm_fturb(rt) result(ft)
!$omp declare target
        real(C_DOUBLE), intent(in) :: rt

        ft = exp(-(0.25d0*rt)**4)
    end function lm_fturb

    ! Equilibrium transition-onset Reynolds number Re_thetat,eq(Tu, lambda)
    ! — the freestream correlation with its capped fixed-point iteration in
    ! the pressure-gradient parameter lambda = theta^2/nu dU/ds. Tu is in
    ! PERCENT and floored at 0.027 (OpenFOAM). Both Tu branches use the
    ! same F_lambda for dU/ds > 0 (exp(-Tu/0.5) = exp(-2 Tu) — OpenFOAM
    ! writes it both ways; 0.5 and 2 are exact binary scalings, so the
    ! merged form is bitwise identical). At lambda = 0 (dU/ds = 0) the
    ! result is closed-form and independent of nu and Us.
    pure real(C_DOUBLE) function lm_rethetat0(tu_in, dusds, nu, us) result(ret0)
!$omp declare target
        real(C_DOUBLE), intent(in) :: tu_in, dusds, nu, us

        real(C_DOUBLE) :: tu, lam, lam0, flam, thetat
        integer :: iter

        tu = max(tu_in, 0.027d0)
        lam = 0.0d0
        thetat = 0.0d0
        do iter = 1, LM_MAX_LAMBDA_ITER
            lam0 = lam
            if (dusds <= 0.0d0) then
                flam = 1.0d0 - (-12.986d0*lam - 123.66d0*lam**2 &
                     - 405.689d0*lam**3)*exp(-(tu/1.5d0)**1.5d0)
            else
                flam = 1.0d0 + 0.275d0*(1.0d0 - exp(-35.0d0*lam))*exp(-2.0d0*tu)
            end if
            if (tu <= 1.3d0) then
                thetat = (1173.51d0 - 589.428d0*tu + 0.2196d0/(tu*tu))*flam*nu/us
            else
                thetat = 331.50d0*(tu - 0.5658d0)**(-0.671d0)*flam*nu/us
            end if
            lam = thetat*thetat/nu*dusds
            lam = max(min(lam, 0.1d0), -0.1d0)
            if (abs(lam - lam0) <= LM_LAMBDA_ERR) exit
        end do
        ret0 = max(thetat*us/nu, 20.0d0)
    end function lm_rethetat0

    ! F_thetat: ~1 inside a boundary layer (freezes the transported
    ! Re_thetat~), -> 0 in the freestream (relaxes it to the local
    ! equilibrium correlation). y/delta = Us^2/(375 Omega nu Re_thetat~)
    ! — the y in delta cancels (OpenFOAM's yBydelta).
    pure real(C_DOUBLE) function lm_fthetat(gam, ret, us, omgmag, wv, y, nu) result(fth)
!$omp declare target
        real(C_DOUBLE), intent(in) :: gam, ret, us, omgmag, wv, y, nu

        real(C_DOUBLE) :: ybydelta, reomega, fwake

        ybydelta = us*us/max(375.0d0*omgmag*nu*ret, LM_DELTAU**2)
        reomega = y*y*wv/nu
        fwake = exp(-(reomega/1.0d5)**2)
        fth = min(max(fwake*exp(-ybydelta**4), &
            1.0d0 - ((gam - 1.0d0/LM_CE2)/(1.0d0 - 1.0d0/LM_CE2))**2), 1.0d0)
    end function lm_fthetat

    ! Unit wall-normal at the constrained (IBM-wall + domain-wall) cells as
    ! the normalized gradient of the RAW dwall field (yeff's floor flattens
    ! exactly the near-wall slope the normal lives on). Per direction the
    ! difference is one-sided AWAY from a solid staggered face and away
    ! from a no-slip physical face (dwall is V-shaped across the immersed
    ! surface and mirrors across a domain wall, so a centred difference
    ! there sees a spurious near-zero slope); centred otherwise. Host-only
    ! init code; halo dwall values are pointwise geometry, hence valid.
    subroutine compute_wall_normals(sst, dns, blk, ibm)
        type(sst_type), intent(inout) :: sst
        type(dns_type), intent(in) :: dns
        type(block_set_type), intent(in) :: blk
        type(ibm_type), intent(in) :: ibm

        integer :: i, j, k, b, nx, ny, nz, d
        real(C_DOUBLE) :: grad(3), gnorm
        logical :: lowbad(3), highbad(3), ibm_on

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))
        ibm_on = dns%ibm_enabled

        !$omp parallel do collapse(2) private(i,j,k,b,d,grad,gnorm,lowbad,highbad)
        do b = 1, int(blk%nBlocks)
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    if (sst%wallcell(i,j,k,b) /= WALL_CELL_WALL .and. &
                        sst%domwall(i,j,k,b) == 0_C_SIGNED_CHAR) cycle

                    lowbad = .false.
                    highbad = .false.
                    if (ibm_on) then
                        lowbad(1) = abs(ibm%coef(i,  j,  k,  VAR_U,b)) > SOLID_FACE_THRESHOLD
                        highbad(1) = abs(ibm%coef(i+1,j,  k,  VAR_U,b)) > SOLID_FACE_THRESHOLD
                        lowbad(2) = abs(ibm%coef(i,  j,  k,  VAR_V,b)) > SOLID_FACE_THRESHOLD
                        highbad(2) = abs(ibm%coef(i,  j+1,k,  VAR_V,b)) > SOLID_FACE_THRESHOLD
                        lowbad(3) = abs(ibm%coef(i,  j,  k,  VAR_W,b)) > SOLID_FACE_THRESHOLD
                        highbad(3) = abs(ibm%coef(i,  j,  k+1,VAR_W,b)) > SOLID_FACE_THRESHOLD
                    end if
                    if (i == 1 .and. blk%physLow(1,b) == FACE_PHYS .and. &
                        sst%facewall(boundary_face_id(1,0)) == 1) lowbad(1) = .true.
                    if (i == nx .and. blk%physHigh(1,b) == FACE_PHYS .and. &
                        sst%facewall(boundary_face_id(1,1)) == 1) highbad(1) = .true.
                    if (j == 1 .and. blk%physLow(2,b) == FACE_PHYS .and. &
                        sst%facewall(boundary_face_id(2,0)) == 1) lowbad(2) = .true.
                    if (j == ny .and. blk%physHigh(2,b) == FACE_PHYS .and. &
                        sst%facewall(boundary_face_id(2,1)) == 1) highbad(2) = .true.
                    if (k == 1 .and. blk%physLow(3,b) == FACE_PHYS .and. &
                        sst%facewall(boundary_face_id(3,0)) == 1) lowbad(3) = .true.
                    if (k == nz .and. blk%physHigh(3,b) == FACE_PHYS .and. &
                        sst%facewall(boundary_face_id(3,1)) == 1) highbad(3) = .true.

                    grad(1) = dwall_slope(sst%dwall(i-1,j,k,b), sst%dwall(i,j,k,b), &
                        sst%dwall(i+1,j,k,b), blk%x(i-1,VAR_P,b), blk%x(i,VAR_P,b), &
                        blk%x(i+1,VAR_P,b), lowbad(1), highbad(1))
                    grad(2) = dwall_slope(sst%dwall(i,j-1,k,b), sst%dwall(i,j,k,b), &
                        sst%dwall(i,j+1,k,b), blk%y(j-1,VAR_P,b), blk%y(j,VAR_P,b), &
                        blk%y(j+1,VAR_P,b), lowbad(2), highbad(2))
                    grad(3) = dwall_slope(sst%dwall(i,j,k-1,b), sst%dwall(i,j,k,b), &
                        sst%dwall(i,j,k+1,b), blk%z(k-1,VAR_P,b), blk%z(k,VAR_P,b), &
                        blk%z(k+1,VAR_P,b), lowbad(3), highbad(3))

                    gnorm = sqrt(grad(1)**2 + grad(2)**2 + grad(3)**2)
                    ! Degenerate gradient (equidistant ridge): leave the
                    ! normal zero -- the production then uses the full
                    ! velocity magnitude, the conservative choice.
                    if (gnorm > 1.0d-12) then
                        do d = 1, 3
                            sst%wnorm(d,i,j,k,b) = grad(d)/gnorm
                        end do
                    end if
                end do
            end do
        end do
        end do
        !$omp end parallel do
    end subroutine compute_wall_normals

    pure real(C_DOUBLE) function dwall_slope(dm, d0, dp, xm, x0, xp, lowbad, highbad) &
            result(slope)
        real(C_DOUBLE), intent(in) :: dm, d0, dp, xm, x0, xp
        logical, intent(in) :: lowbad, highbad

        if (lowbad .and. highbad) then
            slope = 0.0d0
        else if (lowbad) then
            slope = (dp - d0)/(xp - x0)
        else if (highbad) then
            slope = (d0 - dm)/(x0 - xm)
        else
            slope = (dp - dm)/(xp - xm)
        end if
    end function dwall_slope

    ! Host-side constrained-cell values (init/restart time, before the
    ! device maps exist): solid cells carry k = 0 and a benign pinned
    ! omega; IBM wall cells and domain-wall rows carry the pinned omega.
    subroutine set_constrained_cells_host(sst, dns, blk)
        type(sst_type), intent(inout) :: sst
        type(dns_type), intent(in) :: dns
        type(block_set_type), intent(in) :: blk

        integer :: i, j, k, b, nx, ny, nz
        real(C_DOUBLE) :: nu
        logical :: pinned, wallfn

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))
        nu = 1.0d0/dns%re
        wallfn = dns%rans_wall_treatment == 1_C_INT

        do b = 1, int(blk%nBlocks)
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    pinned = sst%wallcell(i,j,k,b) /= WALL_CELL_FLUID .or. &
                             sst%domwall(i,j,k,b) /= 0_C_SIGNED_CHAR
                    if (sst%wallcell(i,j,k,b) == WALL_CELL_SOLID) sst%k(i,j,k,b) = 0.0d0
                    if (pinned) then
                        if (wallfn) then
                            sst%omg(i,j,k,b) = omega_wall_blend(nu, sst%yeff(i,j,k,b), &
                                sst%k(i,j,k,b))
                        else
                            sst%omg(i,j,k,b) = omega_wall_value(nu, sst%yeff(i,j,k,b))
                        end if
                    end if
                end do
            end do
        end do
        end do
    end subroutine set_constrained_cells_host

    ! Device twin of set_constrained_cells_host, run at the top of every
    ! substage BEFORE the transport kernel reads neighbours, so adjacent
    ! cells' diffusion sees the constrained values.
    subroutine rans_set_constrained_cells(sst, dns, blk)
        type(sst_type), intent(inout) :: sst
        type(dns_type), intent(in) :: dns
        type(block_set_type), intent(in) :: blk

        integer :: i, j, k, b, nx, ny, nz, nBlocks
        real(C_DOUBLE) :: nu
        logical :: wallfn

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))
        nBlocks = int(blk%nBlocks)
        nu = 1.0d0/dns%re
        wallfn = dns%rans_wall_treatment == 1_C_INT

        !$omp target teams distribute parallel do collapse(4) &
        !$omp& map(to: nu, wallfn, sst%wallcell, sst%domwall, sst%yeff) &
        !$omp& map(tofrom: sst%k, sst%omg) &
        !$omp& private(i,j,k,b)
        do b = 1, nBlocks
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    if (sst%wallcell(i,j,k,b) == WALL_CELL_SOLID) sst%k(i,j,k,b) = 0.0d0
                    if (sst%wallcell(i,j,k,b) /= WALL_CELL_FLUID .or. &
                        sst%domwall(i,j,k,b) /= 0_C_SIGNED_CHAR) then
                        if (wallfn) then
                            sst%omg(i,j,k,b) = omega_wall_blend(nu, sst%yeff(i,j,k,b), &
                                sst%k(i,j,k,b))
                        else
                            sst%omg(i,j,k,b) = omega_wall_value(nu, sst%yeff(i,j,k,b))
                        end if
                    end if
                end do
            end do
        end do
        end do
        !$omp end target teams distribute parallel do
    end subroutine rans_set_constrained_cells

    ! Cell-centred scalar ghosts at physical domain faces, via the generic
    ! boundary.f90 applicator (this module only supplies the per-scalar mode
    ! tables): k is Dirichlet 0 at no-slip walls (mirror) and zero-gradient
    ! elsewhere; omega is always zero-gradient (its wall value is the pinned
    ! first cell, not a ghost); gamma and Re_thetat~ are zero-gradient
    ! everywhere (kOmegaSSTLM wall condition).
    subroutine rans_apply_scalar_bcs(sst, dns, blk, bc)
        type(sst_type), intent(inout) :: sst
        type(dns_type), intent(in) :: dns
        type(block_set_type), intent(inout) :: blk
        type(boundary_type), intent(in) :: bc

        integer(C_INT) :: kmode(NFACES), copymode(NFACES)
        integer :: f

        do f = 1, NFACES
            kmode(f) = merge(SCALAR_BC_MIRROR, SCALAR_BC_COPY, sst%facewall(f) == 1_C_INT)
        end do
        copymode = SCALAR_BC_COPY

        call apply_scalar_bc(blk, bc, sst%k, kmode)
        call apply_scalar_bc(blk, bc, sst%omg, copymode)
        if (dns%rans_transition) then
            call apply_scalar_bc(blk, bc, sst%gam, copymode)
            call apply_scalar_bc(blk, bc, sst%ret, copymode)
        end if
    end subroutine rans_apply_scalar_bcs

    ! Assemble the SST eddy viscosity into the caller-supplied nut target:
    ! nut = a1 k / max(a1 omega, S F2); zero in IBM wall cells (viscous
    ! limb) and solid cells.
    subroutine rans_assemble_nut(sst, turb, blk, dns, nut)
        type(sst_type), intent(inout) :: sst
        type(turb_type), intent(in) :: turb
        type(block_set_type), intent(inout) :: blk
        type(dns_type), intent(in) :: dns
        real(C_DOUBLE), intent(inout) :: nut(0:,0:,0:,1:)

        integer :: i, j, k, b, nx, ny, nz, nBlocks
        real(C_DOUBLE) :: nu, y, kv, wv, s2, smag, arg2, f2, yplus
        real(C_DOUBLE) :: g11, g12, g13, g21, g22, g23, g31, g32, g33
        logical :: wallfn

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))
        nBlocks = int(blk%nBlocks)
        nu = 1.0d0/dns%re
        wallfn = dns%rans_wall_treatment == 1_C_INT

        !$omp target teams distribute parallel do collapse(4) &
        !$omp& map(to: nu, wallfn, sst%k, sst%omg, sst%yeff, sst%wallcell, sst%domwall, &
        !$omp& blk%q, blk%d1x, blk%d1y, blk%d1z, &
        !$omp& turb%d1xm, turb%d1x0, turb%d1xp, turb%d1ym, turb%d1y0, turb%d1yp, &
        !$omp& turb%d1zm, turb%d1z0, turb%d1zp, &
        !$omp& turb%p_from_u_x, turb%p_from_v_y, turb%p_from_w_z) &
        !$omp& map(tofrom: nut) &
        !$omp& private(i,j,k,b,y,kv,wv,s2,smag,arg2,f2,yplus, &
        !$omp& g11,g12,g13,g21,g22,g23,g31,g32,g33)
        do b = 1, nBlocks
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    ! Constrained cells: solid = 0; IBM wall cells = 0 in
                    ! resolved mode (the Weber viscous limb). T3 wall
                    ! functions OVERWRITE IBM-wall and domain-wall cells
                    ! with the viscous/log blend (nut_wall_value).
                    if (sst%wallcell(i,j,k,b) /= WALL_CELL_FLUID .or. &
                        (wallfn .and. sst%domwall(i,j,k,b) /= 0_C_SIGNED_CHAR)) then
                        nut(i,j,k,b) = 0.0d0
                        if (wallfn .and. sst%wallcell(i,j,k,b) /= WALL_CELL_SOLID) then
                            yplus = WF_CMU25*sqrt(max(sst%k(i,j,k,b), 0.0d0)) &
                                  *sst%yeff(i,j,k,b)/nu
                            nut(i,j,k,b) = nut_wall_value(nu, yplus)
                        end if
                        cycle
                    end if
                    kv = sst%k(i,j,k,b)
                    wv = sst%omg(i,j,k,b)
                    y = sst%yeff(i,j,k,b)

                    call velocity_gradient_tensor(blk, turb, i, j, k, b, &
                        g11, g12, g13, g21, g22, g23, g31, g32, g33)
                    s2 = 2.0d0*(g11*g11 + g22*g22 + g33*g33) &
                       + (g12 + g21)*(g12 + g21) &
                       + (g13 + g31)*(g13 + g31) &
                       + (g23 + g32)*(g23 + g32)
                    smag = sqrt(max(s2, 0.0d0))

                    arg2 = max(2.0d0*sqrt(max(kv, 0.0d0))/(SST_BETA_STAR*wv*y), &
                               500.0d0*nu/(y*y*wv))
                    f2 = tanh(arg2*arg2)
                    nut(i,j,k,b) = SST_A1*max(kv, 0.0d0)/max(SST_A1*wv, smag*f2)
                end do
            end do
        end do
        end do
        !$omp end target teams distribute parallel do
    end subroutine rans_assemble_nut

    ! T3 wall functions: copy the wall-cell nut into the ghost across every
    ! no-slip physical face, so the face-interpolated eddy viscosity the
    ! momentum correction sees at the wall face IS the wall value and the
    ! delivered wall shear is (nu + nut_w) U_1/y_1 -- the whole point of
    ! the log-branch nut. Resolved mode leaves the ghosts alone (nut -> 0
    ! at resolved walls, ghost 0 is consistent).
    subroutine rans_apply_nut_wall_ghosts(sst, blk, bc, nut)
        type(sst_type), intent(in) :: sst
        type(block_set_type), intent(inout) :: blk
        type(boundary_type), intent(in) :: bc
        real(C_DOUBLE), intent(inout) :: nut(0:,0:,0:,1:)

        integer(C_INT) :: mode(NFACES)

        mode = merge(SCALAR_BC_COPY, SCALAR_BC_NONE, sst%facewall == 1_C_INT)
        call apply_scalar_bc(blk, bc, nut, mode)
    end subroutine rans_apply_nut_wall_ghosts

    ! Pre-loop preparation: constrained cells, scalar ghosts/halos, and the
    ! first nut assembly so the momentum predictor starts consistent.
    subroutine rans_prepare(sst, turb, blk, dns, bc, c)
        type(sst_type), intent(inout) :: sst
        type(turb_type), intent(inout) :: turb
        type(block_set_type), intent(inout) :: blk
        type(dns_type), intent(in) :: dns
        type(boundary_type), intent(in) :: bc
        type(comm_type), intent(inout) :: c

        call rans_set_constrained_cells(sst, dns, blk)
        call rans_apply_scalar_bcs(sst, dns, blk, bc)
        call exchange_scalar_halos(c, sst%k, blk)
        call exchange_scalar_halos(c, sst%omg, blk)
        if (dns%rans_transition) then
            call exchange_scalar_halos(c, sst%gam, blk)
            call exchange_scalar_halos(c, sst%ret, blk)
        end if
        call rans_assemble_nut(sst, turb, blk, dns, turb%nut)
        if (dns%rans_wall_treatment == 1_C_INT) &
            call rans_apply_nut_wall_ghosts(sst, blk, bc, turb%nut)
    end subroutine rans_prepare

    ! One RK3 substage of the SST transport: pin, ghosts, halos, the fused
    ! k/omega update into the scratch arrays, copy-back, nut assembly. The
    ! caller exchanges turb%nut afterwards (the consumer chain).
    subroutine rans_substage(sst, turb, blk, dns, ibm, bc, c, dt_alpha, dt_beta)
        type(sst_type), intent(inout) :: sst
        type(turb_type), intent(inout) :: turb
        type(block_set_type), intent(inout) :: blk
        type(dns_type), intent(in) :: dns
        type(ibm_type), intent(in) :: ibm
        type(boundary_type), intent(in) :: bc
        type(comm_type), intent(inout) :: c
        real(C_DOUBLE), intent(in) :: dt_alpha, dt_beta

        call rans_set_constrained_cells(sst, dns, blk)
        call rans_apply_scalar_bcs(sst, dns, blk, bc)
        call exchange_scalar_halos(c, sst%k, blk)
        call exchange_scalar_halos(c, sst%omg, blk)
        if (dns%rans_transition) then
            call exchange_scalar_halos(c, sst%gam, blk)
            call exchange_scalar_halos(c, sst%ret, blk)
        end if
        call rans_transport_kernel(sst, turb, blk, dns, ibm, dt_alpha, dt_beta)
        call rans_copyback(sst, dns, blk)
        call rans_assemble_nut(sst, turb, blk, dns, turb%nut)
        if (dns%rans_wall_treatment == 1_C_INT) &
            call rans_apply_nut_wall_ghosts(sst, blk, bc, turb%nut)
    end subroutine rans_substage

    subroutine rans_copyback(sst, dns, blk)
        type(sst_type), intent(inout) :: sst
        type(dns_type), intent(in) :: dns
        type(block_set_type), intent(in) :: blk

        integer :: i, j, k, b, nx, ny, nz, nBlocks
        logical :: transition

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))
        nBlocks = int(blk%nBlocks)
        transition = dns%rans_transition

        !$omp target teams distribute parallel do collapse(4) &
        !$omp& map(to: transition, sst%ks, sst%omgs, sst%gams, sst%rets) &
        !$omp& map(tofrom: sst%k, sst%omg, sst%gam, sst%ret) &
        !$omp& private(i,j,k,b)
        do b = 1, nBlocks
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    sst%k(i,j,k,b) = sst%ks(i,j,k,b)
                    sst%omg(i,j,k,b) = sst%omgs(i,j,k,b)
                    if (transition) then
                        sst%gam(i,j,k,b) = sst%gams(i,j,k,b)
                        sst%ret(i,j,k,b) = sst%rets(i,j,k,b)
                    end if
                end do
            end do
        end do
        end do
        !$omp end target teams distribute parallel do
    end subroutine rans_copyback

    ! The fused per-substage transport kernel: per cell compute the velocity
    ! gradients, S^2, F1 and the k/omega gradients once, then both scalar
    ! updates (low-storage RK3 like the momentum predictor, point-implicit
    ! sinks, floors). Writes the scratch arrays only.
    subroutine rans_transport_kernel(sst, turb, blk, dns, ibm, dt_alpha, dt_beta)
        type(sst_type), intent(inout) :: sst
        type(turb_type), intent(in) :: turb
        type(block_set_type), intent(inout) :: blk
        type(dns_type), intent(in) :: dns
        type(ibm_type), intent(in) :: ibm
        real(C_DOUBLE), intent(in) :: dt_alpha, dt_beta

        integer :: i, j, k, b, nx, ny, nz, nBlocks
        integer(C_SIGNED_CHAR) :: marker
        logical :: pinned_omega, ibm_wall, wallfn, transition
        real(C_DOUBLE) :: nu, dtsub, y, kv, wv, nutc
        real(C_DOUBLE) :: g11, g12, g13, g21, g22, g23, g31, g32, g33
        real(C_DOUBLE) :: s2, smag, gkx, gky, gkz, gwx, gwy, gwz, gkgw
        real(C_DOUBLE) :: cdkw, arg1, f1, arg2, f2, nutloc, cross
        real(C_DOUBLE) :: sigk, sigw, beta_b, alpha_b, pk
        real(C_DOUBLE) :: uw, ue, vs, vn, wb, wt
        real(C_DOUBLE) :: conv_k, conv_w, diff_k, diff_w
        real(C_DOUBLE) :: fw, fe, dcoef
        real(C_DOUBLE) :: rhsk, rhsw, knew, wnew
        real(C_DOUBLE) :: solid_threshold
        real(C_DOUBLE) :: yplus, nutw, uc, vc, wc, un, ut1, ut2, ut3, magut
        real(C_DOUBLE) :: gv, rv, omgmag, usmag, dusds, rt, rev, tuloc
        real(C_DOUBLE) :: ret0, rethetac, flength, fonset, fturbv, fthetat
        real(C_DOUBLE) :: tcoef, pgam, egam, conv_g, conv_r, diff_g, diff_r
        real(C_DOUBLE) :: diag_r, rhsg, rhsr, gnew, rnew, gameff, dkfac
        logical :: solw, sole, sols, soln, solb, solt

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))
        nBlocks = int(blk%nBlocks)
        nu = 1.0d0/dns%re
        dtsub = dt_alpha + dt_beta
        solid_threshold = SOLID_FACE_THRESHOLD
        wallfn = dns%rans_wall_treatment == 1_C_INT
        transition = dns%rans_transition

        !$omp target teams distribute parallel do collapse(4) &
        !$omp& map(to: nu, dtsub, dt_alpha, dt_beta, solid_threshold, wallfn, transition, &
        !$omp& sst%k, sst%omg, sst%gam, sst%ret, sst%yeff, sst%wallcell, sst%domwall, sst%wnorm, &
        !$omp& blk%q, blk%d1x, blk%d1y, blk%d1z, ibm%coef, &
        !$omp& turb%nut, turb%inv_dx, turb%inv_dy, turb%inv_dz, &
        !$omp& turb%d1xm, turb%d1x0, turb%d1xp, turb%d1ym, turb%d1y0, turb%d1yp, &
        !$omp& turb%d1zm, turb%d1z0, turb%d1zp, &
        !$omp& turb%p_from_u_x, turb%p_from_v_y, turb%p_from_w_z) &
        !$omp& map(tofrom: sst%ks, sst%omgs, sst%koldrhs, sst%omgoldrhs, &
        !$omp& sst%gams, sst%rets, sst%gamoldrhs, sst%retoldrhs) &
        !$omp& private(i,j,k,b,marker,pinned_omega,ibm_wall,y,kv,wv,nutc, &
        !$omp& g11,g12,g13,g21,g22,g23,g31,g32,g33,s2,smag, &
        !$omp& gkx,gky,gkz,gwx,gwy,gwz,gkgw,cdkw,arg1,f1,arg2,f2,nutloc,cross, &
        !$omp& sigk,sigw,beta_b,alpha_b,pk,uw,ue,vs,vn,wb,wt, &
        !$omp& conv_k,conv_w,diff_k,diff_w,fw,fe,dcoef,rhsk,rhsw,knew,wnew, &
        !$omp& yplus,nutw,uc,vc,wc,un,ut1,ut2,ut3,magut, &
        !$omp& gv,rv,omgmag,usmag,dusds,rt,rev,tuloc, &
        !$omp& ret0,rethetac,flength,fonset,fturbv,fthetat, &
        !$omp& tcoef,pgam,egam,conv_g,conv_r,diff_g,diff_r,diag_r, &
        !$omp& rhsg,rhsr,gnew,rnew,gameff,dkfac, &
        !$omp& solw,sole,sols,soln,solb,solt)
        do b = 1, nBlocks
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    marker = sst%wallcell(i,j,k,b)
                    if (marker == WALL_CELL_SOLID) then
                        ! Solid cells keep their benign constrained values.
                        sst%ks(i,j,k,b) = sst%k(i,j,k,b)
                        sst%omgs(i,j,k,b) = sst%omg(i,j,k,b)
                        sst%koldrhs(i,j,k,b) = 0.0d0
                        sst%omgoldrhs(i,j,k,b) = 0.0d0
                        if (transition) then
                            sst%gams(i,j,k,b) = sst%gam(i,j,k,b)
                            sst%rets(i,j,k,b) = sst%ret(i,j,k,b)
                            sst%gamoldrhs(i,j,k,b) = 0.0d0
                            sst%retoldrhs(i,j,k,b) = 0.0d0
                        end if
                        cycle
                    end if
                    ibm_wall = marker == WALL_CELL_WALL
                    pinned_omega = ibm_wall .or. sst%domwall(i,j,k,b) /= 0_C_SIGNED_CHAR

                    kv = sst%k(i,j,k,b)
                    wv = sst%omg(i,j,k,b)
                    y = sst%yeff(i,j,k,b)

                    call velocity_gradient_tensor(blk, turb, i, j, k, b, &
                        g11, g12, g13, g21, g22, g23, g31, g32, g33)
                    s2 = 2.0d0*(g11*g11 + g22*g22 + g33*g33) &
                       + (g12 + g21)*(g12 + g21) &
                       + (g13 + g31)*(g13 + g31) &
                       + (g23 + g32)*(g23 + g32)
                    s2 = max(s2, 0.0d0)
                    smag = sqrt(s2)

                    ! k and omega cell-centred gradients (blend stencils).
                    gkx = turb%d1xm(i,VAR_P,b)*sst%k(i-1,j,k,b) &
                        + turb%d1x0(i,VAR_P,b)*sst%k(i,j,k,b) &
                        + turb%d1xp(i,VAR_P,b)*sst%k(i+1,j,k,b)
                    gky = turb%d1ym(j,VAR_P,b)*sst%k(i,j-1,k,b) &
                        + turb%d1y0(j,VAR_P,b)*sst%k(i,j,k,b) &
                        + turb%d1yp(j,VAR_P,b)*sst%k(i,j+1,k,b)
                    gkz = turb%d1zm(k,VAR_P,b)*sst%k(i,j,k-1,b) &
                        + turb%d1z0(k,VAR_P,b)*sst%k(i,j,k,b) &
                        + turb%d1zp(k,VAR_P,b)*sst%k(i,j,k+1,b)
                    gwx = turb%d1xm(i,VAR_P,b)*sst%omg(i-1,j,k,b) &
                        + turb%d1x0(i,VAR_P,b)*sst%omg(i,j,k,b) &
                        + turb%d1xp(i,VAR_P,b)*sst%omg(i+1,j,k,b)
                    gwy = turb%d1ym(j,VAR_P,b)*sst%omg(i,j-1,k,b) &
                        + turb%d1y0(j,VAR_P,b)*sst%omg(i,j,k,b) &
                        + turb%d1yp(j,VAR_P,b)*sst%omg(i,j+1,k,b)
                    gwz = turb%d1zm(k,VAR_P,b)*sst%omg(i,j,k-1,b) &
                        + turb%d1z0(k,VAR_P,b)*sst%omg(i,j,k,b) &
                        + turb%d1zp(k,VAR_P,b)*sst%omg(i,j,k+1,b)
                    gkgw = gkx*gwx + gky*gwy + gkz*gwz

                    ! Menter blending functions on y_eff.
                    cdkw = max(2.0d0*SST_SIGW2/wv*gkgw, 1.0d-10)
                    arg1 = min(max(sqrt(max(kv, 0.0d0))/(SST_BETA_STAR*wv*y), &
                                   500.0d0*nu/(y*y*wv)), &
                               4.0d0*SST_SIGW2*max(kv, 0.0d0)/(cdkw*y*y))
                    f1 = tanh(arg1**4)
                    ! LM F1 modification (transition only): F1 = max(F1, F3)
                    ! with F3 = exp(-(Ry/120)^8), Ry = y sqrt(k)/nu — keeps
                    ! the inner-set blending active in low-Ry laminar layers
                    ! where the transported k would otherwise flip F1 -> 0.
                    if (transition) &
                        f1 = max(f1, exp(-(y*sqrt(max(kv, 0.0d0))/nu/120.0d0)**8))
                    arg2 = max(2.0d0*sqrt(max(kv, 0.0d0))/(SST_BETA_STAR*wv*y), &
                               500.0d0*nu/(y*y*wv))
                    f2 = tanh(arg2*arg2)
                    nutloc = SST_A1*max(kv, 0.0d0)/max(SST_A1*wv, smag*f2)

                    sigk = f1*SST_SIGK1 + (1.0d0 - f1)*SST_SIGK2
                    sigw = f1*SST_SIGW1 + (1.0d0 - f1)*SST_SIGW2
                    beta_b = f1*SST_BETA1 + (1.0d0 - f1)*SST_BETA2
                    alpha_b = f1*SST_ALPHA1 + (1.0d0 - f1)*SST_ALPHA2

                    ! Production (limited); zero in IBM wall cells (Weber
                    ! viscous limb).
                    pk = min(nutloc*s2, 10.0d0*SST_BETA_STAR*max(kv, 0.0d0)*wv)
                    if (ibm_wall) pk = 0.0d0
                    ! T3 wall functions, log branch: the constrained-cell k
                    ! production comes from the tangential velocity relative
                    ! to the local wall normal at y_eff (Weber Eq. 4.40 /
                    ! OpenFOAM omegaWallFunction G): G = (nu + nut_w)
                    ! * |U_t|/y * C_mu^(1/4) sqrt(k)/(kappa y). The viscous
                    ! branch keeps the resolved-mode rule (IBM wall cells
                    ! produce nothing, domain-wall rows produce normally).
                    if (wallfn .and. pinned_omega) then
                        yplus = WF_CMU25*sqrt(max(kv, 0.0d0))*y/nu
                        if (yplus >= WF_YPLUS_LAM) then
                            nutw = nut_wall_value(nu, yplus)
                            uc = 0.5d0*(blk%q(i,j,k,VAR_U,b) + blk%q(i+1,j,k,VAR_U,b))
                            vc = 0.5d0*(blk%q(i,j,k,VAR_V,b) + blk%q(i,j+1,k,VAR_V,b))
                            wc = 0.5d0*(blk%q(i,j,k,VAR_W,b) + blk%q(i,j,k+1,VAR_W,b))
                            un = uc*sst%wnorm(1,i,j,k,b) + vc*sst%wnorm(2,i,j,k,b) &
                               + wc*sst%wnorm(3,i,j,k,b)
                            ut1 = uc - un*sst%wnorm(1,i,j,k,b)
                            ut2 = vc - un*sst%wnorm(2,i,j,k,b)
                            ut3 = wc - un*sst%wnorm(3,i,j,k,b)
                            magut = sqrt(ut1*ut1 + ut2*ut2 + ut3*ut3)
                            pk = (nu + nutw)*(magut/y) &
                                *WF_CMU25*sqrt(max(kv, 0.0d0))/(WF_KAPPA*y)
                            pk = min(pk, 10.0d0*SST_BETA_STAR*max(kv, 0.0d0)*wv)
                        end if
                    end if

                    ! Solid staggered faces: diffusive fluxes masked.
                    solw = abs(ibm%coef(i,  j,  k,  VAR_U,b)) > solid_threshold
                    sole = abs(ibm%coef(i+1,j,  k,  VAR_U,b)) > solid_threshold
                    sols = abs(ibm%coef(i,  j,  k,  VAR_V,b)) > solid_threshold
                    soln = abs(ibm%coef(i,  j+1,k,  VAR_V,b)) > solid_threshold
                    solb = abs(ibm%coef(i,  j,  k,  VAR_W,b)) > solid_threshold
                    solt = abs(ibm%coef(i,  j,  k+1,VAR_W,b)) > solid_threshold

                    ! First-order upwind convection (see module header).
                    uw = blk%q(i,  j,k,VAR_U,b)
                    ue = blk%q(i+1,j,k,VAR_U,b)
                    vs = blk%q(i,j,  k,VAR_V,b)
                    vn = blk%q(i,j+1,k,VAR_V,b)
                    wb = blk%q(i,j,k,  VAR_W,b)
                    wt = blk%q(i,j,k+1,VAR_W,b)
                    conv_k = (ue*merge(sst%k(i,j,k,b), sst%k(i+1,j,k,b), ue > 0.0d0) &
                            - uw*merge(sst%k(i-1,j,k,b), sst%k(i,j,k,b), uw > 0.0d0)) &
                              *blk%d1x(i,VAR_P,b) &
                           + (vn*merge(sst%k(i,j,k,b), sst%k(i,j+1,k,b), vn > 0.0d0) &
                            - vs*merge(sst%k(i,j-1,k,b), sst%k(i,j,k,b), vs > 0.0d0)) &
                              *blk%d1y(j,VAR_P,b) &
                           + (wt*merge(sst%k(i,j,k,b), sst%k(i,j,k+1,b), wt > 0.0d0) &
                            - wb*merge(sst%k(i,j,k-1,b), sst%k(i,j,k,b), wb > 0.0d0)) &
                              *blk%d1z(k,VAR_P,b)
                    conv_w = (ue*merge(sst%omg(i,j,k,b), sst%omg(i+1,j,k,b), ue > 0.0d0) &
                            - uw*merge(sst%omg(i-1,j,k,b), sst%omg(i,j,k,b), uw > 0.0d0)) &
                              *blk%d1x(i,VAR_P,b) &
                           + (vn*merge(sst%omg(i,j,k,b), sst%omg(i,j+1,k,b), vn > 0.0d0) &
                            - vs*merge(sst%omg(i,j-1,k,b), sst%omg(i,j,k,b), vs > 0.0d0)) &
                              *blk%d1y(j,VAR_P,b) &
                           + (wt*merge(sst%omg(i,j,k,b), sst%omg(i,j,k+1,b), wt > 0.0d0) &
                            - wb*merge(sst%omg(i,j,k-1,b), sst%omg(i,j,k,b), wb > 0.0d0)) &
                              *blk%d1z(k,VAR_P,b)

                    ! Diffusion: face-mean effective diffusivity x central
                    ! gradient, solid faces masked. k first.
                    diff_k = 0.0d0
                    diff_w = 0.0d0

                    dcoef = nu + sigk*0.5d0*(turb%nut(i-1,j,k,b) + turb%nut(i,j,k,b))
                    fw = merge(0.0d0, dcoef*(sst%k(i,j,k,b) - sst%k(i-1,j,k,b)) &
                        *turb%inv_dx(i,VAR_P,b), solw)
                    dcoef = nu + sigk*0.5d0*(turb%nut(i,j,k,b) + turb%nut(i+1,j,k,b))
                    fe = merge(0.0d0, dcoef*(sst%k(i+1,j,k,b) - sst%k(i,j,k,b)) &
                        *turb%inv_dx(i+1,VAR_P,b), sole)
                    diff_k = diff_k + (fe - fw)*blk%d1x(i,VAR_P,b)
                    dcoef = nu + sigk*0.5d0*(turb%nut(i,j-1,k,b) + turb%nut(i,j,k,b))
                    fw = merge(0.0d0, dcoef*(sst%k(i,j,k,b) - sst%k(i,j-1,k,b)) &
                        *turb%inv_dy(j,VAR_P,b), sols)
                    dcoef = nu + sigk*0.5d0*(turb%nut(i,j,k,b) + turb%nut(i,j+1,k,b))
                    fe = merge(0.0d0, dcoef*(sst%k(i,j+1,k,b) - sst%k(i,j,k,b)) &
                        *turb%inv_dy(j+1,VAR_P,b), soln)
                    diff_k = diff_k + (fe - fw)*blk%d1y(j,VAR_P,b)
                    dcoef = nu + sigk*0.5d0*(turb%nut(i,j,k-1,b) + turb%nut(i,j,k,b))
                    fw = merge(0.0d0, dcoef*(sst%k(i,j,k,b) - sst%k(i,j,k-1,b)) &
                        *turb%inv_dz(k,VAR_P,b), solb)
                    dcoef = nu + sigk*0.5d0*(turb%nut(i,j,k,b) + turb%nut(i,j,k+1,b))
                    fe = merge(0.0d0, dcoef*(sst%k(i,j,k+1,b) - sst%k(i,j,k,b)) &
                        *turb%inv_dz(k+1,VAR_P,b), solt)
                    diff_k = diff_k + (fe - fw)*blk%d1z(k,VAR_P,b)

                    dcoef = nu + sigw*0.5d0*(turb%nut(i-1,j,k,b) + turb%nut(i,j,k,b))
                    fw = merge(0.0d0, dcoef*(sst%omg(i,j,k,b) - sst%omg(i-1,j,k,b)) &
                        *turb%inv_dx(i,VAR_P,b), solw)
                    dcoef = nu + sigw*0.5d0*(turb%nut(i,j,k,b) + turb%nut(i+1,j,k,b))
                    fe = merge(0.0d0, dcoef*(sst%omg(i+1,j,k,b) - sst%omg(i,j,k,b)) &
                        *turb%inv_dx(i+1,VAR_P,b), sole)
                    diff_w = diff_w + (fe - fw)*blk%d1x(i,VAR_P,b)
                    dcoef = nu + sigw*0.5d0*(turb%nut(i,j-1,k,b) + turb%nut(i,j,k,b))
                    fw = merge(0.0d0, dcoef*(sst%omg(i,j,k,b) - sst%omg(i,j-1,k,b)) &
                        *turb%inv_dy(j,VAR_P,b), sols)
                    dcoef = nu + sigw*0.5d0*(turb%nut(i,j,k,b) + turb%nut(i,j+1,k,b))
                    fe = merge(0.0d0, dcoef*(sst%omg(i,j+1,k,b) - sst%omg(i,j,k,b)) &
                        *turb%inv_dy(j+1,VAR_P,b), soln)
                    diff_w = diff_w + (fe - fw)*blk%d1y(j,VAR_P,b)
                    dcoef = nu + sigw*0.5d0*(turb%nut(i,j,k-1,b) + turb%nut(i,j,k,b))
                    fw = merge(0.0d0, dcoef*(sst%omg(i,j,k,b) - sst%omg(i,j,k-1,b)) &
                        *turb%inv_dz(k,VAR_P,b), solb)
                    dcoef = nu + sigw*0.5d0*(turb%nut(i,j,k,b) + turb%nut(i,j,k+1,b))
                    fe = merge(0.0d0, dcoef*(sst%omg(i,j,k+1,b) - sst%omg(i,j,k,b)) &
                        *turb%inv_dz(k+1,VAR_P,b), solt)
                    diff_w = diff_w + (fe - fw)*blk%d1z(k,VAR_P,b)

                    ! T4 gamma-Re_thetat transition (Langtry & Menter 2009):
                    ! the two extra scalars share every intermediate above
                    ! (gradients, S, F-blends, face masks). All quantities
                    ! read START-of-substage values (the explicit-RK stance;
                    ! OpenFOAM's segregated sweep uses the freshly solved
                    ! Re_thetat~ in the gamma equation — same fixed point).
                    ! gamma_eff = gamma; the separation-induced gamma_sep
                    ! branch (LM Eq. 18) is a SEPARATE later increment.
                    dkfac = 1.0d0
                    if (transition) then
                        gv = sst%gam(i,j,k,b)
                        rv = sst%ret(i,j,k,b)
                        ! Vorticity magnitude sqrt(2 W:W) and cell-centre
                        ! speed/streamline acceleration dU/ds = u_i u_j
                        ! du_i/dx_j / Us^2 for the correlations.
                        omgmag = sqrt((g12 - g21)**2 + (g13 - g31)**2 + (g23 - g32)**2)
                        uc = 0.5d0*(blk%q(i,j,k,VAR_U,b) + blk%q(i+1,j,k,VAR_U,b))
                        vc = 0.5d0*(blk%q(i,j,k,VAR_V,b) + blk%q(i,j+1,k,VAR_V,b))
                        wc = 0.5d0*(blk%q(i,j,k,VAR_W,b) + blk%q(i,j,k+1,VAR_W,b))
                        usmag = max(sqrt(uc*uc + vc*vc + wc*wc), LM_DELTAU)
                        dusds = (uc*(uc*g11 + vc*g12 + wc*g13) &
                               + vc*(uc*g21 + vc*g22 + wc*g23) &
                               + wc*(uc*g31 + vc*g32 + wc*g33))/(usmag*usmag)

                        rt = max(kv, 0.0d0)/(nu*wv)
                        rev = y*y*smag/nu
                        tuloc = 100.0d0*sqrt(2.0d0/3.0d0*max(kv, 0.0d0))/usmag

                        ! Re_thetat~ relaxation to the local equilibrium
                        ! correlation, coefficient c_thetat/t with
                        ! t = 500 nu/Us^2, gated off inside boundary layers
                        ! by (1 - F_thetat); the -coef*ret part is the
                        ! point-implicit sink (OpenFOAM's fvm::Sp).
                        ret0 = lm_rethetat0(tuloc, dusds, nu, usmag)
                        fthetat = lm_fthetat(gv, rv, usmag, omgmag, wv, y, nu)
                        tcoef = LM_CTHETAT*usmag*usmag/(500.0d0*nu)*(1.0d0 - fthetat)

                        ! gamma sources: production ca1 Flength S
                        ! sqrt(gamma Fonset) (1 - ce1 gamma) and destruction
                        ! ca2 Fturb Omega gamma (ce2 gamma - 1), both in the
                        ! OpenFOAM split +P - ce1 P g and +E - ce2 E g:
                        ! P, E >= 0, so the explicit parts are pure sources
                        ! and the implicit coefficient stays nonnegative on
                        ! BOTH sides of the destruction sign flip at
                        ! gamma = 1/ce2 (Patankar-safe by construction).
                        rethetac = lm_rethetac(rv)
                        flength = lm_flength(rv, y, wv, nu)
                        fonset = lm_fonset(rev, rethetac, rt)
                        fturbv = lm_fturb(rt)
                        pgam = LM_CA1*flength*smag*sqrt(max(gv, 0.0d0)*fonset)
                        egam = LM_CA2*omgmag*fturbv*max(gv, 0.0d0)

                        ! First-order upwind convection (STEP-0 decision in
                        ! the module header) and face-masked diffusion, the
                        ! k/omega pattern; gamma diffuses with nu +
                        ! nut/sigma_f, Re_thetat~ with sigma_thetat (nu + nut).
                        conv_g = (ue*merge(sst%gam(i,j,k,b), sst%gam(i+1,j,k,b), ue > 0.0d0) &
                                - uw*merge(sst%gam(i-1,j,k,b), sst%gam(i,j,k,b), uw > 0.0d0)) &
                                  *blk%d1x(i,VAR_P,b) &
                               + (vn*merge(sst%gam(i,j,k,b), sst%gam(i,j+1,k,b), vn > 0.0d0) &
                                - vs*merge(sst%gam(i,j-1,k,b), sst%gam(i,j,k,b), vs > 0.0d0)) &
                                  *blk%d1y(j,VAR_P,b) &
                               + (wt*merge(sst%gam(i,j,k,b), sst%gam(i,j,k+1,b), wt > 0.0d0) &
                                - wb*merge(sst%gam(i,j,k-1,b), sst%gam(i,j,k,b), wb > 0.0d0)) &
                                  *blk%d1z(k,VAR_P,b)
                        conv_r = (ue*merge(sst%ret(i,j,k,b), sst%ret(i+1,j,k,b), ue > 0.0d0) &
                                - uw*merge(sst%ret(i-1,j,k,b), sst%ret(i,j,k,b), uw > 0.0d0)) &
                                  *blk%d1x(i,VAR_P,b) &
                               + (vn*merge(sst%ret(i,j,k,b), sst%ret(i,j+1,k,b), vn > 0.0d0) &
                                - vs*merge(sst%ret(i,j-1,k,b), sst%ret(i,j,k,b), vs > 0.0d0)) &
                                  *blk%d1y(j,VAR_P,b) &
                               + (wt*merge(sst%ret(i,j,k,b), sst%ret(i,j,k+1,b), wt > 0.0d0) &
                                - wb*merge(sst%ret(i,j,k-1,b), sst%ret(i,j,k,b), wb > 0.0d0)) &
                                  *blk%d1z(k,VAR_P,b)

                        diff_g = 0.0d0
                        diff_r = 0.0d0
                        dcoef = nu + 0.5d0*(turb%nut(i-1,j,k,b) + turb%nut(i,j,k,b))/LM_SIGMAF
                        fw = merge(0.0d0, dcoef*(sst%gam(i,j,k,b) - sst%gam(i-1,j,k,b)) &
                            *turb%inv_dx(i,VAR_P,b), solw)
                        dcoef = nu + 0.5d0*(turb%nut(i,j,k,b) + turb%nut(i+1,j,k,b))/LM_SIGMAF
                        fe = merge(0.0d0, dcoef*(sst%gam(i+1,j,k,b) - sst%gam(i,j,k,b)) &
                            *turb%inv_dx(i+1,VAR_P,b), sole)
                        diff_g = diff_g + (fe - fw)*blk%d1x(i,VAR_P,b)
                        dcoef = nu + 0.5d0*(turb%nut(i,j-1,k,b) + turb%nut(i,j,k,b))/LM_SIGMAF
                        fw = merge(0.0d0, dcoef*(sst%gam(i,j,k,b) - sst%gam(i,j-1,k,b)) &
                            *turb%inv_dy(j,VAR_P,b), sols)
                        dcoef = nu + 0.5d0*(turb%nut(i,j,k,b) + turb%nut(i,j+1,k,b))/LM_SIGMAF
                        fe = merge(0.0d0, dcoef*(sst%gam(i,j+1,k,b) - sst%gam(i,j,k,b)) &
                            *turb%inv_dy(j+1,VAR_P,b), soln)
                        diff_g = diff_g + (fe - fw)*blk%d1y(j,VAR_P,b)
                        dcoef = nu + 0.5d0*(turb%nut(i,j,k-1,b) + turb%nut(i,j,k,b))/LM_SIGMAF
                        fw = merge(0.0d0, dcoef*(sst%gam(i,j,k,b) - sst%gam(i,j,k-1,b)) &
                            *turb%inv_dz(k,VAR_P,b), solb)
                        dcoef = nu + 0.5d0*(turb%nut(i,j,k,b) + turb%nut(i,j,k+1,b))/LM_SIGMAF
                        fe = merge(0.0d0, dcoef*(sst%gam(i,j,k+1,b) - sst%gam(i,j,k,b)) &
                            *turb%inv_dz(k+1,VAR_P,b), solt)
                        diff_g = diff_g + (fe - fw)*blk%d1z(k,VAR_P,b)

                        ! Re_thetat~ diffusion: sigma_thetat = 2 makes its
                        ! diffusivity TWICE the momentum value the Peclet dt
                        ! controller budgets for (k/omega/gamma all have
                        ! sigma <= 1 and live inside that validated margin),
                        ! so a fully explicit treatment is unstable at the
                        ! controller's dt (observed: y-checkerboard growth to
                        ! 1e6 within ~40 steps on the lam30t gate). RK
                        ! hardening in the T2 cross-diffusion spirit: the
                        ! own-cell (diagonal) part of the operator goes into
                        ! the point-implicit denominator (diag_r), the
                        ! neighbour part stays explicit — unconditionally
                        ! stable, positivity-preserving, same steady state
                        ! (and OpenFOAM's fvm::laplacian is fully implicit
                        ! anyway). diff_r accumulates the full operator;
                        ! diag_r the unmasked-face own-cell coefficients.
                        diag_r = 0.0d0
                        dcoef = LM_SIGMATHETAT*(nu + 0.5d0*(turb%nut(i-1,j,k,b) + turb%nut(i,j,k,b)))
                        fw = merge(0.0d0, dcoef*(sst%ret(i,j,k,b) - sst%ret(i-1,j,k,b)) &
                            *turb%inv_dx(i,VAR_P,b), solw)
                        diag_r = diag_r + merge(0.0d0, dcoef*turb%inv_dx(i,VAR_P,b), solw) &
                            *blk%d1x(i,VAR_P,b)
                        dcoef = LM_SIGMATHETAT*(nu + 0.5d0*(turb%nut(i,j,k,b) + turb%nut(i+1,j,k,b)))
                        fe = merge(0.0d0, dcoef*(sst%ret(i+1,j,k,b) - sst%ret(i,j,k,b)) &
                            *turb%inv_dx(i+1,VAR_P,b), sole)
                        diag_r = diag_r + merge(0.0d0, dcoef*turb%inv_dx(i+1,VAR_P,b), sole) &
                            *blk%d1x(i,VAR_P,b)
                        diff_r = diff_r + (fe - fw)*blk%d1x(i,VAR_P,b)
                        dcoef = LM_SIGMATHETAT*(nu + 0.5d0*(turb%nut(i,j-1,k,b) + turb%nut(i,j,k,b)))
                        fw = merge(0.0d0, dcoef*(sst%ret(i,j,k,b) - sst%ret(i,j-1,k,b)) &
                            *turb%inv_dy(j,VAR_P,b), sols)
                        diag_r = diag_r + merge(0.0d0, dcoef*turb%inv_dy(j,VAR_P,b), sols) &
                            *blk%d1y(j,VAR_P,b)
                        dcoef = LM_SIGMATHETAT*(nu + 0.5d0*(turb%nut(i,j,k,b) + turb%nut(i,j+1,k,b)))
                        fe = merge(0.0d0, dcoef*(sst%ret(i,j+1,k,b) - sst%ret(i,j,k,b)) &
                            *turb%inv_dy(j+1,VAR_P,b), soln)
                        diag_r = diag_r + merge(0.0d0, dcoef*turb%inv_dy(j+1,VAR_P,b), soln) &
                            *blk%d1y(j,VAR_P,b)
                        diff_r = diff_r + (fe - fw)*blk%d1y(j,VAR_P,b)
                        dcoef = LM_SIGMATHETAT*(nu + 0.5d0*(turb%nut(i,j,k-1,b) + turb%nut(i,j,k,b)))
                        fw = merge(0.0d0, dcoef*(sst%ret(i,j,k,b) - sst%ret(i,j,k-1,b)) &
                            *turb%inv_dz(k,VAR_P,b), solb)
                        diag_r = diag_r + merge(0.0d0, dcoef*turb%inv_dz(k,VAR_P,b), solb) &
                            *blk%d1z(k,VAR_P,b)
                        dcoef = LM_SIGMATHETAT*(nu + 0.5d0*(turb%nut(i,j,k,b) + turb%nut(i,j,k+1,b)))
                        fe = merge(0.0d0, dcoef*(sst%ret(i,j,k+1,b) - sst%ret(i,j,k,b)) &
                            *turb%inv_dz(k+1,VAR_P,b), solt)
                        diag_r = diag_r + merge(0.0d0, dcoef*turb%inv_dz(k+1,VAR_P,b), solt) &
                            *blk%d1z(k,VAR_P,b)
                        diff_r = diff_r + (fe - fw)*blk%d1z(k,VAR_P,b)

                        rhsg = pgam + egam - conv_g + diff_g
                        gnew = (gv + dt_alpha*rhsg + dt_beta*sst%gamoldrhs(i,j,k,b)) &
                             /(1.0d0 + dtsub*(LM_CE1*pgam + LM_CE2*egam))
                        sst%gams(i,j,k,b) = min(max(gnew, 0.0d0), 1.0d0)
                        sst%gamoldrhs(i,j,k,b) = rhsg

                        ! Explicit rhs carries the diffusion NEIGHBOUR part
                        ! (diff_r + diag_r rv); the diagonal is implicit.
                        rhsr = tcoef*ret0 - conv_r + diff_r + diag_r*rv
                        rnew = (rv + dt_alpha*rhsr + dt_beta*sst%retoldrhs(i,j,k,b)) &
                             /(1.0d0 + dtsub*(tcoef + diag_r))
                        sst%rets(i,j,k,b) = max(rnew, 0.0d0)
                        sst%retoldrhs(i,j,k,b) = rhsr

                        ! Coupling into SST: P_k -> gamma_eff P_k and the k
                        ! destruction scaled by min(max(gamma_eff,0.1),1)
                        ! (dkfac multiplies the point-implicit denominator;
                        ! it is exactly 1.0 when transition is off, which is
                        ! bit-exact).
                        gameff = gv
                        dkfac = min(max(gameff, 0.1d0), 1.0d0)
                        pk = gameff*pk
                    end if

                    ! Explicit RHS; sinks are point-implicit (Patankar). The
                    ! cross-diffusion 2(1-F1) sigma_w2/omega grad k . grad
                    ! omega is SPLIT by sign: its negative part becomes an
                    ! extra implicit sink coefficient, because integrating it
                    ! explicitly can drive omega through zero in one substage
                    ! (steep grad omega against the pinned wall value), and a
                    ! floored omega flips F1 -> 0 through the CD_komega
                    ! branch of arg1, re-enabling the term with its 1/omega
                    ! amplification -- an explosive feedback (observed: omega
                    ! -> 1e150 within steps on the laminar gate).
                    ! The explicit positive part is additionally rate-limited
                    ! to at most doubling omega per substage: at a floored/
                    ! tiny omega the CD_komega branch of arg1 flips F1 -> 0
                    ! and the 1/omega amplification would otherwise cascade.
                    cross = 2.0d0*(1.0d0 - f1)*SST_SIGW2/wv*gkgw
                    rhsk = pk - conv_k + diff_k
                    rhsw = alpha_b*s2 - conv_w + diff_w &
                         + min(max(cross, 0.0d0), wv/max(dtsub, 1.0d-30))

                    knew = (kv + dt_alpha*rhsk + dt_beta*sst%koldrhs(i,j,k,b)) &
                         /(1.0d0 + dtsub*SST_BETA_STAR*wv*dkfac)
                    knew = max(knew, 0.0d0)
                    sst%ks(i,j,k,b) = knew
                    sst%koldrhs(i,j,k,b) = rhsk

                    if (pinned_omega) then
                        sst%omgs(i,j,k,b) = wv
                        sst%omgoldrhs(i,j,k,b) = 0.0d0
                    else
                        wnew = (wv + dt_alpha*rhsw + dt_beta*sst%omgoldrhs(i,j,k,b)) &
                             /(1.0d0 + dtsub*(beta_b*wv + max(-cross, 0.0d0)/wv))
                        sst%omgs(i,j,k,b) = max(wnew, OMEGA_MIN)
                        sst%omgoldrhs(i,j,k,b) = rhsw
                    end if
                end do
            end do
        end do
        end do
        !$omp end target teams distribute parallel do
    end subroutine rans_transport_kernel

end module rans
