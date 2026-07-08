!--------------------------!
!                          !
!      RANS (k-omega       !
!      SST) module         !
!                          !
!--------------------------!
!
! k-omega SST (Menter 2003 constants), resolved-wall mode (IDDES phase T2,
! docs/next_session_iddes.md). The module owns the SST GEOMETRY state
! (phase T1: wall distance + IBM wall cells, below) and the TRANSPORT state
! (k, omega + their low-storage RK3 rhs history). The model is a producer
! of the one cell-centred turb%nut; everything downstream of nut (halo
! exchange, momentum correction, dt limit, io) is the untouched consumer
! chain.
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
!     the block decomposition. Developed-flow gates are source/sink
!     dominated; revisit before the T4 transition fronts.
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
    use :: boundary, only: boundary_type, boundary_face_id, NFACES
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
        ! Per-face no-slip flags for the cell-centred scalar BCs (index =
        ! boundary_face_id).
        integer(C_INT) :: facewall(NFACES) = 0_C_INT
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
    end subroutine init_rans_geometry

    subroutine destroy_rans_geometry(sst)
        type(sst_type), intent(inout) :: sst

        if (allocated(sst%dwall)) deallocate(sst%dwall)
        if (allocated(sst%yeff)) deallocate(sst%yeff)
        if (allocated(sst%wallcell)) deallocate(sst%wallcell)
        if (allocated(sst%k)) deallocate(sst%k, sst%omg, sst%ks, sst%omgs, &
            sst%koldrhs, sst%omgoldrhs, sst%domwall)
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
            !$omp target enter data map(to: sst%koldrhs, sst%omgoldrhs, sst%domwall)
        end if
    end subroutine enter_rans_data

    subroutine exit_rans_data(sst)
        type(sst_type), intent(inout) :: sst

        if (.not. allocated(sst%dwall)) return

        if (allocated(sst%k)) then
            !$omp target exit data map(delete: sst%koldrhs, sst%omgoldrhs, sst%domwall)
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

    ! A domain face is a wall iff its direction is non-periodic and both
    ! tangential velocity components have Dirichlet (no-slip) conditions;
    ! Neumann tangential faces (free slip, symmetry) carry no wall layer.
    logical function domain_face_is_wall(bc, dir, side) result(is_wall)
        type(boundary_type), intent(in) :: bc
        integer, intent(in) :: dir, side

        integer :: face_id, var
        integer, parameter :: DIRICHLET = 0

        is_wall = .false.
        if (bc%isPeriodic(dir)) return

        face_id = boundary_face_id(dir, side)
        is_wall = .true.
        do var = int(VAR_U), int(VAR_W)
            if (var == dir) cycle   ! the normal component does not decide no-slip
            if (bc%faceBcType(var, face_id) /= DIRICHLET) is_wall = .false.
        end do
    end function domain_face_is_wall

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
    subroutine init_rans_transport(sst, dns, blk, bc, has_terminal)
        type(sst_type), intent(inout) :: sst
        type(dns_type), intent(in) :: dns
        type(block_set_type), intent(in) :: blk
        type(boundary_type), intent(in) :: bc
        logical, intent(in) :: has_terminal

        integer :: i, j, k, b, nx, ny, nz, dir, side
        real(C_DOUBLE) :: umag2, uc, vc, wc, tu_frac, nu, kin
        logical :: found_k, found_omg

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
        end if

        call set_constrained_cells_host(sst, dns, blk)
        sst%ks = sst%k
        sst%omgs = sst%omg
        sst%transport_built = .true.
    end subroutine init_rans_transport

    ! Menter viscous-limb wall omega at the effective distance.
    pure real(C_DOUBLE) function omega_wall_value(nu, y) result(w)
!$omp declare target
        real(C_DOUBLE), intent(in) :: nu, y

        w = 6.0d0*nu/(SST_BETA1*y*y)
    end function omega_wall_value

    ! Host-side constrained-cell values (init/restart time, before the
    ! device maps exist): solid cells carry k = 0 and a benign pinned
    ! omega; IBM wall cells and domain-wall rows carry the pinned omega.
    subroutine set_constrained_cells_host(sst, dns, blk)
        type(sst_type), intent(inout) :: sst
        type(dns_type), intent(in) :: dns
        type(block_set_type), intent(in) :: blk

        integer :: i, j, k, b, nx, ny, nz
        real(C_DOUBLE) :: nu
        logical :: pinned

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))
        nu = 1.0d0/dns%re

        do b = 1, int(blk%nBlocks)
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    pinned = sst%wallcell(i,j,k,b) /= WALL_CELL_FLUID .or. &
                             sst%domwall(i,j,k,b) /= 0_C_SIGNED_CHAR
                    if (sst%wallcell(i,j,k,b) == WALL_CELL_SOLID) sst%k(i,j,k,b) = 0.0d0
                    if (pinned) sst%omg(i,j,k,b) = omega_wall_value(nu, sst%yeff(i,j,k,b))
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

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))
        nBlocks = int(blk%nBlocks)
        nu = 1.0d0/dns%re

        !$omp target teams distribute parallel do collapse(4) &
        !$omp& map(to: nu, sst%wallcell, sst%domwall, sst%yeff) &
        !$omp& map(tofrom: sst%k, sst%omg) &
        !$omp& private(i,j,k,b)
        do b = 1, nBlocks
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    if (sst%wallcell(i,j,k,b) == WALL_CELL_SOLID) sst%k(i,j,k,b) = 0.0d0
                    if (sst%wallcell(i,j,k,b) /= WALL_CELL_FLUID .or. &
                        sst%domwall(i,j,k,b) /= 0_C_SIGNED_CHAR) then
                        sst%omg(i,j,k,b) = omega_wall_value(nu, sst%yeff(i,j,k,b))
                    end if
                end do
            end do
        end do
        end do
        !$omp end target teams distribute parallel do
    end subroutine rans_set_constrained_cells

    ! Cell-centred scalar ghosts at physical domain faces: k is Dirichlet 0
    ! at no-slip walls (mirror) and zero-gradient elsewhere; omega is always
    ! zero-gradient (its wall value is the pinned first cell, not a ghost).
    subroutine rans_apply_scalar_bcs(sst, blk, bc)
        type(sst_type), intent(inout) :: sst
        type(block_set_type), intent(inout) :: blk
        type(boundary_type), intent(in) :: bc

        integer :: n, npts, b, i, j, k, face_id, dir, side
        integer :: ghost_idx, interior_idx_dir
        integer :: gi(3), ii(3)
        integer(C_INT) :: local_n(1:3), facewall(NFACES)

        npts = int(bc%nTotal)
        if (npts <= 0) return
        local_n = blk%nb(1:3)
        facewall = sst%facewall

        !$omp target teams distribute parallel do &
        !$omp& map(to: npts, local_n(1:3), facewall(1:NFACES), &
        !$omp& bc%pointFace(1:npts), bc%slot(1:npts), bc%i(1:npts), bc%j(1:npts), bc%k(1:npts)) &
        !$omp& map(tofrom: sst%k, sst%omg) &
        !$omp& private(n,b,i,j,k,face_id,dir,side,ghost_idx,interior_idx_dir,gi,ii)
        do n = 1, npts
            face_id = int(bc%pointFace(n))
            b = int(bc%slot(n))
            dir = (face_id + 1)/2
            side = modulo(face_id - 1, 2)
            i = int(bc%i(n))
            j = int(bc%j(n))
            k = int(bc%k(n))

            if (side == 0) then
                ghost_idx = 0
                interior_idx_dir = 1
            else
                ghost_idx = int(local_n(dir)) + 1
                interior_idx_dir = int(local_n(dir))
            end if
            gi = [i, j, k]
            ii = gi
            gi(dir) = ghost_idx
            ii(dir) = interior_idx_dir

            if (facewall(face_id) == 1_C_INT) then
                sst%k(gi(1),gi(2),gi(3),b) = -sst%k(ii(1),ii(2),ii(3),b)
            else
                sst%k(gi(1),gi(2),gi(3),b) = sst%k(ii(1),ii(2),ii(3),b)
            end if
            sst%omg(gi(1),gi(2),gi(3),b) = sst%omg(ii(1),ii(2),ii(3),b)
        end do
        !$omp end target teams distribute parallel do
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
        real(C_DOUBLE) :: nu, y, kv, wv, s2, smag, arg2, f2
        real(C_DOUBLE) :: g11, g12, g13, g21, g22, g23, g31, g32, g33

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))
        nBlocks = int(blk%nBlocks)
        nu = 1.0d0/dns%re

        !$omp target teams distribute parallel do collapse(4) &
        !$omp& map(to: nu, sst%k, sst%omg, sst%yeff, sst%wallcell, &
        !$omp& blk%q, blk%d1x, blk%d1y, blk%d1z, &
        !$omp& turb%d1xm, turb%d1x0, turb%d1xp, turb%d1ym, turb%d1y0, turb%d1yp, &
        !$omp& turb%d1zm, turb%d1z0, turb%d1zp, &
        !$omp& turb%p_from_u_x, turb%p_from_v_y, turb%p_from_w_z) &
        !$omp& map(tofrom: nut) &
        !$omp& private(i,j,k,b,y,kv,wv,s2,smag,arg2,f2, &
        !$omp& g11,g12,g13,g21,g22,g23,g31,g32,g33)
        do b = 1, nBlocks
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    if (sst%wallcell(i,j,k,b) /= WALL_CELL_FLUID) then
                        nut(i,j,k,b) = 0.0d0
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
        call rans_apply_scalar_bcs(sst, blk, bc)
        call exchange_scalar_halos(c, sst%k, blk)
        call exchange_scalar_halos(c, sst%omg, blk)
        call rans_assemble_nut(sst, turb, blk, dns, turb%nut)
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
        call rans_apply_scalar_bcs(sst, blk, bc)
        call exchange_scalar_halos(c, sst%k, blk)
        call exchange_scalar_halos(c, sst%omg, blk)
        call rans_transport_kernel(sst, turb, blk, dns, ibm, dt_alpha, dt_beta)
        call rans_copyback(sst, blk)
        call rans_assemble_nut(sst, turb, blk, dns, turb%nut)
    end subroutine rans_substage

    subroutine rans_copyback(sst, blk)
        type(sst_type), intent(inout) :: sst
        type(block_set_type), intent(in) :: blk

        integer :: i, j, k, b, nx, ny, nz, nBlocks

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))
        nBlocks = int(blk%nBlocks)

        !$omp target teams distribute parallel do collapse(4) &
        !$omp& map(to: sst%ks, sst%omgs) &
        !$omp& map(tofrom: sst%k, sst%omg) &
        !$omp& private(i,j,k,b)
        do b = 1, nBlocks
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    sst%k(i,j,k,b) = sst%ks(i,j,k,b)
                    sst%omg(i,j,k,b) = sst%omgs(i,j,k,b)
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
        logical :: pinned_omega, ibm_wall
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
        logical :: solw, sole, sols, soln, solb, solt

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))
        nBlocks = int(blk%nBlocks)
        nu = 1.0d0/dns%re
        dtsub = dt_alpha + dt_beta
        solid_threshold = SOLID_FACE_THRESHOLD

        !$omp target teams distribute parallel do collapse(4) &
        !$omp& map(to: nu, dtsub, dt_alpha, dt_beta, solid_threshold, &
        !$omp& sst%k, sst%omg, sst%yeff, sst%wallcell, sst%domwall, &
        !$omp& blk%q, blk%d1x, blk%d1y, blk%d1z, ibm%coef, &
        !$omp& turb%nut, turb%inv_dx, turb%inv_dy, turb%inv_dz, &
        !$omp& turb%d1xm, turb%d1x0, turb%d1xp, turb%d1ym, turb%d1y0, turb%d1yp, &
        !$omp& turb%d1zm, turb%d1z0, turb%d1zp, &
        !$omp& turb%p_from_u_x, turb%p_from_v_y, turb%p_from_w_z) &
        !$omp& map(tofrom: sst%ks, sst%omgs, sst%koldrhs, sst%omgoldrhs) &
        !$omp& private(i,j,k,b,marker,pinned_omega,ibm_wall,y,kv,wv,nutc, &
        !$omp& g11,g12,g13,g21,g22,g23,g31,g32,g33,s2,smag, &
        !$omp& gkx,gky,gkz,gwx,gwy,gwz,gkgw,cdkw,arg1,f1,arg2,f2,nutloc,cross, &
        !$omp& sigk,sigw,beta_b,alpha_b,pk,uw,ue,vs,vn,wb,wt, &
        !$omp& conv_k,conv_w,diff_k,diff_w,fw,fe,dcoef,rhsk,rhsw,knew,wnew, &
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
                         /(1.0d0 + dtsub*SST_BETA_STAR*wv)
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
