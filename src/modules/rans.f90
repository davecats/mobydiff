!--------------------------!
!                          !
!      RANS (k-omega       !
!      SST) module         !
!                          !
!--------------------------!
!
! IDDES phase T1 (docs/next_session_iddes.md): the SST GEOMETRY state only —
! the regularized wall distance and the IBM wall-cell classification the
! k-omega SST wall treatment will consume from phase T2 on. No transport
! arrays and no coupling into the solver exist yet: building this state is
! init-only and nothing downstream reads it, so any run remains bit-exact
! whether or not a [rans] section is configured.
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
    use :: blocks, only: block_set_type
    use :: boundary, only: boundary_type, boundary_face_id
    use :: ibmm, only: ibm_type, isInBody
    use :: walldist, only: walldist_type, build_walldist, destroy_walldist, &
        walldist_distance
    use :: io, only: read_dwall_blocks, write_rans_geometry_file
    use :: comm, only: comm_type, comm_allreduce_sum
    implicit none

    private
    public :: sst_type
    public :: init_rans_geometry, destroy_rans_geometry
    public :: enter_rans_data, exit_rans_data
    public :: write_rans_geometry

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
        ! The transport state (k, omega, ...) arrives in phase T2.
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
        sst%geometry_built = .false.
    end subroutine destroy_rans_geometry

    subroutine enter_rans_data(sst)
        type(sst_type), intent(inout) :: sst

        if (.not. allocated(sst%dwall)) return

        !$omp target enter data map(to: sst)
        !$omp target enter data map(to: sst%dwall, sst%yeff, sst%wallcell)
    end subroutine enter_rans_data

    subroutine exit_rans_data(sst)
        type(sst_type), intent(inout) :: sst

        if (.not. allocated(sst%dwall)) return

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

end module rans
