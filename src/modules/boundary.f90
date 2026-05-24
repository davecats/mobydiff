module boundary
    use, intrinsic :: iso_c_binding
    use :: init, only: dns_type, grid_type, field_type, VAR_U, VAR_V, VAR_W, VAR_P
    implicit none

    integer(C_INT), parameter :: DIR_X = 1_C_INT
    integer(C_INT), parameter :: DIR_Y = 2_C_INT
    integer(C_INT), parameter :: DIR_Z = 3_C_INT
    integer(C_INT), parameter :: SIDE_MIN = 0_C_INT
    integer(C_INT), parameter :: SIDE_MAX = 1_C_INT
    integer, parameter :: NFACES = 6

    type :: boundary_type
        logical(C_BOOL) :: isPeriodic(1:3)
        integer(C_INT) :: nTotal = 0_C_INT

        ! Active physical boundary points. Periodic/MPI halos are handled by comm.f90.
        integer(C_INT), allocatable :: pointFace(:), i(:), j(:), k(:)

        ! Boundary data use the same variable order as field_type: u, v, w, p.
        integer(C_INT) :: faceBcType(VAR_U:VAR_P,1:NFACES) = 0_C_INT
        real(C_DOUBLE) :: faceBcValue(VAR_U:VAR_P,1:NFACES) = 0.0d0
    end type boundary_type

contains

    subroutine init_bc(bc)
        type(boundary_type), intent(inout) :: bc

        call destroy_boundary_faces(bc)
        bc%isPeriodic(1:3) = .true.
        bc%faceBcType = 0_C_INT
        bc%faceBcValue = 0.0d0
        bc%faceBcType(VAR_P,:) = 1_C_INT
    end subroutine init_bc

    subroutine init_boundary_faces(bc, dns)
        type(boundary_type), intent(inout) :: bc
        type(dns_type), intent(in) :: dns
        integer :: nx, ny, nz
        integer :: dir, side, face_id, pos, total

        nx = int(dns%localSize(1,2))
        ny = int(dns%localSize(2,2))
        nz = int(dns%localSize(3,2))

        call destroy_boundary_faces(bc)

        total = 0
        do dir = 1, 3
            do side = 0, 1
                if (is_physical_boundary(bc, dns, dir, side)) then
                    total = total + boundary_face_n(dir, nx, ny, nz)
                end if
            end do
        end do

        bc%nTotal = int(total, C_INT)
        if (total <= 0) return

        allocate(bc%pointFace(total), bc%i(total), bc%j(total), bc%k(total))

        pos = 0
        do dir = 1, 3
            do side = 0, 1
                if (.not. is_physical_boundary(bc, dns, dir, side)) cycle

                face_id = boundary_face_id(dir, side)
                call append_boundary_face_points(bc, face_id, pos, nx, ny, nz)
            end do
        end do
    end subroutine init_boundary_faces

    logical function is_physical_boundary(bc, dns, dir, side)
        type(boundary_type), intent(in) :: bc
        type(dns_type), intent(in) :: dns
        integer, intent(in) :: dir, side

        is_physical_boundary = (.not. bc%isPeriodic(dir)) .and. &
            ((side == SIDE_MIN .and. dns%localSize(dir,0) == 1_C_INT) .or. &
             (side == SIDE_MAX .and. dns%localSize(dir,1) == dns%globalSize(dir)))
    end function is_physical_boundary

    subroutine enter_boundary_data(bc)
        type(boundary_type), intent(inout) :: bc
        integer :: npts

#ifdef USE_OPENMP_OFFLOAD
        npts = int(bc%nTotal)
        !$omp target enter data map(to: bc)
        !$omp target enter data map(to: bc%faceBcType, bc%faceBcValue)
        if (npts > 0) then
            !$omp target enter data map(to: bc%pointFace(1:npts), bc%i(1:npts), bc%j(1:npts), bc%k(1:npts))
        end if
#endif
    end subroutine enter_boundary_data

    subroutine exit_boundary_data(bc)
        type(boundary_type), intent(inout) :: bc
        integer :: npts

#ifdef USE_OPENMP_OFFLOAD
        npts = int(bc%nTotal)
        if (npts > 0) then
            !$omp target exit data map(delete: bc%pointFace(1:npts), bc%i(1:npts), bc%j(1:npts), bc%k(1:npts))
        end if
        !$omp target exit data map(delete: bc%faceBcType, bc%faceBcValue)
        !$omp target exit data map(delete: bc)
#endif
    end subroutine exit_boundary_data

    subroutine destroy_boundary_faces(bc)
        type(boundary_type), intent(inout) :: bc

        if (allocated(bc%pointFace)) deallocate(bc%pointFace)
        if (allocated(bc%i)) deallocate(bc%i)
        if (allocated(bc%j)) deallocate(bc%j)
        if (allocated(bc%k)) deallocate(bc%k)
        bc%nTotal = 0_C_INT
    end subroutine destroy_boundary_faces

    subroutine append_boundary_face_points(bc, face_id, pos, nx, ny, nz)
        type(boundary_type), intent(inout) :: bc
        integer, intent(in) :: face_id, nx, ny, nz
        integer, intent(inout) :: pos
        integer :: i, j, k, dir, side

        dir = boundary_face_dir(face_id)
        side = boundary_face_side(face_id)

        select case (dir)
        case (DIR_X)
            i = merge(1, nx, side == SIDE_MIN)
            do k = 1, nz
                do j = 1, ny
                    pos = pos + 1
                    bc%pointFace(pos) = int(face_id, C_INT)
                    bc%i(pos) = int(i, C_INT)
                    bc%j(pos) = int(j, C_INT)
                    bc%k(pos) = int(k, C_INT)
                end do
            end do
        case (DIR_Y)
            j = merge(1, ny, side == SIDE_MIN)
            do k = 1, nz
                do i = 1, nx
                    pos = pos + 1
                    bc%pointFace(pos) = int(face_id, C_INT)
                    bc%i(pos) = int(i, C_INT)
                    bc%j(pos) = int(j, C_INT)
                    bc%k(pos) = int(k, C_INT)
                end do
            end do
        case (DIR_Z)
            k = merge(1, nz, side == SIDE_MIN)
            do j = 1, ny
                do i = 1, nx
                    pos = pos + 1
                    bc%pointFace(pos) = int(face_id, C_INT)
                    bc%i(pos) = int(i, C_INT)
                    bc%j(pos) = int(j, C_INT)
                    bc%k(pos) = int(k, C_INT)
                end do
            end do
        case default
            error stop "invalid boundary direction"
        end select
    end subroutine append_boundary_face_points

    integer function boundary_face_n(dir, nx, ny, nz) result(n)
        integer, intent(in) :: dir, nx, ny, nz

        select case (dir)
        case (DIR_X)
            n = ny*nz
        case (DIR_Y)
            n = nx*nz
        case (DIR_Z)
            n = nx*ny
        case default
            error stop "invalid boundary direction"
        end select
    end function boundary_face_n

    integer function boundary_face_id(dir, side) result(face_id)
        integer, intent(in) :: dir, side

        face_id = 2*(dir - 1) + side + 1
    end function boundary_face_id

    integer function boundary_face_dir(face_id) result(dir)
        integer, intent(in) :: face_id

        dir = (face_id + 1)/2
    end function boundary_face_dir

    integer function boundary_face_side(face_id) result(side)
        integer, intent(in) :: face_id

        side = modulo(face_id - 1, 2)
    end function boundary_face_side

    subroutine apply_bc(f, dns, g, bc)
        type(field_type), intent(inout) :: f
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        type(boundary_type), intent(in) :: bc
        integer :: n, npts, i, j, k, face_id, var
        integer :: dir, side
        integer :: n_dir, ghost_idx, interior_idx_dir, face_idx_dir, neighbor_idx
        integer(C_INT) :: local_n(1:3)
        integer :: idx(3), interior_idx(3), face_idx(3)
        real(C_DOUBLE) :: dn

        npts = int(bc%nTotal)
        if (npts <= 0) return

        local_n = dns%localSize(1:3,2)

        !$omp target teams distribute parallel do &
        !$omp& map(to: npts, local_n(1:3), g%x, g%y, g%z, &
        !$omp& bc%pointFace(1:npts), bc%i(1:npts), bc%j(1:npts), bc%k(1:npts), &
        !$omp& bc%faceBcType(VAR_U:VAR_P,1:NFACES), bc%faceBcValue(VAR_U:VAR_P,1:NFACES)) &
        !$omp& map(tofrom: f%q) &
        !$omp& private(n,i,j,k,face_id,var,dir,side,n_dir, &
        !$omp& ghost_idx,interior_idx_dir,face_idx_dir,neighbor_idx, &
        !$omp& idx,interior_idx,face_idx,dn)
        do n = 1, npts
            face_id = int(bc%pointFace(n))
            dir = (face_id + 1)/2
            side = modulo(face_id - 1, 2)
            i = int(bc%i(n))
            j = int(bc%j(n))
            k = int(bc%k(n))

            idx = [i, j, k]
            interior_idx = idx
            face_idx = idx

            n_dir = int(local_n(dir))
            if (side == SIDE_MIN) then
                ghost_idx = 0
                interior_idx_dir = 1
                face_idx_dir = 1
                neighbor_idx = 2
            else
                ghost_idx = n_dir + 1
                interior_idx_dir = n_dir
                face_idx_dir = n_dir + 1
                neighbor_idx = n_dir
            end if

            do var = VAR_U, VAR_P
                if (var == dir) then
                    face_idx(dir) = face_idx_dir
                    interior_idx(dir) = neighbor_idx
                    select case (dir)
                    case (DIR_X)
                        dn = g%x(face_idx(1),var) - g%x(interior_idx(1),var)
                    case (DIR_Y)
                        dn = g%y(face_idx(2),var) - g%y(interior_idx(2),var)
                    case (DIR_Z)
                        dn = g%z(face_idx(3),var) - g%z(interior_idx(3),var)
                    end select
                    if (bc%faceBcType(var,face_id) == 0_C_INT) then
                        f%q(face_idx(1),face_idx(2),face_idx(3),var) = bc%faceBcValue(var,face_id)
                    else
                        f%q(face_idx(1),face_idx(2),face_idx(3),var) = &
                            f%q(interior_idx(1),interior_idx(2),interior_idx(3),var) &
                          + dn*bc%faceBcValue(var,face_id)
                    end if
                else
                    face_idx(dir) = ghost_idx
                    interior_idx(dir) = interior_idx_dir
                    select case (dir)
                    case (DIR_X)
                        dn = g%x(face_idx(1),var) - g%x(interior_idx(1),var)
                    case (DIR_Y)
                        dn = g%y(face_idx(2),var) - g%y(interior_idx(2),var)
                    case (DIR_Z)
                        dn = g%z(face_idx(3),var) - g%z(interior_idx(3),var)
                    end select
                    if (bc%faceBcType(var,face_id) == 0_C_INT) then
                        f%q(face_idx(1),face_idx(2),face_idx(3),var) = &
                            2.0d0*bc%faceBcValue(var,face_id) &
                          - f%q(interior_idx(1),interior_idx(2),interior_idx(3),var)
                    else
                        f%q(face_idx(1),face_idx(2),face_idx(3),var) = &
                            f%q(interior_idx(1),interior_idx(2),interior_idx(3),var) &
                          + dn*bc%faceBcValue(var,face_id)
                    end if
                end if
            end do
        end do
        !$omp end target teams distribute parallel do
    end subroutine apply_bc

end module boundary
