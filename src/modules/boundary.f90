module boundary
    use, intrinsic :: iso_c_binding
    use :: init, only: dns_type, VAR_U, VAR_V, VAR_W, VAR_P
    use :: blocks, only: block_set_type, FACE_PHYS
    implicit none

    integer(C_INT), parameter :: DIR_X = 1_C_INT
    integer(C_INT), parameter :: DIR_Y = 2_C_INT
    integer(C_INT), parameter :: DIR_Z = 3_C_INT
    integer(C_INT), parameter :: SIDE_MIN = 0_C_INT
    integer(C_INT), parameter :: SIDE_MAX = 1_C_INT
    integer, parameter :: NFACES = 6

    ! Domain-face patch types ([boundary] <dir>_<side>_patch = wall | patch).
    ! Meaningful on non-periodic faces only. UNSET (no declaration) falls back
    ! to the historical inference in domain_face_is_wall, so existing inis run
    ! unchanged.
    integer(C_INT), parameter :: PATCH_UNSET = 0_C_INT
    integer(C_INT), parameter :: PATCH_GENERIC = 1_C_INT
    integer(C_INT), parameter :: PATCH_WALL = 2_C_INT

    ! Per-face ghost modes for the generic cell-centred scalar BC applicator
    ! (apply_scalar_bc). MIRROR gives a zero face value (Dirichlet 0), COPY a
    ! zero normal gradient, VALUE a prescribed face value; NONE leaves the
    ! ghost untouched.
    integer(C_INT), parameter :: SCALAR_BC_NONE = 0_C_INT
    integer(C_INT), parameter :: SCALAR_BC_COPY = 1_C_INT
    integer(C_INT), parameter :: SCALAR_BC_MIRROR = 2_C_INT
    integer(C_INT), parameter :: SCALAR_BC_VALUE = 3_C_INT

    type :: boundary_type
        logical(C_BOOL) :: isPeriodic(1:3)
        integer(C_INT) :: nTotal = 0_C_INT

        ! Active physical boundary points, block-local with the owning block
        ! slot. Periodic and block/MPI halos are handled by comm.f90.
        integer(C_INT), allocatable :: pointFace(:), slot(:), i(:), j(:), k(:)

        ! Face defaults seed pointwise values; apply_bc uses pointBcValue only.
        integer(C_INT) :: faceBcType(VAR_U:VAR_P,1:NFACES) = 0_C_INT
        real(C_DOUBLE) :: faceBcDefaultValue(VAR_U:VAR_P,1:NFACES) = 0.0d0
        real(C_DOUBLE), allocatable :: pointBcValue(:,:)

        ! Declared patch type per face (index = boundary_face_id).
        integer(C_INT) :: facePatchType(1:NFACES) = PATCH_UNSET
    end type boundary_type

contains

    subroutine init_bc(bc)
        type(boundary_type), intent(inout) :: bc

        call destroy_boundary_faces(bc)
        bc%isPeriodic(1:3) = .true.
        bc%faceBcType = 0_C_INT
        bc%faceBcDefaultValue = 0.0d0
        bc%faceBcType(VAR_P,:) = 1_C_INT
        bc%facePatchType = PATCH_UNSET
    end subroutine init_bc

    ! A domain face is a wall iff its direction is non-periodic and the face
    ! is a wall patch: an explicit [boundary] <dir>_<side>_patch declaration
    ! wins; when absent, the historical inference applies (Dirichlet on both
    ! tangential velocity components = no-slip; Neumann tangential faces --
    ! free slip, symmetry -- carry no wall layer). NOTE the inference reads a
    ! Dirichlet velocity INLET as a wall; declare such a face `patch`.
    logical function domain_face_is_wall(bc, dir, side) result(is_wall)
        type(boundary_type), intent(in) :: bc
        integer, intent(in) :: dir, side

        integer :: face_id, var
        integer, parameter :: DIRICHLET = 0

        is_wall = .false.
        if (bc%isPeriodic(dir)) return

        face_id = boundary_face_id(dir, side)
        select case (bc%facePatchType(face_id))
        case (PATCH_WALL)
            is_wall = .true.
        case (PATCH_GENERIC)
            is_wall = .false.
        case default
            is_wall = .true.
            do var = int(VAR_U), int(VAR_W)
                if (var == dir) cycle   ! the normal component does not decide no-slip
                if (bc%faceBcType(var, face_id) /= DIRICHLET) is_wall = .false.
            end do
        end select
    end function domain_face_is_wall

    ! [boundary] <dir>_<side>_patch is meaningful on non-periodic faces only;
    ! declaring one on a periodic direction is a config error (checked here,
    ! after both the patch keys and the periodic_* flags are final).
    subroutine validate_patch_types(bc)
        type(boundary_type), intent(in) :: bc
        integer :: dir, side

        do dir = 1, 3
            if (.not. bc%isPeriodic(dir)) cycle
            do side = 0, 1
                if (bc%facePatchType(boundary_face_id(dir, side)) /= PATCH_UNSET) then
                    error stop "[boundary] patch type declared on a periodic direction"
                end if
            end do
        end do
    end subroutine validate_patch_types

    subroutine init_boundary_faces(bc, blk)
        type(boundary_type), intent(inout) :: bc
        type(block_set_type), intent(in) :: blk
        integer :: nx, ny, nz
        integer :: b, dir, side, face_id, pos, total

        call validate_patch_types(bc)

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))

        call destroy_boundary_faces(bc)

        total = 0
        do b = 1, int(blk%nBlocks)
            do dir = 1, 3
                do side = 0, 1
                    if (block_face_is_physical(blk, b, dir, side)) then
                        total = total + boundary_face_n(dir, nx, ny, nz)
                    end if
                end do
            end do
        end do

        bc%nTotal = int(total, C_INT)
        if (total <= 0) return

        allocate(bc%pointFace(total), bc%slot(total), bc%i(total), bc%j(total), bc%k(total))
        allocate(bc%pointBcValue(VAR_U:VAR_P,total))

        pos = 0
        do b = 1, int(blk%nBlocks)
            do dir = 1, 3
                do side = 0, 1
                    if (.not. block_face_is_physical(blk, b, dir, side)) cycle

                    face_id = boundary_face_id(dir, side)
                    call append_boundary_face_points(bc, face_id, b, pos, nx, ny, nz)
                end do
            end do
        end do

        call update_boundary_values(bc)
    end subroutine init_boundary_faces

    logical function block_face_is_physical(blk, b, dir, side)
        type(block_set_type), intent(in) :: blk
        integer, intent(in) :: b, dir, side

        ! FACE_PHYS only: FACE_CLOSED faces get zeroed halos, not boundary
        ! conditions.
        if (side == SIDE_MIN) then
            block_face_is_physical = blk%physLow(dir,b) == FACE_PHYS
        else
            block_face_is_physical = blk%physHigh(dir,b) == FACE_PHYS
        end if
    end function block_face_is_physical

    subroutine enter_boundary_data(bc)
        type(boundary_type), intent(inout) :: bc
        integer :: npts

#ifdef USE_OPENMP_OFFLOAD
        npts = int(bc%nTotal)
        !$omp target enter data map(to: bc)
        !$omp target enter data map(to: bc%faceBcType)
        if (npts > 0) then
            !$omp target enter data map(to: bc%pointFace(1:npts), bc%slot(1:npts), &
            !$omp& bc%i(1:npts), bc%j(1:npts), bc%k(1:npts))
            !$omp target enter data map(to: bc%pointBcValue(VAR_U:VAR_P,1:npts))
        end if
#endif
    end subroutine enter_boundary_data

    subroutine exit_boundary_data(bc)
        type(boundary_type), intent(inout) :: bc
        integer :: npts

#ifdef USE_OPENMP_OFFLOAD
        npts = int(bc%nTotal)
        if (npts > 0) then
            !$omp target exit data map(delete: bc%pointBcValue(VAR_U:VAR_P,1:npts))
            !$omp target exit data map(delete: bc%pointFace(1:npts), bc%slot(1:npts), &
            !$omp& bc%i(1:npts), bc%j(1:npts), bc%k(1:npts))
        end if
        !$omp target exit data map(delete: bc%faceBcType)
        !$omp target exit data map(delete: bc)
#endif
    end subroutine exit_boundary_data

    subroutine destroy_boundary_faces(bc)
        type(boundary_type), intent(inout) :: bc

        if (allocated(bc%pointFace)) deallocate(bc%pointFace)
        if (allocated(bc%slot)) deallocate(bc%slot)
        if (allocated(bc%i)) deallocate(bc%i)
        if (allocated(bc%j)) deallocate(bc%j)
        if (allocated(bc%k)) deallocate(bc%k)
        if (allocated(bc%pointBcValue)) deallocate(bc%pointBcValue)
        bc%nTotal = 0_C_INT
    end subroutine destroy_boundary_faces

    subroutine update_boundary_values(bc)
        type(boundary_type), intent(inout) :: bc
        integer :: n, var, face_id, npts

        npts = int(bc%nTotal)
        if (npts <= 0) return

        do n = 1, npts
            face_id = int(bc%pointFace(n))
            do var = VAR_U, VAR_P
                bc%pointBcValue(var,n) = bc%faceBcDefaultValue(var,face_id)
            end do
        end do
    end subroutine update_boundary_values

    subroutine append_boundary_face_points(bc, face_id, slot, pos, nx, ny, nz)
        type(boundary_type), intent(inout) :: bc
        integer, intent(in) :: face_id, slot, nx, ny, nz
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
                    bc%slot(pos) = int(slot, C_INT)
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
                    bc%slot(pos) = int(slot, C_INT)
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
                    bc%slot(pos) = int(slot, C_INT)
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

    subroutine apply_bc(blk, bc)
        type(block_set_type), intent(inout) :: blk
        type(boundary_type), intent(in) :: bc
        integer :: n, npts, b, i, j, k, face_id, var
        integer :: dir, side
        integer :: n_dir, ghost_idx, interior_idx_dir, face_idx_dir, neighbor_idx
        integer(C_INT) :: local_n(1:3)
        integer :: idx(3), interior_idx(3), face_idx(3)
        real(C_DOUBLE) :: dn, bc_value

        npts = int(bc%nTotal)
        if (npts <= 0) return

        local_n = blk%nb(1:3)

        !$omp target teams distribute parallel do &
        !$omp& map(to: npts, local_n(1:3), blk%x, blk%y, blk%z, &
        !$omp& bc%pointFace(1:npts), bc%slot(1:npts), bc%i(1:npts), bc%j(1:npts), bc%k(1:npts), &
        !$omp& bc%faceBcType(VAR_U:VAR_P,1:NFACES), bc%pointBcValue(VAR_U:VAR_P,1:npts)) &
        !$omp& map(tofrom: blk%q) &
        !$omp& private(n,b,i,j,k,face_id,var,dir,side,n_dir, &
        !$omp& ghost_idx,interior_idx_dir,face_idx_dir,neighbor_idx, &
        !$omp& idx,interior_idx,face_idx,dn,bc_value)
        do n = 1, npts
            ! Each entry is one active physical boundary point of one block.
            face_id = int(bc%pointFace(n))
            b = int(bc%slot(n))
            dir = (face_id + 1)/2
            side = modulo(face_id - 1, 2)
            i = int(bc%i(n))
            j = int(bc%j(n))
            k = int(bc%k(n))

            idx = [i, j, k]
            interior_idx = idx
            face_idx = idx

            n_dir = int(local_n(dir))
            ! Pick the boundary/ghost/interior indices along the wall-normal direction.
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
                bc_value = bc%pointBcValue(var,n)
                if (var == dir) then
                    ! Normal velocity lives on the boundary face itself in the staggered layout.
                    face_idx(dir) = face_idx_dir
                    interior_idx(dir) = neighbor_idx
                    select case (dir)
                    case (DIR_X)
                        dn = blk%x(face_idx(1),var,b) - blk%x(interior_idx(1),var,b)
                    case (DIR_Y)
                        dn = blk%y(face_idx(2),var,b) - blk%y(interior_idx(2),var,b)
                    case (DIR_Z)
                        dn = blk%z(face_idx(3),var,b) - blk%z(interior_idx(3),var,b)
                    end select
                    if (bc%faceBcType(var,face_id) == 0_C_INT) then
                        blk%q(face_idx(1),face_idx(2),face_idx(3),var,b) = bc_value
                    else
                        ! Neumann data are stored as the normal derivative.
                        blk%q(face_idx(1),face_idx(2),face_idx(3),var,b) = &
                            blk%q(interior_idx(1),interior_idx(2),interior_idx(3),var,b) &
                          + dn*bc_value
                    end if
                else
                    ! Tangential velocities and pressure use the ghost layer next to the boundary.
                    face_idx(dir) = ghost_idx
                    interior_idx(dir) = interior_idx_dir
                    select case (dir)
                    case (DIR_X)
                        dn = blk%x(face_idx(1),var,b) - blk%x(interior_idx(1),var,b)
                    case (DIR_Y)
                        dn = blk%y(face_idx(2),var,b) - blk%y(interior_idx(2),var,b)
                    case (DIR_Z)
                        dn = blk%z(face_idx(3),var,b) - blk%z(interior_idx(3),var,b)
                    end select
                    if (bc%faceBcType(var,face_id) == 0_C_INT) then
                        ! Dirichlet ghost value chosen so the boundary midpoint has bc_value.
                        blk%q(face_idx(1),face_idx(2),face_idx(3),var,b) = &
                            2.0d0*bc_value &
                          - blk%q(interior_idx(1),interior_idx(2),interior_idx(3),var,b)
                    else
                        blk%q(face_idx(1),face_idx(2),face_idx(3),var,b) = &
                            blk%q(interior_idx(1),interior_idx(2),interior_idx(3),var,b) &
                          + dn*bc_value
                    end if
                end if
            end do
        end do
        !$omp end target teams distribute parallel do
    end subroutine apply_bc

    ! Generic cell-centred scalar ghost update at physical domain faces:
    ! one SCALAR_BC_* mode per face over the bc point lists (mechanics only
    ! -- the caller supplies the per-scalar mode table). The optional value
    ! array feeds the VALUE mode (face value via the ghost-mirror identity),
    ! the hook a scalar inlet needs.
    subroutine apply_scalar_bc(blk, bc, s, mode, value)
        type(block_set_type), intent(in) :: blk
        type(boundary_type), intent(in) :: bc
        real(C_DOUBLE), intent(inout) :: s(0:,0:,0:,1:)
        integer(C_INT), intent(in) :: mode(NFACES)
        real(C_DOUBLE), intent(in), optional :: value(NFACES)

        integer :: n, npts, b, i, j, k, face_id, dir, side, m
        integer :: ghost_idx, interior_idx_dir
        integer :: gi(3), ii(3)
        integer(C_INT) :: local_n(1:3), mode_l(NFACES)
        real(C_DOUBLE) :: value_l(NFACES)

        npts = int(bc%nTotal)
        if (npts <= 0) return
        if (all(mode == SCALAR_BC_NONE)) return
        local_n = blk%nb(1:3)
        mode_l = mode
        value_l = 0.0d0
        if (present(value)) value_l = value

        !$omp target teams distribute parallel do &
        !$omp& map(to: npts, local_n(1:3), mode_l(1:NFACES), value_l(1:NFACES), &
        !$omp& bc%pointFace(1:npts), bc%slot(1:npts), bc%i(1:npts), bc%j(1:npts), bc%k(1:npts)) &
        !$omp& map(tofrom: s) &
        !$omp& private(n,b,i,j,k,face_id,dir,side,m,ghost_idx,interior_idx_dir,gi,ii)
        do n = 1, npts
            face_id = int(bc%pointFace(n))
            m = int(mode_l(face_id))
            if (m == SCALAR_BC_NONE) cycle
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

            select case (m)
            case (SCALAR_BC_COPY)
                s(gi(1),gi(2),gi(3),b) = s(ii(1),ii(2),ii(3),b)
            case (SCALAR_BC_MIRROR)
                s(gi(1),gi(2),gi(3),b) = -s(ii(1),ii(2),ii(3),b)
            case (SCALAR_BC_VALUE)
                ! Ghost chosen so the face midpoint carries the value.
                s(gi(1),gi(2),gi(3),b) = 2.0d0*value_l(face_id) &
                                       - s(ii(1),ii(2),ii(3),b)
            end select
        end do
        !$omp end target teams distribute parallel do
    end subroutine apply_scalar_bc

end module boundary
