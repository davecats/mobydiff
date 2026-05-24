module comm
    use, intrinsic :: iso_c_binding
    use :: mpi_f08
    use :: init, only: dns_type, field_type, NVAR
    use :: boundary, only: boundary_type
#ifdef USE_OPENMP_OFFLOAD
    use omp_lib
#endif
    implicit none

    private

    ! Diagonal halos are needed by staggered cross-fluxes at processor edges.
    integer, parameter :: MAX_NEIGHBORS = 26

    type, public :: comm_type
        logical :: initialized = .false.
        logical :: exchangeActive = .false.

        type(MPI_Comm) :: cart_comm = MPI_COMM_NULL
        integer :: world_rank = 0
        integer :: world_size = 1
        integer :: cart_rank = 0
        integer :: cart_size = 1
        integer :: local_rank = 0
        logical :: has_terminal = .true.

        integer :: dims(3) = [1, 1, 1]
        integer :: coords(3) = [0, 0, 0]
        logical :: periodic(3) = [.false., .false., .false.]
        logical :: physicalLow(3) = [.false., .false., .false.]
        logical :: physicalHigh(3) = [.false., .false., .false.]

        integer :: nNeighbors = 0
        integer :: neighborRank(MAX_NEIGHBORS) = MPI_PROC_NULL
        integer :: offset(3,MAX_NEIGHBORS) = 0
        integer :: sendLo(3,MAX_NEIGHBORS) = 0
        integer :: sendHi(3,MAX_NEIGHBORS) = 0
        integer :: recvLo(3,MAX_NEIGHBORS) = 0
        integer :: recvHi(3,MAX_NEIGHBORS) = 0
        integer :: nPoints(MAX_NEIGHBORS) = 0

        integer :: maxBufferCount = 0
        real(C_DOUBLE), allocatable :: sendbuf(:,:)
        real(C_DOUBLE), allocatable :: recvbuf(:,:)

        type(MPI_Request) :: request(2*MAX_NEIGHBORS) = MPI_REQUEST_NULL
        integer(C_INT) :: activeVars(NVAR) = 0_C_INT
        integer :: nActiveVars = 0
        integer :: activeCount(MAX_NEIGHBORS) = 0
    end type comm_type

    public :: comm_init_world, comm_init, comm_finalize
    public :: comm_allreduce_max
    public :: start_halo_exchange, finish_halo_exchange, exchange_halos

contains

    subroutine comm_init_world(c)
        type(comm_type), intent(inout) :: c

        integer :: ierr
        logical :: mpi_is_initialized

        call MPI_Initialized(mpi_is_initialized, ierr)
        if (.not. mpi_is_initialized) call MPI_Init(ierr)

        call MPI_Comm_rank(MPI_COMM_WORLD, c%world_rank, ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, c%world_size, ierr)
        c%cart_rank = c%world_rank
        c%cart_size = c%world_size
        c%has_terminal = (c%world_rank == 0)
    end subroutine comm_init_world

    subroutine comm_init(c, dns, bc, dims_in)
        type(comm_type), intent(inout) :: c
        type(dns_type), intent(inout) :: dns
        type(boundary_type), intent(in) :: bc
        integer, intent(in), optional :: dims_in(3)

        type(MPI_Comm) :: local_comm
        integer :: ierr, dir
        integer :: local_n(3)

        call comm_init_world(c)

        c%periodic = bc%isPeriodic
        c%dims = [0, 0, 0]
        if (present(dims_in)) then
            c%dims = dims_in
        end if

        if (any(c%dims < 0)) then
            error stop "MPI Cartesian dimensions must be non-negative"
        end if
        if (any(c%dims == 0)) then
            call MPI_Dims_create(c%world_size, 3, c%dims, ierr)
            if (ierr /= MPI_SUCCESS) error stop "MPI_Dims_create failed"
        end if

        if (product(c%dims) /= c%world_size) then
            error stop "MPI Cartesian dimensions do not match the number of ranks"
        end if

        call MPI_Cart_create(MPI_COMM_WORLD, 3, c%dims, c%periodic, .true., c%cart_comm, ierr)
        call MPI_Comm_rank(c%cart_comm, c%cart_rank, ierr)
        call MPI_Comm_size(c%cart_comm, c%cart_size, ierr)
        call MPI_Cart_coords(c%cart_comm, c%cart_rank, 3, c%coords, ierr)
        if (c%has_terminal) then
            write(*,'(A,1X,I0,1X,I0,1X,I0)') "MPI Cartesian dimensions:", c%dims
        end if
        c%physicalLow = (.not. c%periodic) .and. (c%coords == 0)
        c%physicalHigh = (.not. c%periodic) .and. (c%coords == c%dims - 1)

        do dir = 1, 3
            call local_range(int(dns%globalSize(dir)), c%dims(dir), c%coords(dir), &
                             dns%localSize(dir,0), dns%localSize(dir,1))
            dns%localSize(dir,2) = dns%localSize(dir,1) - dns%localSize(dir,0) + 1_C_INT
            local_n(dir) = int(dns%localSize(dir,2))
        end do

        if (any(local_n <= 0)) error stop "MPI decomposition produced an empty local block"

        call build_neighbors(c, local_n)
        call ensure_buffer_capacity(c, max(1, max_neighbor_points(c)))

        call MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, local_comm, ierr)
        call MPI_Comm_rank(local_comm, c%local_rank, ierr)
        call MPI_Comm_free(local_comm, ierr)

#ifdef USE_OPENMP_OFFLOAD
        if (omp_get_num_devices() > 0) then
            call omp_set_default_device(mod(c%local_rank, omp_get_num_devices()))
        end if
#endif

        c%initialized = .true.
    end subroutine comm_init

    subroutine comm_finalize(c)
        type(comm_type), intent(inout) :: c

        integer :: ierr

        if (c%exchangeActive) error stop "cannot finalize MPI while a halo exchange is active"

        if (allocated(c%sendbuf)) deallocate(c%sendbuf)
        if (allocated(c%recvbuf)) deallocate(c%recvbuf)

        if (c%cart_comm /= MPI_COMM_NULL) then
            call MPI_Comm_free(c%cart_comm, ierr)
        end if

        call MPI_Finalize(ierr)

        c%initialized = .false.
    end subroutine comm_finalize

    subroutine comm_allreduce_max(c, values)
        type(comm_type), intent(in) :: c
        real(C_DOUBLE), intent(inout) :: values(:)

        integer :: ierr

        call MPI_Allreduce(MPI_IN_PLACE, values, size(values), MPI_DOUBLE_PRECISION, MPI_MAX, c%cart_comm, ierr)
    end subroutine comm_allreduce_max

    subroutine start_halo_exchange(c, f, vars)
        type(comm_type), intent(inout) :: c
        type(field_type), intent(inout) :: f
        integer(C_INT), intent(in) :: vars(:)

        integer :: ierr, n, count, recvOffset(3)

        call require_ready(c)
        if (c%exchangeActive) error stop "halo exchange already active"

        call set_active_vars(c, vars)
        if (c%nNeighbors == 0 .or. c%nActiveVars == 0) return

        count = max_neighbor_points(c) * c%nActiveVars
        call ensure_buffer_capacity(c, count)

#ifdef USE_OPENMP_OFFLOAD
        call update_send_boxes_from_device(c, f)
#endif

        c%request = MPI_REQUEST_NULL
        do n = 1, c%nNeighbors
            c%activeCount(n) = c%nPoints(n) * c%nActiveVars
            recvOffset = -c%offset(:,n)
            call MPI_Irecv(c%recvbuf(1,n), c%activeCount(n), MPI_DOUBLE_PRECISION, &
                c%neighborRank(n), halo_tag(recvOffset), c%cart_comm, c%request(n), ierr)
        end do

        do n = 1, c%nNeighbors
            call pack_q_box(f, c%sendLo(:,n), c%sendHi(:,n), &
                            c%activeVars(1:c%nActiveVars), c%sendbuf(:,n))
            call MPI_Isend(c%sendbuf(1,n), c%activeCount(n), MPI_DOUBLE_PRECISION, &
                c%neighborRank(n), halo_tag(c%offset(:,n)), c%cart_comm, &
                c%request(c%nNeighbors+n), ierr)
        end do

        c%exchangeActive = .true.
    end subroutine start_halo_exchange

    subroutine finish_halo_exchange(c, f)
        type(comm_type), intent(inout) :: c
        type(field_type), intent(inout) :: f

        integer :: ierr, n, nRequest

        if (.not. c%exchangeActive) return

        nRequest = 2*c%nNeighbors
        call MPI_Waitall(nRequest, c%request(1:nRequest), MPI_STATUSES_IGNORE, ierr)

        do n = 1, c%nNeighbors
            call unpack_q_box(f, c%recvLo(:,n), c%recvHi(:,n), &
                              c%activeVars(1:c%nActiveVars), c%recvbuf(:,n))
        end do

#ifdef USE_OPENMP_OFFLOAD
        call update_recv_boxes_to_device(c, f)
#endif

        c%request = MPI_REQUEST_NULL
        c%activeCount = 0
        c%activeVars = 0_C_INT
        c%nActiveVars = 0
        c%exchangeActive = .false.
    end subroutine finish_halo_exchange

    subroutine exchange_halos(c, f, vars)
        type(comm_type), intent(inout) :: c
        type(field_type), intent(inout) :: f
        integer(C_INT), intent(in) :: vars(:)

        call start_halo_exchange(c, f, vars)
        call finish_halo_exchange(c, f)
    end subroutine exchange_halos

    subroutine build_neighbors(c, local_n)
        type(comm_type), intent(inout) :: c
        integer, intent(in) :: local_n(3)

        integer :: ox, oy, oz, ierr, n
        integer :: off(3), neighborCoords(3)
        logical :: valid

        c%nNeighbors = 0
        do oz = -1, 1
            do oy = -1, 1
                do ox = -1, 1
                    off = [ox, oy, oz]
                    if (all(off == 0)) cycle

                    call get_neighbor_coords(c, off, neighborCoords, valid)
                    if (.not. valid) cycle

                    n = c%nNeighbors + 1
                    if (n > MAX_NEIGHBORS) error stop "too many MPI halo neighbors"

                    c%offset(:,n) = off
                    call MPI_Cart_rank(c%cart_comm, neighborCoords, c%neighborRank(n), ierr)
                    call set_neighbor_boxes(local_n, off, c%physicalLow, c%physicalHigh, &
                                            c%sendLo(:,n), c%sendHi(:,n), &
                                            c%recvLo(:,n), c%recvHi(:,n))
                    c%nPoints(n) = box_point_count(c%sendLo(:,n), c%sendHi(:,n))
                    c%nNeighbors = n
                end do
            end do
        end do
    end subroutine build_neighbors

    subroutine get_neighbor_coords(c, off, neighborCoords, valid)
        type(comm_type), intent(in) :: c
        integer, intent(in) :: off(3)
        integer, intent(out) :: neighborCoords(3)
        logical, intent(out) :: valid

        integer :: dir, coord

        valid = .true.
        do dir = 1, 3
            coord = c%coords(dir) + off(dir)
            if (coord < 0 .or. coord >= c%dims(dir)) then
                if (.not. c%periodic(dir)) then
                    valid = .false.
                    return
                end if
                coord = modulo(coord, c%dims(dir))
            end if
            neighborCoords(dir) = coord
        end do
    end subroutine get_neighbor_coords

    subroutine set_neighbor_boxes(local_n, off, physicalLow, physicalHigh, sendLo, sendHi, recvLo, recvHi)
        integer, intent(in) :: local_n(3), off(3)
        logical, intent(in) :: physicalLow(3), physicalHigh(3)
        integer, intent(out) :: sendLo(3), sendHi(3), recvLo(3), recvHi(3)

        integer :: dir

        do dir = 1, 3
            select case (off(dir))
            case (-1)
                sendLo(dir) = 1
                sendHi(dir) = 1
                recvLo(dir) = 0
                recvHi(dir) = 0
            case (0)
                sendLo(dir) = merge(0, 1, physicalLow(dir))
                sendHi(dir) = merge(local_n(dir)+1, local_n(dir), physicalHigh(dir))
                recvLo(dir) = sendLo(dir)
                recvHi(dir) = sendHi(dir)
            case (1)
                sendLo(dir) = local_n(dir)
                sendHi(dir) = local_n(dir)
                recvLo(dir) = local_n(dir) + 1
                recvHi(dir) = local_n(dir) + 1
            case default
                error stop "invalid MPI neighbor offset"
            end select
        end do
    end subroutine set_neighbor_boxes

    integer function box_point_count(lo, hi) result(n)
        integer, intent(in) :: lo(3), hi(3)

        n = (hi(1)-lo(1)+1) * (hi(2)-lo(2)+1) * (hi(3)-lo(3)+1)
    end function box_point_count

    integer function max_neighbor_points(c) result(n)
        type(comm_type), intent(in) :: c

        if (c%nNeighbors == 0) then
            n = 1
        else
            n = maxval(c%nPoints(1:c%nNeighbors))
        end if
    end function max_neighbor_points

    integer function halo_tag(off) result(tag)
        integer, intent(in) :: off(3)

        tag = 1000 + (off(1)+1)*9 + (off(2)+1)*3 + (off(3)+1)
    end function halo_tag

    subroutine set_active_vars(c, vars)
        type(comm_type), intent(inout) :: c
        integer(C_INT), intent(in) :: vars(:)

        integer :: n

        if (size(vars) > NVAR) error stop "too many halo variables requested"

        c%activeVars = 0_C_INT
        c%nActiveVars = size(vars)
        do n = 1, c%nActiveVars
            if (vars(n) < 1_C_INT .or. vars(n) > NVAR) error stop "invalid halo variable requested"
            c%activeVars(n) = vars(n)
        end do
    end subroutine set_active_vars

    subroutine pack_q_box(f, lo, hi, vars, buf)
        type(field_type), intent(in) :: f
        integer, intent(in) :: lo(3), hi(3)
        integer(C_INT), intent(in) :: vars(:)
        real(C_DOUBLE), intent(out) :: buf(:)

        integer :: i, j, k, n, nv, var

        n = 0
        do nv = 1, size(vars)
            var = int(vars(nv))
            do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                    do i = lo(1), hi(1)
                        n = n + 1
                        buf(n) = f%q(i,j,k,var)
                    end do
                end do
            end do
        end do
    end subroutine pack_q_box

    subroutine unpack_q_box(f, lo, hi, vars, buf)
        type(field_type), intent(inout) :: f
        integer, intent(in) :: lo(3), hi(3)
        integer(C_INT), intent(in) :: vars(:)
        real(C_DOUBLE), intent(in) :: buf(:)

        integer :: i, j, k, n, nv, var

        n = 0
        do nv = 1, size(vars)
            var = int(vars(nv))
            do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                    do i = lo(1), hi(1)
                        n = n + 1
                        f%q(i,j,k,var) = buf(n)
                    end do
                end do
            end do
        end do
    end subroutine unpack_q_box

#ifdef USE_OPENMP_OFFLOAD
    subroutine update_send_boxes_from_device(c, f)
        type(comm_type), intent(in) :: c
        type(field_type), intent(inout) :: f

        integer :: n, nv

        do n = 1, c%nNeighbors
            do nv = 1, c%nActiveVars
                call update_q_box_from_device(f, c%sendLo(:,n), c%sendHi(:,n), c%activeVars(nv))
            end do
        end do
    end subroutine update_send_boxes_from_device

    subroutine update_recv_boxes_to_device(c, f)
        type(comm_type), intent(in) :: c
        type(field_type), intent(inout) :: f

        integer :: n, nv

        do n = 1, c%nNeighbors
            do nv = 1, c%nActiveVars
                call update_q_box_to_device(f, c%recvLo(:,n), c%recvHi(:,n), c%activeVars(nv))
            end do
        end do
    end subroutine update_recv_boxes_to_device

    subroutine update_q_box_from_device(f, lo, hi, var)
        type(field_type), intent(inout) :: f
        integer, intent(in) :: lo(3), hi(3)
        integer(C_INT), intent(in) :: var

        !$omp target update from(f%q(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),var))
    end subroutine update_q_box_from_device

    subroutine update_q_box_to_device(f, lo, hi, var)
        type(field_type), intent(inout) :: f
        integer, intent(in) :: lo(3), hi(3)
        integer(C_INT), intent(in) :: var

        !$omp target update to(f%q(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),var))
    end subroutine update_q_box_to_device
#endif

    subroutine local_range(n_global, nproc_dir, coord, first, last)
        integer, intent(in) :: n_global, nproc_dir, coord
        integer(C_INT), intent(out) :: first, last

        first = int((coord*n_global)/nproc_dir + 1, C_INT)
        last = int(((coord+1)*n_global)/nproc_dir, C_INT)
    end subroutine local_range

    subroutine ensure_buffer_capacity(c, count)
        type(comm_type), intent(inout) :: c
        integer, intent(in) :: count

        integer :: capacity

        capacity = max(1, count)
        if (allocated(c%sendbuf) .and. capacity <= c%maxBufferCount) return

        if (allocated(c%sendbuf)) deallocate(c%sendbuf)
        if (allocated(c%recvbuf)) deallocate(c%recvbuf)

        c%maxBufferCount = capacity
        allocate(c%sendbuf(c%maxBufferCount, MAX_NEIGHBORS))
        allocate(c%recvbuf(c%maxBufferCount, MAX_NEIGHBORS))
    end subroutine ensure_buffer_capacity

    subroutine require_ready(c)
        type(comm_type), intent(in) :: c

        if (.not. c%initialized) error stop "comm_type has not been initialized"
        if (.not. allocated(c%sendbuf)) error stop "comm_type buffers have not been allocated"
    end subroutine require_ready

end module comm
