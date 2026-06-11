module comm
    use, intrinsic :: iso_c_binding
    use :: mpi_f08
    use :: init, only: dns_type, NVAR
    use :: blocks, only: block_set_type
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

        integer :: dims(3) = [0, 0, 0]
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
        integer :: totalActiveCount = 0
        integer :: bufferOffset(MAX_NEIGHBORS) = 0
        real(C_DOUBLE), allocatable :: sendbuf(:,:)
        real(C_DOUBLE), allocatable :: recvbuf(:,:)

        type(MPI_Request) :: request(2*MAX_NEIGHBORS) = MPI_REQUEST_NULL
        integer(C_INT) :: activeVars(NVAR) = 0_C_INT
        integer :: nActiveVars = 0
        integer :: activeCount(MAX_NEIGHBORS) = 0
    end type comm_type

    public :: comm_init_world, comm_init, comm_finalize
    public :: comm_allreduce_max, comm_allreduce_sum
    public :: start_halo_exchange, finish_halo_exchange, exchange_halos, exchange_scalar_halos

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

    subroutine comm_init(c, dns, bc)
        type(comm_type), intent(inout) :: c
        type(dns_type), intent(inout) :: dns
        type(boundary_type), intent(in) :: bc

        type(MPI_Comm) :: local_comm
        integer :: ierr, dir
        integer :: local_n(3)

        call comm_init_world(c)

        c%periodic = bc%isPeriodic

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

        call MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, local_comm, ierr)
        call MPI_Comm_rank(local_comm, c%local_rank, ierr)
        call MPI_Comm_free(local_comm, ierr)

#ifdef USE_OPENMP_OFFLOAD
        if (omp_get_num_devices() > 0) then
            call omp_set_default_device(mod(c%local_rank, omp_get_num_devices()))
        end if
#endif

        call ensure_buffer_capacity(c, max(1, max_neighbor_points(c)))

        c%initialized = .true.
    end subroutine comm_init

    subroutine comm_finalize(c)
        type(comm_type), intent(inout) :: c

        integer :: ierr

        if (c%exchangeActive) error stop "cannot finalize MPI while a halo exchange is active"

        if (allocated(c%sendbuf)) then
#ifdef USE_OPENMP_OFFLOAD
            !$omp target exit data map(delete: c%sendbuf, c%recvbuf)
#endif
            deallocate(c%sendbuf)
            deallocate(c%recvbuf)
        end if

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

    subroutine comm_allreduce_sum(c, values)
        type(comm_type), intent(in) :: c
        real(C_DOUBLE), intent(inout) :: values(:)

        integer :: ierr

        call MPI_Allreduce(MPI_IN_PLACE, values, size(values), MPI_DOUBLE_PRECISION, MPI_SUM, c%cart_comm, ierr)
    end subroutine comm_allreduce_sum

    subroutine start_halo_exchange(c, blk, vars)
        type(comm_type), intent(inout) :: c
        type(block_set_type), intent(inout) :: blk
        integer(C_INT), intent(in) :: vars(:)

        integer :: ierr, n, count, recvOffset(3)

        call require_ready(c)
        if (c%exchangeActive) error stop "halo exchange already active"

        call set_active_vars(c, vars)
        if (c%nNeighbors == 0 .or. c%nActiveVars == 0) return

        count = max_neighbor_points(c) * c%nActiveVars
        call ensure_buffer_capacity(c, count)

        c%request = MPI_REQUEST_NULL
        c%bufferOffset = 0
        c%totalActiveCount = 0
        do n = 1, c%nNeighbors
            c%activeCount(n) = c%nPoints(n) * c%nActiveVars
            c%bufferOffset(n) = c%totalActiveCount
            c%totalActiveCount = c%totalActiveCount + c%activeCount(n)
        end do

        call pack_q_boxes(c, blk)

#ifdef USE_OPENMP_OFFLOAD
        !$omp target data use_device_addr(c%sendbuf, c%recvbuf)
#endif
        do n = 1, c%nNeighbors
            recvOffset = -c%offset(:,n)
            call MPI_Irecv(c%recvbuf(1,n), c%activeCount(n), MPI_DOUBLE_PRECISION, &
                c%neighborRank(n), halo_tag(recvOffset), c%cart_comm, c%request(n), ierr)
        end do

        do n = 1, c%nNeighbors
            call MPI_Isend(c%sendbuf(1,n), c%activeCount(n), MPI_DOUBLE_PRECISION, &
                c%neighborRank(n), halo_tag(c%offset(:,n)), c%cart_comm, &
                c%request(c%nNeighbors+n), ierr)
        end do
#ifdef USE_OPENMP_OFFLOAD
        !$omp end target data
#endif

        c%exchangeActive = .true.
    end subroutine start_halo_exchange

    subroutine finish_halo_exchange(c, blk)
        type(comm_type), intent(inout) :: c
        type(block_set_type), intent(inout) :: blk

        integer :: ierr, n, nRequest

        if (.not. c%exchangeActive) return

        nRequest = 2*c%nNeighbors
        call MPI_Waitall(nRequest, c%request(1:nRequest), MPI_STATUSES_IGNORE, ierr)

        call unpack_q_boxes(c, blk)

        c%request = MPI_REQUEST_NULL
        c%activeCount = 0
        c%bufferOffset = 0
        c%totalActiveCount = 0
        c%activeVars = 0_C_INT
        c%nActiveVars = 0
        c%exchangeActive = .false.
    end subroutine finish_halo_exchange

    subroutine exchange_halos(c, blk, vars)
        type(comm_type), intent(inout) :: c
        type(block_set_type), intent(inout) :: blk
        integer(C_INT), intent(in) :: vars(:)

        call require_ready(c)
        if (c%exchangeActive) error stop "halo exchange already active"

        call set_active_vars(c, vars)
        if (c%nNeighbors == 0 .or. c%nActiveVars == 0) return

        call start_halo_exchange(c, blk, vars)
        call finish_halo_exchange(c, blk)
    end subroutine exchange_halos

    subroutine exchange_scalar_halos(c, scalar)
        type(comm_type), intent(inout) :: c
        real(C_DOUBLE), intent(inout) :: scalar(0:,0:,0:)

        integer :: ierr, n, nRequest, recvOffset(3)

        call require_ready(c)
        if (c%exchangeActive) error stop "halo exchange already active"
        if (c%nNeighbors == 0) return

        call ensure_buffer_capacity(c, max(1, max_neighbor_points(c)))

        c%request = MPI_REQUEST_NULL
        c%bufferOffset = 0
        c%totalActiveCount = 0
        do n = 1, c%nNeighbors
            c%activeCount(n) = c%nPoints(n)
            c%bufferOffset(n) = c%totalActiveCount
            c%totalActiveCount = c%totalActiveCount + c%activeCount(n)
        end do

        call pack_scalar_boxes(c, scalar)

#ifdef USE_OPENMP_OFFLOAD
        !$omp target data use_device_addr(c%sendbuf, c%recvbuf)
#endif
        do n = 1, c%nNeighbors
            recvOffset = -c%offset(:,n)
            call MPI_Irecv(c%recvbuf(1,n), c%activeCount(n), MPI_DOUBLE_PRECISION, &
                c%neighborRank(n), halo_tag(recvOffset), c%cart_comm, c%request(n), ierr)
        end do

        do n = 1, c%nNeighbors
            call MPI_Isend(c%sendbuf(1,n), c%activeCount(n), MPI_DOUBLE_PRECISION, &
                c%neighborRank(n), halo_tag(c%offset(:,n)), c%cart_comm, &
                c%request(c%nNeighbors+n), ierr)
        end do
#ifdef USE_OPENMP_OFFLOAD
        !$omp end target data
#endif

        c%exchangeActive = .true.
        nRequest = 2*c%nNeighbors
        call MPI_Waitall(nRequest, c%request(1:nRequest), MPI_STATUSES_IGNORE, ierr)
        call unpack_scalar_boxes(c, scalar)

        c%request = MPI_REQUEST_NULL
        call clear_active_counts(c)
        c%exchangeActive = .false.
    end subroutine exchange_scalar_halos

    subroutine prepare_active_counts(c)
        type(comm_type), intent(inout) :: c

        integer :: n

        c%bufferOffset = 0
        c%totalActiveCount = 0
        do n = 1, c%nNeighbors
            c%activeCount(n) = c%nPoints(n) * c%nActiveVars
            c%bufferOffset(n) = c%totalActiveCount
            c%totalActiveCount = c%totalActiveCount + c%activeCount(n)
        end do
    end subroutine prepare_active_counts

    subroutine clear_active_counts(c)
        type(comm_type), intent(inout) :: c

        c%activeCount = 0
        c%bufferOffset = 0
        c%totalActiveCount = 0
        c%activeVars = 0_C_INT
        c%nActiveVars = 0
    end subroutine clear_active_counts

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
                                            c%sendLo(:,n), c%sendHi(:,n), c%recvLo(:,n), c%recvHi(:,n))
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

    ! The send/recv boxes span this rank's box, which in Phase 0 is exactly
    ! block 1. Phase 1 replaces them with per-block-pair exchange entries.
    subroutine pack_q_boxes(c, blk)
        type(comm_type), intent(inout) :: c
        type(block_set_type), intent(in) :: blk

        integer :: pAll, p, q, n, nv, var
        integer :: i, j, k, ni, nj

#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do &
        !$omp& map(to: c%nNeighbors, c%totalActiveCount, c%bufferOffset, c%nPoints, &
        !$omp& c%activeCount, c%sendLo, c%sendHi, c%activeVars, blk%q) &
        !$omp& map(tofrom: c%sendbuf) &
        !$omp& private(pAll,p,q,n,nv,var,i,j,k,ni,nj)
#endif
        do pAll = 1, c%totalActiveCount
            n = 1
            do while (n < c%nNeighbors .and. pAll > c%bufferOffset(n) + c%activeCount(n))
                n = n + 1
            end do

            p = pAll - c%bufferOffset(n)
            ni = c%sendHi(1,n) - c%sendLo(1,n) + 1
            nj = c%sendHi(2,n) - c%sendLo(2,n) + 1
            nv = (p - 1) / c%nPoints(n) + 1
            q = modulo(p - 1, c%nPoints(n))
            i = c%sendLo(1,n) + modulo(q, ni)
            j = c%sendLo(2,n) + modulo(q / ni, nj)
            k = c%sendLo(3,n) + q / (ni*nj)
            var = int(c%activeVars(nv))
            c%sendbuf(p,n) = blk%q(i,j,k,var,1)
        end do
#ifdef USE_OPENMP_OFFLOAD
        !$omp end target teams distribute parallel do
#endif
    end subroutine pack_q_boxes

    subroutine unpack_q_boxes(c, blk)
        type(comm_type), intent(in) :: c
        type(block_set_type), intent(inout) :: blk

        integer :: pAll, p, q, n, nv, var
        integer :: i, j, k, ni, nj

#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do &
        !$omp& map(to: c%nNeighbors, c%totalActiveCount, c%bufferOffset, c%nPoints, &
        !$omp& c%activeCount, c%recvLo, c%recvHi, c%activeVars, c%recvbuf) &
        !$omp& map(tofrom: blk%q) &
        !$omp& private(pAll,p,q,n,nv,var,i,j,k,ni,nj)
#endif
        do pAll = 1, c%totalActiveCount
            n = 1
            do while (n < c%nNeighbors .and. pAll > c%bufferOffset(n) + c%activeCount(n))
                n = n + 1
            end do

            p = pAll - c%bufferOffset(n)
            ni = c%recvHi(1,n) - c%recvLo(1,n) + 1
            nj = c%recvHi(2,n) - c%recvLo(2,n) + 1
            nv = (p - 1) / c%nPoints(n) + 1
            q = modulo(p - 1, c%nPoints(n))
            i = c%recvLo(1,n) + modulo(q, ni)
            j = c%recvLo(2,n) + modulo(q / ni, nj)
            k = c%recvLo(3,n) + q / (ni*nj)
            var = int(c%activeVars(nv))
            blk%q(i,j,k,var,1) = c%recvbuf(p,n)
        end do
#ifdef USE_OPENMP_OFFLOAD
        !$omp end target teams distribute parallel do
#endif
    end subroutine unpack_q_boxes

    subroutine pack_scalar_boxes(c, scalar)
        type(comm_type), intent(inout) :: c
        real(C_DOUBLE), intent(in) :: scalar(0:,0:,0:)

        integer :: pAll, p, q, n
        integer :: i, j, k, ni, nj

#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do &
        !$omp& map(to: c%nNeighbors, c%totalActiveCount, c%bufferOffset, c%nPoints, &
        !$omp& c%activeCount, c%sendLo, c%sendHi, scalar) &
        !$omp& map(tofrom: c%sendbuf) &
        !$omp& private(pAll,p,q,n,i,j,k,ni,nj)
#endif
        do pAll = 1, c%totalActiveCount
            n = 1
            do while (n < c%nNeighbors .and. pAll > c%bufferOffset(n) + c%activeCount(n))
                n = n + 1
            end do

            p = pAll - c%bufferOffset(n)
            ni = c%sendHi(1,n) - c%sendLo(1,n) + 1
            nj = c%sendHi(2,n) - c%sendLo(2,n) + 1
            q = p - 1
            i = c%sendLo(1,n) + modulo(q, ni)
            j = c%sendLo(2,n) + modulo(q / ni, nj)
            k = c%sendLo(3,n) + q / (ni*nj)
            c%sendbuf(p,n) = scalar(i,j,k)
        end do
#ifdef USE_OPENMP_OFFLOAD
        !$omp end target teams distribute parallel do
#endif
    end subroutine pack_scalar_boxes

    subroutine unpack_scalar_boxes(c, scalar)
        type(comm_type), intent(in) :: c
        real(C_DOUBLE), intent(inout) :: scalar(0:,0:,0:)

        integer :: pAll, p, q, n
        integer :: i, j, k, ni, nj

#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do &
        !$omp& map(to: c%nNeighbors, c%totalActiveCount, c%bufferOffset, c%nPoints, &
        !$omp& c%activeCount, c%recvLo, c%recvHi, c%recvbuf) &
        !$omp& map(tofrom: scalar) &
        !$omp& private(pAll,p,q,n,i,j,k,ni,nj)
#endif
        do pAll = 1, c%totalActiveCount
            n = 1
            do while (n < c%nNeighbors .and. pAll > c%bufferOffset(n) + c%activeCount(n))
                n = n + 1
            end do

            p = pAll - c%bufferOffset(n)
            ni = c%recvHi(1,n) - c%recvLo(1,n) + 1
            nj = c%recvHi(2,n) - c%recvLo(2,n) + 1
            q = p - 1
            i = c%recvLo(1,n) + modulo(q, ni)
            j = c%recvLo(2,n) + modulo(q / ni, nj)
            k = c%recvLo(3,n) + q / (ni*nj)
            scalar(i,j,k) = c%recvbuf(p,n)
        end do
#ifdef USE_OPENMP_OFFLOAD
        !$omp end target teams distribute parallel do
#endif
    end subroutine unpack_scalar_boxes

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

        if (allocated(c%sendbuf)) then
#ifdef USE_OPENMP_OFFLOAD
            !$omp target exit data map(delete: c%sendbuf, c%recvbuf)
#endif
            deallocate(c%sendbuf)
            deallocate(c%recvbuf)
        end if

        c%maxBufferCount = capacity
        allocate(c%sendbuf(c%maxBufferCount, MAX_NEIGHBORS))
        allocate(c%recvbuf(c%maxBufferCount, MAX_NEIGHBORS))
#ifdef USE_OPENMP_OFFLOAD
        !$omp target enter data map(alloc: c%sendbuf, c%recvbuf)
#endif
    end subroutine ensure_buffer_capacity

    subroutine require_ready(c)
        type(comm_type), intent(in) :: c

        if (.not. c%initialized) error stop "comm_type has not been initialized"
        if (.not. allocated(c%sendbuf)) error stop "comm_type buffers have not been allocated"
    end subroutine require_ready

end module comm
