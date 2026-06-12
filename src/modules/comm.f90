module comm
    use, intrinsic :: iso_c_binding
    use :: mpi_f08
    use :: init, only: dns_type, NVAR
    use :: blocks, only: block_set_type, DIST_ZORDER, zorder_owner, zorder_start, zorder_count, &
        leaf_at, level_cells
    use :: boundary, only: boundary_type
#ifdef USE_OPENMP_OFFLOAD
    use omp_lib
#endif
    implicit none

    private

    integer, parameter :: HALO_TAG = 1000

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

        ! Block-pair halo exchange entries, built once by
        ! init_block_exchange. One entry is the transfer of one block
        ! face/edge/corner region from its owner block to a neighbour
        ! block's halo. Same-rank entries are executed as a single flat
        ! device copy kernel; off-rank entries are grouped into one message
        ! per peer rank, with a canonical entry order (receiver's blocks in
        ! slot order, fixed direction order) that both ends derive
        ! independently, so the wire format needs no negotiation.
        integer :: nLocal = 0
        integer :: nLocalPts = 0
        integer, allocatable :: lSrcSlot(:), lDstSlot(:)   ! (nLocal)
        integer, allocatable :: lSrcLo(:,:), lDstLo(:,:), lExt(:,:) ! (3,nLocal)
        integer, allocatable :: lOff(:)                    ! (0:nLocal) point prefix

        integer :: nPeers = 0
        integer, allocatable :: peerRank(:)
        integer, allocatable :: peerSendOff(:)             ! (0:nPeers) point prefix per peer
        integer, allocatable :: peerRecvOff(:)
        integer :: nSend = 0, nRecv = 0
        integer, allocatable :: sSlot(:), sPeer(:)         ! (nSend)
        integer, allocatable :: sLo(:,:), sExt(:,:)        ! (3,nSend)
        integer, allocatable :: sOff(:)                    ! (0:nSend) point prefix, peer-major
        integer, allocatable :: rSlot(:), rPeer(:)
        integer, allocatable :: rLo(:,:), rExt(:,:)
        integer, allocatable :: rOff(:)

        integer :: maxBufferCount = 0
        real(C_DOUBLE), allocatable :: sendbuf(:,:)        ! (maxBufferCount, nPeers)
        real(C_DOUBLE), allocatable :: recvbuf(:,:)

        type(MPI_Request), allocatable :: request(:)       ! (2*nPeers)
        integer(C_INT) :: activeVars(NVAR) = 0_C_INT
        integer :: nActiveVars = 0
    end type comm_type

    public :: comm_init_world, comm_init, comm_finalize
    public :: comm_allreduce_max, comm_allreduce_sum
    public :: init_block_exchange
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
        do dir = 1, 3
            call local_range(int(dns%globalSize(dir)), c%dims(dir), c%coords(dir), &
                             dns%localSize(dir,0), dns%localSize(dir,1))
            dns%localSize(dir,2) = dns%localSize(dir,1) - dns%localSize(dir,0) + 1_C_INT
            local_n(dir) = int(dns%localSize(dir,2))
        end do

        if (any(local_n <= 0)) error stop "MPI decomposition produced an empty local block"

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

    ! Enumerate the block-pair exchange entries (Section 5 of the strategy
    ! document). For each (destination block, 26-direction) pair the source
    ! is the neighbour block on the global block lattice; box extents follow
    ! the per-component rule of the old rank exchange, including the
    ! tangential extension into physical-boundary halos that fills
    ! wall-adjacent corners from apply_bc-set values.
    subroutine init_block_exchange(c, blk, dns)
        type(comm_type), intent(inout) :: c
        type(block_set_type), intent(in) :: blk
        type(dns_type), intent(in) :: dns

        integer :: off(3,26)
        integer :: nb(3), gn(3)
        integer :: b, d, p, pass, owner, slot
        integer :: to(3), srcLo(3), dstLo(3), ext(3)
        integer :: peerCoords(3), peerFirst(3), peerLast(3)
        integer :: peerBlocks, peerStart, pb, dorigin(3), dlevel
        integer :: nLocal, nSend, nRecv, pts, maxCount, ierr
        logical :: haveNeighbor

        call build_direction_table(off)
        nb = int(blk%nb)
        gn = int(dns%globalSize)

        call free_block_exchange(c)
        call collect_peers(c, blk, dns, off)

        ! Pass 1 counts, pass 2 fills; identical loop structure keeps the
        ! enumeration canonical.
        do pass = 1, 2
            nLocal = 0
            c%nLocalPts = 0

            ! Local entries: my blocks' halos served by my own blocks
            ! (including the periodic wrap inside a single rank).
            do b = 1, int(blk%nBlocks)
                do d = 1, 26
                    call neighbor_block(c, blk, dns, b, off(:,d), haveNeighbor, owner, slot, to)
                    if (.not. haveNeighbor) cycle
                    if (owner /= c%cart_rank) cycle
                    call entry_boxes(c, blk, dns, int(blk%level(b)), int(blk%origin(:,b)), &
                        off(:,d), nb, srcLo, dstLo, ext)
                    nLocal = nLocal + 1
                    pts = ext(1)*ext(2)*ext(3)
                    if (pass == 2) then
                        c%lSrcSlot(nLocal) = slot
                        c%lDstSlot(nLocal) = b
                        c%lSrcLo(:,nLocal) = srcLo
                        c%lDstLo(:,nLocal) = dstLo
                        c%lExt(:,nLocal) = ext
                        c%lOff(nLocal) = c%lOff(nLocal-1) + pts
                    end if
                    c%nLocalPts = c%nLocalPts + pts
                end do
            end do
            c%nLocal = nLocal

            ! Recv entries: my blocks in slot order, fixed direction order,
            ! grouped peer-major. Send entries mirror this by enumerating the
            ! PEER's blocks in the peer's slot order, so both ends of each
            ! message agree on the entry sequence without communication.
            nRecv = 0
            nSend = 0
            c%peerSendOff(0) = 0
            c%peerRecvOff(0) = 0
            do p = 1, c%nPeers
                c%peerRecvOff(p) = c%peerRecvOff(p-1)
                do b = 1, int(blk%nBlocks)
                    do d = 1, 26
                        call neighbor_block(c, blk, dns, b, off(:,d), haveNeighbor, owner, slot, to)
                        if (.not. haveNeighbor) cycle
                        if (owner /= c%peerRank(p)) cycle
                        call entry_boxes(c, blk, dns, int(blk%level(b)), int(blk%origin(:,b)), &
                            off(:,d), nb, srcLo, dstLo, ext)
                        nRecv = nRecv + 1
                        pts = ext(1)*ext(2)*ext(3)
                        if (pass == 2) then
                            c%rSlot(nRecv) = b
                            c%rPeer(nRecv) = p
                            c%rLo(:,nRecv) = dstLo
                            c%rExt(:,nRecv) = ext
                            c%rOff(nRecv) = c%rOff(nRecv-1) + pts
                        end if
                        c%peerRecvOff(p) = c%peerRecvOff(p) + pts
                    end do
                end do

                ! Send side of the message me -> peer p: walk the peer's
                ! blocks exactly as the peer walks its own slot order above.
                if (blk%distMode == DIST_ZORDER) then
                    peerStart = int(zorder_start(blk%nBlocksGlobal, int(c%cart_size, C_INT), &
                        int(c%peerRank(p), C_INT)))
                    peerBlocks = int(zorder_count(blk%nBlocksGlobal, int(c%cart_size, C_INT), &
                        int(c%peerRank(p), C_INT)))
                else
                    call MPI_Cart_coords(c%cart_comm, c%peerRank(p), 3, peerCoords, ierr)
                    call rank_box(dns, c%dims, peerCoords, peerFirst, peerLast)
                    peerStart = 0
                    peerBlocks = 1
                end if
                c%peerSendOff(p) = c%peerSendOff(p-1)
                do pb = 1, peerBlocks
                    if (blk%distMode == DIST_ZORDER) then
                        dorigin = int(blk%leafCoord(:,peerStart + pb))*nb
                        dlevel = int(blk%leafLevel(peerStart + pb))
                    else
                        dorigin = peerFirst - 1
                        dlevel = 0
                    end if
                    do d = 1, 26
                        call neighbor_origin(c, dns, dlevel, dorigin, off(:,d), nb, haveNeighbor, to)
                        if (.not. haveNeighbor) cycle
                        if (.not. origin_is_mine(blk, dns, dlevel, to)) cycle
                        call entry_boxes(c, blk, dns, dlevel, dorigin, off(:,d), nb, srcLo, dstLo, ext)
                        nSend = nSend + 1
                        pts = ext(1)*ext(2)*ext(3)
                        if (pass == 2) then
                            c%sSlot(nSend) = my_slot_of(blk, dns, dlevel, to)
                            c%sPeer(nSend) = p
                            c%sLo(:,nSend) = srcLo
                            c%sExt(:,nSend) = ext
                            c%sOff(nSend) = c%sOff(nSend-1) + pts
                        end if
                        c%peerSendOff(p) = c%peerSendOff(p) + pts
                    end do
                end do
            end do
            c%nRecv = nRecv
            c%nSend = nSend

            if (pass == 1) then
                allocate(c%lSrcSlot(max(1,nLocal)), c%lDstSlot(max(1,nLocal)))
                allocate(c%lSrcLo(3,max(1,nLocal)), c%lDstLo(3,max(1,nLocal)), c%lExt(3,max(1,nLocal)))
                allocate(c%lOff(0:max(1,nLocal)))
                allocate(c%sSlot(max(1,nSend)), c%sPeer(max(1,nSend)))
                allocate(c%sLo(3,max(1,nSend)), c%sExt(3,max(1,nSend)))
                allocate(c%sOff(0:max(1,nSend)))
                allocate(c%rSlot(max(1,nRecv)), c%rPeer(max(1,nRecv)))
                allocate(c%rLo(3,max(1,nRecv)), c%rExt(3,max(1,nRecv)))
                allocate(c%rOff(0:max(1,nRecv)))
                c%lOff = 0
                c%sOff = 0
                c%rOff = 0
                c%peerSendOff = 0
                c%peerRecvOff = 0
            end if
        end do

        maxCount = 1
        do p = 1, c%nPeers
            maxCount = max(maxCount, (c%peerSendOff(p) - c%peerSendOff(p-1))*int(NVAR))
            maxCount = max(maxCount, (c%peerRecvOff(p) - c%peerRecvOff(p-1))*int(NVAR))
        end do
        c%maxBufferCount = maxCount
        allocate(c%sendbuf(c%maxBufferCount, max(1, c%nPeers)))
        allocate(c%recvbuf(c%maxBufferCount, max(1, c%nPeers)))
        allocate(c%request(2*max(1, c%nPeers)))
        c%request = MPI_REQUEST_NULL

#ifdef USE_OPENMP_OFFLOAD
        !$omp target enter data map(to: &
        !$omp& c%lSrcSlot, c%lDstSlot, c%lSrcLo, c%lDstLo, c%lExt, c%lOff, &
        !$omp& c%sSlot, c%sPeer, c%sLo, c%sExt, c%sOff, &
        !$omp& c%rSlot, c%rPeer, c%rLo, c%rExt, c%rOff)
        !$omp target enter data map(alloc: c%sendbuf, c%recvbuf)
#endif
    end subroutine init_block_exchange

    subroutine free_block_exchange(c)
        type(comm_type), intent(inout) :: c

        if (allocated(c%sendbuf)) then
#ifdef USE_OPENMP_OFFLOAD
            !$omp target exit data map(delete: c%sendbuf, c%recvbuf)
            !$omp target exit data map(delete: &
            !$omp& c%lSrcSlot, c%lDstSlot, c%lSrcLo, c%lDstLo, c%lExt, c%lOff, &
            !$omp& c%sSlot, c%sPeer, c%sLo, c%sExt, c%sOff, &
            !$omp& c%rSlot, c%rPeer, c%rLo, c%rExt, c%rOff)
#endif
            deallocate(c%sendbuf, c%recvbuf)
            deallocate(c%lSrcSlot, c%lDstSlot, c%lSrcLo, c%lDstLo, c%lExt, c%lOff)
            deallocate(c%sSlot, c%sPeer, c%sLo, c%sExt, c%sOff)
            deallocate(c%rSlot, c%rPeer, c%rLo, c%rExt, c%rOff)
            deallocate(c%request)
        end if
        if (allocated(c%peerRank)) deallocate(c%peerRank)
        if (allocated(c%peerSendOff)) deallocate(c%peerSendOff)
        if (allocated(c%peerRecvOff)) deallocate(c%peerRecvOff)
        c%nLocal = 0
        c%nLocalPts = 0
        c%nPeers = 0
        c%nSend = 0
        c%nRecv = 0
        c%maxBufferCount = 0
    end subroutine free_block_exchange

    subroutine build_direction_table(off)
        integer, intent(out) :: off(3,26)

        integer :: ox, oy, oz, d

        d = 0
        do oz = -1, 1
            do oy = -1, 1
                do ox = -1, 1
                    if (ox == 0 .and. oy == 0 .and. oz == 0) cycle
                    d = d + 1
                    off(:,d) = [ox, oy, oz]
                end do
            end do
        end do
    end subroutine build_direction_table

    ! Neighbour of one of MY blocks: global origin (level-l cells), owner
    ! rank, local slot on the owner.
    subroutine neighbor_block(c, blk, dns, b, off, haveNeighbor, owner, slot, to)
        type(comm_type), intent(in) :: c
        type(block_set_type), intent(in) :: blk
        type(dns_type), intent(in) :: dns
        integer, intent(in) :: b, off(3)
        logical, intent(out) :: haveNeighbor
        integer, intent(out) :: owner, slot, to(3)

        integer :: dorigin(3)

        dorigin = int(blk%origin(:,b))
        call neighbor_origin(c, dns, int(blk%level(b)), dorigin, off, int(blk%nb), haveNeighbor, to)
        owner = -1
        slot = -1
        if (.not. haveNeighbor) return
        call owner_of_origin(c, blk, dns, int(blk%level(b)), to, owner, slot)
        ! Removed (solid-buried) neighbour: a FACE_CLOSED face, no entry.
        if (owner < 0) haveNeighbor = .false.
    end subroutine neighbor_block

    ! Level-l cell origin of the neighbour block in direction off, with
    ! periodic wrap; haveNeighbor is false outside non-periodic boundaries.
    subroutine neighbor_origin(c, dns, level, dorigin, off, nb, haveNeighbor, to)
        type(comm_type), intent(in) :: c
        type(dns_type), intent(in) :: dns
        integer, intent(in) :: level, dorigin(3), off(3), nb(3)
        logical, intent(out) :: haveNeighbor
        integer, intent(out) :: to(3)

        integer :: d, gnl

        haveNeighbor = .true.
        to = 0
        do d = 1, 3
            gnl = int(level_cells(dns, d, int(level, C_INT)))
            to(d) = dorigin(d) + off(d)*nb(d)
            if (to(d) < 0 .or. to(d) >= gnl) then
                if (.not. c%periodic(d)) then
                    haveNeighbor = .false.
                    return
                end if
                to(d) = modulo(to(d), gnl)
            end if
        end do
    end subroutine neighbor_origin

    ! Owner rank and owner-local slot of the leaf at level-l origin `to`.
    subroutine owner_of_origin(c, blk, dns, level, to, owner, slot)
        type(comm_type), intent(in) :: c
        type(block_set_type), intent(in) :: blk
        type(dns_type), intent(in) :: dns
        integer, intent(in) :: level, to(3)
        integer, intent(out) :: owner, slot

        integer :: d, r, id, ierr, sx, sy, sz
        integer(C_INT) :: first, last
        integer :: ownerCoords(3), bcoord(3)

        if (blk%distMode == DIST_ZORDER) then
            bcoord = to/int(blk%nb)
            id = int(leaf_at(blk, level, bcoord))
            if (id < 0) then
                ! A coarser or finer occupant is a 2:1 interface; transfer
                ! operators arrive in Phase 3b/3c.
                if (int(leaf_at(blk, level - 1, bcoord/2)) >= 0) then
                    error stop "2:1 level interfaces are not supported yet (Phase 3b)"
                end if
                do sz = 0, 1
                    do sy = 0, 1
                        do sx = 0, 1
                            if (int(leaf_at(blk, level + 1, 2*bcoord + [sx, sy, sz])) >= 0) then
                                error stop "2:1 level interfaces are not supported yet (Phase 3b)"
                            end if
                        end do
                    end do
                end do
                owner = -1
                slot = -1
                return
            end if
            owner = int(zorder_owner(int(id, C_INT), blk%nBlocksGlobal, int(c%cart_size, C_INT)))
            slot = id - int(zorder_start(blk%nBlocksGlobal, int(c%cart_size, C_INT), &
                int(owner, C_INT))) + 1
            return
        end if

        ! Rank-box mode: one block per rank, owner from the Cartesian ranges.
        do d = 1, 3
            ownerCoords(d) = -1
            do r = 0, c%dims(d) - 1
                call local_range(int(dns%globalSize(d)), c%dims(d), r, first, last)
                if (to(d) >= int(first) - 1 .and. to(d) <= int(last) - 1) then
                    ownerCoords(d) = r
                    exit
                end if
            end do
            if (ownerCoords(d) < 0) error stop "block neighbour outside the rank decomposition"
        end do

        call MPI_Cart_rank(c%cart_comm, ownerCoords, owner, ierr)
        slot = 1
    end subroutine owner_of_origin

    logical function origin_is_mine(blk, dns, level, to)
        type(block_set_type), intent(in) :: blk
        type(dns_type), intent(in) :: dns
        integer, intent(in) :: level, to(3)

        integer :: id

        if (blk%distMode == DIST_ZORDER) then
            id = int(leaf_at(blk, level, to/int(blk%nb)))
            origin_is_mine = id >= int(blk%idStart) .and. &
                             id < int(blk%idStart) + int(blk%nBlocks)
        else
            origin_is_mine = all(to >= int(dns%localSize(1:3,0)) - 1) .and. &
                             all(to <= int(dns%localSize(1:3,1)) - 1)
        end if
    end function origin_is_mine

    integer function my_slot_of(blk, dns, level, to) result(slot)
        type(block_set_type), intent(in) :: blk
        type(dns_type), intent(in) :: dns
        integer, intent(in) :: level, to(3)

        if (blk%distMode == DIST_ZORDER) then
            slot = int(leaf_at(blk, level, to/int(blk%nb))) - int(blk%idStart) + 1
        else
            slot = 1
        end if
        associate(unused => dns)
        end associate
    end function my_slot_of

    ! Source/destination boxes for the entry (dst block B, direction off):
    ! src is in the neighbour block's local frame, dst in B's. The
    ! tangential range in dim d extends into the halo on side s exactly
    ! when the combined neighbour at offset off + s*e_d is absent (outside
    ! a non-periodic boundary or removed): the corner data then cannot
    ! come from an edge/corner entry, and the source block's own halo on
    ! that side is an apply_bc wall halo or a zeroed closed halo, never
    ! exchange-written, so the copy is race-free. With nothing removed
    ! this reduces to the old physical-wall rule.
    subroutine entry_boxes(c, blk, dns, level, dorigin, off, nb, srcLo, dstLo, ext)
        type(comm_type), intent(in) :: c
        type(block_set_type), intent(in) :: blk
        type(dns_type), intent(in) :: dns
        integer, intent(in) :: level, dorigin(3), off(3), nb(3)
        integer, intent(out) :: srcLo(3), dstLo(3), ext(3)

        integer :: d, lo, hi

        do d = 1, 3
            select case (off(d))
            case (1)
                srcLo(d) = 1
                dstLo(d) = nb(d) + 1
                ext(d) = 1
            case (-1)
                srcLo(d) = nb(d)
                dstLo(d) = 0
                ext(d) = 1
            case default
                lo = merge(0, 1, &
                    .not. combined_neighbor_exists(c, blk, dns, level, dorigin, off, nb, d, -1))
                hi = merge(nb(d) + 1, nb(d), &
                    .not. combined_neighbor_exists(c, blk, dns, level, dorigin, off, nb, d, +1))
                srcLo(d) = lo
                dstLo(d) = lo
                ext(d) = hi - lo + 1
            end select
        end do
    end subroutine entry_boxes

    ! Does the block at lattice offset off + side*e_d from `dorigin` exist
    ! (inside the domain or periodic image, and not removed)?
    logical function combined_neighbor_exists(c, blk, dns, level, dorigin, off, nb, d, side) result(exists)
        type(comm_type), intent(in) :: c
        type(block_set_type), intent(in) :: blk
        type(dns_type), intent(in) :: dns
        integer, intent(in) :: level, dorigin(3), off(3), nb(3), d, side

        integer :: offc(3), to(3)

        offc = off
        offc(d) = side
        call neighbor_origin(c, dns, level, dorigin, offc, nb, exists, to)
        if (.not. exists) return
        if (blk%distMode == DIST_ZORDER) then
            exists = leaf_at(blk, level, to/nb) >= 0_C_INT
        end if
    end function combined_neighbor_exists

    subroutine collect_peers(c, blk, dns, off)
        type(comm_type), intent(inout) :: c
        type(block_set_type), intent(in) :: blk
        type(dns_type), intent(in) :: dns
        integer, intent(in) :: off(3,26)

        integer :: b, d, owner, slot, to(3), p
        integer, allocatable :: found(:)
        logical :: haveNeighbor, known

        allocate(found(max(1, c%cart_size)))
        c%nPeers = 0
        do b = 1, int(blk%nBlocks)
            do d = 1, 26
                call neighbor_block(c, blk, dns, b, off(:,d), haveNeighbor, owner, slot, to)
                if (.not. haveNeighbor) cycle
                if (owner == c%cart_rank) cycle
                known = .false.
                do p = 1, c%nPeers
                    if (found(p) == owner) then
                        known = .true.
                        exit
                    end if
                end do
                if (.not. known) then
                    c%nPeers = c%nPeers + 1
                    found(c%nPeers) = owner
                end if
            end do
        end do

        allocate(c%peerRank(max(1, c%nPeers)))
        allocate(c%peerSendOff(0:max(1, c%nPeers)))
        allocate(c%peerRecvOff(0:max(1, c%nPeers)))
        c%peerRank(1:c%nPeers) = found(1:c%nPeers)
        c%peerSendOff = 0
        c%peerRecvOff = 0
        deallocate(found)

        call sort_peers(c)
    end subroutine collect_peers

    subroutine sort_peers(c)
        type(comm_type), intent(inout) :: c

        integer :: i, j, tmp

        do i = 2, c%nPeers
            tmp = c%peerRank(i)
            j = i - 1
            do while (j >= 1)
                if (c%peerRank(j) <= tmp) exit
                c%peerRank(j+1) = c%peerRank(j)
                j = j - 1
            end do
            c%peerRank(j+1) = tmp
        end do
    end subroutine sort_peers

    subroutine comm_finalize(c)
        type(comm_type), intent(inout) :: c

        integer :: ierr

        if (c%exchangeActive) error stop "cannot finalize MPI while a halo exchange is active"

        call free_block_exchange(c)

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

        integer :: ierr, p, nv

        call require_ready(c)
        if (c%exchangeActive) error stop "halo exchange already active"

        call set_active_vars(c, vars)
        nv = c%nActiveVars
        if (nv == 0) return

        c%request = MPI_REQUEST_NULL
        if (c%nPeers > 0) then
            call pack_entries(c, blk)
#ifdef USE_OPENMP_OFFLOAD
            !$omp target data use_device_addr(c%sendbuf, c%recvbuf)
#endif
            do p = 1, c%nPeers
                call MPI_Irecv(c%recvbuf(1,p), (c%peerRecvOff(p) - c%peerRecvOff(p-1))*nv, &
                    MPI_DOUBLE_PRECISION, c%peerRank(p), HALO_TAG, c%cart_comm, c%request(p), ierr)
            end do
            do p = 1, c%nPeers
                call MPI_Isend(c%sendbuf(1,p), (c%peerSendOff(p) - c%peerSendOff(p-1))*nv, &
                    MPI_DOUBLE_PRECISION, c%peerRank(p), HALO_TAG, c%cart_comm, &
                    c%request(c%nPeers+p), ierr)
            end do
#ifdef USE_OPENMP_OFFLOAD
            !$omp end target data
#endif
        end if

        ! Same-rank block-pair copies overlap with the messages in flight.
        call copy_local_entries(c, blk)

        c%exchangeActive = .true.
    end subroutine start_halo_exchange

    subroutine finish_halo_exchange(c, blk)
        type(comm_type), intent(inout) :: c
        type(block_set_type), intent(inout) :: blk

        integer :: ierr, nRequest

        if (.not. c%exchangeActive) return

        if (c%nPeers > 0) then
            nRequest = 2*c%nPeers
            call MPI_Waitall(nRequest, c%request(1:nRequest), MPI_STATUSES_IGNORE, ierr)
            call unpack_entries(c, blk)
        end if

        c%request = MPI_REQUEST_NULL
        c%activeVars = 0_C_INT
        c%nActiveVars = 0
        c%exchangeActive = .false.
    end subroutine finish_halo_exchange

    subroutine exchange_halos(c, blk, vars)
        type(comm_type), intent(inout) :: c
        type(block_set_type), intent(inout) :: blk
        integer(C_INT), intent(in) :: vars(:)

        call start_halo_exchange(c, blk, vars)
        call finish_halo_exchange(c, blk)
    end subroutine exchange_halos

    ! Per-block scalar halos (e.g. les%nut) on the same exchange entries.
    subroutine exchange_scalar_halos(c, scalar)
        type(comm_type), intent(inout) :: c
        real(C_DOUBLE), intent(inout) :: scalar(0:,0:,0:,1:)

        integer :: ierr, p, nRequest

        call require_ready(c)
        if (c%exchangeActive) error stop "halo exchange already active"

        if (c%nPeers > 0) then
            call pack_scalar_entries(c, scalar)
            c%request = MPI_REQUEST_NULL
#ifdef USE_OPENMP_OFFLOAD
            !$omp target data use_device_addr(c%sendbuf, c%recvbuf)
#endif
            do p = 1, c%nPeers
                call MPI_Irecv(c%recvbuf(1,p), c%peerRecvOff(p) - c%peerRecvOff(p-1), &
                    MPI_DOUBLE_PRECISION, c%peerRank(p), HALO_TAG, c%cart_comm, c%request(p), ierr)
            end do
            do p = 1, c%nPeers
                call MPI_Isend(c%sendbuf(1,p), c%peerSendOff(p) - c%peerSendOff(p-1), &
                    MPI_DOUBLE_PRECISION, c%peerRank(p), HALO_TAG, c%cart_comm, &
                    c%request(c%nPeers+p), ierr)
            end do
#ifdef USE_OPENMP_OFFLOAD
            !$omp end target data
#endif
        end if

        call copy_local_scalar_entries(c, scalar)

        if (c%nPeers > 0) then
            nRequest = 2*c%nPeers
            call MPI_Waitall(nRequest, c%request(1:nRequest), MPI_STATUSES_IGNORE, ierr)
            call unpack_scalar_entries(c, scalar)
            c%request = MPI_REQUEST_NULL
        end if
    end subroutine exchange_scalar_halos

    subroutine copy_local_entries(c, blk)
        type(comm_type), intent(inout) :: c
        type(block_set_type), intent(inout) :: blk

        integer :: p, gp, v, e, pt, ni, nj
        integer :: si, sj, sk, di, dj, dk, var, nv, totalItems

        nv = c%nActiveVars
        totalItems = c%nLocalPts*nv
        if (totalItems == 0) return

#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do &
        !$omp& map(to: totalItems, nv, c%nLocal, c%lOff, c%lSrcSlot, c%lDstSlot, &
        !$omp& c%lSrcLo, c%lDstLo, c%lExt, c%activeVars) &
        !$omp& map(tofrom: blk%q) &
        !$omp& private(p,gp,v,e,pt,ni,nj,si,sj,sk,di,dj,dk,var)
#endif
        do p = 1, totalItems
            gp = (p - 1)/nv
            v = p - 1 - gp*nv
            e = find_entry(c%lOff, c%nLocal, gp)
            pt = gp - c%lOff(e-1)
            ni = c%lExt(1,e)
            nj = c%lExt(2,e)
            si = c%lSrcLo(1,e) + modulo(pt, ni)
            sj = c%lSrcLo(2,e) + modulo(pt/ni, nj)
            sk = c%lSrcLo(3,e) + pt/(ni*nj)
            di = c%lDstLo(1,e) + modulo(pt, ni)
            dj = c%lDstLo(2,e) + modulo(pt/ni, nj)
            dk = c%lDstLo(3,e) + pt/(ni*nj)
            var = int(c%activeVars(v+1))
            blk%q(di,dj,dk,var,c%lDstSlot(e)) = blk%q(si,sj,sk,var,c%lSrcSlot(e))
        end do
#ifdef USE_OPENMP_OFFLOAD
        !$omp end target teams distribute parallel do
#endif
    end subroutine copy_local_entries

    subroutine copy_local_scalar_entries(c, scalar)
        type(comm_type), intent(inout) :: c
        real(C_DOUBLE), intent(inout) :: scalar(0:,0:,0:,1:)

        integer :: p, e, pt, ni, nj
        integer :: si, sj, sk, di, dj, dk, totalItems

        totalItems = c%nLocalPts
        if (totalItems == 0) return

#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do &
        !$omp& map(to: totalItems, c%nLocal, c%lOff, c%lSrcSlot, c%lDstSlot, &
        !$omp& c%lSrcLo, c%lDstLo, c%lExt) &
        !$omp& map(tofrom: scalar) &
        !$omp& private(p,e,pt,ni,nj,si,sj,sk,di,dj,dk)
#endif
        do p = 1, totalItems
            e = find_entry(c%lOff, c%nLocal, p - 1)
            pt = p - 1 - c%lOff(e-1)
            ni = c%lExt(1,e)
            nj = c%lExt(2,e)
            si = c%lSrcLo(1,e) + modulo(pt, ni)
            sj = c%lSrcLo(2,e) + modulo(pt/ni, nj)
            sk = c%lSrcLo(3,e) + pt/(ni*nj)
            di = c%lDstLo(1,e) + modulo(pt, ni)
            dj = c%lDstLo(2,e) + modulo(pt/ni, nj)
            dk = c%lDstLo(3,e) + pt/(ni*nj)
            scalar(di,dj,dk,c%lDstSlot(e)) = scalar(si,sj,sk,c%lSrcSlot(e))
        end do
#ifdef USE_OPENMP_OFFLOAD
        !$omp end target teams distribute parallel do
#endif
    end subroutine copy_local_scalar_entries

    subroutine pack_entries(c, blk)
        type(comm_type), intent(inout) :: c
        type(block_set_type), intent(in) :: blk

        integer :: p, gp, v, e, pt, ni, nj
        integer :: i, j, k, var, peer, pos, nv, totalItems

        nv = c%nActiveVars
        totalItems = c%peerSendOff(c%nPeers)*nv
        if (totalItems == 0) return

#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do &
        !$omp& map(to: totalItems, nv, c%nSend, c%sOff, c%sSlot, c%sPeer, &
        !$omp& c%sLo, c%sExt, c%peerSendOff, c%activeVars, blk%q) &
        !$omp& map(tofrom: c%sendbuf) &
        !$omp& private(p,gp,v,e,pt,ni,nj,i,j,k,var,peer,pos)
#endif
        do p = 1, totalItems
            gp = (p - 1)/nv
            v = p - 1 - gp*nv
            e = find_entry(c%sOff, c%nSend, gp)
            pt = gp - c%sOff(e-1)
            ni = c%sExt(1,e)
            nj = c%sExt(2,e)
            i = c%sLo(1,e) + modulo(pt, ni)
            j = c%sLo(2,e) + modulo(pt/ni, nj)
            k = c%sLo(3,e) + pt/(ni*nj)
            var = int(c%activeVars(v+1))
            peer = c%sPeer(e)
            pos = (gp - c%peerSendOff(peer-1))*nv + v + 1
            c%sendbuf(pos,peer) = blk%q(i,j,k,var,c%sSlot(e))
        end do
#ifdef USE_OPENMP_OFFLOAD
        !$omp end target teams distribute parallel do
#endif
    end subroutine pack_entries

    subroutine unpack_entries(c, blk)
        type(comm_type), intent(in) :: c
        type(block_set_type), intent(inout) :: blk

        integer :: p, gp, v, e, pt, ni, nj
        integer :: i, j, k, var, peer, pos, nv, totalItems

        nv = c%nActiveVars
        totalItems = c%peerRecvOff(c%nPeers)*nv
        if (totalItems == 0) return

#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do &
        !$omp& map(to: totalItems, nv, c%nRecv, c%rOff, c%rSlot, c%rPeer, &
        !$omp& c%rLo, c%rExt, c%peerRecvOff, c%activeVars, c%recvbuf) &
        !$omp& map(tofrom: blk%q) &
        !$omp& private(p,gp,v,e,pt,ni,nj,i,j,k,var,peer,pos)
#endif
        do p = 1, totalItems
            gp = (p - 1)/nv
            v = p - 1 - gp*nv
            e = find_entry(c%rOff, c%nRecv, gp)
            pt = gp - c%rOff(e-1)
            ni = c%rExt(1,e)
            nj = c%rExt(2,e)
            i = c%rLo(1,e) + modulo(pt, ni)
            j = c%rLo(2,e) + modulo(pt/ni, nj)
            k = c%rLo(3,e) + pt/(ni*nj)
            var = int(c%activeVars(v+1))
            peer = c%rPeer(e)
            pos = (gp - c%peerRecvOff(peer-1))*nv + v + 1
            blk%q(i,j,k,var,c%rSlot(e)) = c%recvbuf(pos,peer)
        end do
#ifdef USE_OPENMP_OFFLOAD
        !$omp end target teams distribute parallel do
#endif
    end subroutine unpack_entries

    subroutine pack_scalar_entries(c, scalar)
        type(comm_type), intent(inout) :: c
        real(C_DOUBLE), intent(in) :: scalar(0:,0:,0:,1:)

        integer :: p, e, pt, ni, nj
        integer :: i, j, k, peer, pos, totalItems

        totalItems = c%peerSendOff(c%nPeers)
        if (totalItems == 0) return

#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do &
        !$omp& map(to: totalItems, c%nSend, c%sOff, c%sSlot, c%sPeer, &
        !$omp& c%sLo, c%sExt, c%peerSendOff, scalar) &
        !$omp& map(tofrom: c%sendbuf) &
        !$omp& private(p,e,pt,ni,nj,i,j,k,peer,pos)
#endif
        do p = 1, totalItems
            e = find_entry(c%sOff, c%nSend, p - 1)
            pt = p - 1 - c%sOff(e-1)
            ni = c%sExt(1,e)
            nj = c%sExt(2,e)
            i = c%sLo(1,e) + modulo(pt, ni)
            j = c%sLo(2,e) + modulo(pt/ni, nj)
            k = c%sLo(3,e) + pt/(ni*nj)
            peer = c%sPeer(e)
            pos = p - c%peerSendOff(peer-1)
            c%sendbuf(pos,peer) = scalar(i,j,k,c%sSlot(e))
        end do
#ifdef USE_OPENMP_OFFLOAD
        !$omp end target teams distribute parallel do
#endif
    end subroutine pack_scalar_entries

    subroutine unpack_scalar_entries(c, scalar)
        type(comm_type), intent(in) :: c
        real(C_DOUBLE), intent(inout) :: scalar(0:,0:,0:,1:)

        integer :: p, e, pt, ni, nj
        integer :: i, j, k, peer, pos, totalItems

        totalItems = c%peerRecvOff(c%nPeers)
        if (totalItems == 0) return

#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do &
        !$omp& map(to: totalItems, c%nRecv, c%rOff, c%rSlot, c%rPeer, &
        !$omp& c%rLo, c%rExt, c%peerRecvOff, c%recvbuf) &
        !$omp& map(tofrom: scalar) &
        !$omp& private(p,e,pt,ni,nj,i,j,k,peer,pos)
#endif
        do p = 1, totalItems
            e = find_entry(c%rOff, c%nRecv, p - 1)
            pt = p - 1 - c%rOff(e-1)
            ni = c%rExt(1,e)
            nj = c%rExt(2,e)
            i = c%rLo(1,e) + modulo(pt, ni)
            j = c%rLo(2,e) + modulo(pt/ni, nj)
            k = c%rLo(3,e) + pt/(ni*nj)
            peer = c%rPeer(e)
            pos = p - c%peerRecvOff(peer-1)
            scalar(i,j,k,c%rSlot(e)) = c%recvbuf(pos,peer)
        end do
#ifdef USE_OPENMP_OFFLOAD
        !$omp end target teams distribute parallel do
#endif
    end subroutine unpack_scalar_entries

    ! Largest e in 1..n with off(e-1) <= gp (off is a 0:n point prefix).
    pure integer function find_entry(off, n, gp) result(e)
!$omp declare target
        integer, intent(in) :: off(0:), n, gp

        integer :: lo, hi, mid

        lo = 1
        hi = n
        do while (lo < hi)
            mid = (lo + hi)/2
            if (gp >= off(mid)) then
                lo = mid + 1
            else
                hi = mid
            end if
        end do
        e = lo
    end function find_entry

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

    subroutine rank_box(dns, dims, coords, first, last)
        type(dns_type), intent(in) :: dns
        integer, intent(in) :: dims(3), coords(3)
        integer, intent(out) :: first(3), last(3)

        integer :: d
        integer(C_INT) :: f, l

        do d = 1, 3
            call local_range(int(dns%globalSize(d)), dims(d), coords(d), f, l)
            first(d) = int(f)
            last(d) = int(l)
        end do
    end subroutine rank_box

    subroutine local_range(n_global, nproc_dir, coord, first, last)
        integer, intent(in) :: n_global, nproc_dir, coord
        integer(C_INT), intent(out) :: first, last

        first = int((coord*n_global)/nproc_dir + 1, C_INT)
        last = int(((coord+1)*n_global)/nproc_dir, C_INT)
    end subroutine local_range

    subroutine require_ready(c)
        type(comm_type), intent(in) :: c

        if (.not. c%initialized) error stop "comm_type has not been initialized"
        if (.not. allocated(c%sendbuf)) error stop "comm_type exchange entries have not been built"
    end subroutine require_ready

end module comm
