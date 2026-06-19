module comm
    use, intrinsic :: iso_c_binding
    use :: mpi_f08
    use :: init, only: dns_type, NVAR, VAR_P
    use :: blocks, only: block_set_type, DIST_ZORDER, zorder_owner, zorder_start, zorder_count, &
        leaf_at, level_cells
    use :: boundary, only: boundary_type
#ifdef USE_OPENMP_OFFLOAD
    use omp_lib
#endif
    implicit none

    private

    integer, parameter :: HALO_TAG = 1000

    ! Exchange entry operations. RESTRICT averages the fine samples
    ! covering each coarse destination point (8 cell-centred, 4 with one
    ! face-staggered dimension, ...); PROLONG injects the covering coarse
    ! value. Sampling happens on the SOURCE side (pack/local copy), so the
    ! wire always carries destination-point values.
    integer, parameter :: OP_COPY = 0
    integer, parameter :: OP_RESTRICT = 1
    integer, parameter :: OP_PROLONG = 2

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
        ! block's halo, executed as a single weighted gather: per
        ! dimension, the source rows for destination index i are
        ! base..base+cnt-1 with base = ishft(ga*i + gb, -gs), a form that
        ! covers same-level copies, restrictions and prolongations alike
        ! (entry_gather_map). Entries are ordered same-level COPY first
        ! (per peer), with prefix counts, so a copy-only exchange is a
        ! prefix of the full one. Same-rank entries run as one flat device
        ! kernel; off-rank entries form one message per peer rank, in a
        ! canonical order (receiver's blocks in slot order, fixed
        ! direction order, copies first) both ends derive independently,
        ! so the wire format needs no negotiation.
        integer :: nLocal = 0
        integer :: nLocalPts = 0
        integer :: nLocalCopyPts = 0
        integer, allocatable :: lSrcSlot(:), lDstSlot(:)   ! (nLocal)
        integer, allocatable :: lDstLo(:,:), lExt(:,:)     ! (3,nLocal)
        ! Per-dim affine gather map (3,nLocal): destination index d reads source
        ! (GA*d + GB) >> GS, averaging GC source cells (GC>1 only for RESTRICT).
        integer, allocatable :: lGA(:,:), lGB(:,:), lGS(:,:), lGC(:,:)
        integer, allocatable :: lLin(:,:)                  ! (3,nLocal) 1 = PROLONG tangential dim (linear when linProlong)
        integer, allocatable :: lFaceNrm(:)                ! 2:1 interface owned-face normal dir (0 = none)
        integer, allocatable :: lOff(:)                    ! (0:nLocal) point prefix
        integer, allocatable :: lEntryOf(:)                ! (0:nLocalPts-1) owning entry per point

        integer :: nPeers = 0
        integer, allocatable :: peerRank(:)
        integer, allocatable :: peerSendOff(:)             ! (0:nPeers) point prefix per peer
        integer, allocatable :: peerRecvOff(:)
        integer, allocatable :: peerSendCopyOff(:)         ! (0:nPeers) same-level copy prefix
        integer, allocatable :: peerRecvCopyOff(:)
        integer :: nSend = 0, nRecv = 0
        integer, allocatable :: sSlot(:), sPeer(:)         ! (nSend)
        integer, allocatable :: sExt(:,:)                  ! (3,nSend)
        integer, allocatable :: sGA(:,:), sGB(:,:), sGS(:,:), sGC(:,:)
        integer, allocatable :: sLin(:,:)                  ! (3,nSend) 1 = PROLONG tangential dim
        integer, allocatable :: sDstLo(:,:)                ! dst box lo (gather indexing)
        integer, allocatable :: sFaceNrm(:)                ! 2:1 interface owned-face normal dir (0 = none)
        integer, allocatable :: sOff(:)                    ! (0:nSend) point prefix, peer-major
        integer, allocatable :: sEntryOf(:)                ! (0:nSendPts-1) owning send entry per point
        integer, allocatable :: rSlot(:), rPeer(:)
        integer, allocatable :: rLo(:,:), rExt(:,:)
        integer, allocatable :: rFaceNrm(:)                ! 2:1 interface owned-face normal dir (0 = none)
        integer, allocatable :: rOff(:)
        integer, allocatable :: rEntryOf(:)                ! (0:nRecvPts-1) owning recv entry per point

        integer :: maxBufferCount = 0
        real(C_DOUBLE), allocatable :: sendbuf(:,:)        ! (maxBufferCount, nPeers)
        real(C_DOUBLE), allocatable :: recvbuf(:,:)

        type(MPI_Request), allocatable :: request(:)       ! (2*nPeers)
        integer(C_INT) :: activeVars(NVAR) = 0_C_INT
        integer :: nActiveVars = 0
        ! When true, the exchange touches only the same-level COPY prefix
        ! of the entry lists, leaving interface ghosts and face copies
        ! frozen (the messages shrink to the per-peer copy prefixes).
        ! Used by the per-colour exchanges inside the pressure projection.
        logical :: copyOnly = .false.
        ! When true, VAR_P is gathered only on the cross-level (interp) entries,
        ! not on the same-level COPY prefix. The intermediate composite-projection
        ! exchanges use this: the sweep reads a neighbour pressure only at 2:1
        ! interfaces, so a same-level halo pressure is never read between sweeps.
        logical :: pSkipCopy = .false.
        ! When true, the coarse->fine (PROLONG) velocity halos are filled by
        ! linear interpolation instead of piecewise-constant injection. Used
        ! ONLY for the projection's final exchange (which feeds the next
        ! substage's momentum predictor): injection inside the relaxation keeps
        ! it contracting, while the momentum sees an O(h^2)-accurate halo. E3.
        logical :: linProlong = .false.
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
        integer :: b, d, p, pass, cand, ncand
        integer :: owner(4), slot(4), opc(4), tqc(3,4)
        integer :: srcLo(3), dstLo(3), ext(3)
        integer :: peerCoords(3), peerFirst(3), peerLast(3)
        integer :: peerBlocks, peerStart, pb, dorigin(3), dlevel
        integer :: nLocal, nSend, nRecv, pts, maxCount, ierr, round, ent, gpt

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
            ! (including the periodic wrap inside a single rank). Round 1
            ! emits the same-level copies, round 2 the cross-level
            ! entries, so the copy-only view is a prefix.
            do round = 1, 2
                do b = 1, int(blk%nBlocks)
                    do d = 1, 26
                        call resolve_neighbors(c, blk, dns, int(blk%level(b)), int(blk%origin(:,b)), &
                            off(:,d), ncand, owner, slot, opc, tqc)
                        do cand = 1, ncand
                            if (owner(cand) /= c%cart_rank) cycle
                            if ((opc(cand) == OP_COPY) .neqv. (round == 1)) cycle
                            call candidate_boxes(c, blk, dns, int(blk%level(b)), int(blk%origin(:,b)), &
                                off(:,d), nb, opc(cand), tqc(:,cand), srcLo, dstLo, ext)
                            nLocal = nLocal + 1
                            pts = ext(1)*ext(2)*ext(3)
                            if (pass == 2) then
                                c%lSrcSlot(nLocal) = slot(cand)
                                c%lDstSlot(nLocal) = b
                                c%lDstLo(:,nLocal) = dstLo
                                c%lExt(:,nLocal) = ext
                                c%lFaceNrm(nLocal) = face_normal(opc(cand), off(:,d))
                                call entry_gather_map(opc(cand), off(:,d), tqc(:,cand), nb, &
                                    srcLo, dstLo, c%lGA(:,nLocal), c%lGB(:,nLocal), &
                                    c%lGS(:,nLocal), c%lGC(:,nLocal), c%lLin(:,nLocal))
                                c%lOff(nLocal) = c%lOff(nLocal-1) + pts
                            end if
                            c%nLocalPts = c%nLocalPts + pts
                        end do
                    end do
                end do
                if (round == 1) c%nLocalCopyPts = c%nLocalPts
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
            c%peerSendCopyOff(0) = 0
            c%peerRecvCopyOff(0) = 0
            do p = 1, c%nPeers
                c%peerRecvOff(p) = c%peerRecvOff(p-1)
                do round = 1, 2
                    do b = 1, int(blk%nBlocks)
                        do d = 1, 26
                            call resolve_neighbors(c, blk, dns, int(blk%level(b)), int(blk%origin(:,b)), &
                                off(:,d), ncand, owner, slot, opc, tqc)
                            do cand = 1, ncand
                                if (owner(cand) /= c%peerRank(p)) cycle
                                if ((opc(cand) == OP_COPY) .neqv. (round == 1)) cycle
                                call candidate_boxes(c, blk, dns, int(blk%level(b)), int(blk%origin(:,b)), &
                                    off(:,d), nb, opc(cand), tqc(:,cand), srcLo, dstLo, ext)
                                nRecv = nRecv + 1
                                pts = ext(1)*ext(2)*ext(3)
                                if (pass == 2) then
                                    c%rSlot(nRecv) = b
                                    c%rPeer(nRecv) = p
                                    c%rLo(:,nRecv) = dstLo
                                    c%rExt(:,nRecv) = ext
                                    c%rFaceNrm(nRecv) = face_normal(opc(cand), off(:,d))
                                    c%rOff(nRecv) = c%rOff(nRecv-1) + pts
                                end if
                                c%peerRecvOff(p) = c%peerRecvOff(p) + pts
                            end do
                        end do
                    end do
                    if (round == 1) c%peerRecvCopyOff(p) = c%peerRecvCopyOff(p-1) &
                        + (c%peerRecvOff(p) - c%peerRecvOff(p-1))
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
                do round = 1, 2
                    do pb = 1, peerBlocks
                        if (blk%distMode == DIST_ZORDER) then
                            dorigin = int(blk%leafCoord(:,peerStart + pb))*nb
                            dlevel = int(blk%leafLevel(peerStart + pb))
                        else
                            dorigin = peerFirst - 1
                            dlevel = 0
                        end if
                        do d = 1, 26
                            call resolve_neighbors(c, blk, dns, dlevel, dorigin, &
                                off(:,d), ncand, owner, slot, opc, tqc)
                            do cand = 1, ncand
                                if (owner(cand) /= c%cart_rank) cycle
                                if ((opc(cand) == OP_COPY) .neqv. (round == 1)) cycle
                                call candidate_boxes(c, blk, dns, dlevel, dorigin, &
                                    off(:,d), nb, opc(cand), tqc(:,cand), srcLo, dstLo, ext)
                                nSend = nSend + 1
                                pts = ext(1)*ext(2)*ext(3)
                                if (pass == 2) then
                                    c%sSlot(nSend) = slot(cand)
                                    c%sPeer(nSend) = p
                                    c%sDstLo(:,nSend) = dstLo
                                    c%sExt(:,nSend) = ext
                                    c%sFaceNrm(nSend) = face_normal(opc(cand), off(:,d))
                                    call entry_gather_map(opc(cand), off(:,d), tqc(:,cand), nb, &
                                        srcLo, dstLo, c%sGA(:,nSend), c%sGB(:,nSend), &
                                        c%sGS(:,nSend), c%sGC(:,nSend), c%sLin(:,nSend))
                                    c%sOff(nSend) = c%sOff(nSend-1) + pts
                                end if
                                c%peerSendOff(p) = c%peerSendOff(p) + pts
                            end do
                        end do
                    end do
                    if (round == 1) c%peerSendCopyOff(p) = c%peerSendCopyOff(p-1) &
                        + (c%peerSendOff(p) - c%peerSendOff(p-1))
                end do
            end do
            c%nRecv = nRecv
            c%nSend = nSend

            if (pass == 1) then
                allocate(c%lSrcSlot(max(1,nLocal)), c%lDstSlot(max(1,nLocal)))
                allocate(c%lDstLo(3,max(1,nLocal)), c%lExt(3,max(1,nLocal)))
                allocate(c%lGA(3,max(1,nLocal)), c%lGB(3,max(1,nLocal)))
                allocate(c%lGS(3,max(1,nLocal)), c%lGC(3,max(1,nLocal)))
                allocate(c%lLin(3,max(1,nLocal)))
                allocate(c%lFaceNrm(max(1,nLocal)))
                allocate(c%lOff(0:max(1,nLocal)))
                allocate(c%sSlot(max(1,nSend)), c%sPeer(max(1,nSend)))
                allocate(c%sExt(3,max(1,nSend)))
                allocate(c%sGA(3,max(1,nSend)), c%sGB(3,max(1,nSend)))
                allocate(c%sGS(3,max(1,nSend)), c%sGC(3,max(1,nSend)))
                allocate(c%sLin(3,max(1,nSend)))
                allocate(c%sDstLo(3,max(1,nSend)))
                allocate(c%sFaceNrm(max(1,nSend)))
                allocate(c%sOff(0:max(1,nSend)))
                allocate(c%rSlot(max(1,nRecv)), c%rPeer(max(1,nRecv)))
                allocate(c%rLo(3,max(1,nRecv)), c%rExt(3,max(1,nRecv)))
                allocate(c%rFaceNrm(max(1,nRecv)))
                allocate(c%rOff(0:max(1,nRecv)))
                c%lOff = 0
                c%sOff = 0
                c%rOff = 0
                c%peerSendOff = 0
                c%peerRecvOff = 0
                c%peerSendCopyOff = 0
                c%peerRecvCopyOff = 0
            end if
        end do

        ! Per-point -> owning-entry lookups: the copy/pack kernels read these
        ! directly instead of binary-searching the *Off prefixes per point.
        allocate(c%lEntryOf(0:max(0, c%nLocalPts-1)))
        do ent = 1, c%nLocal
            do gpt = c%lOff(ent-1), c%lOff(ent)-1
                c%lEntryOf(gpt) = ent
            end do
        end do
        allocate(c%sEntryOf(0:max(0, c%sOff(c%nSend)-1)))
        do ent = 1, c%nSend
            do gpt = c%sOff(ent-1), c%sOff(ent)-1
                c%sEntryOf(gpt) = ent
            end do
        end do
        allocate(c%rEntryOf(0:max(0, c%rOff(c%nRecv)-1)))
        do ent = 1, c%nRecv
            do gpt = c%rOff(ent-1), c%rOff(ent)-1
                c%rEntryOf(gpt) = ent
            end do
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
        !$omp& c%lSrcSlot, c%lDstSlot, c%lDstLo, c%lExt, c%lOff, c%lEntryOf, &
        !$omp& c%lGA, c%lGB, c%lGS, c%lGC, c%lLin, c%lFaceNrm, &
        !$omp& c%sSlot, c%sPeer, c%sExt, c%sOff, c%sEntryOf, &
        !$omp& c%sGA, c%sGB, c%sGS, c%sGC, c%sLin, c%sDstLo, c%sFaceNrm, &
        !$omp& c%peerSendOff, c%peerRecvOff, c%peerSendCopyOff, c%peerRecvCopyOff, &
        !$omp& c%rSlot, c%rPeer, c%rLo, c%rExt, c%rFaceNrm, c%rOff, c%rEntryOf)
        !$omp target enter data map(alloc: c%sendbuf, c%recvbuf)
#endif
    end subroutine init_block_exchange

    subroutine free_block_exchange(c)
        type(comm_type), intent(inout) :: c

        if (allocated(c%sendbuf)) then
#ifdef USE_OPENMP_OFFLOAD
            !$omp target exit data map(delete: c%sendbuf, c%recvbuf)
            !$omp target exit data map(delete: &
            !$omp& c%lSrcSlot, c%lDstSlot, c%lDstLo, c%lExt, c%lOff, c%lEntryOf, &
            !$omp& c%lGA, c%lGB, c%lGS, c%lGC, c%lLin, c%lFaceNrm, &
            !$omp& c%sSlot, c%sPeer, c%sExt, c%sOff, c%sEntryOf, &
            !$omp& c%sGA, c%sGB, c%sGS, c%sGC, c%sLin, c%sDstLo, c%sFaceNrm, &
            !$omp& c%peerSendOff, c%peerRecvOff, c%peerSendCopyOff, c%peerRecvCopyOff, &
            !$omp& c%rSlot, c%rPeer, c%rLo, c%rExt, c%rFaceNrm, c%rOff, c%rEntryOf)
#endif
            deallocate(c%sendbuf, c%recvbuf)
            deallocate(c%lSrcSlot, c%lDstSlot, c%lDstLo, c%lExt, c%lOff, c%lEntryOf)
            deallocate(c%lGA, c%lGB, c%lGS, c%lGC, c%lLin, c%lFaceNrm)
            deallocate(c%sSlot, c%sPeer, c%sExt, c%sOff, c%sEntryOf)
            deallocate(c%sGA, c%sGB, c%sGS, c%sGC, c%sLin, c%sDstLo, c%sFaceNrm)
            deallocate(c%rSlot, c%rPeer, c%rLo, c%rExt, c%rFaceNrm, c%rOff, c%rEntryOf)
            deallocate(c%request)
        end if
        if (allocated(c%peerRank)) deallocate(c%peerRank)
        if (allocated(c%peerSendOff)) deallocate(c%peerSendOff)
        if (allocated(c%peerRecvOff)) deallocate(c%peerRecvOff)
        if (allocated(c%peerSendCopyOff)) deallocate(c%peerSendCopyOff)
        if (allocated(c%peerRecvCopyOff)) deallocate(c%peerRecvCopyOff)
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

    ! Neighbours of the level-`level` block at lattice cell origin
    ! `dorigin` across direction off: one same-level or coarser leaf
    ! (COPY/PROLONG), or up to 2^(tangential dims) finer leaves
    ! (RESTRICT), in fixed child order so the enumeration is canonical on
    ! both ends. tq holds the fine quarter (RESTRICT: child parity;
    ! PROLONG: the destination block's own parity).
    subroutine resolve_neighbors(c, blk, dns, level, dorigin, off, n, owner, slot, op, tq)
        type(comm_type), intent(in) :: c
        type(block_set_type), intent(in) :: blk
        type(dns_type), intent(in) :: dns
        integer, intent(in) :: level, dorigin(3), off(3)
        integer, intent(out) :: n, owner(4), slot(4), op(4), tq(3,4)

        integer :: to(3), cl(3), cc(3), sub(3), id, d
        integer :: sx, sy, sz
        logical :: haveNeighbor

        n = 0
        call neighbor_origin(c, dns, level, dorigin, off, int(blk%nb), haveNeighbor, to)
        if (.not. haveNeighbor) return

        if (blk%distMode /= DIST_ZORDER) then
            call owner_of_origin(c, blk, dns, level, to, owner(1), slot(1))
            if (owner(1) >= 0) then
                n = 1
                op(1) = OP_COPY
                tq(:,1) = 0
            end if
            return
        end if

        cl = to/int(blk%nb)
        id = int(leaf_at(blk, level, cl))
        if (id >= 0) then
            n = 1
            op(1) = OP_COPY
            tq(:,1) = 0
            call id_owner_slot(c, blk, id, owner(1), slot(1))
            return
        end if

        ! Coarser occupant: this block is the fine side of a 2:1 interface.
        id = int(leaf_at(blk, level - 1, cl/2))
        if (id >= 0) then
            n = 1
            op(1) = OP_PROLONG
            tq(:,1) = modulo(dorigin/int(blk%nb), 2)
            call id_owner_slot(c, blk, id, owner(1), slot(1))
            return
        end if

        ! Finer occupants: the children adjacent to this block across off.
        do sz = 0, 1
            do sy = 0, 1
                do sx = 0, 1
                    sub = [sx, sy, sz]
                    do d = 1, 3
                        if (off(d) == 1 .and. sub(d) /= 0) sub(d) = -9
                        if (off(d) == -1 .and. sub(d) /= 1) sub(d) = -9
                    end do
                    if (any(sub == -9)) cycle
                    cc = 2*cl + sub
                    id = int(leaf_at(blk, level + 1, cc))
                    if (id < 0) cycle
                    n = n + 1
                    op(n) = OP_RESTRICT
                    tq(:,n) = sub
                    call id_owner_slot(c, blk, id, owner(n), slot(n))
                end do
            end do
        end do
    end subroutine resolve_neighbors

    subroutine id_owner_slot(c, blk, id, owner, slot)
        type(comm_type), intent(in) :: c
        type(block_set_type), intent(in) :: blk
        integer, intent(in) :: id
        integer, intent(out) :: owner, slot

        owner = int(zorder_owner(int(id, C_INT), blk%nBlocksGlobal, int(c%cart_size, C_INT)))
        slot = id - int(zorder_start(blk%nBlocksGlobal, int(c%cart_size, C_INT), &
            int(owner, C_INT))) + 1
    end subroutine id_owner_slot

    subroutine candidate_boxes(c, blk, dns, level, dorigin, off, nb, op, tq, srcLo, dstLo, ext)
        type(comm_type), intent(in) :: c
        type(block_set_type), intent(in) :: blk
        type(dns_type), intent(in) :: dns
        integer, intent(in) :: level, dorigin(3), off(3), nb(3), op, tq(3)
        integer, intent(out) :: srcLo(3), dstLo(3), ext(3)

        integer :: d

        if (op == OP_COPY) then
            call entry_boxes(c, blk, dns, level, dorigin, off, nb, srcLo, dstLo, ext)
        else
            call interface_boxes(op, off, nb, srcLo, dstLo, ext)
            if (op == OP_RESTRICT) then
                do d = 1, 3
                    if (off(d) == 0) dstLo(d) = tq(d)*nb(d)/2 + 1
                end do
            end if
        end if
    end subroutine candidate_boxes

    ! Normal dimension of a 2:1 interface FACE entry that owns the shared
    ! face on the destination side (doc 6a, uniform-B). Returns the normal
    ! direction d (1/2/3) for the two fine-low orientations and 0 otherwise:
    !   PROLONG off=+d : the fine block computes its own top face v(nb+1),
    !                    so the prolong fills only the nb+2 stencil halo for
    !                    the normal component (and nb+1 for tangential/p).
    !   RESTRICT off=-d: the coarse block's interface face v(1) is the
    !                    restriction of the four fine v(nb+1) (index-1 write).
    ! Edges/corners (|off|>1) and the already-correct fine-high orientations
    ! (PROLONG off=-d, RESTRICT off=+d) return 0 - handled as before.
    integer function face_normal(op, off)
        integer, intent(in) :: op, off(3)
        integer :: d
        face_normal = 0
        if (sum(abs(off)) /= 1) return
        do d = 1, 3
            if (op == OP_PROLONG .and. off(d) == 1) face_normal = d
            if (op == OP_RESTRICT .and. off(d) == -1) face_normal = d
        end do
    end function face_normal

    ! Destination box and per-dim source bases for a 2:1 interface entry.
    ! Tangential source indices follow from tq and the destination index
    ! inside the kernels; the normal-dimension base is stored in srcLo.
    subroutine interface_boxes(op, off, nb, srcLo, dstLo, ext)
        integer, intent(in) :: op, off(3), nb(3)
        integer, intent(out) :: srcLo(3), dstLo(3), ext(3)

        integer :: d

        do d = 1, 3
            select case (off(d))
            case (1)
                ! A fine block's high face toward a coarser neighbour
                ! (PROLONG) fills TWO halo layers so its redundant top-face
                ! momentum reaches the nb+2 stencil cell (doc 6a); both layers
                ! inject the same covering coarse value (the gather base is
                ! di-independent for PROLONG). A coarse block restricting from
                ! finer neighbours above keeps one layer (it never computes its
                ! own top face).
                dstLo(d) = nb(d) + 1
                ext(d) = merge(2, 1, op == OP_PROLONG)
                srcLo(d) = 1
            case (-1)
                ! A coarse block restricting from finer neighbours below
                ! (RESTRICT face) writes TWO layers: index 0 (the halo below
                ! the interface, tangential/pressure) and index 1 (its own
                ! interface face v(1) = restriction of the four fine v(nb+1),
                ! normal component). The gather slope ga=2 maps index 0 -> the
                ! halo rows and index 1 -> the fine top face nb+1.
                dstLo(d) = 0
                ext(d) = merge(2, 1, op == OP_RESTRICT .and. sum(abs(off)) == 1)
                srcLo(d) = merge(nb(d) - 1, nb(d), op == OP_RESTRICT)
            case default
                if (op == OP_RESTRICT) then
                    dstLo(d) = 1          ! quarter offset added via tq in the builder
                    ext(d) = nb(d)/2
                else
                    dstLo(d) = 1
                    ext(d) = nb(d)
                end if
                srcLo(d) = 1
            end select
        end do
    end subroutine interface_boxes

    ! Per-dimension affine gather map of one exchange entry: the source
    ! rows for destination index i are base..base+cnt-1 with
    ! base = ishft(ga*i + gb, -gs); variables staggered along a dimension
    ! sample 1 row there instead of cnt (the matching face). One form
    ! covers every transfer: same-level copies map index for index
    ! (ga=1), restrictions gather the two covering fine rows (ga=2
    ! tangentially, the srcLo row pair across the face), prolongations
    ! read the covering coarse row (halving shift tangentially, the
    ! tq-dependent covering row across the face).
    subroutine entry_gather_map(op, off, tq, nb, srcLo, dstLo, ga, gb, gs, gc, lin)
        integer, intent(in) :: op, off(3), tq(3), nb(3), srcLo(3), dstLo(3)
        integer, intent(out) :: ga(3), gb(3), gs(3), gc(3)
        integer, intent(out) :: lin(3)

        integer :: d

        do d = 1, 3
            gs(d) = 0
            gc(d) = 1
            ! Tangential dims of a PROLONG transfer can use linear interpolation
            ! (gated at runtime by comm%linProlong; injection otherwise). E3.
            lin(d) = merge(1, 0, op == OP_PROLONG .and. off(d) == 0)
            if (op == OP_COPY) then
                ga(d) = 1
                gb(d) = srcLo(d) - dstLo(d)
            else if (off(d) /= 0) then
                ga(d) = 0
                if (op == OP_RESTRICT) then
                    gb(d) = srcLo(d)
                    gc(d) = 2
                    ! Restrict interface FACE (coarse looking down at fine):
                    ! ga=2 sends dst index 0 to the halo rows (srcLo, +1) and
                    ! dst index 1 to the fine top face nb+1, so the coarse owns
                    ! v(1) = average of the four fine v(nb+1) (doc 6a).
                    if (face_normal(op, off) == d) ga(d) = 2
                else
                    ! Coarse row covering the halo layer: depends on the
                    ! fine block's parity tq because across edge/corner
                    ! directions the coarse neighbour spans past the
                    ! shared boundary (on faces 2:1 smoothing pins tq to
                    ! the touching side and this reduces to 1 or nb).
                    gb(d) = nb(d) - tq(d)*nb(d)/2
                    if (off(d) == 1) gb(d) = nb(d)/2 + 1 - tq(d)*nb(d)/2
                end if
            else
                if (op == OP_RESTRICT) then
                    ga(d) = 2
                    gb(d) = -tq(d)*nb(d) - 1
                    gc(d) = 2
                else
                    ga(d) = 1
                    gb(d) = tq(d)*nb(d) + 1
                    gs(d) = 1
                end if
            end if
        end do
    end subroutine entry_gather_map

    ! Ghost-interpolation weight for the pressure on a PROLONG face entry.
    ! Plain injection puts the coarse cell value at the fine halo centre,
    ! so the face-normal pressure gradient at the interface is evaluated
    ! over a gap 1.5x larger than the fine metric assumes - a systematic
    ! overdriving of the interface velocity that feeds back through the
    ! projection and blows up. Blending the coarse value with the first
    ! interior fine cell places the ghost where the fine stencil expects
    ! it (uniform 2:1: ghost = (2*coarse + fine)/3). Edges and corners
    ! keep plain injection: only face halos enter the pressure gradient.
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
                ! Fill TWO high-side halo layers (nb+1, nb+2) of the normal
                ! direction: the redundant top-face momentum (doc 6a) reaches
                ! v(nb+2). The gather (ga=1, gb=srcLo-dstLo=-nb) maps halo
                ! nb+1,nb+2 to neighbour interior 1,2. The low side stays one
                ! layer (the top-face stencil never reads below v(0)).
                srcLo(d) = 1
                dstLo(d) = nb(d) + 1
                ext(d) = 2
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
            ! Occupied at the same, coarser, or finer level all count: the
            ! corner data then arrives through an edge/corner entry (COPY,
            ! RESTRICT or PROLONG), so the face entry must not extend.
            exists = occupied_any_level(blk, level, to/nb)
        end if
    end function combined_neighbor_exists

    logical function occupied_any_level(blk, level, cl) result(occ)
        type(block_set_type), intent(in) :: blk
        integer, intent(in) :: level, cl(3)

        integer :: sx, sy, sz

        occ = leaf_at(blk, level, cl) >= 0_C_INT
        if (occ) return
        occ = leaf_at(blk, level - 1, cl/2) >= 0_C_INT
        if (occ) return
        do sz = 0, 1
            do sy = 0, 1
                do sx = 0, 1
                    if (leaf_at(blk, level + 1, 2*cl + [sx, sy, sz]) >= 0_C_INT) then
                        occ = .true.
                        return
                    end if
                end do
            end do
        end do
    end function occupied_any_level

    subroutine collect_peers(c, blk, dns, off)
        type(comm_type), intent(inout) :: c
        type(block_set_type), intent(in) :: blk
        type(dns_type), intent(in) :: dns
        integer, intent(in) :: off(3,26)

        integer :: b, d, p, cand, ncand
        integer :: owner(4), slot(4), opc(4), tqc(3,4)
        integer, allocatable :: found(:)
        logical :: known

        allocate(found(max(1, c%cart_size)))
        c%nPeers = 0
        do b = 1, int(blk%nBlocks)
            do d = 1, 26
                call resolve_neighbors(c, blk, dns, int(blk%level(b)), int(blk%origin(:,b)), &
                    off(:,d), ncand, owner, slot, opc, tqc)
                do cand = 1, ncand
                    if (owner(cand) == c%cart_rank) cycle
                    known = .false.
                    do p = 1, c%nPeers
                        if (found(p) == owner(cand)) then
                            known = .true.
                            exit
                        end if
                    end do
                    if (.not. known) then
                        c%nPeers = c%nPeers + 1
                        found(c%nPeers) = owner(cand)
                    end if
                end do
            end do
        end do

        allocate(c%peerRank(max(1, c%nPeers)))
        allocate(c%peerSendOff(0:max(1, c%nPeers)))
        allocate(c%peerRecvOff(0:max(1, c%nPeers)))
        allocate(c%peerSendCopyOff(0:max(1, c%nPeers)))
        allocate(c%peerRecvCopyOff(0:max(1, c%nPeers)))
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

        integer :: ierr, p, nv, nRecvPts, nSendPts

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
            ! A copy-only exchange is a prefix of the full message on both
            ! ends (entries are ordered same-level first per peer).
            do p = 1, c%nPeers
                nRecvPts = c%peerRecvOff(p) - c%peerRecvOff(p-1)
                if (c%copyOnly) nRecvPts = c%peerRecvCopyOff(p) - c%peerRecvCopyOff(p-1)
                call MPI_Irecv(c%recvbuf(1,p), nRecvPts*nv, &
                    MPI_DOUBLE_PRECISION, c%peerRank(p), HALO_TAG, c%cart_comm, c%request(p), ierr)
            end do
            do p = 1, c%nPeers
                nSendPts = c%peerSendOff(p) - c%peerSendOff(p-1)
                if (c%copyOnly) nSendPts = c%peerSendCopyOff(p) - c%peerSendCopyOff(p-1)
                call MPI_Isend(c%sendbuf(1,p), nSendPts*nv, &
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

    subroutine exchange_halos(c, blk, vars, interp, p_interface_only, linear_prolong)
        type(comm_type), intent(inout) :: c
        type(block_set_type), intent(inout) :: blk
        integer(C_INT), intent(in) :: vars(:)
        logical, intent(in), optional :: interp
        logical, intent(in), optional :: p_interface_only
        logical, intent(in), optional :: linear_prolong

        c%copyOnly = .false.
        if (present(interp)) c%copyOnly = .not. interp
        c%pSkipCopy = .false.
        if (present(p_interface_only)) c%pSkipCopy = p_interface_only
        c%linProlong = .false.
        if (present(linear_prolong)) c%linProlong = linear_prolong
        call start_halo_exchange(c, blk, vars)
        call finish_halo_exchange(c, blk)
        c%copyOnly = .false.
        c%pSkipCopy = .false.
        c%linProlong = .false.
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
        integer :: di, dj, dk, var, nv, totalItems, nd, layerN, lin

        nv = c%nActiveVars
        totalItems = merge(c%nLocalCopyPts, c%nLocalPts, c%copyOnly)*nv
        if (totalItems == 0) return

#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do &
        !$omp& map(to: totalItems, nv, c%nLocal, c%nLocalCopyPts, c%pSkipCopy, c%linProlong, &
        !$omp& c%lOff, c%lEntryOf, c%lSrcSlot, c%lDstSlot, &
        !$omp& c%lDstLo, c%lExt, c%lGA, c%lGB, c%lGS, c%lGC, c%lLin, c%lFaceNrm, &
        !$omp& c%activeVars, blk%nb) &
        !$omp& map(tofrom: blk%q) &
        !$omp& private(p,gp,v,e,pt,ni,nj,di,dj,dk,var,nd,layerN,lin)
#endif
        do p = 1, totalItems
            gp = (p - 1)/nv
            v = p - 1 - gp*nv
            var = int(c%activeVars(v+1))
            ! Skip pressure on the same-level COPY prefix when requested: its halo
            ! is never read between projection sweeps (only interface halos are).
            if (c%pSkipCopy .and. var == VAR_P .and. gp < c%nLocalCopyPts) cycle
            e = c%lEntryOf(gp)
            pt = gp - c%lOff(e-1)
            ni = c%lExt(1,e)
            nj = c%lExt(2,e)
            di = c%lDstLo(1,e) + modulo(pt, ni)
            dj = c%lDstLo(2,e) + modulo(pt/ni, nj)
            dk = c%lDstLo(3,e) + pt/(ni*nj)
            ! 2:1 interface owned face. PROLONG (dstLo=nb+1>0): the fine owns
            ! its near face, so the normal component is written only to the
            ! deep stencil layer nb+2 and the others to the nb+1 halo.
            ! RESTRICT (dstLo=0): the normal component writes BOTH the coarse's
            ! -side halo (index 0, needed by its own momentum stencil) and its
            ! interface face (index 1 = average of the fine faces); the others
            ! write only the index-0 halo.
            ! Above-block-owns: the block ABOVE a 2:1 shared face owns it at its
            ! interior v(1); the block below holds v(nb+1) as a slaved halo. The
            ! normal component is therefore filled like any halo -- PROLONG injects
            ! the covering coarse v(1) (both halo layers nb+1, nb+2); RESTRICT
            ! writes ONLY the v(0) halo (layerN 0), never the coarse interior v(1).
            nd = c%lFaceNrm(e)
            if (nd /= 0) then
                layerN = merge(di, merge(dj, dk, nd == 2), nd == 1) - c%lDstLo(nd,e)
                if (var == nd) then
                    if (c%lDstLo(nd,e) == 0 .and. layerN /= 0) cycle
                else if (layerN /= 0) then
                    cycle
                end if
            end if
            ! Linear PROLONG only when enabled (final exchange) and only for the
            ! velocity components; pressure and all non-final exchanges inject.
            lin = merge(1, 0, c%linProlong .and. var /= VAR_P)
            blk%q(di,dj,dk,var,c%lDstSlot(e)) = gather_point(blk%q, c%lSrcSlot(e), &
                c%lGA(:,e), c%lGB(:,e), c%lGS(:,e), c%lGC(:,e), c%lLin(:,e), &
                lin, di, dj, dk, var, blk%nb)
        end do
#ifdef USE_OPENMP_OFFLOAD
        !$omp end target teams distribute parallel do
#endif
    end subroutine copy_local_entries

    subroutine copy_local_scalar_entries(c, scalar)
        type(comm_type), intent(inout) :: c
        real(C_DOUBLE), intent(inout) :: scalar(0:,0:,0:,1:)

        integer :: p, e, pt, ni, nj
        integer :: di, dj, dk, totalItems, nd, layerN
        integer :: b1, b2, b3, c1, c2, c3, s1, s2, s3
        real(C_DOUBLE) :: val

        totalItems = c%nLocalPts
        if (totalItems == 0) return

#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do &
        !$omp& map(to: totalItems, c%nLocal, c%lOff, c%lEntryOf, c%lSrcSlot, c%lDstSlot, &
        !$omp& c%lDstLo, c%lExt, c%lGA, c%lGB, c%lGS, c%lGC, c%lFaceNrm) &
        !$omp& map(tofrom: scalar) &
        !$omp& private(p,e,pt,ni,nj,di,dj,dk,nd,layerN,b1,b2,b3,c1,c2,c3,s1,s2,s3,val)
#endif
        do p = 1, totalItems
            e = c%lEntryOf(p - 1)
            pt = p - 1 - c%lOff(e-1)
            ni = c%lExt(1,e)
            nj = c%lExt(2,e)
            di = c%lDstLo(1,e) + modulo(pt, ni)
            dj = c%lDstLo(2,e) + modulo(pt/ni, nj)
            dk = c%lDstLo(3,e) + pt/(ni*nj)
            ! A cell-centred scalar (LES nut) at a 2:1 owned face writes only
            ! the shallow halo layer; the deep layer is the velocity normal
            ! component's (PROLONG nb+2 / RESTRICT coarse interior v(1)).
            nd = c%lFaceNrm(e)
            if (nd /= 0) then
                layerN = merge(di, merge(dj, dk, nd == 2), nd == 1) - c%lDstLo(nd,e)
                if (layerN /= 0) cycle
            end if
            b1 = ishft(c%lGA(1,e)*di + c%lGB(1,e), -c%lGS(1,e))
            b2 = ishft(c%lGA(2,e)*dj + c%lGB(2,e), -c%lGS(2,e))
            b3 = ishft(c%lGA(3,e)*dk + c%lGB(3,e), -c%lGS(3,e))
            c1 = c%lGC(1,e)
            c2 = c%lGC(2,e)
            c3 = c%lGC(3,e)
            val = 0.0d0
            do s3 = 0, c3 - 1
                do s2 = 0, c2 - 1
                    do s1 = 0, c1 - 1
                        val = val + scalar(b1+s1, b2+s2, b3+s3, c%lSrcSlot(e))
                    end do
                end do
            end do
            scalar(di,dj,dk,c%lDstSlot(e)) = val/real(c1*c2*c3, C_DOUBLE)
        end do
#ifdef USE_OPENMP_OFFLOAD
        !$omp end target teams distribute parallel do
#endif
    end subroutine copy_local_scalar_entries

    subroutine pack_entries(c, blk)
        type(comm_type), intent(inout) :: c
        type(block_set_type), intent(in) :: blk

        integer :: p, gp, v, e, pt, ni, nj
        integer :: di, dj, dk, var, peer, pos, nv, totalItems, copyOnly, lin
        real(C_DOUBLE) :: val

        nv = c%nActiveVars
        totalItems = merge(c%peerSendCopyOff(c%nPeers), c%peerSendOff(c%nPeers), c%copyOnly)*nv
        if (totalItems == 0) return
        copyOnly = merge(1, 0, c%copyOnly)

#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do &
        !$omp& map(to: totalItems, nv, copyOnly, c%pSkipCopy, c%linProlong, c%nPeers, c%nSend, &
        !$omp& c%sOff, c%sEntryOf, c%sSlot, c%sPeer, &
        !$omp& c%sDstLo, c%sExt, c%sGA, c%sGB, c%sGS, c%sGC, c%sLin, &
        !$omp& c%peerSendOff, c%peerSendCopyOff, c%activeVars, blk%q, blk%nb) &
        !$omp& map(tofrom: c%sendbuf) &
        !$omp& private(p,gp,v,e,pt,ni,nj,di,dj,dk,var,peer,pos,lin,val)
#endif
        do p = 1, totalItems
            gp = (p - 1)/nv
            v = p - 1 - gp*nv
            ! A copy-only point index maps into the full enumeration via
            ! the per-peer copy prefixes.
            if (copyOnly == 1) then
                peer = find_entry(c%peerSendCopyOff, c%nPeers, gp)
                gp = c%peerSendOff(peer-1) + (gp - c%peerSendCopyOff(peer-1))
            end if
            e = c%sEntryOf(gp)
            peer = c%sPeer(e)
            var = int(c%activeVars(v+1))
            ! Skip pressure on the same-level COPY prefix of this peer (see
            ! copy_local_entries): leaves the buffer slot stale, which the
            ! receiver's unpack skips identically.
            if (c%pSkipCopy .and. var == VAR_P .and. &
                gp - c%peerSendOff(peer-1) < c%peerSendCopyOff(peer) - c%peerSendCopyOff(peer-1)) cycle
            pt = gp - c%sOff(e-1)
            ni = c%sExt(1,e)
            nj = c%sExt(2,e)
            di = c%sDstLo(1,e) + modulo(pt, ni)
            dj = c%sDstLo(2,e) + modulo(pt/ni, nj)
            dk = c%sDstLo(3,e) + pt/(ni*nj)
            lin = merge(1, 0, c%linProlong .and. var /= VAR_P)
            val = gather_point(blk%q, c%sSlot(e), c%sGA(:,e), c%sGB(:,e), c%sGS(:,e), &
                c%sGC(:,e), c%sLin(:,e), lin, di, dj, dk, var, blk%nb)
            pos = (gp - c%peerSendOff(peer-1))*nv + v + 1
            c%sendbuf(pos,peer) = val
        end do
#ifdef USE_OPENMP_OFFLOAD
        !$omp end target teams distribute parallel do
#endif
    end subroutine pack_entries

    subroutine unpack_entries(c, blk)
        type(comm_type), intent(in) :: c
        type(block_set_type), intent(inout) :: blk

        integer :: p, gp, v, e, pt, ni, nj
        integer :: i, j, k, var, peer, pos, nv, totalItems, copyOnly, nd, layerN
        real(C_DOUBLE) :: val

        nv = c%nActiveVars
        totalItems = merge(c%peerRecvCopyOff(c%nPeers), c%peerRecvOff(c%nPeers), c%copyOnly)*nv
        if (totalItems == 0) return
        copyOnly = merge(1, 0, c%copyOnly)

#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do &
        !$omp& map(to: totalItems, nv, copyOnly, c%pSkipCopy, c%nPeers, c%nRecv, &
        !$omp& c%rOff, c%rEntryOf, c%rSlot, c%rPeer, &
        !$omp& c%rLo, c%rExt, c%rFaceNrm, &
        !$omp& c%peerRecvOff, c%peerRecvCopyOff, c%activeVars, c%recvbuf) &
        !$omp& map(tofrom: blk%q) &
        !$omp& private(p,gp,v,e,pt,ni,nj,i,j,k,var,peer,pos,nd,layerN,val)
#endif
        do p = 1, totalItems
            gp = (p - 1)/nv
            v = p - 1 - gp*nv
            if (copyOnly == 1) then
                peer = find_entry(c%peerRecvCopyOff, c%nPeers, gp)
                gp = c%peerRecvOff(peer-1) + (gp - c%peerRecvCopyOff(peer-1))
            end if
            e = c%rEntryOf(gp)
            peer = c%rPeer(e)
            var = int(c%activeVars(v+1))
            ! Skip pressure on the same-level COPY prefix (mirror of the sender's
            ! pack_entries skip): the buffer slot was left stale, never read here.
            if (c%pSkipCopy .and. var == VAR_P .and. &
                gp - c%peerRecvOff(peer-1) < c%peerRecvCopyOff(peer) - c%peerRecvCopyOff(peer-1)) cycle
            pt = gp - c%rOff(e-1)
            ni = c%rExt(1,e)
            nj = c%rExt(2,e)
            i = c%rLo(1,e) + modulo(pt, ni)
            j = c%rLo(2,e) + modulo(pt/ni, nj)
            k = c%rLo(3,e) + pt/(ni*nj)
            ! Above-block-owns (mirror of copy_local_entries): the normal
            ! component is a slaved halo on the BELOW block -- PROLONG injects the
            ! covering coarse v(1) into both halo layers; RESTRICT writes only the
            ! v(0) halo, never the coarse interior v(1) (it owns it).
            nd = c%rFaceNrm(e)
            if (nd /= 0) then
                layerN = merge(i, merge(j, k, nd == 2), nd == 1) - c%rLo(nd,e)
                if (var == nd) then
                    if (c%rLo(nd,e) == 0 .and. layerN /= 0) cycle
                else if (layerN /= 0) then
                    cycle
                end if
            end if
            pos = (gp - c%peerRecvOff(peer-1))*nv + v + 1
            val = c%recvbuf(pos,peer)
            blk%q(i,j,k,var,c%rSlot(e)) = val
        end do
#ifdef USE_OPENMP_OFFLOAD
        !$omp end target teams distribute parallel do
#endif
    end subroutine unpack_entries

    subroutine pack_scalar_entries(c, scalar)
        type(comm_type), intent(inout) :: c
        real(C_DOUBLE), intent(in) :: scalar(0:,0:,0:,1:)

        integer :: p, e, pt, ni, nj
        integer :: di, dj, dk, peer, pos, totalItems, nd, layerN
        integer :: b1, b2, b3, c1, c2, c3, s1, s2, s3
        real(C_DOUBLE) :: val

        totalItems = c%peerSendOff(c%nPeers)
        if (totalItems == 0) return

#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do &
        !$omp& map(to: totalItems, c%nSend, c%sOff, c%sEntryOf, c%sSlot, c%sPeer, &
        !$omp& c%sDstLo, c%sExt, c%sGA, c%sGB, c%sGS, c%sGC, c%sFaceNrm, c%peerSendOff, scalar) &
        !$omp& map(tofrom: c%sendbuf) &
        !$omp& private(p,e,pt,ni,nj,di,dj,dk,peer,pos,nd,layerN,b1,b2,b3,c1,c2,c3,s1,s2,s3,val)
#endif
        do p = 1, totalItems
            e = c%sEntryOf(p - 1)
            pt = p - 1 - c%sOff(e-1)
            ni = c%sExt(1,e)
            nj = c%sExt(2,e)
            di = c%sDstLo(1,e) + modulo(pt, ni)
            dj = c%sDstLo(2,e) + modulo(pt/ni, nj)
            dk = c%sDstLo(3,e) + pt/(ni*nj)
            ! Cell-centred scalar packs only the shallow halo layer; the deep
            ! layer's gather would read nut(nb+2), outside its 0:nb+1 bound.
            nd = c%sFaceNrm(e)
            if (nd /= 0) then
                layerN = merge(di, merge(dj, dk, nd == 2), nd == 1) - c%sDstLo(nd,e)
                if (layerN /= 0) cycle
            end if
            b1 = ishft(c%sGA(1,e)*di + c%sGB(1,e), -c%sGS(1,e))
            b2 = ishft(c%sGA(2,e)*dj + c%sGB(2,e), -c%sGS(2,e))
            b3 = ishft(c%sGA(3,e)*dk + c%sGB(3,e), -c%sGS(3,e))
            c1 = c%sGC(1,e)
            c2 = c%sGC(2,e)
            c3 = c%sGC(3,e)
            val = 0.0d0
            do s3 = 0, c3 - 1
                do s2 = 0, c2 - 1
                    do s1 = 0, c1 - 1
                        val = val + scalar(b1+s1, b2+s2, b3+s3, c%sSlot(e))
                    end do
                end do
            end do
            peer = c%sPeer(e)
            pos = p - c%peerSendOff(peer-1)
            c%sendbuf(pos,peer) = val/real(c1*c2*c3, C_DOUBLE)
        end do
#ifdef USE_OPENMP_OFFLOAD
        !$omp end target teams distribute parallel do
#endif
    end subroutine pack_scalar_entries

    subroutine unpack_scalar_entries(c, scalar)
        type(comm_type), intent(in) :: c
        real(C_DOUBLE), intent(inout) :: scalar(0:,0:,0:,1:)

        integer :: p, e, pt, ni, nj
        integer :: i, j, k, peer, pos, totalItems, nd, layerN

        totalItems = c%peerRecvOff(c%nPeers)
        if (totalItems == 0) return

#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do &
        !$omp& map(to: totalItems, c%nRecv, c%rOff, c%rEntryOf, c%rSlot, c%rPeer, &
        !$omp& c%rLo, c%rExt, c%rFaceNrm, c%peerRecvOff, c%recvbuf) &
        !$omp& map(tofrom: scalar) &
        !$omp& private(p,e,pt,ni,nj,i,j,k,peer,pos,nd,layerN)
#endif
        do p = 1, totalItems
            e = c%rEntryOf(p - 1)
            pt = p - 1 - c%rOff(e-1)
            ni = c%rExt(1,e)
            nj = c%rExt(2,e)
            i = c%rLo(1,e) + modulo(pt, ni)
            j = c%rLo(2,e) + modulo(pt/ni, nj)
            k = c%rLo(3,e) + pt/(ni*nj)
            ! Cell-centred scalar: shallow halo layer only (see copy_local).
            nd = c%rFaceNrm(e)
            if (nd /= 0) then
                layerN = merge(i, merge(j, k, nd == 2), nd == 1) - c%rLo(nd,e)
                if (layerN /= 0) cycle
            end if
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

    ! Source taps (up to 2) for one dimension of the weighted gather. lin=0
    ! reproduces the plain affine map (inject / gc-cell equal average, with the
    ! face-staggered c=1 override) bit-exactly. lin=1 (PROLONG tangential dim,
    ! only when linProlong is set) uses linear interpolation of the coarse
    ! field: face-linear for the variable's own staggered dim (fine face
    ! coincident -> 1 tap; midpoint -> 1/2,1/2), cell-linear otherwise (3/4
    ! covering + 1/4 adjacent coarse cell). E3.
    pure subroutine gather_taps(ga, gb, gs, gc, lin, dst, var, dim, nbDim, idx, w, n)
!$omp declare target
        integer, intent(in) :: ga, gb, gs, gc, lin, dst, var, dim, nbDim
        integer, intent(out) :: idx(0:1), n
        real(C_DOUBLE), intent(out) :: w(0:1)
        integer :: raw, base, par, cd, hi

        raw = ga*dst + gb
        base = ishft(raw, -gs)
        idx(0) = base; idx(1) = base; w(0) = 1.0d0; w(1) = 0.0d0; n = 1
        if (lin == 1) then
            par = raw - 2*ishft(raw, -1)
            if (var == dim) then
                if (par /= 0) then
                    idx(1) = base + 1; w(0) = 0.5d0; w(1) = 0.5d0; n = 2
                end if
            else
                idx(1) = base + (2*par - 1); w(0) = 0.75d0; w(1) = 0.25d0; n = 2
            end if
            ! Keep the adjacent tap inside the source block's interior cells
            ! 1..nb. A halo there (e.g. the staggered face nb+1) is never read by
            ! injection, so it can be rank-dependently stale during the exchange;
            ! reading it would break exact rank independence. Fall back to
            ! injection at the prolong region's edges instead.
            hi = nbDim
            if (idx(1) < 1 .or. idx(1) > hi) idx(1) = idx(0)
        else
            cd = merge(1, gc, var == dim)
            n = cd; w(0) = 1.0d0/real(cd, C_DOUBLE)
            if (cd == 2) then; idx(1) = base + 1; w(1) = 0.5d0; end if
        end if
    end subroutine gather_taps

    ! Weighted gather of one destination point (di,dj,dk) of variable var from
    ! source slot. Shared by the local-copy and off-rank-pack kernels so the
    ! arithmetic is compiled once and both produce bit-identical results (the
    ! same source inlined into two kernels can otherwise reassociate the
    ! weighted sum differently, breaking exact rank independence). linFlag=1
    ! enables linear PROLONG on the tangential dims flagged by lin(:).
    pure real(C_DOUBLE) function gather_point(q, slot, ga, gb, gs, gc, lin, &
            linFlag, di, dj, dk, var, nb) result(val)
!$omp declare target
        real(C_DOUBLE), intent(in) :: q(0:,0:,0:,1:,1:)
        integer, intent(in) :: slot, ga(3), gb(3), gs(3), gc(3), lin(3)
        integer, intent(in) :: linFlag, di, dj, dk, var, nb(3)
        integer :: i1(0:1), i2(0:1), i3(0:1), n1, n2, n3, a1, a2, a3
        real(C_DOUBLE) :: w1(0:1), w2(0:1), w3(0:1)

        call gather_taps(ga(1), gb(1), gs(1), gc(1), linFlag*lin(1), di, var, 1, nb(1), i1, w1, n1)
        call gather_taps(ga(2), gb(2), gs(2), gc(2), linFlag*lin(2), dj, var, 2, nb(2), i2, w2, n2)
        call gather_taps(ga(3), gb(3), gs(3), gc(3), linFlag*lin(3), dk, var, 3, nb(3), i3, w3, n3)
        val = 0.0d0
        do a3 = 0, n3 - 1
            do a2 = 0, n2 - 1
                do a1 = 0, n1 - 1
                    val = val + w1(a1)*w2(a2)*w3(a3)*q(i1(a1), i2(a2), i3(a3), var, slot)
                end do
            end do
        end do
    end function gather_point

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
