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
        integer, allocatable :: lGA(:,:), lGB(:,:), lGS(:,:), lGC(:,:) ! gather map (3,nLocal)
        integer, allocatable :: lDir(:,:)                  ! direction (dst-completion adjacency)
        real(C_DOUBLE), allocatable :: lWp(:), lWpDst(:)   ! pressure dst-completion weights
        integer, allocatable :: lOff(:)                    ! (0:nLocal) point prefix
        ! Per-point -> entry index, precomputed once so the exchange kernels do an
        ! O(1) lookup instead of a per-point binary search over lOff/sOff/rOff
        ! (find_entry). The exchange is the dominant refined-grid cost and runs many
        ! times per step; the search was ~half the per-point work.
        integer, allocatable :: lPointEntry(:)             ! (0:nLocalPts-1)
        integer, allocatable :: sPointEntry(:)             ! (0:totalSendPts-1)
        integer, allocatable :: rPointEntry(:)             ! (0:totalRecvPts-1)
        ! Fine-owned interface normal velocity: dimension d of the normal
        ! velocity component this entry must NOT overwrite (a PROLONG onto a
        ! fine block's high interface face, off(d)==+1), else 0. The fine
        ! block owns that face (predictor + jacobi reconstruction); letting
        ! the coarse value prolong over it differences a coarse-resolution
        ! flux into the fine divergence (an O(1) interface inconsistency).
        integer, allocatable :: lNrm(:)                    ! (nLocal)
        ! Signed normal dim (+d/-d) of an interface RESTRICT pure-face entry,
        ! else 0. Used only by the projection's phi scalar exchange
        ! (phiIfaceRow): the coarse interface-correction ghost must be the
        ! restrict of the FINE INTERFACE ROW (the single fine row touching the
        ! face), not the cell-centred 2-normal-row restrict -- otherwise the
        ! coarse and fine interface corrections are not exactly mean-equal and
        ! the projection leaks a little mass. The averaging happens on the
        ! source side, so only the local-copy and send (pack) paths apply it.
        integer, allocatable :: lPhiN(:)                   ! (nLocal)

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
        integer, allocatable :: sDstLo(:,:)                ! dst box lo (gather indexing)
        integer, allocatable :: sOff(:)                    ! (0:nSend) point prefix, peer-major
        integer, allocatable :: sPhiN(:)                   ! (nSend) see lPhiN
        integer, allocatable :: rSlot(:), rPeer(:)
        integer, allocatable :: rLo(:,:), rExt(:,:)
        integer, allocatable :: rDir(:,:)
        real(C_DOUBLE), allocatable :: rWp(:), rWpDst(:)
        integer, allocatable :: rOff(:)
        integer, allocatable :: rNrm(:)                    ! (nRecv) see lNrm

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
        ! When true, the interface normal-velocity skip (lNrm/rNrm) is
        ! DISABLED, so the cross-level PROLONG (orientation A) and RESTRICT
        ! (orientation B) write the shared face again -- the conservation
        ! SYNC. The post-predictor exchange sets it so the two stored copies
        ! of each 2:1 face start the projection mean-consistent
        ! (avg(fine)=coarse); the projection then OWNS the face (default,
        ! skip) and the composite stencil's mean-equal corrections keep it
        ! conservative throughout.
        logical :: syncFace = .false.
        ! When true, the scalar (phi) exchange restricts the FINE INTERFACE
        ! ROW into the coarse interface-correction ghost (see lPhiN). Off for
        ! cell-centred scalars (e.g. les%nut).
        logical :: phiIfaceRow = .false.
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
        integer :: nLocal, nSend, nRecv, pts, maxCount, ierr, e, round

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
                                c%lDir(:,nLocal) = off(:,d)
                                call entry_gather_map(opc(cand), off(:,d), tqc(:,cand), nb, &
                                    srcLo, dstLo, c%lGA(:,nLocal), c%lGB(:,nLocal), &
                                    c%lGS(:,nLocal), c%lGC(:,nLocal))
                                c%lWp(nLocal) = entry_blend(blk, dns, int(blk%level(b)), &
                                    int(blk%origin(:,b)), off(:,d), opc(cand))
                                c%lWpDst(nLocal) = 1.0d0 - c%lWp(nLocal)
                                c%lNrm(nLocal) = interface_normal_dim(opc(cand), off(:,d))
                                c%lPhiN(nLocal) = iface_restrict_normal(opc(cand), off(:,d))
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
                                    c%rDir(:,nRecv) = off(:,d)
                                    c%rWp(nRecv) = entry_blend(blk, dns, int(blk%level(b)), &
                                        int(blk%origin(:,b)), off(:,d), opc(cand))
                                    c%rWpDst(nRecv) = 1.0d0 - c%rWp(nRecv)
                                    c%rNrm(nRecv) = interface_normal_dim(opc(cand), off(:,d))
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
                                    call entry_gather_map(opc(cand), off(:,d), tqc(:,cand), nb, &
                                        srcLo, dstLo, c%sGA(:,nSend), c%sGB(:,nSend), &
                                        c%sGS(:,nSend), c%sGC(:,nSend))
                                    c%sPhiN(nSend) = iface_restrict_normal(opc(cand), off(:,d))
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
                allocate(c%lDir(3,max(1,nLocal)))
                allocate(c%lWp(max(1,nLocal)), c%lWpDst(max(1,nLocal)))
                allocate(c%lNrm(max(1,nLocal)))
                allocate(c%lPhiN(max(1,nLocal)))
                allocate(c%lOff(0:max(1,nLocal)))
                allocate(c%sSlot(max(1,nSend)), c%sPeer(max(1,nSend)))
                allocate(c%sExt(3,max(1,nSend)))
                allocate(c%sGA(3,max(1,nSend)), c%sGB(3,max(1,nSend)))
                allocate(c%sGS(3,max(1,nSend)), c%sGC(3,max(1,nSend)))
                allocate(c%sDstLo(3,max(1,nSend)))
                allocate(c%sPhiN(max(1,nSend)))
                allocate(c%sOff(0:max(1,nSend)))
                allocate(c%rSlot(max(1,nRecv)), c%rPeer(max(1,nRecv)))
                allocate(c%rLo(3,max(1,nRecv)), c%rExt(3,max(1,nRecv)))
                allocate(c%rDir(3,max(1,nRecv)))
                allocate(c%rWp(max(1,nRecv)), c%rWpDst(max(1,nRecv)))
                allocate(c%rNrm(max(1,nRecv)))
                allocate(c%rOff(0:max(1,nRecv)))
                c%lWp = 1.0d0
                c%lWpDst = 0.0d0
                c%rWp = 1.0d0
                c%rWpDst = 0.0d0
                c%lNrm = 0
                c%rNrm = 0
                c%lPhiN = 0
                c%sPhiN = 0
                c%lOff = 0
                c%sOff = 0
                c%rOff = 0
                c%peerSendOff = 0
                c%peerRecvOff = 0
                c%peerSendCopyOff = 0
                c%peerRecvCopyOff = 0
            end if
        end do

        ! Precompute the per-point -> entry maps (one-time; replaces the per-point
        ! find_entry binary search inside the exchange kernels).
        allocate(c%lPointEntry(0:max(0, c%nLocalPts - 1)))
        allocate(c%sPointEntry(0:max(0, c%peerSendOff(c%nPeers) - 1)))
        allocate(c%rPointEntry(0:max(0, c%peerRecvOff(c%nPeers) - 1)))
        do e = 1, c%nLocal
            do p = c%lOff(e-1), c%lOff(e) - 1
                c%lPointEntry(p) = e
            end do
        end do
        do e = 1, c%nSend
            do p = c%sOff(e-1), c%sOff(e) - 1
                c%sPointEntry(p) = e
            end do
        end do
        do e = 1, c%nRecv
            do p = c%rOff(e-1), c%rOff(e) - 1
                c%rPointEntry(p) = e
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
        !$omp& c%lPointEntry, c%sPointEntry, c%rPointEntry, &
        !$omp& c%lSrcSlot, c%lDstSlot, c%lDstLo, c%lExt, c%lOff, &
        !$omp& c%lGA, c%lGB, c%lGS, c%lGC, c%lDir, c%lWp, c%lWpDst, c%lNrm, c%lPhiN, &
        !$omp& c%sSlot, c%sPeer, c%sExt, c%sOff, &
        !$omp& c%sGA, c%sGB, c%sGS, c%sGC, c%sDstLo, c%sPhiN, &
        !$omp& c%peerSendOff, c%peerRecvOff, c%peerSendCopyOff, c%peerRecvCopyOff, &
        !$omp& c%rSlot, c%rPeer, c%rLo, c%rExt, c%rDir, c%rWp, c%rWpDst, c%rOff, c%rNrm)
        !$omp target enter data map(alloc: c%sendbuf, c%recvbuf)
#endif
    end subroutine init_block_exchange

    subroutine free_block_exchange(c)
        type(comm_type), intent(inout) :: c

        if (allocated(c%sendbuf)) then
#ifdef USE_OPENMP_OFFLOAD
            !$omp target exit data map(delete: c%sendbuf, c%recvbuf)
            !$omp target exit data map(delete: &
            !$omp& c%lPointEntry, c%sPointEntry, c%rPointEntry, &
            !$omp& c%lSrcSlot, c%lDstSlot, c%lDstLo, c%lExt, c%lOff, &
            !$omp& c%lGA, c%lGB, c%lGS, c%lGC, c%lDir, c%lWp, c%lWpDst, c%lNrm, c%lPhiN, &
            !$omp& c%sSlot, c%sPeer, c%sExt, c%sOff, &
            !$omp& c%sGA, c%sGB, c%sGS, c%sGC, c%sDstLo, c%sPhiN, &
            !$omp& c%peerSendOff, c%peerRecvOff, c%peerSendCopyOff, c%peerRecvCopyOff, &
            !$omp& c%rSlot, c%rPeer, c%rLo, c%rExt, c%rDir, c%rWp, c%rWpDst, c%rOff, c%rNrm)
#endif
            deallocate(c%sendbuf, c%recvbuf)
            deallocate(c%lPointEntry, c%sPointEntry, c%rPointEntry)
            deallocate(c%lSrcSlot, c%lDstSlot, c%lDstLo, c%lExt, c%lOff)
            deallocate(c%lGA, c%lGB, c%lGS, c%lGC, c%lDir, c%lWp, c%lWpDst, c%lNrm, c%lPhiN)
            deallocate(c%sSlot, c%sPeer, c%sExt, c%sOff)
            deallocate(c%sGA, c%sGB, c%sGS, c%sGC, c%sDstLo, c%sPhiN)
            deallocate(c%rSlot, c%rPeer, c%rLo, c%rExt, c%rDir, c%rWp, c%rWpDst, c%rOff, c%rNrm)
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
                dstLo(d) = nb(d) + 1
                ext(d) = 1
                srcLo(d) = 1
            case (-1)
                dstLo(d) = 0
                ext(d) = 1
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
    subroutine entry_gather_map(op, off, tq, nb, srcLo, dstLo, ga, gb, gs, gc)
        integer, intent(in) :: op, off(3), tq(3), nb(3), srcLo(3), dstLo(3)
        integer, intent(out) :: ga(3), gb(3), gs(3), gc(3)

        integer :: d

        do d = 1, 3
            gs(d) = 0
            gc(d) = 1
            if (op == OP_COPY) then
                ga(d) = 1
                gb(d) = srcLo(d) - dstLo(d)
            else if (off(d) /= 0) then
                ga(d) = 0
                if (op == OP_RESTRICT) then
                    gb(d) = srcLo(d)
                    gc(d) = 2
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

    ! Normal-velocity dimension a cross-level face entry must leave to its own
    ! block (owned), or 0. At a 2:1 face the block on the LOW-coordinate side
    ! holds the shared normal face as its HIGH halo q(nb+1) in dim d; that face
    ! is its own cell's interface normal velocity (read by its divergence) and
    ! is reconstructed inside the projection from its pressure plus the
    ! cross-level ghost pressure (ifGrad). Letting the neighbour clobber it
    ! differences a foreign-resolution flux into the divergence -- the O(1)
    ! interface inconsistency. This happens for BOTH orientations of a pure
    ! face (off(d)==+1): the FINE block's halo is fed by PROLONG (orientation
    ! A, fine below a coarse -- the fine cell would read a coarse-resolution
    ! flux), the COARSE block's halo by RESTRICT (orientation B, coarse below a
    ! fine -- the coarse cell would read the averaged fine flux, slaved, so its
    ! pressure cannot drive it divergence-free). Both are skipped; the
    ! interface corrections stay mean-equal (avg(prolong)=coarse, avg(fine)=
    ! restrict), so conservation is preserved if the two stored copies start
    ! consistent. off(d)==-1 entries write the LOW halo q(0), below the owned
    ! interface face q(1), a genuine momentum halo, so they are kept.
    pure integer function interface_normal_dim(op, off) result(dnorm)
        integer, intent(in) :: op, off(3)
        integer :: d
        dnorm = 0
        if (op /= OP_PROLONG .and. op /= OP_RESTRICT) return
        if (sum(abs(off)) /= 1) return
        do d = 1, 3
            if (off(d) == 1) dnorm = d
        end do
    end function interface_normal_dim

    ! Signed normal dim (off(d)*d) of an interface RESTRICT pure-face entry,
    ! else 0. The sign tells the phi exchange which of the two cell-centred
    ! source rows is the one touching the interface: for off=+1 the lower row
    ! (base, no shift), for off=-1 the upper row (base + 1).
    pure integer function iface_restrict_normal(op, off) result(sd)
        integer, intent(in) :: op, off(3)
        integer :: d
        sd = 0
        if (op /= OP_RESTRICT) return
        if (sum(abs(off)) /= 1) return
        do d = 1, 3
            if (off(d) /= 0) sd = off(d)*d
        end do
    end function iface_restrict_normal

    ! Ghost-interpolation weight for the pressure on a PROLONG face entry.
    ! Plain injection puts the coarse cell value at the fine halo centre,
    ! so the face-normal pressure gradient at the interface is evaluated
    ! over a gap 1.5x larger than the fine metric assumes - a systematic
    ! overdriving of the interface velocity that feeds back through the
    ! projection and blows up. Blending the coarse value with the first
    ! interior fine cell places the ghost where the fine stencil expects
    ! it (uniform 2:1: ghost = (2*coarse + fine)/3). Edges and corners
    ! keep plain injection: only face halos enter the pressure gradient.
    function entry_blend(blk, dns, level, dorigin, off, op) result(w)
        type(block_set_type), intent(in) :: blk
        type(dns_type), intent(in) :: dns
        integer, intent(in) :: level, dorigin(3), off(3), op
        real(C_DOUBLE) :: w

        integer :: d, g, gnl
        real(C_DOUBLE) :: aHalf, bHalf, cHalf

        w = 1.0d0
        if (op /= OP_PROLONG) return
        if (sum(abs(off)) /= 1) return
        do d = 1, 3
            if (off(d) == 0) cycle
            gnl = int(level_cells(dns, d, int(level, C_INT)))
            if (off(d) == -1) then
                g = modulo(dorigin(d) - 1, gnl)
            else
                g = modulo(dorigin(d) + int(blk%nb(d)), gnl)
            end if
            ! Half-widths from the level node lines: fine halo cell (bHalf),
            ! its covering coarse cell (aHalf), first interior cell (cHalf).
            ! All three centres sit on the face normal; the face is a shared
            ! node of both lines, so distances reduce to half-widths.
            bHalf = 0.5d0*level_cell_width(blk, d, level, g)
            aHalf = 0.5d0*level_cell_width(blk, d, level - 1, g/2)
            if (off(d) == -1) then
                cHalf = 0.5d0*level_cell_width(blk, d, level, dorigin(d))
            else
                cHalf = 0.5d0*level_cell_width(blk, d, level, dorigin(d) + int(blk%nb(d)) - 1)
            end if
            w = 1.0d0 - (aHalf - bHalf)/(aHalf + cHalf)
        end do
    end function entry_blend

    ! Width of 0-based cell g of the level-l node line in direction d.
    function level_cell_width(blk, d, level, g) result(width)
        type(block_set_type), intent(in) :: blk
        integer, intent(in) :: d, level, g
        real(C_DOUBLE) :: width

        select case (d)
        case (1)
            width = blk%lineX(g + 1, level + 1) - blk%lineX(g, level + 1)
        case (2)
            width = blk%lineY(g + 1, level + 1) - blk%lineY(g, level + 1)
        case default
            width = blk%lineZ(g + 1, level + 1) - blk%lineZ(g, level + 1)
        end select
    end function level_cell_width

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

    subroutine exchange_halos(c, blk, vars, interp, syncface)
        type(comm_type), intent(inout) :: c
        type(block_set_type), intent(inout) :: blk
        integer(C_INT), intent(in) :: vars(:)
        logical, intent(in), optional :: interp, syncface

        c%copyOnly = .false.
        if (present(interp)) c%copyOnly = .not. interp
        c%syncFace = .false.
        if (present(syncface)) c%syncFace = syncface
        call start_halo_exchange(c, blk, vars)
        call finish_halo_exchange(c, blk)
        c%copyOnly = .false.
        c%syncFace = .false.
    end subroutine exchange_halos

    ! Per-block scalar halos (e.g. les%nut, the projection's phi) on the same
    ! exchange entries. ifaceRow restricts the fine INTERFACE ROW (not the
    ! cell-centred 2-row average) into the coarse interface-correction ghost --
    ! used for phi so the interface corrections are conservative; off for
    ! cell-centred scalars like les%nut.
    subroutine exchange_scalar_halos(c, scalar, ifaceRow)
        type(comm_type), intent(inout) :: c
        real(C_DOUBLE), intent(inout) :: scalar(0:,0:,0:,1:)
        logical, intent(in), optional :: ifaceRow

        integer :: ierr, p, nRequest

        call require_ready(c)
        if (c%exchangeActive) error stop "halo exchange already active"

        c%phiIfaceRow = .false.
        if (present(ifaceRow)) c%phiIfaceRow = ifaceRow

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
        c%phiIfaceRow = .false.
    end subroutine exchange_scalar_halos

    subroutine copy_local_entries(c, blk)
        type(comm_type), intent(inout) :: c
        type(block_set_type), intent(inout) :: blk

        integer :: p, gp, v, e, pt, ni, nj
        integer :: di, dj, dk, var, nv, totalItems, sf
        integer :: b1, b2, b3, c1, c2, c3, s1, s2, s3
        real(C_DOUBLE) :: val
        logical :: doBlend

        nv = c%nActiveVars
        totalItems = merge(c%nLocalCopyPts, c%nLocalPts, c%copyOnly)*nv
        if (totalItems == 0) return
        sf = merge(1, 0, c%syncFace)

#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do &
        !$omp& map(to: totalItems, nv, sf, c%nLocal, c%lOff, c%lPointEntry, c%lSrcSlot, c%lDstSlot, &
        !$omp& c%lDstLo, c%lExt, c%lGA, c%lGB, c%lGS, c%lGC, c%lDir, c%lWp, c%lWpDst, &
        !$omp& c%lNrm, c%activeVars) &
        !$omp& map(tofrom: blk%q) &
        !$omp& private(p,gp,v,e,pt,ni,nj,di,dj,dk,var,b1,b2,b3,c1,c2,c3,s1,s2,s3,val,doBlend)
#endif
        do p = 1, totalItems
            gp = (p - 1)/nv
            v = p - 1 - gp*nv
            e = c%lPointEntry(gp)
            pt = gp - c%lOff(e-1)
            ni = c%lExt(1,e)
            nj = c%lExt(2,e)
            di = c%lDstLo(1,e) + modulo(pt, ni)
            dj = c%lDstLo(2,e) + modulo(pt/ni, nj)
            dk = c%lDstLo(3,e) + pt/(ni*nj)
            var = int(c%activeVars(v+1))
            b1 = ishft(c%lGA(1,e)*di + c%lGB(1,e), -c%lGS(1,e))
            b2 = ishft(c%lGA(2,e)*dj + c%lGB(2,e), -c%lGS(2,e))
            b3 = ishft(c%lGA(3,e)*dk + c%lGB(3,e), -c%lGS(3,e))
            c1 = merge(1, c%lGC(1,e), var == 1)
            c2 = merge(1, c%lGC(2,e), var == 2)
            c3 = merge(1, c%lGC(3,e), var == 3)
            val = 0.0d0
            do s3 = 0, c3 - 1
                do s2 = 0, c2 - 1
                    do s1 = 0, c1 - 1
                        val = val + blk%q(b1+s1, b2+s2, b3+s3, var, c%lSrcSlot(e))
                    end do
                end do
            end do
            val = val/real(c1*c2*c3, C_DOUBLE)
            ! 2:1 face ghost blend: place the prolonged value where the fine
            ! stencil expects it (ghost = (2*coarse + fine)/3 for uniform 2:1).
            ! Plain injection puts the coarse value at the fine halo centre, so
            ! BOTH the pressure gradient (projection) and the momentum predictor's
            ! advection / diffusion read it at the wrong location -> an O(1)
            ! interface truncation. Applies to the pressure AND the velocity
            ! components TANGENTIAL to this face (lDir(var)==0): a cell-centred
            ! halo whose blend geometry matches the pressure's. These halos feed
            ! only the momentum stencil, never the divergence, so conservation is
            ! untouched. The component NORMAL to the face (lDir(var)/=0) is
            ! face-staggered with a different geometry and is left to the
            ! face-staggered increment (its interface face is also fine-owned).
            ! The blend weight lives partly on the destination's first interior
            ! cell (entries write only halo cells, so that read is never a dst).
            doBlend = c%lWpDst(e) /= 0.0d0
            if (doBlend .and. var <= 3) doBlend = c%lDir(var,e) == 0
            if (doBlend) then
                val = c%lWp(e)*val + c%lWpDst(e) &
                    *blk%q(di-c%lDir(1,e), dj-c%lDir(2,e), dk-c%lDir(3,e), var, c%lDstSlot(e))
            end if
            ! Leave the owned interface normal velocity face untouched, unless
            ! this is the conservation SYNC (sf), which writes it.
            if (sf == 1 .or. var /= c%lNrm(e)) blk%q(di,dj,dk,var,c%lDstSlot(e)) = val
        end do
#ifdef USE_OPENMP_OFFLOAD
        !$omp end target teams distribute parallel do
#endif
    end subroutine copy_local_entries

    subroutine copy_local_scalar_entries(c, scalar)
        type(comm_type), intent(inout) :: c
        real(C_DOUBLE), intent(inout) :: scalar(0:,0:,0:,1:)

        integer :: p, e, pt, ni, nj
        integer :: di, dj, dk, totalItems, sfr, pn
        integer :: b1, b2, b3, c1, c2, c3, s1, s2, s3
        real(C_DOUBLE) :: val

        totalItems = c%nLocalPts
        if (totalItems == 0) return
        sfr = merge(1, 0, c%phiIfaceRow)

#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do &
        !$omp& map(to: totalItems, sfr, c%nLocal, c%lOff, c%lPointEntry, c%lSrcSlot, c%lDstSlot, &
        !$omp& c%lDstLo, c%lExt, c%lGA, c%lGB, c%lGS, c%lGC, c%lPhiN) &
        !$omp& map(tofrom: scalar) &
        !$omp& private(p,e,pt,ni,nj,di,dj,dk,b1,b2,b3,c1,c2,c3,s1,s2,s3,val,pn)
#endif
        do p = 1, totalItems
            e = c%lPointEntry(p - 1)
            pt = p - 1 - c%lOff(e-1)
            ni = c%lExt(1,e)
            nj = c%lExt(2,e)
            di = c%lDstLo(1,e) + modulo(pt, ni)
            dj = c%lDstLo(2,e) + modulo(pt/ni, nj)
            dk = c%lDstLo(3,e) + pt/(ni*nj)
            b1 = ishft(c%lGA(1,e)*di + c%lGB(1,e), -c%lGS(1,e))
            b2 = ishft(c%lGA(2,e)*dj + c%lGB(2,e), -c%lGS(2,e))
            b3 = ishft(c%lGA(3,e)*dk + c%lGB(3,e), -c%lGS(3,e))
            c1 = c%lGC(1,e)
            c2 = c%lGC(2,e)
            c3 = c%lGC(3,e)
            ! phi interface-row restrict: read only the fine row touching the
            ! face along the normal dim (base + 1 for off=-1, base for off=+1).
            pn = c%lPhiN(e)
            if (sfr == 1 .and. pn /= 0) then
                select case (abs(pn))
                case (1); b1 = b1 + merge(1, 0, pn < 0); c1 = 1
                case (2); b2 = b2 + merge(1, 0, pn < 0); c2 = 1
                case (3); b3 = b3 + merge(1, 0, pn < 0); c3 = 1
                end select
            end if
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
        integer :: di, dj, dk, var, peer, pos, nv, totalItems, copyOnly
        integer :: b1, b2, b3, c1, c2, c3, s1, s2, s3
        real(C_DOUBLE) :: val

        nv = c%nActiveVars
        totalItems = merge(c%peerSendCopyOff(c%nPeers), c%peerSendOff(c%nPeers), c%copyOnly)*nv
        if (totalItems == 0) return
        copyOnly = merge(1, 0, c%copyOnly)

#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do &
        !$omp& map(to: totalItems, nv, copyOnly, c%nPeers, c%nSend, c%sOff, c%sPointEntry, c%sSlot, c%sPeer, &
        !$omp& c%sDstLo, c%sExt, c%sGA, c%sGB, c%sGS, c%sGC, &
        !$omp& c%peerSendOff, c%peerSendCopyOff, c%activeVars, blk%q) &
        !$omp& map(tofrom: c%sendbuf) &
        !$omp& private(p,gp,v,e,pt,ni,nj,di,dj,dk,var,peer,pos,b1,b2,b3,c1,c2,c3,s1,s2,s3,val)
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
            e = c%sPointEntry(gp)
            pt = gp - c%sOff(e-1)
            ni = c%sExt(1,e)
            nj = c%sExt(2,e)
            di = c%sDstLo(1,e) + modulo(pt, ni)
            dj = c%sDstLo(2,e) + modulo(pt/ni, nj)
            dk = c%sDstLo(3,e) + pt/(ni*nj)
            var = int(c%activeVars(v+1))
            b1 = ishft(c%sGA(1,e)*di + c%sGB(1,e), -c%sGS(1,e))
            b2 = ishft(c%sGA(2,e)*dj + c%sGB(2,e), -c%sGS(2,e))
            b3 = ishft(c%sGA(3,e)*dk + c%sGB(3,e), -c%sGS(3,e))
            c1 = merge(1, c%sGC(1,e), var == 1)
            c2 = merge(1, c%sGC(2,e), var == 2)
            c3 = merge(1, c%sGC(3,e), var == 3)
            val = 0.0d0
            do s3 = 0, c3 - 1
                do s2 = 0, c2 - 1
                    do s1 = 0, c1 - 1
                        val = val + blk%q(b1+s1, b2+s2, b3+s3, var, c%sSlot(e))
                    end do
                end do
            end do
            peer = c%sPeer(e)
            pos = (gp - c%peerSendOff(peer-1))*nv + v + 1
            c%sendbuf(pos,peer) = val/real(c1*c2*c3, C_DOUBLE)
        end do
#ifdef USE_OPENMP_OFFLOAD
        !$omp end target teams distribute parallel do
#endif
    end subroutine pack_entries

    subroutine unpack_entries(c, blk)
        type(comm_type), intent(in) :: c
        type(block_set_type), intent(inout) :: blk

        integer :: p, gp, v, e, pt, ni, nj
        integer :: i, j, k, var, peer, pos, nv, totalItems, copyOnly, sf
        real(C_DOUBLE) :: val
        logical :: doBlend

        nv = c%nActiveVars
        totalItems = merge(c%peerRecvCopyOff(c%nPeers), c%peerRecvOff(c%nPeers), c%copyOnly)*nv
        if (totalItems == 0) return
        copyOnly = merge(1, 0, c%copyOnly)
        sf = merge(1, 0, c%syncFace)

#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do &
        !$omp& map(to: totalItems, nv, copyOnly, sf, c%nPeers, c%nRecv, c%rOff, c%rPointEntry, c%rSlot, c%rPeer, &
        !$omp& c%rLo, c%rExt, c%rDir, c%rWp, c%rWpDst, c%rNrm, &
        !$omp& c%peerRecvOff, c%peerRecvCopyOff, c%activeVars, c%recvbuf) &
        !$omp& map(tofrom: blk%q) &
        !$omp& private(p,gp,v,e,pt,ni,nj,i,j,k,var,peer,pos,val,doBlend)
#endif
        do p = 1, totalItems
            gp = (p - 1)/nv
            v = p - 1 - gp*nv
            if (copyOnly == 1) then
                peer = find_entry(c%peerRecvCopyOff, c%nPeers, gp)
                gp = c%peerRecvOff(peer-1) + (gp - c%peerRecvCopyOff(peer-1))
            end if
            e = c%rPointEntry(gp)
            pt = gp - c%rOff(e-1)
            ni = c%rExt(1,e)
            nj = c%rExt(2,e)
            i = c%rLo(1,e) + modulo(pt, ni)
            j = c%rLo(2,e) + modulo(pt/ni, nj)
            k = c%rLo(3,e) + pt/(ni*nj)
            var = int(c%activeVars(v+1))
            peer = c%rPeer(e)
            pos = (gp - c%peerRecvOff(peer-1))*nv + v + 1
            val = c%recvbuf(pos,peer)
            ! 2:1 face ghost blend (see copy_local_entries): pressure + the
            ! velocity components TANGENTIAL to this face (rDir(var)==0) placed at
            ! the fine halo centre; the normal component is left to the
            ! face-staggered increment.
            doBlend = c%rWpDst(e) /= 0.0d0
            if (doBlend .and. var <= 3) doBlend = c%rDir(var,e) == 0
            if (doBlend) then
                val = c%rWp(e)*val + c%rWpDst(e) &
                    *blk%q(i-c%rDir(1,e), j-c%rDir(2,e), k-c%rDir(3,e), var, c%rSlot(e))
            end if
            ! Leave the owned interface normal velocity face untouched, unless
            ! this is the conservation SYNC (sf), which writes it.
            if (sf == 1 .or. var /= c%rNrm(e)) blk%q(i,j,k,var,c%rSlot(e)) = val
        end do
#ifdef USE_OPENMP_OFFLOAD
        !$omp end target teams distribute parallel do
#endif
    end subroutine unpack_entries

    subroutine pack_scalar_entries(c, scalar)
        type(comm_type), intent(inout) :: c
        real(C_DOUBLE), intent(in) :: scalar(0:,0:,0:,1:)

        integer :: p, e, pt, ni, nj
        integer :: di, dj, dk, peer, pos, totalItems, sfr, pn
        integer :: b1, b2, b3, c1, c2, c3, s1, s2, s3
        real(C_DOUBLE) :: val

        totalItems = c%peerSendOff(c%nPeers)
        if (totalItems == 0) return
        sfr = merge(1, 0, c%phiIfaceRow)

#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do &
        !$omp& map(to: totalItems, sfr, c%nSend, c%sOff, c%sPointEntry, c%sSlot, c%sPeer, &
        !$omp& c%sDstLo, c%sExt, c%sGA, c%sGB, c%sGS, c%sGC, c%sPhiN, c%peerSendOff, scalar) &
        !$omp& map(tofrom: c%sendbuf) &
        !$omp& private(p,e,pt,ni,nj,di,dj,dk,peer,pos,b1,b2,b3,c1,c2,c3,s1,s2,s3,val,pn)
#endif
        do p = 1, totalItems
            e = c%sPointEntry(p - 1)
            pt = p - 1 - c%sOff(e-1)
            ni = c%sExt(1,e)
            nj = c%sExt(2,e)
            di = c%sDstLo(1,e) + modulo(pt, ni)
            dj = c%sDstLo(2,e) + modulo(pt/ni, nj)
            dk = c%sDstLo(3,e) + pt/(ni*nj)
            b1 = ishft(c%sGA(1,e)*di + c%sGB(1,e), -c%sGS(1,e))
            b2 = ishft(c%sGA(2,e)*dj + c%sGB(2,e), -c%sGS(2,e))
            b3 = ishft(c%sGA(3,e)*dk + c%sGB(3,e), -c%sGS(3,e))
            c1 = c%sGC(1,e)
            c2 = c%sGC(2,e)
            c3 = c%sGC(3,e)
            pn = c%sPhiN(e)
            if (sfr == 1 .and. pn /= 0) then
                select case (abs(pn))
                case (1); b1 = b1 + merge(1, 0, pn < 0); c1 = 1
                case (2); b2 = b2 + merge(1, 0, pn < 0); c2 = 1
                case (3); b3 = b3 + merge(1, 0, pn < 0); c3 = 1
                end select
            end if
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
        integer :: i, j, k, peer, pos, totalItems

        totalItems = c%peerRecvOff(c%nPeers)
        if (totalItems == 0) return

#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do &
        !$omp& map(to: totalItems, c%nRecv, c%rOff, c%rPointEntry, c%rSlot, c%rPeer, &
        !$omp& c%rLo, c%rExt, c%peerRecvOff, c%recvbuf) &
        !$omp& map(tofrom: scalar) &
        !$omp& private(p,e,pt,ni,nj,i,j,k,peer,pos)
#endif
        do p = 1, totalItems
            e = c%rPointEntry(p - 1)
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
