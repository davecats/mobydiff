!--------------------------!
!                          !
!     Block-structured     !
!       mesh module        !
!                          !
!--------------------------!
!
! Equal-size Cartesian blocks following the Building-Cube Method
! (Nakahashi & Kim, AIAA 2004-434) in the flat, tree-free formulation of
! Jansson et al. (IJHPCA 33(4), 2019). Design notes and the phased plan
! live in docs/block_refinement_strategy.md.
!
! Phase 0 (this file): a block_set_type holding ONE block per MPI rank
! whose metric content is bit-identical to the existing grid_type, plus
! field storage shaped like field_type with a trailing block index.
! Later phases add many blocks per rank, 2:1 level interfaces, and the
! removal of blocks buried inside the immersed boundary.
!
! GPU layout rationale: every per-block quantity is one flat contiguous
! array with the block slot as the LAST index. Kernels then add a single
! outer block loop (collapsed with k,j,i) and the device mapping stays a
! handful of `target enter data` clauses. There are no per-block
! allocations and no indirection inside stencils.

module blocks
    use, intrinsic :: iso_c_binding
    use, intrinsic :: iso_fortran_env, only: int64
    use :: init, only: dns_type, grid_type, NVAR, NVEL, VAR_U, VAR_V, VAR_W, slice_grid_direction
    implicit none

    private
    public :: block_set_type
    public :: init_block_set, destroy_block_set
    public :: enter_block_data, exit_block_data
    public :: subdivide_node_line
    public :: zorder_owner, zorder_start, zorder_count
    public :: zero_closed_halos, face_kind
    public :: DIST_RANKBOX, DIST_ZORDER
    public :: FACE_OPEN, FACE_PHYS, FACE_CLOSED

    ! Block ownership: one block per rank box (default), or the global
    ! Z-order lattice split linearly over the ranks ([blocks] nb).
    integer(C_INT), parameter :: DIST_RANKBOX = 0_C_INT
    integer(C_INT), parameter :: DIST_ZORDER  = 1_C_INT

    ! Per-block face kinds held in physLow/physHigh. Kernels treat any
    ! nonzero kind as a no-flux face; only FACE_PHYS faces receive
    ! boundary conditions, and FACE_CLOSED halos are zeroed once at init.
    integer(C_INT), parameter :: FACE_OPEN   = 0_C_INT
    integer(C_INT), parameter :: FACE_PHYS   = 1_C_INT
    integer(C_INT), parameter :: FACE_CLOSED = 2_C_INT

    type :: block_set_type
        ! All blocks of a set share the same cell count nb(1:3): the cubic
        ! [blocks] nb when configured, otherwise the rank-local box (one
        ! block per rank, the Phase-0 layout).
        integer(C_INT) :: nb(1:3) = 0_C_INT
        integer(C_INT) :: nBlocks = 0_C_INT
        integer(C_INT) :: nLevels = 1_C_INT

        ! Distribution of the global block table over the ranks.
        !   DIST_RANKBOX: one block per rank box; globalId == rank.
        !   DIST_ZORDER:  global lattice gnbt = globalSize/nb numbered along
        !                 a Z-order (Morton) curve; rank p owns the
        !                 consecutive ids [zorder_start(p), +zorder_count(p)).
        ! Local slot s holds global id idStart + s - 1 in both modes.
        integer(C_INT) :: distMode = DIST_RANKBOX
        integer(C_INT) :: gnbt(1:3) = 0_C_INT
        integer(C_INT) :: nBlocksGlobal = 0_C_INT
        integer(C_INT) :: idStart = 0_C_INT
        ! Z-order lookup tables, identical on every rank (host only).
        integer(C_INT), allocatable :: zidOf(:,:,:)  ! lattice coords -> id
        integer(C_INT), allocatable :: zcoord(:,:)   ! (3, id+1) -> lattice coords

        ! Per-block metadata (small, host + device).
        !   level:    refinement level; 0 = base grid, +1 per cell bisection
        !   origin:   zero-based cell origin of the block in the level-l
        !             global index space
        !   globalId: position in the global (Z-ordered) block table
        !   physLow/physHigh: per-direction face kind (FACE_OPEN /
        !             FACE_PHYS on a non-periodic global boundary /
        !             FACE_CLOSED towards a removed solid block); the
        !             per-block face masks used by the momentum starts,
        !             the SOR sweep and the boundary-condition point list
        integer(C_INT), allocatable :: level(:)        ! (nBlocks)
        integer(C_INT), allocatable :: origin(:,:)     ! (3,nBlocks)
        integer(C_INT), allocatable :: globalId(:)     ! (nBlocks)
        integer(C_INT), allocatable :: physLow(:,:)    ! (3,nBlocks)
        integer(C_INT), allocatable :: physHigh(:,:)   ! (3,nBlocks)

        ! Per-block staggered coordinates and finite-difference metrics,
        ! sliced from the (level-specific) global node lines. Shapes mirror
        ! grid_type with one trailing block index.
        real(C_DOUBLE), allocatable :: x(:,:,:), y(:,:,:), z(:,:,:)       ! (-1:nb+2,NVAR,nBlocks)
        real(C_DOUBLE), allocatable :: d1x(:,:,:), d1y(:,:,:), d1z(:,:,:) ! (0:nb+1,NVAR,nBlocks)
        real(C_DOUBLE), allocatable :: lapXm(:,:,:), lapX0(:,:,:), lapXp(:,:,:)
        real(C_DOUBLE), allocatable :: lapYm(:,:,:), lapY0(:,:,:), lapYp(:,:,:)
        real(C_DOUBLE), allocatable :: lapZm(:,:,:), lapZ0(:,:,:), lapZp(:,:,:)

        ! Flow state with one halo cell per side (second-order stencils).
        real(C_DOUBLE), allocatable :: q(:,:,:,:,:)      ! (0:nb+1,...,NVAR,nBlocks)
        real(C_DOUBLE), allocatable :: qs(:,:,:,:,:)     ! (0:nb+1,...,NVEL,nBlocks)
        real(C_DOUBLE), allocatable :: oldrhs(:,:,:,:,:) ! (1:nb,...,NVEL,nBlocks)
    end type block_set_type

contains

    ! Build this rank's blocks. With [blocks] nb the global grid is a
    ! uniform block lattice numbered along a Z-order curve and split
    ! linearly over the ranks (each owns floor((N+P-p-1)/P) consecutive
    ! ids); without it every rank owns its Cartesian box as one block.
    ! Coordinates and metrics are sliced from the global node lines via
    ! slice_grid_direction either way.
    subroutine init_block_set(blk, dns, g, periodic, nranks, myrank, active)
        type(block_set_type), intent(inout) :: blk
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        logical(C_BOOL), intent(in) :: periodic(1:3)
        integer(C_INT), intent(in) :: nranks, myrank
        ! Optional per-lattice-cell keep flags (x-fastest raster order, 1 =
        ! keep). Blocks buried inside the immersed boundary are dropped from
        ! the global table; only meaningful with [blocks] nb.
        integer(C_INT), intent(in), optional :: active(:)

        integer :: nx, ny, nz, b, d, id, nRemoved

        call destroy_block_set(blk)

        if (dns%block_nb > 0_C_INT) then
            blk%distMode = DIST_ZORDER
            blk%nb = dns%block_nb
            do d = 1, 3
                if (mod(dns%globalSize(d), dns%block_nb) /= 0_C_INT) then
                    print *, "block size nb =", dns%block_nb, "does not divide the grid", &
                        dns%globalSize(1:3), "in direction", d
                    error stop "[blocks] nb must divide the global grid in every direction"
                end if
            end do
            blk%gnbt = dns%globalSize(1:3)/blk%nb
            call build_zorder_table(blk, active)
            nRemoved = int(product(blk%gnbt)) - int(blk%nBlocksGlobal)
            if (nRemoved > 0 .and. myrank == 0_C_INT) then
                print *, "removed", nRemoved, "of", product(blk%gnbt), &
                    "blocks buried inside the immersed boundary"
            end if
            blk%idStart = zorder_start(blk%nBlocksGlobal, nranks, myrank)
            blk%nBlocks = zorder_count(blk%nBlocksGlobal, nranks, myrank)
        else
            blk%distMode = DIST_RANKBOX
            blk%nb = dns%localSize(1:3,2)
            blk%nBlocksGlobal = nranks
            blk%idStart = myrank
            blk%nBlocks = 1_C_INT
        end if
        blk%nLevels = 1_C_INT
        if (blk%nBlocks < 1_C_INT) error stop "rank owns no blocks; use fewer ranks or smaller nb"

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))

        allocate(blk%level(blk%nBlocks))
        allocate(blk%origin(3,blk%nBlocks))
        allocate(blk%globalId(blk%nBlocks))
        allocate(blk%physLow(3,blk%nBlocks))
        allocate(blk%physHigh(3,blk%nBlocks))
        blk%level = 0_C_INT

        do b = 1, int(blk%nBlocks)
            id = int(blk%idStart) + b - 1
            blk%globalId(b) = int(id, C_INT)
            if (blk%distMode == DIST_ZORDER) then
                blk%origin(:,b) = blk%zcoord(:,id+1)*blk%nb
            else
                blk%origin(:,b) = dns%localSize(1:3,0) - 1_C_INT
            end if
            do d = 1, 3
                blk%physLow(d,b) = face_kind(blk, dns, periodic, blk%origin(:,b), d, -1)
                blk%physHigh(d,b) = face_kind(blk, dns, periodic, blk%origin(:,b), d, +1)
            end do
        end do

        allocate(blk%x(-1:nx+2,NVAR,blk%nBlocks), blk%d1x(0:nx+1,NVAR,blk%nBlocks))
        allocate(blk%y(-1:ny+2,NVAR,blk%nBlocks), blk%d1y(0:ny+1,NVAR,blk%nBlocks))
        allocate(blk%z(-1:nz+2,NVAR,blk%nBlocks), blk%d1z(0:nz+1,NVAR,blk%nBlocks))
        allocate(blk%lapXm(0:nx+1,NVAR,blk%nBlocks), blk%lapX0(0:nx+1,NVAR,blk%nBlocks), &
                 blk%lapXp(0:nx+1,NVAR,blk%nBlocks))
        allocate(blk%lapYm(0:ny+1,NVAR,blk%nBlocks), blk%lapY0(0:ny+1,NVAR,blk%nBlocks), &
                 blk%lapYp(0:ny+1,NVAR,blk%nBlocks))
        allocate(blk%lapZm(0:nz+1,NVAR,blk%nBlocks), blk%lapZ0(0:nz+1,NVAR,blk%nBlocks), &
                 blk%lapZp(0:nz+1,NVAR,blk%nBlocks))

        ! Slice per-block metrics from the global node lines. From Phase 3 on
        ! the line passed here is the block-level line obtained by midpoint
        ! subdivision (subdivide_node_line), which preserves stretched base
        ! grids and makes a coarse cell the exact union of its children.
        do b = 1, int(blk%nBlocks)
            call slice_grid_direction(g%xNode, blk%x(:,:,b), blk%d1x(:,:,b), &
                blk%lapXm(:,:,b), blk%lapX0(:,:,b), blk%lapXp(:,:,b), &
                dns%globalSize(1), blk%origin(1,b) + 1_C_INT, nx, dns%leng(1), periodic(1), 1)
            call slice_grid_direction(g%yNode, blk%y(:,:,b), blk%d1y(:,:,b), &
                blk%lapYm(:,:,b), blk%lapY0(:,:,b), blk%lapYp(:,:,b), &
                dns%globalSize(2), blk%origin(2,b) + 1_C_INT, ny, dns%leng(2), periodic(2), 2)
            call slice_grid_direction(g%zNode, blk%z(:,:,b), blk%d1z(:,:,b), &
                blk%lapZm(:,:,b), blk%lapZ0(:,:,b), blk%lapZp(:,:,b), &
                dns%globalSize(3), blk%origin(3,b) + 1_C_INT, nz, dns%leng(3), periodic(3), 3)
        end do

        allocate(blk%q(0:nx+1,0:ny+1,0:nz+1,NVAR,blk%nBlocks))
        allocate(blk%qs(0:nx+1,0:ny+1,0:nz+1,NVEL,blk%nBlocks))
        allocate(blk%oldrhs(1:nx,1:ny,1:nz,NVEL,blk%nBlocks))
        blk%q = 0.0d0
        blk%qs = 0.0d0
        blk%oldrhs = 0.0d0
    end subroutine init_block_set

    subroutine destroy_block_set(blk)
        type(block_set_type), intent(inout) :: blk

        if (allocated(blk%level)) deallocate(blk%level)
        if (allocated(blk%origin)) deallocate(blk%origin)
        if (allocated(blk%globalId)) deallocate(blk%globalId)
        if (allocated(blk%physLow)) deallocate(blk%physLow)
        if (allocated(blk%physHigh)) deallocate(blk%physHigh)
        if (allocated(blk%x)) deallocate(blk%x)
        if (allocated(blk%y)) deallocate(blk%y)
        if (allocated(blk%z)) deallocate(blk%z)
        if (allocated(blk%d1x)) deallocate(blk%d1x)
        if (allocated(blk%d1y)) deallocate(blk%d1y)
        if (allocated(blk%d1z)) deallocate(blk%d1z)
        if (allocated(blk%lapXm)) deallocate(blk%lapXm)
        if (allocated(blk%lapX0)) deallocate(blk%lapX0)
        if (allocated(blk%lapXp)) deallocate(blk%lapXp)
        if (allocated(blk%lapYm)) deallocate(blk%lapYm)
        if (allocated(blk%lapY0)) deallocate(blk%lapY0)
        if (allocated(blk%lapYp)) deallocate(blk%lapYp)
        if (allocated(blk%lapZm)) deallocate(blk%lapZm)
        if (allocated(blk%lapZ0)) deallocate(blk%lapZ0)
        if (allocated(blk%lapZp)) deallocate(blk%lapZp)
        if (allocated(blk%q)) deallocate(blk%q)
        if (allocated(blk%qs)) deallocate(blk%qs)
        if (allocated(blk%oldrhs)) deallocate(blk%oldrhs)

        if (allocated(blk%zidOf)) deallocate(blk%zidOf)
        if (allocated(blk%zcoord)) deallocate(blk%zcoord)

        blk%nb = 0_C_INT
        blk%gnbt = 0_C_INT
        blk%nBlocks = 0_C_INT
        blk%nBlocksGlobal = 0_C_INT
        blk%idStart = 0_C_INT
        blk%distMode = DIST_RANKBOX
        blk%nLevels = 1_C_INT
    end subroutine destroy_block_set

    ! Device mapping, mirroring enter_grid_data/enter_field_data in
    ! gpu_runtime: map the container once, then its flat member arrays.
    subroutine enter_block_data(blk)
        type(block_set_type), intent(inout) :: blk

#ifdef USE_OPENMP_OFFLOAD
        !$omp target enter data map(to: blk)
        !$omp target enter data map(to: &
        !$omp& blk%origin, blk%physLow, blk%physHigh, &
        !$omp& blk%x, blk%y, blk%z, blk%d1x, blk%d1y, blk%d1z, &
        !$omp& blk%lapXm, blk%lapX0, blk%lapXp, &
        !$omp& blk%lapYm, blk%lapY0, blk%lapYp, &
        !$omp& blk%lapZm, blk%lapZ0, blk%lapZp, &
        !$omp& blk%q, blk%qs, blk%oldrhs)
#endif
    end subroutine enter_block_data

    subroutine exit_block_data(blk)
        type(block_set_type), intent(inout) :: blk

#ifdef USE_OPENMP_OFFLOAD
        !$omp target exit data map(delete: &
        !$omp& blk%origin, blk%physLow, blk%physHigh, &
        !$omp& blk%x, blk%y, blk%z, blk%d1x, blk%d1y, blk%d1z, &
        !$omp& blk%lapXm, blk%lapX0, blk%lapXp, &
        !$omp& blk%lapYm, blk%lapY0, blk%lapYp, &
        !$omp& blk%lapZm, blk%lapZ0, blk%lapZp, &
        !$omp& blk%q, blk%qs, blk%oldrhs)
        !$omp target exit data map(delete: blk)
#endif
    end subroutine exit_block_data

    ! Number the global block lattice along a Z-order (Morton) curve:
    ! sort lattice cells by their bit-interleaved coordinates, then assign
    ! compact ids to the surviving (active) blocks only; removed lattice
    ! cells get zidOf = -1. Identical on every rank, so the distribution
    ! needs no communication.
    subroutine build_zorder_table(blk, active)
        type(block_set_type), intent(inout) :: blk
        integer(C_INT), intent(in), optional :: active(:)

        integer :: n, gx, gy, gz, i, id, raster
        integer(int64), allocatable :: keys(:)
        integer, allocatable :: order(:)
        integer(C_INT), allocatable :: coordTmp(:,:)
        logical :: keep

        n = int(product(blk%gnbt))
        if (present(active)) then
            if (size(active) /= n) error stop "block active mask does not match the lattice"
        end if
        allocate(blk%zidOf(0:int(blk%gnbt(1))-1, 0:int(blk%gnbt(2))-1, 0:int(blk%gnbt(3))-1))
        allocate(coordTmp(3,n))
        allocate(keys(n), order(n))
        blk%zidOf = -1_C_INT

        i = 0
        do gz = 0, int(blk%gnbt(3)) - 1
            do gy = 0, int(blk%gnbt(2)) - 1
                do gx = 0, int(blk%gnbt(1)) - 1
                    i = i + 1
                    keys(i) = morton_key(gx, gy, gz)
                    coordTmp(:,i) = int([gx, gy, gz], C_INT)
                end do
            end do
        end do

        call heapsort_index(keys, order)

        ! Walk the curve and hand out ids to survivors; zcoord(:,id+1) are
        ! the lattice coords of Z-id `id`.
        allocate(blk%zcoord(3,n))
        id = 0
        do i = 1, n
            keep = .true.
            if (present(active)) then
                raster = 1 + int(coordTmp(1,order(i))) &
                    + int(blk%gnbt(1))*(int(coordTmp(2,order(i))) &
                    + int(blk%gnbt(2))*int(coordTmp(3,order(i))))
                keep = active(raster) /= 0_C_INT
            end if
            if (.not. keep) cycle
            id = id + 1
            blk%zcoord(:,id) = coordTmp(:,order(i))
            blk%zidOf(coordTmp(1,order(i)), coordTmp(2,order(i)), coordTmp(3,order(i))) = &
                int(id - 1, C_INT)
        end do
        blk%nBlocksGlobal = int(id, C_INT)

        deallocate(keys, order, coordTmp)
    end subroutine build_zorder_table

    ! Face kind of the block at `origin` in direction d, side -1/+1: a
    ! physical wall outside non-periodic boundaries, closed towards a
    ! removed block, open otherwise.
    integer(C_INT) function face_kind(blk, dns, periodic, origin, d, side) result(fk)
        type(block_set_type), intent(in) :: blk
        type(dns_type), intent(in) :: dns
        logical(C_BOOL), intent(in) :: periodic(1:3)
        integer(C_INT), intent(in) :: origin(3)
        integer, intent(in) :: d, side

        integer :: to(3)

        to = int(origin)
        to(d) = to(d) + side*int(blk%nb(d))
        if (to(d) < 0 .or. to(d) >= int(dns%globalSize(d))) then
            if (.not. periodic(d)) then
                fk = FACE_PHYS
                return
            end if
            to(d) = modulo(to(d), int(dns%globalSize(d)))
        end if
        fk = FACE_OPEN
        if (blk%distMode == DIST_ZORDER) then
            if (blk%zidOf(to(1)/int(blk%nb(1)), to(2)/int(blk%nb(2)), &
                          to(3)/int(blk%nb(3))) < 0_C_INT) fk = FACE_CLOSED
        end if
    end function face_kind

    ! FACE_CLOSED faces are exact zero-flux faces: the halo layer and the
    ! pinned interface velocity are zeroed once here and never written
    ! again (momentum skips the face, the sweep corrections are masked,
    ! and no exchange entries point at removed blocks).
    subroutine zero_closed_halos(blk)
        type(block_set_type), intent(inout) :: blk

        integer :: b, nx, ny, nz, nClosed

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))

        nClosed = count(blk%physLow == FACE_CLOSED) + count(blk%physHigh == FACE_CLOSED)
        if (nClosed == 0) return

        do b = 1, int(blk%nBlocks)
            if (blk%physLow(1,b) == FACE_CLOSED) then
                blk%q(0,:,:,:,b) = 0.0d0
                blk%q(1,:,:,VAR_U,b) = 0.0d0
            end if
            if (blk%physHigh(1,b) == FACE_CLOSED) blk%q(nx+1,:,:,:,b) = 0.0d0
            if (blk%physLow(2,b) == FACE_CLOSED) then
                blk%q(:,0,:,:,b) = 0.0d0
                blk%q(:,1,:,VAR_V,b) = 0.0d0
            end if
            if (blk%physHigh(2,b) == FACE_CLOSED) blk%q(:,ny+1,:,:,b) = 0.0d0
            if (blk%physLow(3,b) == FACE_CLOSED) then
                blk%q(:,:,0,:,b) = 0.0d0
                blk%q(:,:,1,VAR_W,b) = 0.0d0
            end if
            if (blk%physHigh(3,b) == FACE_CLOSED) blk%q(:,:,nz+1,:,b) = 0.0d0
        end do

#ifdef USE_OPENMP_OFFLOAD
        !$omp target update to(blk%q)
#endif
    end subroutine zero_closed_halos

    pure integer(int64) function morton_key(gx, gy, gz) result(key)
        integer, intent(in) :: gx, gy, gz

        integer :: bit

        key = 0_int64
        do bit = 0, 20
            key = ior(key, ishft(iand(int(gx, int64), ishft(1_int64, bit)), 2*bit))
            key = ior(key, ishft(iand(int(gy, int64), ishft(1_int64, bit)), 2*bit + 1))
            key = ior(key, ishft(iand(int(gz, int64), ishft(1_int64, bit)), 2*bit + 2))
        end do
    end function morton_key

    subroutine heapsort_index(keys, order)
        integer(int64), intent(in) :: keys(:)
        integer, intent(out) :: order(:)

        integer :: n, i, tmp

        n = size(keys)
        do i = 1, n
            order(i) = i
        end do

        do i = n/2, 1, -1
            call sift_down(keys, order, i, n)
        end do
        do i = n, 2, -1
            tmp = order(1)
            order(1) = order(i)
            order(i) = tmp
            call sift_down(keys, order, 1, i - 1)
        end do
    end subroutine heapsort_index

    subroutine sift_down(keys, order, start, last)
        integer(int64), intent(in) :: keys(:)
        integer, intent(inout) :: order(:)
        integer, intent(in) :: start, last

        integer :: root, child, tmp

        root = start
        do while (2*root <= last)
            child = 2*root
            if (child < last) then
                if (keys(order(child)) < keys(order(child+1))) child = child + 1
            end if
            if (keys(order(root)) >= keys(order(child))) return
            tmp = order(root)
            order(root) = order(child)
            order(child) = tmp
            root = child
        end do
    end subroutine sift_down

    ! Linear distribution of N Z-ordered blocks over P ranks: rank p owns
    ! floor((N + P - p - 1)/P) consecutive ids.
    pure integer(C_INT) function zorder_count(n, p_total, p) result(cnt)
        integer(C_INT), intent(in) :: n, p_total, p

        cnt = int((int(n) + int(p_total) - int(p) - 1)/int(p_total), C_INT)
    end function zorder_count

    pure integer(C_INT) function zorder_start(n, p_total, p) result(start)
        integer(C_INT), intent(in) :: n, p_total, p

        integer :: q, r

        q = int(n)/int(p_total)
        r = mod(int(n), int(p_total))
        start = int(int(p)*q + min(int(p), r), C_INT)
    end function zorder_start

    pure integer(C_INT) function zorder_owner(id, n, p_total) result(owner)
        integer(C_INT), intent(in) :: id, n, p_total

        integer :: q, r, split

        q = int(n)/int(p_total)
        r = mod(int(n), int(p_total))
        split = r*(q + 1)
        if (int(id) < split) then
            owner = int(int(id)/(q + 1), C_INT)
        else
            owner = int(r + (int(id) - split)/max(q, 1), C_INT)
        end if
    end function zorder_owner

    ! Midpoint subdivision of a node line: each coarse cell is split into two
    ! halves, so a level-l cell is exactly the union of its children. This is
    ! what makes fine-to-coarse face averaging conservative on stretched base
    ! grids (used from Phase 3 on).
    subroutine subdivide_node_line(coarse, fine)
        real(C_DOUBLE), intent(in) :: coarse(0:)
        real(C_DOUBLE), intent(out) :: fine(0:)

        integer :: i, nc

        nc = ubound(coarse, 1)
        if (ubound(fine, 1) /= 2*nc) error stop "subdivide_node_line: fine line must have 2*n cells"

        do i = 0, nc - 1
            fine(2*i) = coarse(i)
            fine(2*i+1) = 0.5d0*(coarse(i) + coarse(i+1))
        end do
        fine(2*nc) = coarse(nc)
    end subroutine subdivide_node_line

end module blocks
