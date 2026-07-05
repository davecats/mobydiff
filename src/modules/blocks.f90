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
    public :: zero_closed_halos, face_kind, leaf_at, level_cells
    public :: level_cell_width, occupied_any_level
    public :: DIST_RANKBOX, DIST_ZORDER
    public :: FACE_OPEN, FACE_PHYS, FACE_CLOSED, FACE_COARSE, FACE_FINE

    ! Block ownership: one block per rank box (default), or the global
    ! Z-order lattice split linearly over the ranks ([blocks] nb).
    integer(C_INT), parameter :: DIST_RANKBOX = 0_C_INT
    integer(C_INT), parameter :: DIST_ZORDER  = 1_C_INT

    ! Per-block face kinds held in physLow/physHigh. The momentum starts
    ! and the SOR sweep treat PHYS and CLOSED as no-flux faces; only
    ! FACE_PHYS faces receive boundary conditions; FACE_CLOSED halos are
    ! zeroed once at init; COARSE/FINE mark 2:1 level interfaces.
    integer(C_INT), parameter :: FACE_OPEN   = 0_C_INT
    integer(C_INT), parameter :: FACE_PHYS   = 1_C_INT
    integer(C_INT), parameter :: FACE_CLOSED = 2_C_INT
    integer(C_INT), parameter :: FACE_COARSE = 3_C_INT
    integer(C_INT), parameter :: FACE_FINE   = 4_C_INT

    type :: block_set_type
        ! All blocks of a set share the same cell count nb(1:3): the cubic
        ! [blocks] nb when configured, otherwise the rank-local box (one
        ! block per rank, the Phase-0 layout).
        integer(C_INT) :: nb(1:3) = 0_C_INT
        integer(C_INT) :: nBlocks = 0_C_INT
        integer(C_INT) :: nLevels = 1_C_INT

        ! Distribution of the global block table over the ranks.
        !   DIST_RANKBOX: one block per rank box; globalId == rank.
        !   DIST_ZORDER:  global lattice nTiles = globalSize/nb numbered along
        !                 a Z-order (Morton) curve; rank p owns the
        !                 consecutive ids [zorder_start(p), +zorder_count(p)).
        ! Local slot s holds global id idStart + s - 1 in both modes.
        integer(C_INT) :: distMode = DIST_RANKBOX
        integer(C_INT) :: nTiles(1:3) = 0_C_INT        ! root (level-0) lattice
        integer(C_INT) :: nBlocksGlobal = 0_C_INT
        integer(C_INT) :: idStart = 0_C_INT

        ! Global leaf table, identical on every rank (host only). Leaves of
        ! the refinement forest are numbered along the Z-order curve of the
        ! finest lattice (coord << (nLevels-1-level) per dimension).
        integer(C_INT), allocatable :: leafLevel(:)  ! (nBlocksGlobal), id+1
        integer(C_INT), allocatable :: leafCoord(:,:)! (3, id+1) level-l lattice coords
        ! Per-level lattice -> leaf id (or -1), x-fastest, level l starting
        ! at levelOffset(l)+1.
        integer(C_INT), allocatable :: lidOf(:)
        integer(C_INT), allocatable :: levelOffset(:)     ! (0:nLevels)

        ! Per-level global node lines (level l in column l+1; the level-l
        ! line has globalSize*2^l + 1 nodes, allocated at the finest length).
        real(C_DOUBLE), allocatable :: lineX(:,:), lineY(:,:), lineZ(:,:)

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
    subroutine init_block_set(blk, dns, g, periodic, nranks, myrank, active, touch, buried)
        type(block_set_type), intent(inout) :: blk
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        logical(C_BOOL), intent(in) :: periodic(1:3)
        integer(C_INT), intent(in) :: nranks, myrank
        ! Optional per-lattice-cell keep flags (x-fastest raster order, 1 =
        ! keep). Blocks buried inside the immersed boundary are dropped from
        ! the global table; only meaningful with [blocks] nb.
        integer(C_INT), intent(in), optional :: active(:)
        ! Optional per-level geometry masks (raster, level) for
        ! geometry-driven refinement: touch = dilated block straddles the
        ! immersed surface, buried = dilated block fully solid.
        integer(C_INT), intent(in), optional :: touch(:,:), buried(:,:)

        integer :: d, nRemoved

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
            blk%nTiles = dns%globalSize(1:3)/blk%nb
            blk%nLevels = 1_C_INT
            if (refine_box_set(dns) .or. dns%block_refine_body) then
                blk%nLevels = 1_C_INT + max(dns%block_refine_levels, 0_C_INT)
            end if
            call build_level_lines(blk, dns, g)
            call build_leaf_table(blk, dns, periodic, active, myrank, touch, buried)
            nRemoved = int(product(blk%nTiles)) - count_level0_leaves(blk)
            if (nRemoved > 0 .and. myrank == 0_C_INT .and. blk%nLevels == 1_C_INT) then
                print *, "removed", nRemoved, "of", product(blk%nTiles), &
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
            blk%nLevels = 1_C_INT
        end if
        if (blk%nBlocks < 1_C_INT) error stop "rank owns no blocks; use fewer ranks or smaller nb"

        call build_block_metadata(blk, dns, periodic)
        call build_block_metrics(blk, dns, g, periodic)
        call allocate_flow_state(blk)
    end subroutine init_block_set

    ! Allocate and fill this rank's per-block metadata: global id, refinement
    ! level, cell origin (level-l index space) and the six face kinds.
    subroutine build_block_metadata(blk, dns, periodic)
        type(block_set_type), intent(inout) :: blk
        type(dns_type), intent(in) :: dns
        logical(C_BOOL), intent(in) :: periodic(1:3)

        integer :: b, d, id

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
                blk%level(b) = blk%leafLevel(id+1)
                blk%origin(:,b) = blk%leafCoord(:,id+1)*blk%nb
            else
                blk%origin(:,b) = dns%localSize(1:3,0) - 1_C_INT
            end if
            do d = 1, 3
                blk%physLow(d,b) = face_kind(blk, dns, periodic, blk%origin(:,b), &
                    blk%level(b), d, -1)
                blk%physHigh(d,b) = face_kind(blk, dns, periodic, blk%origin(:,b), &
                    blk%level(b), d, +1)
            end do
        end do
    end subroutine build_block_metadata

    ! Allocate this rank's per-block staggered coordinates and finite-difference
    ! metrics, then slice them from the (level-specific) global node lines:
    ! level l is the level-0 line after l midpoint subdivisions, so a coarse
    ! cell is the exact union of its children even on stretched grids. lcol is
    ! the node-line column for the block's level (level 0 in column 1).
    subroutine build_block_metrics(blk, dns, g, periodic)
        type(block_set_type), intent(inout) :: blk
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        logical(C_BOOL), intent(in) :: periodic(1:3)

        integer :: b, nx, ny, nz, lcol

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))

        allocate(blk%x(-1:nx+2,NVAR,blk%nBlocks), blk%d1x(0:nx+1,NVAR,blk%nBlocks))
        allocate(blk%y(-1:ny+2,NVAR,blk%nBlocks), blk%d1y(0:ny+1,NVAR,blk%nBlocks))
        allocate(blk%z(-1:nz+2,NVAR,blk%nBlocks), blk%d1z(0:nz+1,NVAR,blk%nBlocks))
        allocate(blk%lapXm(0:nx+1,NVAR,blk%nBlocks), blk%lapX0(0:nx+1,NVAR,blk%nBlocks), &
                 blk%lapXp(0:nx+1,NVAR,blk%nBlocks))
        allocate(blk%lapYm(0:ny+1,NVAR,blk%nBlocks), blk%lapY0(0:ny+1,NVAR,blk%nBlocks), &
                 blk%lapYp(0:ny+1,NVAR,blk%nBlocks))
        allocate(blk%lapZm(0:nz+1,NVAR,blk%nBlocks), blk%lapZ0(0:nz+1,NVAR,blk%nBlocks), &
                 blk%lapZp(0:nz+1,NVAR,blk%nBlocks))

        do b = 1, int(blk%nBlocks)
            if (blk%distMode == DIST_ZORDER) then
                lcol = int(blk%level(b)) + 1
                call slice_grid_direction(blk%lineX(:,lcol), blk%x(:,:,b), blk%d1x(:,:,b), &
                    blk%lapXm(:,:,b), blk%lapX0(:,:,b), blk%lapXp(:,:,b), &
                    level_cells(dns, 1, blk%level(b)), blk%origin(1,b) + 1_C_INT, nx, &
                    dns%leng(1), periodic(1), 1)
                call slice_grid_direction(blk%lineY(:,lcol), blk%y(:,:,b), blk%d1y(:,:,b), &
                    blk%lapYm(:,:,b), blk%lapY0(:,:,b), blk%lapYp(:,:,b), &
                    level_cells(dns, 2, blk%level(b)), blk%origin(2,b) + 1_C_INT, ny, &
                    dns%leng(2), periodic(2), 2)
                call slice_grid_direction(blk%lineZ(:,lcol), blk%z(:,:,b), blk%d1z(:,:,b), &
                    blk%lapZm(:,:,b), blk%lapZ0(:,:,b), blk%lapZp(:,:,b), &
                    level_cells(dns, 3, blk%level(b)), blk%origin(3,b) + 1_C_INT, nz, &
                    dns%leng(3), periodic(3), 3)
            else
                call slice_grid_direction(g%xNode, blk%x(:,:,b), blk%d1x(:,:,b), &
                    blk%lapXm(:,:,b), blk%lapX0(:,:,b), blk%lapXp(:,:,b), &
                    dns%globalSize(1), blk%origin(1,b) + 1_C_INT, nx, dns%leng(1), periodic(1), 1)
                call slice_grid_direction(g%yNode, blk%y(:,:,b), blk%d1y(:,:,b), &
                    blk%lapYm(:,:,b), blk%lapY0(:,:,b), blk%lapYp(:,:,b), &
                    dns%globalSize(2), blk%origin(2,b) + 1_C_INT, ny, dns%leng(2), periodic(2), 2)
                call slice_grid_direction(g%zNode, blk%z(:,:,b), blk%d1z(:,:,b), &
                    blk%lapZm(:,:,b), blk%lapZ0(:,:,b), blk%lapZp(:,:,b), &
                    dns%globalSize(3), blk%origin(3,b) + 1_C_INT, nz, dns%leng(3), periodic(3), 3)
            end if
        end do
    end subroutine build_block_metrics

    ! Allocate the flow state (one halo cell per side) and zero it.
    subroutine allocate_flow_state(blk)
        type(block_set_type), intent(inout) :: blk

        integer :: nx, ny, nz

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))

        allocate(blk%q(0:nx+1,0:ny+1,0:nz+1,NVAR,blk%nBlocks))
        allocate(blk%qs(0:nx+1,0:ny+1,0:nz+1,NVEL,blk%nBlocks))
        allocate(blk%oldrhs(1:nx,1:ny,1:nz,NVEL,blk%nBlocks))
        blk%q = 0.0d0
        blk%qs = 0.0d0
        blk%oldrhs = 0.0d0
    end subroutine allocate_flow_state

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

        if (allocated(blk%leafLevel)) deallocate(blk%leafLevel)
        if (allocated(blk%leafCoord)) deallocate(blk%leafCoord)
        if (allocated(blk%lidOf)) deallocate(blk%lidOf)
        if (allocated(blk%levelOffset)) deallocate(blk%levelOffset)
        if (allocated(blk%lineX)) deallocate(blk%lineX)
        if (allocated(blk%lineY)) deallocate(blk%lineY)
        if (allocated(blk%lineZ)) deallocate(blk%lineZ)

        blk%nb = 0_C_INT
        blk%nTiles = 0_C_INT
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

    logical function refine_box_set(dns)
        type(dns_type), intent(in) :: dns

        refine_box_set = dns%block_refine_nboxes > 0_C_INT
    end function refine_box_set

    ! Cells of the level-l global grid in direction d.
    pure integer(C_INT) function level_cells(dns, d, level) result(n)
        type(dns_type), intent(in) :: dns
        integer, intent(in) :: d
        integer(C_INT), intent(in) :: level

        n = dns%globalSize(d)*int(2**int(level), C_INT)
    end function level_cells

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

    ! Per-level node lines: level 0 is the configured grid, level l+1 the
    ! midpoint subdivision of level l (subdivide_node_line).
    subroutine build_level_lines(blk, dns, g)
        type(block_set_type), intent(inout) :: blk
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g

        integer :: l, nf(3)

        nf = int(dns%globalSize)*2**(int(blk%nLevels) - 1)
        allocate(blk%lineX(0:nf(1), blk%nLevels))
        allocate(blk%lineY(0:nf(2), blk%nLevels))
        allocate(blk%lineZ(0:nf(3), blk%nLevels))
        blk%lineX(0:int(dns%globalSize(1)),1) = g%xNode
        blk%lineY(0:int(dns%globalSize(2)),1) = g%yNode
        blk%lineZ(0:int(dns%globalSize(3)),1) = g%zNode
        do l = 2, int(blk%nLevels)
            call subdivide_node_line(blk%lineX(0:int(dns%globalSize(1))*2**(l-2), l-1), &
                                     blk%lineX(0:int(dns%globalSize(1))*2**(l-1), l))
            call subdivide_node_line(blk%lineY(0:int(dns%globalSize(2))*2**(l-2), l-1), &
                                     blk%lineY(0:int(dns%globalSize(2))*2**(l-1), l))
            call subdivide_node_line(blk%lineZ(0:int(dns%globalSize(3))*2**(l-2), l-1), &
                                     blk%lineZ(0:int(dns%globalSize(3))*2**(l-1), l))
        end do
    end subroutine build_level_lines

    pure integer function lid_index(blk, level, c) result(idx)
        type(block_set_type), intent(in) :: blk
        integer, intent(in) :: level, c(3)

        idx = int(blk%levelOffset(level)) + level_raster(blk, level, c)
    end function lid_index

    ! Leaf id occupying the level-`level` lattice cell c, or -1.
    integer(C_INT) function leaf_at(blk, level, c) result(id)
        type(block_set_type), intent(in) :: blk
        integer, intent(in) :: level, c(3)

        id = -1_C_INT
        if (level < 0 .or. level >= int(blk%nLevels)) return
        if (any(c < 0) .or. any(c >= int(blk%nTiles)*2**level)) return
        id = blk%lidOf(lid_index(blk, level, c))
    end function leaf_at

    ! Build the global leaf table: root tiling minus removed blocks, box
    ! refinement rounds, 2:1 smoothing over the 26-neighbourhood, then
    ! leaf ids along the Z-order curve of the finest lattice. Identical on
    ! every rank, so the distribution needs no communication.
    subroutine build_leaf_table(blk, dns, periodic, active, myrank, touch, buried)
        type(block_set_type), intent(inout) :: blk
        type(dns_type), intent(in) :: dns
        logical(C_BOOL), intent(in) :: periodic(1:3)
        integer(C_INT), intent(in), optional :: active(:)
        integer(C_INT), intent(in) :: myrank
        integer(C_INT), intent(in), optional :: touch(:,:), buried(:,:)

        ! lidOf markers during construction: -1 none/covered, -2 split,
        ! 0 leaf placeholder (real ids assigned at the end).
        integer, parameter :: M_NONE = -1, M_SPLIT = -2, M_LEAF = 0

        integer :: l, c(3), cn(3), cc(3), gx, gy, gz, i, n, id
        integer :: ox, oy, oz, sx, sy, sz, lmax, round, box
        integer(int64), allocatable :: keys(:)
        integer, allocatable :: order(:), tmpLevel(:), tmpCoord(:,:)
        logical :: changed, hit
        real(C_DOUBLE) :: lo(3), hi(3)

        lmax = int(blk%nLevels) - 1
        allocate(blk%levelOffset(0:int(blk%nLevels)))
        blk%levelOffset(0) = 0_C_INT
        do l = 0, lmax
            blk%levelOffset(l+1) = blk%levelOffset(l) + int(product(int(blk%nTiles)*2**l), C_INT)
        end do
        allocate(blk%lidOf(int(blk%levelOffset(int(blk%nLevels)))))
        blk%lidOf = int(M_NONE, C_INT)

        ! Root tiling, minus blocks removed inside the immersed boundary.
        i = 0
        do gz = 0, int(blk%nTiles(3)) - 1
            do gy = 0, int(blk%nTiles(2)) - 1
                do gx = 0, int(blk%nTiles(1)) - 1
                    i = i + 1
                    if (present(active)) then
                        if (active(i) == 0_C_INT) cycle
                    end if
                    blk%lidOf(lid_index(blk, 0, [gx, gy, gz])) = int(M_LEAF, C_INT)
                end do
            end do
        end do

        ! Refinement rounds: split leaves whose physical region intersects
        ! the [blocks] refine box, and (refine_body) leaves that straddle
        ! the immersed surface or 26-neighbour one that does (the >= 1
        ! block buffer of strategy doc Section 4).
        do round = 1, lmax
            l = round - 1
            do gz = 0, int(blk%nTiles(3))*2**l - 1
                do gy = 0, int(blk%nTiles(2))*2**l - 1
                    do gx = 0, int(blk%nTiles(1))*2**l - 1
                        c = [gx, gy, gz]
                        if (blk%lidOf(lid_index(blk, l, c)) /= int(M_LEAF, C_INT)) cycle
                        hit = .false.
                        if (refine_box_set(dns)) then
                            lo(1) = blk%lineX(c(1)*int(blk%nb(1)), l+1)
                            hi(1) = blk%lineX((c(1)+1)*int(blk%nb(1)), l+1)
                            lo(2) = blk%lineY(c(2)*int(blk%nb(2)), l+1)
                            hi(2) = blk%lineY((c(2)+1)*int(blk%nb(2)), l+1)
                            lo(3) = blk%lineZ(c(3)*int(blk%nb(3)), l+1)
                            hi(3) = blk%lineZ((c(3)+1)*int(blk%nb(3)), l+1)
                            do box = 1, int(dns%block_refine_nboxes)
                                hit = hit .or. &
                                  (lo(1) < dns%block_refine_box(2,box) .and. hi(1) > dns%block_refine_box(1,box) &
                              .and. lo(2) < dns%block_refine_box(4,box) .and. hi(2) > dns%block_refine_box(3,box) &
                              .and. lo(3) < dns%block_refine_box(6,box) .and. hi(3) > dns%block_refine_box(5,box))
                            end do
                        end if
                        if (.not. hit .and. present(touch)) then
                            hit = touch(level_raster(blk, l, c), l+1) /= 0_C_INT
                            if (.not. hit) then
                                buffer: do oz = -1, 1
                                do oy = -1, 1
                                do ox = -1, 1
                                    if (ox == 0 .and. oy == 0 .and. oz == 0) cycle
                                    cn = c + [ox, oy, oz]
                                    if (.not. wrap_lattice(blk, periodic, l, cn)) cycle
                                    if (touch(level_raster(blk, l, cn), l+1) /= 0_C_INT) then
                                        hit = .true.
                                        exit buffer
                                    end if
                                end do
                                end do
                                end do buffer
                            end if
                        end if
                        if (hit) call split_leaf(blk, l, c)
                    end do
                end do
            end do
        end do

        ! 2:1 smoothing: a leaf whose 26-neighbour cell is split into
        ! grandchildren (level >= l+2 below it) must split too. Iterate to
        ! a fixed point.
        changed = .true.
        do while (changed)
            changed = .false.
            do l = 0, lmax - 2
                do gz = 0, int(blk%nTiles(3))*2**l - 1
                    do gy = 0, int(blk%nTiles(2))*2**l - 1
                        do gx = 0, int(blk%nTiles(1))*2**l - 1
                            c = [gx, gy, gz]
                            if (blk%lidOf(lid_index(blk, l, c)) /= int(M_LEAF, C_INT)) cycle
                            outer: do oz = -1, 1
                            do oy = -1, 1
                            do ox = -1, 1
                                if (ox == 0 .and. oy == 0 .and. oz == 0) cycle
                                cn = c + [ox, oy, oz]
                                if (.not. wrap_lattice(blk, periodic, l, cn)) cycle
                                if (blk%lidOf(lid_index(blk, l, cn)) /= int(M_SPLIT, C_INT)) cycle
                                ! any split child => neighbour reaches l+2
                                do sz = 0, 1
                                do sy = 0, 1
                                do sx = 0, 1
                                    cc = 2*cn + [sx, sy, sz]
                                    if (blk%lidOf(lid_index(blk, l+1, cc)) == int(M_SPLIT, C_INT)) then
                                        call split_leaf(blk, l, c)
                                        changed = .true.
                                        exit outer
                                    end if
                                end do
                                end do
                                end do
                            end do
                            end do
                            end do outer
                        end do
                    end do
                end do
            end do
        end do

        ! Collect leaves and number them along the finest-lattice Z-order.
        n = 0
        do l = 0, lmax
            n = n + count(blk%lidOf(int(blk%levelOffset(l))+1:int(blk%levelOffset(l+1))) == int(M_LEAF, C_INT))
        end do
        allocate(tmpLevel(n), tmpCoord(3,n), keys(n), order(n))
        i = 0
        do l = 0, lmax
            do gz = 0, int(blk%nTiles(3))*2**l - 1
                do gy = 0, int(blk%nTiles(2))*2**l - 1
                    do gx = 0, int(blk%nTiles(1))*2**l - 1
                        c = [gx, gy, gz]
                        if (blk%lidOf(lid_index(blk, l, c)) /= int(M_LEAF, C_INT)) cycle
                        if (present(buried)) then
                            ! Removal at every level (Phase 2 generalized).
                            if (buried(level_raster(blk, l, c), l+1) /= 0_C_INT) then
                                blk%lidOf(lid_index(blk, l, c)) = int(M_NONE, C_INT)
                                cycle
                            end if
                        end if
                        i = i + 1
                        tmpLevel(i) = l
                        tmpCoord(:,i) = c
                        keys(i) = morton_key(gx*2**(lmax-l), gy*2**(lmax-l), gz*2**(lmax-l))
                    end do
                end do
            end do
        end do
        ! Buried leaves were dropped during collection; sort the survivors.
        n = i
        call heapsort_index(keys(1:n), order(1:n))
        blk%nBlocksGlobal = int(n, C_INT)
        allocate(blk%leafLevel(max(1,n)), blk%leafCoord(3,max(1,n)))
        do i = 1, n
            id = i - 1
            blk%leafLevel(i) = int(tmpLevel(order(i)), C_INT)
            blk%leafCoord(:,i) = int(tmpCoord(:,order(i)), C_INT)
            blk%lidOf(lid_index(blk, tmpLevel(order(i)), tmpCoord(:,order(i)))) = int(id, C_INT)
        end do
        ! Clear construction markers so lookups see only leaf ids and -1.
        where (blk%lidOf == int(M_SPLIT, C_INT)) blk%lidOf = -1_C_INT

        if (lmax > 0 .and. myrank == 0_C_INT) then
            print *, "block refinement:", n, "leaves,", &
                count(blk%leafLevel(1:n) > 0_C_INT), "refined"
        end if

        deallocate(tmpLevel, tmpCoord, keys, order)
    end subroutine build_leaf_table

    pure integer function level_raster(blk, l, c) result(r)
        type(block_set_type), intent(in) :: blk
        integer, intent(in) :: l, c(3)

        integer :: nl(3)

        nl = int(blk%nTiles)*2**l
        r = 1 + c(1) + nl(1)*(c(2) + nl(2)*c(3))
    end function level_raster

    subroutine split_leaf(blk, l, c)
        type(block_set_type), intent(inout) :: blk
        integer, intent(in) :: l, c(3)

        integer :: sx, sy, sz

        blk%lidOf(lid_index(blk, l, c)) = -2_C_INT
        do sz = 0, 1
            do sy = 0, 1
                do sx = 0, 1
                    blk%lidOf(lid_index(blk, l+1, 2*c + [sx, sy, sz])) = 0_C_INT
                end do
            end do
        end do
    end subroutine split_leaf

    ! Wrap level-l lattice coords in periodic directions; .false. outside
    ! non-periodic boundaries.
    logical function wrap_lattice(blk, periodic, l, c)
        type(block_set_type), intent(in) :: blk
        logical(C_BOOL), intent(in) :: periodic(1:3)
        integer, intent(in) :: l
        integer, intent(inout) :: c(3)

        integer :: d, nl

        wrap_lattice = .true.
        do d = 1, 3
            nl = int(blk%nTiles(d))*2**l
            if (c(d) < 0 .or. c(d) >= nl) then
                if (.not. periodic(d)) then
                    wrap_lattice = .false.
                    return
                end if
                c(d) = modulo(c(d), nl)
            end if
        end do
    end function wrap_lattice

    integer(C_INT) function count_level0_leaves(blk) result(n)
        type(block_set_type), intent(in) :: blk

        n = int(count(blk%lidOf(1:int(blk%levelOffset(1))) >= 0_C_INT), C_INT)
    end function count_level0_leaves

    ! Face kind of the level-`level` block at `origin` (level-l cells) in
    ! direction d, side -1/+1: physical wall outside non-periodic
    ! boundaries; open towards a same-level leaf; COARSE/FINE towards a
    ! coarser/finer leaf (2:1 interface); closed towards a removed block.
    integer(C_INT) function face_kind(blk, dns, periodic, origin, level, d, side) result(fk)
        type(block_set_type), intent(in) :: blk
        type(dns_type), intent(in) :: dns
        logical(C_BOOL), intent(in) :: periodic(1:3)
        integer(C_INT), intent(in) :: origin(3), level
        integer, intent(in) :: d, side

        integer :: to(3), cl(3), cc(3), l, sx, sy, sz

        if (blk%distMode /= DIST_ZORDER) then
            ! Rank-box mode: one block per rank, walls only.
            to = int(origin)
            to(d) = to(d) + side*int(blk%nb(d))
            fk = FACE_OPEN
            if ((to(d) < 0 .or. to(d) >= int(dns%globalSize(d))) .and. .not. periodic(d)) fk = FACE_PHYS
            return
        end if

        l = int(level)
        to = int(origin)
        to(d) = to(d) + side*int(blk%nb(d))
        if (to(d) < 0 .or. to(d) >= int(level_cells(dns, d, level))) then
            if (.not. periodic(d)) then
                fk = FACE_PHYS
                return
            end if
            to(d) = modulo(to(d), int(level_cells(dns, d, level)))
        end if

        cl = to/int(blk%nb)
        if (leaf_at(blk, l, cl) >= 0_C_INT) then
            fk = FACE_OPEN
            return
        end if
        if (leaf_at(blk, l - 1, cl/2) >= 0_C_INT) then
            fk = FACE_COARSE
            return
        end if
        do sz = 0, 1
            do sy = 0, 1
                do sx = 0, 1
                    cc = 2*cl + [sx, sy, sz]
                    if (leaf_at(blk, l + 1, cc) >= 0_C_INT) then
                        fk = FACE_FINE
                        return
                    end if
                end do
            end do
        end do
        fk = FACE_CLOSED
    end function face_kind

    ! True if the level-l lattice cell cl is occupied by a leaf at the same,
    ! the coarser (cl/2), or any of the finer (2*cl+child) levels -- i.e. a
    ! non-wall neighbour exists at some level (the 2:1-aware occupancy test).
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
