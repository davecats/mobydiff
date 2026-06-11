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
    use :: init, only: dns_type, grid_type, NVAR, NVEL, slice_grid_direction
    implicit none

    private
    public :: block_set_type
    public :: init_block_set, destroy_block_set
    public :: enter_block_data, exit_block_data
    public :: subdivide_node_line

    type :: block_set_type
        ! All blocks of a set share the same cell count nb(1:3).
        ! Phase 0: nb = rank-local box (one block); later phases: a cube.
        integer(C_INT) :: nb(1:3) = 0_C_INT
        integer(C_INT) :: nBlocks = 0_C_INT
        integer(C_INT) :: nLevels = 1_C_INT

        ! Per-block metadata (small, host-resident).
        !   level:    refinement level; 0 = base grid, +1 per cell bisection
        !   origin:   zero-based cell origin of the block in the level-l
        !             global index space
        !   globalId: position in the global (Z-ordered) block table
        ! Neighbour adjacency and per-face kinds (same/coarse/fine/physical/
        ! closed) are introduced in Phase 1 together with the block-pair
        ! halo exchange; see the strategy document.
        integer(C_INT), allocatable :: level(:)        ! (nBlocks)
        integer(C_INT), allocatable :: origin(:,:)     ! (3,nBlocks)
        integer(C_INT), allocatable :: globalId(:)     ! (nBlocks)

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

    ! Phase 0 constructor: one block spanning this rank's box of the global
    ! grid, with coordinates and metrics sliced from the global node lines
    ! by slice_grid_direction.
    subroutine init_block_set(blk, dns, g, periodic, global_id)
        type(block_set_type), intent(inout) :: blk
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        logical(C_BOOL), intent(in) :: periodic(1:3)
        integer(C_INT), intent(in), optional :: global_id

        integer :: nx, ny, nz, b

        call destroy_block_set(blk)

        nx = int(dns%localSize(1,2))
        ny = int(dns%localSize(2,2))
        nz = int(dns%localSize(3,2))

        blk%nb = dns%localSize(1:3,2)
        blk%nBlocks = 1_C_INT
        blk%nLevels = 1_C_INT

        allocate(blk%level(blk%nBlocks))
        allocate(blk%origin(3,blk%nBlocks))
        allocate(blk%globalId(blk%nBlocks))
        blk%level = 0_C_INT
        blk%origin(:,1) = dns%localSize(1:3,0) - 1_C_INT
        blk%globalId = 0_C_INT
        if (present(global_id)) blk%globalId(1) = global_id

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

        blk%nb = 0_C_INT
        blk%nBlocks = 0_C_INT
        blk%nLevels = 1_C_INT
    end subroutine destroy_block_set

    ! Device mapping, mirroring enter_grid_data/enter_field_data in
    ! gpu_runtime: map the container once, then its flat member arrays.
    subroutine enter_block_data(blk)
        type(block_set_type), intent(inout) :: blk

#ifdef USE_OPENMP_OFFLOAD
        !$omp target enter data map(to: blk)
        !$omp target enter data map(to: &
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
        !$omp& blk%x, blk%y, blk%z, blk%d1x, blk%d1y, blk%d1z, &
        !$omp& blk%lapXm, blk%lapX0, blk%lapXp, &
        !$omp& blk%lapYm, blk%lapY0, blk%lapYp, &
        !$omp& blk%lapZm, blk%lapZ0, blk%lapZp, &
        !$omp& blk%q, blk%qs, blk%oldrhs)
        !$omp target exit data map(delete: blk)
#endif
    end subroutine exit_block_data

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
