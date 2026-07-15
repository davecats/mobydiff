! Gate driver for the leaf-table builder (2D refinement phase R2D-0,
! gate b — docs/next_session_refine2d.md): builds the PRODUCTION leaf
! table (blocks.f90 build_leaf_table) for a box-refined case and prints
! it row-by-row in the mobygeom blocks-dataset convention (origin in
! level-l cells for refined directions, global cells for fixed ones;
! rows ordered along the canonical Morton curve of the refine_dims
! mode). validation compares this against mobygeom's Python mirror
! (build_leaf_table_py). Run:
!   mpirun -n 1 build_cpu/leaftable_test nx ny nz lx ly lz nb levels \
!       dims x0 x1 y0 y1 z0 z1 [px py pz]
! dims = xyz | xz; px py pz = 1/0 periodicity flags (default periodic).
program test_leaftable
    use, intrinsic :: iso_c_binding
    use :: init, only: dns_type, grid_type, init_grid
    use :: blocks, only: block_set_type, init_block_set
    implicit none

    type(dns_type) :: dns
    type(grid_type) :: g
    type(block_set_type) :: blk
    logical(C_BOOL) :: periodic(1:3)
    character(len=64) :: arg
    integer :: i, nargs, p

    nargs = command_argument_count()
    if (nargs < 15) then
        print *, "usage: leaftable_test nx ny nz lx ly lz nb levels dims" &
            // " x0 x1 y0 y1 z0 z1 [px py pz]"
        error stop 1
    end if

    do i = 1, 3
        call get_command_argument(i, arg)
        read(arg, *) dns%globalSize(i)
        call get_command_argument(3 + i, arg)
        read(arg, *) dns%leng(i)
    end do
    call get_command_argument(7, arg)
    read(arg, *) dns%block_nb
    call get_command_argument(8, arg)
    read(arg, *) dns%block_refine_levels
    call get_command_argument(9, arg)
    select case (trim(arg))
    case ("xyz")
        dns%block_refine_mask = 1_C_INT
    case ("xz")
        dns%block_refine_mask = [1_C_INT, 0_C_INT, 1_C_INT]
    case default
        error stop "dims must be xyz or xz"
    end select
    do i = 1, 6
        call get_command_argument(9 + i, arg)
        read(arg, *) dns%block_refine_box(i, 1)
    end do
    dns%block_refine_nboxes = 1_C_INT

    periodic = .true._C_BOOL
    do i = 1, 3
        if (nargs >= 15 + i) then
            call get_command_argument(15 + i, arg)
            read(arg, *) p
            periodic(i) = p /= 0
        end if
    end do

    call init_grid(g, dns, periodic)
    call init_block_set(blk, dns, g, periodic, 1_C_INT, 0_C_INT)

    print '(A,1X,I0)', "leaftable", blk%nBlocksGlobal
    do i = 1, int(blk%nBlocksGlobal)
        print '(I0,1X,I0,1X,I0,1X,I0,1X,I0)', i - 1, &
            blk%leafCoord(1,i)*blk%nb(1), &
            blk%leafCoord(2,i)*blk%nb(2), &
            blk%leafCoord(3,i)*blk%nb(3), &
            blk%leafLevel(i)
    end do

end program test_leaftable
