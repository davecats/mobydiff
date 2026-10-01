! Gate driver for the minimum-surface key-bit order (blocks.f90
! min_surface_key_order, docs/next_session_after_step9.md item 2) and the
! partition audit that goes with it (partition_surface).
!
! Part A: the bit order the rule gives for the layouts it was decided on,
! against the sequences of the independent Python mirror
! (tools/partition_analysis.py --order minsurface) and the two the handout
! states (base_jacobi at 8 and 16 ranks):
!
!   4096 x 176 x 192, blocks 2048 x 88 x 96, z periodic    x y z
!   the same in blocks 1024 x 88 x 96                       x x y z
!   the same in blocks 64 x 44 x 48 (rect_jacobi)           x x x x x y x z z y
!   channel 256 x 128 x 256, nb 32, x and z periodic        x x y z z x y z
!   64 x 1024 x 8 in blocks 16 x 256 x 8, no periodicity    y y x x
!   non-periodic cube 64^3, nb 16: the plain interleave     x y z x y z
!   two levels double every direction's tile count          (one more bit each)
!
! Part B: a real block set in the shape of base_jacobi scaled down
! (128 x 8 x 16 in 2 x 2 x 2 blocks, z periodic), built by the PRODUCTION
! builder, split over 8 ranks on two nodes: the builder must order it x-major
! (the node boundary on the 4 x 8 x-faces: 4 * 32 = 128 cells, where the old
! z-major order cut 2 * 4 * 64 * 4 = 2048), and partition_surface must count
! every interface once -- 128 (x) + 2048 (y) + 2048 (z, both periodic
! interfaces) = 4224 across ranks.
!
! Run:  mpirun -n 1 build_cpu/keyorder_test      (exit 0 = PASS)
program test_keyorder
    use, intrinsic :: iso_c_binding
    use, intrinsic :: iso_fortran_env, only: int64
    use :: init, only: dns_type, grid_type, init_grid
    use :: blocks, only: block_set_type, init_block_set, KEY_MAX_BITS, &
        min_surface_key_order, partition_surface
    implicit none

    logical(C_BOOL), parameter :: T = .true._C_BOOL, F = .false._C_BOOL
    integer :: nfail

    nfail = 0
    call check_order("base_jacobi, 8 ranks ", [4096, 176, 192], [2048, 88, 96], 1, [F, F, T], "xyz")
    call check_order("base_jacobi, 16 ranks", [4096, 176, 192], [1024, 88, 96], 1, [F, F, T], "xxyz")
    call check_order("rect_jacobi          ", [4096, 176, 192], [64, 44, 48], 1, [F, F, T], "xxxxxyxzzy")
    call check_order("channel nb 32        ", [256, 128, 256], [32, 32, 32], 1, [T, F, T], "xxyzzxyz")
    call check_order("non-cubic blocks     ", [64, 1024, 8], [16, 256, 8], 1, [F, F, F], "yyxx")
    call check_order("non-periodic cube    ", [64, 64, 64], [16, 16, 16], 1, [F, F, F], "xyzxyz")
    call check_order("the cube, two levels ", [64, 64, 64], [16, 16, 16], 2, [F, F, F], "xyzxyzxyz")

    call check_partition()

    if (nfail > 0) then
        print '(A,I0,A)', "key order: ", nfail, " FAILURE(S)"
        error stop 1
    end if
    print '(A)', "key order: PASS"

contains

    subroutine check_order(label, grid, nb, nLevels, periodic, expected)
        character(len=*), intent(in) :: label, expected
        integer, intent(in) :: grid(3), nb(3), nLevels
        logical(C_BOOL), intent(in) :: periodic(3)

        integer(C_INT) :: order(KEY_MAX_BITS)
        integer :: n, i
        character(len=KEY_MAX_BITS) :: got

        call min_surface_key_order(int(grid, C_INT), int(nb, C_INT), int(nLevels, C_INT), &
            periodic, order, n)
        got = ""
        do i = 1, n
            got(i:i) = "xyz"(order(i):order(i))
        end do
        if (trim(got) == expected) then
            print '(3A)', "  ", label, "  " // trim(got)
        else
            print '(5A)', "  ", label, "  " // trim(got), "   FAIL: expected ", expected
            nfail = nfail + 1
        end if
    end subroutine check_order

    subroutine check_partition()
        type(dns_type) :: dns
        type(grid_type) :: g
        type(block_set_type) :: blk
        logical(C_BOOL) :: periodic(3)
        integer(int64) :: cellsRank, cellsNode
        integer :: i

        dns%globalSize(1:3) = [128_C_INT, 8_C_INT, 16_C_INT]
        dns%leng(1:3) = [16.0d0, 1.0d0, 2.0d0]
        dns%block_nb = [64_C_INT, 4_C_INT, 8_C_INT]
        periodic = [F, F, T]
        call init_grid(g, dns, periodic)
        call init_block_set(blk, dns, g, periodic, 1_C_INT, 0_C_INT)

        ! x-major: the first four leaves are the x = 0 half.
        do i = 1, int(blk%nBlocksGlobal)
            if ((blk%leafCoord(1, i) == 0_C_INT) .neqv. (i <= 4)) then
                print '(A,I0)', "  FAIL: the builder's order is not x-major at leaf ", i - 1
                nfail = nfail + 1
                exit
            end if
        end do

        call partition_surface(blk, periodic, 8_C_INT, [0, 0, 0, 0, 1, 1, 1, 1], cellsRank, cellsNode)
        print '(A,I0,A,I0)', "  partition 128 x 8 x 16 in 2 x 2 x 2, 8 ranks on 2 nodes: across ranks ", &
            cellsRank, ", across nodes ", cellsNode
        if (cellsRank /= 4224_int64 .or. cellsNode /= 128_int64) then
            print '(A)', "  FAIL: expected 4224 across ranks and 128 across nodes"
            nfail = nfail + 1
        end if
    end subroutine check_partition
end program test_keyorder
