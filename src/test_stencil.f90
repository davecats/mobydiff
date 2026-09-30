! Gate driver for the viscous stencil (numerics review F1, step 8, 2026-09-30):
! slice_grid_direction must build a FLUX-FORM second-derivative operator for
! every variable on every node line, i.e. for any values u(0:n+1)
!
!   sum_i  W_i * (lapM_i u_{i-1} + lap0_i u_i + lapP_i u_{i+1})
!        =  (u_{n+1}-u_n)/h_{n+1/2}  -  (u_1-u_0)/h_{1/2}          (telescoping)
!   W_i lapP_i = W_{i+1} lapM_{i+1}                                (symmetry)
!
! with W_i = 1/d1_i the width of the variable's OWN control volume and h the
! DOF-to-DOF distances. Both are algebraic identities of the coefficients,
! so they hold to round-off or not at all: the Taylor three-point weights
! 2/(hm(hm+hp)) satisfy them only where W_i = (hm+hp)/2 -- a component's
! face-staggered direction, or a uniform line. Every stretched line in a
! cell-centred direction failed both before step 8; the informational
! column prints how far the Taylor weights sit from the flux weights.
! Run:  mpirun -n 1 build_cpu/stencil_test      (exit 0 = PASS)
program test_stencil
    use, intrinsic :: iso_c_binding
    use :: init, only: slice_grid_direction, NVAR, VAR_U, VAR_V, VAR_W, VAR_P
    implicit none

    integer, parameter :: NLINES = 4
    character(len=12), parameter :: lineName(NLINES) = &
        [character(len=12) :: "uniform", "geometric", "cosine", "tanh"]
    integer :: il, dir, var, nfail
    real(C_DOUBLE) :: cons, symm, taylor

    nfail = 0
    print '(A)', "line         dir var   conservation   symmetry   max|taylor/flux-1|"
    do il = 1, NLINES
        do dir = 1, 3
            do var = VAR_U, VAR_P
                call check_line(il, dir, var, cons, symm, taylor)
                print '(A12,I4,I4,2ES14.2,ES12.2)', lineName(il), dir, var, cons, symm, taylor
                if (cons > 1.0d-12 .or. symm > 1.0d-12) nfail = nfail + 1
            end do
        end do
    end do
    if (nfail == 0) then
        print '(A)', "stencil_test: PASS (flux form on every line, direction and variable)"
    else
        print '(A,I0,A)', "stencil_test: FAIL (", nfail, " line/dir/var combinations)"
        error stop 1
    end if

contains

    subroutine check_line(il, dir, var, cons, symm, taylor)
        integer, intent(in) :: il, dir, var
        real(C_DOUBLE), intent(out) :: cons, symm, taylor
        integer, parameter :: n = 48
        real(C_DOUBLE), parameter :: pi = 3.1415926535897932384626433832795d0
        real(C_DOUBLE), parameter :: length = 2.0d0
        real(C_DOUBLE) :: node(0:n), coord(-1:n+2,NVAR), d1(0:n+1,NVAR)
        real(C_DOUBLE) :: lapM(0:n+1,NVAR), lap0(0:n+1,NVAR), lapP(0:n+1,NVAR)
        real(C_DOUBLE) :: u(0:n+1), w(1:n), s, r, hm, hp, lsum, ftop, fbot, scale, tm, tp
        logical(C_BOOL) :: periodic
        integer :: i

        ! The node lines: the identities must hold on every one of them.
        do i = 0, n
            s = real(i, C_DOUBLE)/real(n, C_DOUBLE)
            select case (il)
            case (1); node(i) = length*s
            case (2); r = 1.10d0; node(i) = length*(r**i - 1.0d0)/(r**n - 1.0d0)
            case (3); node(i) = 0.5d0*length*(1.0d0 - cos(pi*s))
            case (4); node(i) = 0.5d0*length*(1.0d0 + tanh(2.0d0*(2.0d0*s - 1.0d0))/tanh(2.0d0))
            end select
        end do
        periodic = (il == 1)        ! the uniform line doubles as the periodic control
        coord = 0.0d0; d1 = 0.0d0; lapM = 0.0d0; lap0 = 0.0d0; lapP = 0.0d0
        call slice_grid_direction(node, coord, d1, lapM, lap0, lapP, int(n, C_INT), 1_C_INT, n, &
                                  length, periodic, dir)

        ! Arbitrary (deterministic) values, ghosts included: the identities
        ! are properties of the coefficients, not of any particular field.
        do i = 0, n+1
            u(i) = sin(1.7d0*i + 0.3d0) + 0.25d0*cos(0.61d0*i*i)
        end do
        do i = 1, n
            w(i) = 1.0d0/d1(i,var)
        end do

        lsum = 0.0d0; scale = 0.0d0
        do i = 1, n
            lsum = lsum + w(i)*(lapM(i,var)*u(i-1) + lap0(i,var)*u(i) + lapP(i,var)*u(i+1))
            scale = scale + abs(w(i)*lapM(i,var)*u(i-1)) + abs(w(i)*lapP(i,var)*u(i+1))
        end do
        ftop = (u(n+1) - u(n))/(coord(n+1,var) - coord(n,var))
        fbot = (u(1) - u(0))/(coord(1,var) - coord(0,var))
        cons = abs(lsum - (ftop - fbot))/scale

        symm = 0.0d0; taylor = 0.0d0
        do i = 1, n-1
            symm = max(symm, abs(w(i)*lapP(i,var) - w(i+1)*lapM(i+1,var)) / abs(w(i)*lapP(i,var)))
        end do
        do i = 1, n
            hm = coord(i,var) - coord(i-1,var)
            hp = coord(i+1,var) - coord(i,var)
            tm = 2.0d0/(hm*(hm + hp)); tp = 2.0d0/(hp*(hm + hp))
            taylor = max(taylor, abs(tm/lapM(i,var) - 1.0d0), abs(tp/lapP(i,var) - 1.0d0))
        end do
    end subroutine check_line

end program test_stencil
