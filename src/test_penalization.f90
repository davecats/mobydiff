! Gate driver for the penalization integrating factors (numerics review F7;
! the rational AMPHIBIOUS form since 2026-10-01). With x = lambda*dt and
! P3 = 1 + x + x^2/2 + x^3/6 (the third-order Taylor polynomial of e^x),
! ibm.f90 provides
!
!   penal_state_factor(x)     = 1/P3                      (for e^{-x})
!   penal_incr_factor(x)      = (1 + x/2 + x^2/6)/P3      (for (1 - e^{-x})/x)
!   penal_state_minus_incr(x) = -x (1/2 + x/6)/P3
!
! Checked against reference values from exact rational arithmetic (Python
! fractions on the same doubles, 1/6 as the double the code multiplies with)
! and against what the scheme relies on:
!
!   state + x*incr = 1            (the steady fixed point lambda u = R is exact)
!   state - incr = minus_incr     (the predictor's state correction)
!   incr(0) = state(0) = 1, minus_incr(0) = 0    EXACTLY (no body: untouched)
!   solid cell (x ~ 1e28): 1 - state = 1 EXACTLY (a Dirichlet scalar holds the
!       body value to the last bit), incr = 1/x and minus_incr = -1/x to round-off
!   positive and monotonically decreasing in x >= 0
!   third-order agreement with the exponential: |state - e^{-x}| ~ x^4/24
!
! Run:  mpirun -n 1 build_cpu/penalization_test      (exit 0 = PASS)
program test_penalization
    use, intrinsic :: iso_c_binding
    use :: ibmm, only: penal_incr_factor, penal_state_factor, penal_state_minus_incr
    implicit none

    integer, parameter :: N = 13
    real(C_DOUBLE), parameter :: xs(N) = [0.0d0, 1.0d-12, 1.0d-6, 1.0d-3, 0.05d0, &
        0.1d0, 0.5d0, 1.0d0, 4.3d0, 10.0d0, 100.0d0, 1.0d6, 1.0d28]
    real(C_DOUBLE), parameter :: incrRef(N) = [1.0d0, 0.9999999999995d0, &
        0.9999995000001667d0, 0.9995001665834166d0, 0.9754067497671469d0, &
        0.9515910119137385d0, 0.7848101265822784d0, 0.625d0, 0.22419158517061705d0, &
        0.09956076134699854d0, 0.00999994178182545d0, 1.0d-06, 1.0000000000000001d-28]
    real(C_DOUBLE), parameter :: stateRef(N) = [1.0d0, 0.999999999999d0, &
        0.9999990000005d0, 0.9990004998334165d0, 0.9512296625116426d0, &
        0.9048408988086262d0, 0.6075949367088608d0, 0.375d0, 0.03597618376634668d0, &
        0.004392386530014641d0, 5.821817454973094d-06, 5.9999820000180005d-18, &
        6.000000000000001d-84]
    real(C_DOUBLE), parameter :: diffRef(N) = [0.0d0, -4.999999999996667d-13, &
        -4.9999966666675d-07, -0.0004996667500000139d0, -0.024177087255504253d0, &
        -0.04675011310511235d0, -0.17721518987341772d0, -0.25d0, -0.18821540140427037d0, &
        -0.0951683748169839d0, -0.009994119964370477d0, -9.99999999994d-07, &
        -1.0000000000000001d-28]
    real(C_DOUBLE), parameter :: TOL = 4.0d-15
    integer :: i, nfail
    real(C_DOUBLE) :: x, mu, a, d, eIncr, eState, eDiff, eSum, eSub, prevMu, prevA, ratio

    nfail = 0
    print '(A)', "         x        incr rel.err   state rel.err    diff rel.err  |state+x*incr-1|  |state-incr-diff|"
    do i = 1, N
        x = xs(i)
        mu = penal_incr_factor(x)
        a = penal_state_factor(x)
        d = penal_state_minus_incr(x)
        eIncr = abs(mu - incrRef(i))/incrRef(i)
        eState = abs(a - stateRef(i))/stateRef(i)
        eDiff = abs(d - diffRef(i))
        if (diffRef(i) /= 0.0d0) eDiff = eDiff/abs(diffRef(i))
        eSum = abs(a + x*mu - 1.0d0)
        ! state - incr cancels at small x; the absolute error is the statement.
        eSub = abs((a - mu) - d)
        print '(ES10.2,5ES16.3)', x, eIncr, eState, eDiff, eSum, eSub
        if (eIncr > TOL .or. eState > TOL .or. eDiff > TOL .or. eSum > TOL .or. eSub > TOL) &
            nfail = nfail + 1
    end do

    ! The exact limits the scheme is built on.
    if (penal_incr_factor(0.0d0) /= 1.0d0 .or. penal_state_factor(0.0d0) /= 1.0d0 &
        .or. penal_state_minus_incr(0.0d0) /= 0.0d0) then
        print '(A)', "FAIL: x = 0 is not exactly (1, 1, 0)"; nfail = nfail + 1
    end if
    if (1.0d0 - penal_state_factor(1.0d28) /= 1.0d0 .or. 1.0d0 - penal_state_factor(1.0d20) /= 1.0d0) then
        print '(A)', "FAIL: 1 - state is not exactly 1 in a solid cell"; nfail = nfail + 1
    end if

    ! Positive and monotonically decreasing over the whole range a
    ! coefficient can take (a graded cut cell to a solid one).
    prevMu = 1.0d0
    prevA = 1.0d0
    x = 1.0d-8
    do while (x < 1.0d30)
        mu = penal_incr_factor(x)
        a = penal_state_factor(x)
        if (.not. (mu > 0.0d0 .and. a > 0.0d0 .and. mu <= prevMu .and. a <= prevA .and. a <= mu)) then
            print '(A,ES12.4)', "FAIL: not positive / monotone at x = ", x; nfail = nfail + 1
        end if
        prevMu = mu
        prevA = a
        x = x*1.25d0
    end do

    ! Third-order agreement with the exponential: state - e^{-x} -> x^4/24.
    do i = 1, 3
        x = 0.2d0/2.0d0**i
        ratio = (penal_state_factor(x) - exp(-x))/(x**4/24.0d0)
        print '(A,F8.5,A,F8.5)', "   x = ", x, "   (state - e^-x)/(x^4/24) = ", ratio
        if (abs(ratio - 1.0d0) > 2.0d0*x) then
            print '(A)', "FAIL: the state factor is not the third-order one"; nfail = nfail + 1
        end if
    end do

    if (nfail > 0) then
        print '(A,I0,A)', "penalization factors: ", nfail, " FAILURE(S)"
        error stop 1
    end if
    print '(A)', "penalization factors: PASS"
end program test_penalization
