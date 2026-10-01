! Gate driver for the exact penalization factors (numerics review F7, step 9,
! 2026-10-01): ibm.f90 penal_incr_factor(x) = (1 - e^{-x})/x and
! penal_state_factor(x) = e^{-x}, x = lambda*dt.
!
! Checked against reference values from Python (-math.expm1(-x)/x and
! math.exp(-x)) across the series/closed-form switch at x = 0.1, and against
! the identities the scheme relies on:
!
!   state + x*incr = 1            (the two factors integrate du/dt = R - lambda u)
!   incr(0) = state(0) = 1        EXACTLY (no body: the predictor is untouched)
!   state(x) = 0, incr(x) = 1/x   EXACTLY for a solid cell (x ~ 1e28)
!   incr = 1/(x + B), state = B/(x + B),  B = x/(e^x - 1)   (the review's form)
!
! Run:  mpirun -n 1 build_cpu/penalization_test      (exit 0 = PASS)
program test_penalization
    use, intrinsic :: iso_c_binding
    use :: ibmm, only: penal_incr_factor, penal_state_factor
    implicit none

    integer, parameter :: N = 14
    real(C_DOUBLE), parameter :: xs(N) = [0.0d0, 1.0d-12, 1.0d-6, 1.0d-3, 0.05d0, &
        0.0999d0, 0.1d0, 0.5d0, 1.0d0, 10.0d0, 100.0d0, 700.0d0, 800.0d0, 1.0d28]
    real(C_DOUBLE), parameter :: incrRef(N) = [1.0d0, 0.9999999999995d0, &
        0.9999995000001668d0, 0.9995001666250084d0, 0.9754115099857198d0, &
        0.9516726095885778d0, 0.9516258196404043d0, 0.7869386805747332d0, &
        0.6321205588285577d0, 0.09999546000702375d0, 0.01d0, &
        0.0014285714285714286d0, 0.00125d0, 1.0000000000000001d-28]
    real(C_DOUBLE), parameter :: stateRef(N) = [1.0d0, 0.999999999999d0, &
        0.9999990000005d0, 0.999000499833375d0, 0.951229424500714d0, &
        0.9049279063021011d0, 0.9048374180359595d0, 0.6065306597126334d0, &
        0.36787944117144233d0, 4.5399929762484854d-05, 3.720075976020836d-44, &
        9.85967654375977d-305, 0.0d0, 0.0d0]
    real(C_DOUBLE), parameter :: TOL = 4.0d-15
    integer :: i, nfail
    real(C_DOUBLE) :: x, mu, a, eIncr, eState, eSum, eB, bfac

    nfail = 0
    print '(A)', "         x        incr rel.err   state rel.err   |state+x*incr-1|   |incr-1/(x+B)|"
    do i = 1, N
        x = xs(i)
        mu = penal_incr_factor(x)
        a = penal_state_factor(x)
        eIncr = abs(mu - incrRef(i))/incrRef(i)
        eState = abs(a - stateRef(i))
        if (stateRef(i) > 0.0d0) eState = eState/stateRef(i)
        eSum = abs(a + x*mu - 1.0d0)
        ! The review's B form, where e^x does not overflow.
        eB = 0.0d0
        if (x > 0.0d0 .and. x < 700.0d0) then
            bfac = x/(exp(x) - 1.0d0)
            if (x < 1.0d-3) bfac = 1.0d0 - 0.5d0*x + x*x/12.0d0
            eB = abs(mu - 1.0d0/(x + bfac))/mu
        end if
        print '(ES10.2,4ES16.3)', x, eIncr, eState, eSum, eB
        if (eIncr > TOL .or. eState > TOL .or. eSum > TOL .or. eB > 1.0d-11) nfail = nfail + 1
    end do

    ! The exact limits the scheme is built on.
    if (penal_incr_factor(0.0d0) /= 1.0d0 .or. penal_state_factor(0.0d0) /= 1.0d0) then
        print '(A)', "FAIL: x = 0 is not exactly (1, 1)"; nfail = nfail + 1
    end if
    if (penal_state_factor(1.0d28) /= 0.0d0 .or. penal_incr_factor(1.0d28) /= 1.0d0/1.0d28) then
        print '(A)', "FAIL: the solid limit is not exactly (0, 1/x)"; nfail = nfail + 1
    end if

    if (nfail > 0) then
        print '(A,I0,A)', "penalization factors: ", nfail, " FAILURE(S)"
        error stop 1
    end if
    print '(A)', "penalization factors: PASS"
end program test_penalization
