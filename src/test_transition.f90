! Unit test for the T4 gamma-Re_thetat transition correlations (IDDES
! phase T4, docs/next_session_iddes.md): the pure functions in rans.f90,
! transcribed from OpenFOAM kOmegaSSTLM.C / Langtry & Menter (2009), are
! checked HOST-SIDE against tabulated values from an independent Python
! transcription of the same source (scratch lm_ref.py) before they ever
! run in a device kernel. Every piecewise branch is hit at least once.
! Run: mpirun -n 1 build_cpu/transition_test
program test_transition
    use, intrinsic :: iso_c_binding
    use :: rans, only: lm_rethetac, lm_flength, lm_fonset, lm_fturb, &
        lm_rethetat0, lm_fthetat
    implicit none

    integer :: nfail
    real(C_DOUBLE), parameter :: NU180 = 1.0d0/180.0d0

    nfail = 0

    ! Re_thetac: below/above the 1870 branch switch.
    call check("rethetac(20)",   lm_rethetac(20.0d0),   15.939214191200003d0)
    call check("rethetac(100)",  lm_rethetac(100.0d0),  89.2430055d0)
    call check("rethetac(300)",  lm_rethetac(300.0d0),  238.91404150000002d0)
    call check("rethetac(500)",  lm_rethetac(500.0d0),  361.1966375d0)
    call check("rethetac(800)",  lm_rethetac(800.0d0),  535.322594d0)
    call check("rethetac(1200)", lm_rethetac(1200.0d0), 802.8054099999999d0)
    call check("rethetac(1870)", lm_rethetac(1870.0d0), 1278.0731150689503d0)
    call check("rethetac(2500)", lm_rethetac(2500.0d0), 1603.23d0)

    ! F_length: all four fit branches with the sublayer blend off
    ! (y^2 omega/(200 nu) large -> Fsublayer ~ 0), then the blend itself.
    call check("flength(100)",  lm_flength(100.0d0,  1.0d0, 1.0d3, NU180), 37.300529999999995d0)
    call check("flength(399)",  lm_flength(399.0d0,  1.0d0, 1.0d3, NU180), 13.955228033000001d0)
    call check("flength(400)",  lm_flength(400.0d0,  1.0d0, 1.0d3, NU180), 13.84000000000006d0)
    call check("flength(595)",  lm_flength(595.0d0,  1.0d0, 1.0d3, NU180), 0.5002013687499982d0)
    call check("flength(596)",  lm_flength(596.0d0,  1.0d0, 1.0d3, NU180), 0.5d0)
    call check("flength(1199)", lm_flength(1199.0d0, 1.0d0, 1.0d3, NU180), 0.31910000000000005d0)
    call check("flength(1200)", lm_flength(1200.0d0, 1.0d0, 1.0d3, NU180), 0.3188d0)
    call check("flength(2000)", lm_flength(2000.0d0, 1.0d0, 1.0d3, NU180), 0.3188d0)
    call check("flength sublayer strong", &
        lm_flength(300.0d0, 0.05d0, 5000.0d0, NU180), 24.30977d0)
    call check("flength sublayer mid", &
        lm_flength(300.0d0, 0.02d0, 5000.0d0, NU180), 24.92426052179898d0)

    ! F_onset: sub-critical (0), the min/max limiter branches, the R_T
    ! shutoff, and the cap at 2.
    call check("fonset sub-crit", lm_fonset(100.0d0, 200.0d0, 1.0d0),  0.0d0)
    call check("fonset onset",    lm_fonset(500.0d0, 200.0d0, 1.0d0),  0.7529061143557976d0)
    call check("fonset high-RT",  lm_fonset(500.0d0, 200.0d0, 10.0d0), 1.6889061143557975d0)
    call check("fonset low",      lm_fonset(50.0d0,  300.0d0, 0.1d0),  0.0d0)
    call check("fonset capped",   lm_fonset(900.0d0, 200.0d0, 0.5d0),  1.008d0)

    call check("fturb(0.5)", lm_fturb(0.5d0), 0.9997558891748972d0)
    call check("fturb(2)",   lm_fturb(2.0d0), 0.9394130628134758d0)
    call check("fturb(8)",   lm_fturb(8.0d0), 1.1253517471925912d-7)

    ! Re_thetat,eq: both Tu branches, the Tu floor (0.027), the zero-
    ! pressure-gradient closed form, and converging lambda fixed points on
    ! both dU/ds signs (including the clamp at -0.1).
    call check("rethetat0 Tu=1",     lm_rethetat0(1.0d0, 0.0d0, 1.0d0, 1.0d0), 584.3016d0)
    call check("rethetat0 Tu=5",     lm_rethetat0(5.0d0, 0.0d0, 1.0d0, 1.0d0), 122.03108746641541d0)
    call check("rethetat0 Tu floor", lm_rethetat0(0.01d0, 0.0d0, 1.0d0, 1.0d0), 1458.8300119012347d0)
    call check("rethetat0 Tu=0.5",   lm_rethetat0(0.5d0, 0.0d0, 1.0d0, 1.0d0), 879.6744000000001d0)
    call check("rethetat0 Tu=1.3",   lm_rethetat0(1.3d0, 0.0d0, 1.0d0, 1.0d0), 407.3835408284024d0)
    call check("rethetat0 dUsds>0",  lm_rethetat0(3.0d0, 1.0d0, 1.0d-6, 1.0d0), 182.57541193398401d0)
    call check("rethetat0 dUsds<0",  lm_rethetat0(1.0d0, -1.0d0, 1.0d-3, 1.0d0), 425.74137334377053d0)
    call check("rethetat0 lam interior", &
        lm_rethetat0(5.0d0, 2.0d-4, 1.0d-2, 2.0d0), 122.03143698996632d0)

    ! F_thetat: the Fwake branch, the gamma branch (= 1 at gamma = 1/ce2),
    ! the freestream limit, and a mixed interior value.
    call check("fthetat BL edge",   lm_fthetat(1.0d0, 300.0d0, 1.0d0, 1.0d0, &
        50.0d0, 0.5d0, NU180), 0.9994938781163593d0)
    call check("fthetat gamma=1/ce2", lm_fthetat(0.02d0, 300.0d0, 1.0d0, 1.0d0, &
        50.0d0, 0.5d0, NU180), 1.0d0)
    call check("fthetat freestream", lm_fthetat(1.0d0, 300.0d0, 10.0d0, 0.1d0, &
        5.0d0, 0.05d0, NU180), 0.0014249764374710692d0)
    call check("fthetat mixed",     lm_fthetat(0.5d0, 100.0d0, 0.5d0, 20.0d0, &
        2000.0d0, 0.9d0, NU180), 0.7600999583506872d0)

    if (nfail > 0) then
        print '(A,I0,A)', "transition_test: ", nfail, " FAILURES"
        error stop
    end if
    print *, "transition_test: ALL PASS"

contains

    ! Relative comparison; the reference values come from Python (C libm),
    ! so allow a few-ulp cross-runtime spread on the transcendentals.
    subroutine check(name, got, want)
        character(len=*), intent(in) :: name
        real(C_DOUBLE), intent(in) :: got, want

        real(C_DOUBLE), parameter :: TOL = 1.0d-12

        if (abs(got - want) > TOL*max(abs(want), 1.0d0)) then
            print '(A,A,ES24.16,A,ES24.16)', name, ": got ", got, " want ", want
            nfail = nfail + 1
        end if
    end subroutine check

end program test_transition
