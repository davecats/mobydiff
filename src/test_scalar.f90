! Unit test for the passive-scalar turbulent-Prandtl correlation (increment
! S2, docs/next_session_scalar.md Section 3): the pure declare-target
! function prt_kays in scalar.f90 is checked HOST-SIDE against tabulated
! values from an independent mpmath transcription of Kays-Crawford
! (scratch kc_ref.py, 50-1000 digits so the large-Pe_t cancellation is
! resolved exactly) before it ever runs in a device kernel.
!
! Every branch is hit: the Pe_t = 0 guard, the direct expression, the
! small-x series, both sides of the x = 1/2 crossover (a continuity check
! on the switch itself), and the two limits Pe_t -> 0 (Pr_t = 2 Prt_inf)
! and Pe_t -> infinity (Pr_t -> Prt_inf), the latter at Pe_t = 1e300 where
! the direct expression would overflow a^2 and cancel to nothing.
! Run: mpirun -n 1 build_cpu/scalar_test
!
! Increment S5a adds the THERMAL WALL FUNCTION's two correlations to the same
! driver, against the same kind of independent mpmath transcription:
! Jayatilleke's P, the thermal sublayer thickness y+_T (a root, so the
! bisection is checked against mpmath's), and the wall-cell eddy diffusivity
! itself -- including that it is EXACTLY zero on the conduction branch and
! continuous across the switch.
!
! Increment C3 adds the CONJUGATE CELL CAPACITY's geometry: the fluid volume
! fraction of a box cut by a plane (plane_box_fraction / scalar_volume_
! fraction, scalar.f90). It is checked two ways -- against hand-derivable
! exact values in each of the three degeneracy regimes (the plane parallel to
! two axes, to one, to none), and against BRUTE FORCE on a 300^3 midpoint
! lattice, which is the only check that cannot share a mistake with the
! closed form. The symmetry f(phi) + f(-phi) = 1 and the clipping to 0/1
! outside the cut band are checked too.
program test_scalar
    use, intrinsic :: iso_c_binding
    use :: scalar, only: prt_kays, jayatilleke_p, thermal_yplus, wall_diffusivity, &
        plane_box_fraction, scalar_volume_fraction
    implicit none

    integer :: nfail
    ! The branch switch sits at x = 1/(C Pe_t sqrt(Prt_inf)) = 1/2, i.e.
    ! Pe_t = 2/(0.3 sqrt(0.85)) = 7.2310152606218724 at the default Prt_inf.
    real(C_DOUBLE), parameter :: XOVER = 7.2310152606218724d0
    real(C_DOUBLE) :: lo, hi

    nfail = 0

    ! --- the Pe_t = 0 guard and the molecular limit Pr_t = 2 Prt_inf ------
    call check("prt(0, 0.85)",     prt_kays(0.0d0, 0.85d0),     1.7d0)
    call check("prt(0, 0.5)",      prt_kays(0.0d0, 0.5d0),      1.0d0)
    call check("prt(-1, 0.85)",    prt_kays(-1.0d0, 0.85d0),    1.7d0)   ! clamped
    call check("prt(1e-8, 0.85)",  prt_kays(1.0d-8, 0.85d0),    1.69999999059606473d0)

    ! --- the direct branch (x >= 1/2) ------------------------------------
    call check("prt(0.01, 0.85)",  prt_kays(0.01d0, 0.85d0),    1.69067352347742831d0)
    call check("prt(0.1, 0.85)",   prt_kays(0.1d0, 0.85d0),     1.61322883057133966d0)
    call check("prt(0.5, 0.85)",   prt_kays(0.5d0, 0.85d0),     1.37277890392146789d0)
    call check("prt(1, 0.85)",     prt_kays(1.0d0, 0.85d0),     1.21057713978955374d0)
    call check("prt(2, 0.85)",     prt_kays(2.0d0, 0.85d0),     1.06601325382850231d0)
    call check("prt(5, 0.85)",     prt_kays(5.0d0, 0.85d0),     0.946060549363522067d0)
    call check("prt(7, 0.85)",     prt_kays(7.0d0, 0.85d0),     0.91993751432939629d0)

    ! --- the crossover: last direct value vs first series value ----------
    lo = XOVER*(1.0d0 - 1.0d-9)
    hi = XOVER*(1.0d0 + 1.0d-9)
    call check("prt(xover-, 0.85)", prt_kays(lo, 0.85d0), 0.917805012427862318d0)
    call check("prt(xover+, 0.85)", prt_kays(hi, 0.85d0), 0.917805012298422177d0)
    ! The two branches must join smoothly: over a 2e-9 relative step in
    ! Pe_t the correlation may move only by its own slope (~1.3e-10).
    if (abs(prt_kays(hi, 0.85d0) - prt_kays(lo, 0.85d0)) > 1.0d-9) then
        print *, "crossover discontinuity: ", prt_kays(lo, 0.85d0), prt_kays(hi, 0.85d0)
        nfail = nfail + 1
    end if

    ! --- the series branch (x < 1/2) and the Pe_t -> infinity limit ------
    call check("prt(7.5, 0.85)",   prt_kays(7.5d0, 0.85d0),     0.915479835517276634d0)
    call check("prt(10, 0.85)",    prt_kays(10.0d0, 0.85d0),    0.899645234090461921d0)
    call check("prt(50, 0.85)",    prt_kays(50.0d0, 0.85d0),    0.860181922113475087d0)
    call check("prt(100, 0.85)",   prt_kays(100.0d0, 0.85d0),   0.855106500461899404d0)
    call check("prt(1e3, 0.85)",   prt_kays(1.0d3, 0.85d0),     0.850512042556171432d0)
    call check("prt(1e6, 0.85)",   prt_kays(1.0d6, 0.85d0),     0.850000512196759973d0)
    call check("prt(1e12, 0.85)",  prt_kays(1.0d12, 0.85d0),    0.850000000000512197d0)
    ! Pe_t = 1e300: a^2 would overflow and 1 - exp(-x) would round to 0 in
    ! the direct expression; the series gives the exact asymptote.
    call check("prt(1e300, 0.85)", prt_kays(1.0d300, 0.85d0),   0.85d0)

    ! --- other Prt_inf values (the key is per-scalar) ---------------------
    call check("prt(0.5, 0.5)",    prt_kays(0.5d0, 0.5d0),      0.840594780459121946d0)
    call check("prt(5, 0.5)",      prt_kays(5.0d0, 0.5d0),      0.572158470598205533d0)
    call check("prt(100, 0.5)",    prt_kays(100.0d0, 0.5d0),    0.50391289145265619d0)
    call check("prt(0.5, 1)",      prt_kays(0.5d0, 1.0d0),      1.59355278054902505d0)
    call check("prt(5, 1)",        prt_kays(5.0d0, 1.0d0),      1.10474225016094663d0)
    call check("prt(100, 1)",      prt_kays(100.0d0, 1.0d0),    1.00554008968553513d0)
    call check("prt(3, 0.4)",      prt_kays(3.0d0, 0.4d0),      0.499294742094104203d0)
    call check("prt(30, 1.2)",     prt_kays(30.0d0, 1.2d0),     1.22011344680583436d0)

    ! --- monotonicity: Pr_t falls from 2 Prt_inf to Prt_inf --------------
    call check_monotone()

    ! === S5a: the thermal wall function ==================================
    ! Jayatilleke's P(Pr/Pr_t): zero at Pr = Pr_t (the Reynolds analogy),
    ! negative below it, and strongly positive for high-Pr fluids.
    call check("P(0.71/0.85)", jayatilleke_p(0.71d0/0.85d0), -1.49146084477208598d0)
    call check("P(1)",         jayatilleke_p(1.0d0),          0.0d0)
    call check("P(2)",         jayatilleke_p(2.0d0),          8.03917714490408820d0)
    call check("P(0.1)",       jayatilleke_p(0.1d0),         -9.72250491069673163d0)
    call check("P(7)",         jayatilleke_p(7.0d0),         38.6626559386263292d0)
    call check("P(100)",       jayatilleke_p(100.0d0),      322.297542629418095d0)

    ! y+_T, the thermal sublayer thickness: the root of Pr y+ = Pr_t
    ! [ln(E y+)/kappa + P] beyond the minimum of the difference. Air
    ! (Pr < Pr_t) has a slightly THICKER thermal sublayer than the momentum
    ! y+_lam = 11.530; oil (Pr = 7) a much thinner one; a liquid metal
    ! (Pr = 0.025) is conductive far into the log layer.
    call check("ypt(0.71,0.85)",  thermal_yplus(0.71d0, 0.85d0, jayatilleke_p(0.71d0/0.85d0)), &
        12.1776453329441076d0)
    call check("ypt(1,0.85)",     thermal_yplus(1.0d0, 0.85d0, jayatilleke_p(1.0d0/0.85d0)), &
        11.0047461578706663d0)
    call check("ypt(7,0.9)",      thermal_yplus(7.0d0, 0.9d0, jayatilleke_p(7.0d0/0.9d0)), &
        6.81462087400107665d0)
    call check("ypt(0.025,0.85)", thermal_yplus(0.025d0, 0.85d0, jayatilleke_p(0.025d0/0.85d0)), &
        284.245504876139746d0)
    call check("ypt(0.71,0.9)",   thermal_yplus(0.71d0, 0.9d0, jayatilleke_p(0.71d0/0.9d0)), &
        12.4010219438244231d0)

    ! The wall-cell eddy diffusivity nu(y+/theta+ - 1/Pr), Re = 180.
    call check("alpha(5)",   alpha_air(5.0d0),   0.0d0)      ! conduction branch
    call check("alpha(12)",  alpha_air(12.0d0),  0.0d0)      ! still below y+_T
    call check("alpha(30)",  alpha_air(30.0d0),  0.00802520995676793475d0)
    call check("alpha(45)",  alpha_air(45.0d0),  0.0141902851232290114d0)
    call check("alpha(100)", alpha_air(100.0d0), 0.0348731015462603677d0)
    call check("alpha(30,Pr=7)", wall_diffusivity(30.0d0, 7.0d0, 0.9d0, 180.0d0, &
        jayatilleke_p(7.0d0/0.9d0), thermal_yplus(7.0d0, 0.9d0, jayatilleke_p(7.0d0/0.9d0))), &
        0.00247715805096350751d0)
    ! Continuity across the switch: theta+ = Pr y+ AT y+_T by construction,
    ! so the log branch must start from zero there too.
    if (abs(alpha_air(12.1776453329441076d0*(1.0d0 + 1.0d-12))) > 1.0d-13) then
        print *, "thermal wall function discontinuous at y+_T: ", &
            alpha_air(12.1776453329441076d0*(1.0d0 + 1.0d-12))
        nfail = nfail + 1
    end if

    ! === C3: the fluid volume fraction of a cut cell =====================
    call check_fraction_exact()
    call check_fraction_symmetry()
    call check_fraction_brute()

    if (nfail > 0) then
        print '(A,I0,A)', "scalar_test: ", nfail, " FAILURES"
        error stop
    end if
    print *, "scalar_test: ALL PASS"

contains

    ! The S5a wall diffusivity at Pr = 0.71, Pr_t = 0.85, Re = 180 (air in
    ! the gate channels), with its two constants formed the way init_scalar
    ! forms them.
    real(C_DOUBLE) function alpha_air(yplus) result(a)
        real(C_DOUBLE), intent(in) :: yplus

        real(C_DOUBLE) :: p

        p = jayatilleke_p(0.71d0/0.85d0)
        a = wall_diffusivity(yplus, 0.71d0, 0.85d0, 180.0d0, p, &
            thermal_yplus(0.71d0, 0.85d0, p))
    end function alpha_air

    ! Relative comparison; the reference values come from mpmath (its own
    ! exp/sqrt), so allow a few-ulp cross-runtime spread.
    subroutine check(name, got, want)
        character(len=*), intent(in) :: name
        real(C_DOUBLE), intent(in) :: got, want

        real(C_DOUBLE), parameter :: TOL = 1.0d-13

        if (abs(got - want) > TOL*max(abs(want), 1.0d0)) then
            print '(A,A,ES24.16,A,ES24.16)', name, ": got ", got, " want ", want
            nfail = nfail + 1
        end if
    end subroutine check

    ! Pr_t(Pe_t) must decrease monotonically from 2 Prt_inf to Prt_inf --
    ! any branch mismatch or series divergence shows up as a bump.
    subroutine check_monotone()
        integer :: i
        real(C_DOUBLE) :: pet, prev, p

        prev = prt_kays(0.0d0, 0.85d0)
        do i = 1, 600
            pet = 10.0d0**(-4.0d0 + 0.02d0*real(i, C_DOUBLE))
            p = prt_kays(pet, 0.85d0)
            if (p > prev + 1.0d-15 .or. p < 0.85d0) then
                print '(A,ES12.4,A,ES24.16,A,ES24.16)', &
                    "monotonicity broken at Pe_t = ", pet, ": ", p, " after ", prev
                nfail = nfail + 1
                return
            end if
            prev = p
        end do
        if (abs(prev - 0.85d0) > 1.0d-8) then
            print *, "Pe_t -> infinity limit not reached: ", prev
            nfail = nfail + 1
        end if
    end subroutine check_monotone

    !------------------------------------------------------------------
    ! C3: the fluid volume fraction of a cut cell
    !------------------------------------------------------------------

    ! Exact values, one per degeneracy regime, each derivable by hand.
    subroutine check_fraction_exact()
        real(C_DOUBLE), parameter :: R2 = 1.4142135623730951d0
        real(C_DOUBLE), parameter :: R3 = 1.7320508075688772d0

        ! (a) the plane parallel to two axes -- a grid-aligned wall, which is
        ! the common case and where the closed form IS clip(1/2 + phi/h).
        call check("vfrac 1D +0.25",  scalar_volume_fraction( 0.25d0, &
            0.0d0, 1.0d0, 0.0d0, 1.0d0, 1.0d0, 1.0d0), 0.75d0)
        call check("vfrac 1D -0.25",  scalar_volume_fraction(-0.25d0, &
            0.0d0, 1.0d0, 0.0d0, 1.0d0, 1.0d0, 1.0d0), 0.25d0)
        ! ...and a normal that is only ALMOST axis-aligned must give the same
        ! answer: this is the degenerate-direction threshold doing its job,
        ! where the full 3D form would divide a cancelled numerator by zero.
        call check("vfrac 1D wobble", scalar_volume_fraction( 0.25d0, &
            1.0d-14, 1.0d0, -1.0d-13, 1.0d0, 1.0d0, 1.0d0), 0.75d0)
        ! non-cubic cell: h_y = 2 halves the fraction's sensitivity to phi
        call check("vfrac 1D h=2",    scalar_volume_fraction( 0.25d0, &
            0.0d0, 1.0d0, 0.0d0, 1.0d0, 2.0d0, 1.0d0), 0.625d0)

        ! (b) parallel to ONE axis: a_1 = a_2 = 1/sqrt(2), the corner wedge is
        ! a triangle, f = s^2/(2 a1 a2) = s^2 with s = 0.3.
        call check("vfrac 2D corner", scalar_volume_fraction(0.3d0 - 0.5d0*R2, &
            1.0d0, 1.0d0, 0.0d0, 1.0d0, 1.0d0, 1.0d0), 0.09d0)
        call check("vfrac 2D middle", scalar_volume_fraction(0.0d0, &
            1.0d0, 1.0d0, 0.0d0, 1.0d0, 1.0d0, 1.0d0), 0.5d0)

        ! (c) the generic plane: a_i = 1/sqrt(3), the corner wedge is a
        ! tetrahedron, f = s^3/(6 a1 a2 a3) = s^3 sqrt(3)/2 with s = 0.2.
        call check("vfrac 3D corner", scalar_volume_fraction(0.2d0 - 0.5d0*R3, &
            1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0), &
            0.2d0**3*0.5d0*R3)
        call check("vfrac 3D middle", scalar_volume_fraction(0.0d0, &
            1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0), 0.5d0)

        ! Clipping: a cell the plane misses is one material, exactly.
        call check("vfrac all fluid", scalar_volume_fraction(3.0d0, &
            1.0d0, 2.0d0, -0.5d0, 1.0d0, 1.0d0, 1.0d0), 1.0d0)
        call check("vfrac all solid", scalar_volume_fraction(-3.0d0, &
            1.0d0, 2.0d0, -0.5d0, 1.0d0, 1.0d0, 1.0d0), 0.0d0)
        ! No usable normal (a medial axis): the pointwise marker, as the rest
        ! of the scheme does there.
        call check("vfrac no normal", scalar_volume_fraction(0.1d0, &
            0.0d0, 1.0d-3, 0.0d0, 1.0d0, 1.0d0, 1.0d0), 1.0d0)
    end subroutine check_fraction_exact

    ! f(phi) + f(-phi) = 1 for every plane: the two materials fill the cell.
    subroutine check_fraction_symmetry()
        integer :: i
        real(C_DOUBLE) :: nx, ny, nz, ph, f1, f2

        do i = 1, 40
            nx = sin(0.7d0*i); ny = cos(1.3d0*i); nz = sin(2.1d0*i + 0.4d0)
            ph = 0.02d0*(i - 20)
            f1 = scalar_volume_fraction( ph, nx, ny, nz, 1.0d0, 0.5d0, 2.0d0)
            f2 = scalar_volume_fraction(-ph, nx, ny, nz, 1.0d0, 0.5d0, 2.0d0)
            if (abs(f1 + f2 - 1.0d0) > 1.0d-14) then
                print '(A,I3,A,ES24.16)', "vfrac symmetry broken at i = ", i, &
                    ": f(phi) + f(-phi) - 1 = ", f1 + f2 - 1.0d0
                nfail = nfail + 1
                return
            end if
        end do
    end subroutine check_fraction_symmetry

    ! BRUTE FORCE: count midpoints of an N^3 lattice on the fluid side. The
    ! quadrature error of a plane cut is O(1/N), so this cannot confirm the
    ! last digits -- it confirms the FORM, which is what an independent check
    ! is for. (The face-area fraction of C2 was validated the same way.)
    subroutine check_fraction_brute()
        integer, parameter :: N = 300
        integer :: c, i, j, k, hit
        real(C_DOUBLE) :: nx, ny, nz, gn, ph, h(3), f, fb, x, y, z, worst

        worst = 0.0d0
        do c = 1, 8
            nx = sin(1.1d0*c); ny = cos(0.6d0*c + 0.2d0); nz = sin(2.3d0*c + 1.0d0)
            gn = sqrt(nx*nx + ny*ny + nz*nz)
            nx = nx/gn; ny = ny/gn; nz = nz/gn
            ph = 0.15d0*(c - 4)
            h = [1.0d0, 0.6d0, 1.7d0]
            if (c > 4) h = [1.0d0, 1.0d0, 1.0d0]
            f = scalar_volume_fraction(ph, nx, ny, nz, h(1), h(2), h(3))
            hit = 0
            do k = 1, N
                z = ((k - 0.5d0)/N - 0.5d0)*h(3)
                do j = 1, N
                    y = ((j - 0.5d0)/N - 0.5d0)*h(2)
                    do i = 1, N
                        x = ((i - 0.5d0)/N - 0.5d0)*h(1)
                        if (ph + nx*x + ny*y + nz*z > 0.0d0) hit = hit + 1
                    end do
                end do
            end do
            fb = real(hit, C_DOUBLE)/real(N, C_DOUBLE)**3
            worst = max(worst, abs(f - fb))
            if (abs(f - fb) > 3.0d-3) then
                print '(A,I2,A,ES14.6,A,ES14.6)', "vfrac brute force case ", c, &
                    ": closed form ", f, " vs sampled ", fb
                nfail = nfail + 1
            end if
        end do
        print '(A,ES11.4,A,I0,A)', " vfrac vs brute force: max deviation ", worst, &
            " (", N, "^3 midpoints; the quadrature error is O(1/N))"
    end subroutine check_fraction_brute

end program test_scalar
