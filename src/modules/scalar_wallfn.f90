! The PLAIN scalar family's face diffusivity (split out of scalar.f90 on
! 2026-10-21): molecular + eddy (prt_kays, the Kays-Crawford Pr_t) + the
! S5a thermal wall function (Jayatilleke), delivered as a wall-cell eddy
! diffusivity exactly as rans.f90 delivers the wall shear. ONE definition
! shared by the transport kernel (scalar.f90) and the statistics
! (scalar_stats.f90), so a wall row's flux IS the flux the kernel applied.
module scalar_wallfn
    use, intrinsic :: iso_c_binding
    use :: init, only: WF_KAPPA, WF_E
    use :: scalar_types, only: SC_PRT_KAYS
    implicit none
    private
    public :: eddy_diffusivity, jayatilleke_p, plain_face_diffusivity, prt_kays, thermal_yplus
    public :: wall_diffusivity, wall_face_diffusivity

contains

    ! Kays-Crawford turbulent Prandtl number (Kays 1994; Weigand, Ferguson &
    ! Crawford 1997), a pointwise function of the turbulent Peclet number
    ! Pe_t = (nu_t/nu) Pr:
    !
    !   1/Pr_t = 1/(2 Prt_inf) + C Pe_t/sqrt(Prt_inf)
    !            - (C Pe_t)^2 [1 - exp(-1/(C Pe_t sqrt(Prt_inf)))],   C = 0.3
    !
    ! Limits: Pe_t -> 0 gives Pr_t = 2 Prt_inf (the molecular limit, inert --
    ! it multiplies nu_t = 0); Pe_t -> infinity gives Pr_t -> Prt_inf, which
    ! is the per-scalar [scalar.N] prt (0.85 by default), so `constant` and
    ! `kays` share one configuration key.
    !
    ! NUMERICS: written with x = 1/(C Pe_t sqrt(Prt_inf)), the bracket is
    ! (C Pe_t)^2 [1 - exp(-x)] = a/b - 1/(2b^2) + ... and the leading a/b
    ! cancels the second term of the sum exactly -- catastrophically so once
    ! Pe_t is large. The small-x branch therefore evaluates the analytically
    ! equivalent series
    !
    !   1/Pr_t = (1/Prt_inf) [ 1 - x/3! + x^2/4! - x^3/5! + ... ]
    !
    ! which is also where the Pe_t -> infinity limit is exact by inspection.
    ! The two branches agree to round-off at the x = 1/2 crossover.
    real(C_DOUBLE) function prt_kays(pet, prtinf) result(prt)
        !$omp declare target
        real(C_DOUBLE), intent(in) :: pet      ! turbulent Peclet number
        real(C_DOUBLE), intent(in) :: prtinf   ! Pr_t at Pe_t -> infinity

        real(C_DOUBLE), parameter :: CKC = 0.3d0
        real(C_DOUBLE) :: a, b, x, u, sum, invprt
        integer :: m

        b = sqrt(prtinf)
        a = CKC*max(pet, 0.0d0)
        if (a*b < 1.0d-300) then
            prt = 2.0d0*prtinf              ! the molecular limit
            return
        end if
        x = 1.0d0/(a*b)

        if (x < 0.5d0) then
            sum = 1.0d0
            u = -x/6.0d0                    ! the m = 1 term, -x/3!
            do m = 1, 16                    ! machine precision for x < 1/2
                sum = sum + u
                u = -u*x/real(m + 3, C_DOUBLE)
            end do
            invprt = sum/prtinf
        else
            invprt = 0.5d0/prtinf + a/b - a*a*(1.0d0 - exp(-x))
        end if

        prt = 1.0d0/invprt
    end function prt_kays

    ! Face eddy diffusivity nu_t/Pr_t for one scalar. Called per face inside
    ! the transport kernel, so it is a declare-target pure function.
    real(C_DOUBLE) function eddy_diffusivity(nutFace, pr, prt, prtModel, re) result(d)
        !$omp declare target
        real(C_DOUBLE), intent(in) :: nutFace, pr, prt, re
        integer(C_INT), intent(in) :: prtModel

        real(C_DOUBLE) :: nt

        nt = max(nutFace, 0.0d0)
        if (prtModel == SC_PRT_KAYS) then
            ! Pe_t = (nu_t/nu) Pr = nu_t Re Pr.
            d = nt/prt_kays(nt*re*pr, prt)
        else
            d = nt/prt
        end if
    end function eddy_diffusivity

    ! Jayatilleke's P (Jayatilleke 1969; the OpenFOAM Psmooth form): the
    ! conductive sublayer's extra resistance relative to the momentum one.
    ! P = 0 at Pr = Pr_t (the Reynolds analogy), positive for Pr > Pr_t.
    pure real(C_DOUBLE) function jayatilleke_p(prat) result(p)
        real(C_DOUBLE), intent(in) :: prat        ! Pr/Pr_t

        p = 9.24d0*(prat**0.75d0 - 1.0d0)*(1.0d0 + 0.28d0*exp(-0.007d0*prat))
    end function jayatilleke_p

    ! The thermal sublayer thickness y+_T: where the conductive branch
    ! Pr y+ meets the log branch. f(y) = Pr y - Pr_t[ln(E y)/kappa + P] has
    ! its minimum at y* = Pr_t/(kappa Pr) and, whenever a conductive
    ! sublayer exists at all, exactly one root beyond it -- the physical
    ! crossing (the root BELOW y* is spurious: there the log branch is the
    ! larger of the two). Bisection, host-side, once per scalar at init;
    ! OpenFOAM iterates the same fixed point.
    !
    ! Degenerate case f(y*) >= 0 (Pr well below Pr_t -- liquid metals): the
    ! log branch never rises above the conductive one, i.e. the thermal
    ! sublayer has no extent. y+_T = y* is then the switch, which keeps the
    ! log branch positive (it is increasing beyond y*) and reduces to
    ! conduction below it.
    real(C_DOUBLE) function thermal_yplus(pr, prt, p) result(ypt)
        real(C_DOUBLE), intent(in) :: pr, prt, p

        real(C_DOUBLE) :: ylo, yhi, ym, fm
        integer :: it

        ylo = prt/(WF_KAPPA*pr)                   ! the minimum of f
        if (pr*ylo - prt*(log(WF_E*ylo)/WF_KAPPA + p) >= 0.0d0) then
            ypt = ylo
            return
        end if
        yhi = ylo
        do it = 1, 200                            ! bracket the root above
            yhi = 2.0d0*yhi
            if (pr*yhi - prt*(log(WF_E*yhi)/WF_KAPPA + p) > 0.0d0) exit
        end do
        do it = 1, 200                            ! ~1e-60 of the bracket
            ym = 0.5d0*(ylo + yhi)
            fm = pr*ym - prt*(log(WF_E*ym)/WF_KAPPA + p)
            if (fm < 0.0d0) then
                ylo = ym
            else
                yhi = ym
            end if
        end do
        ypt = 0.5d0*(ylo + yhi)
    end function thermal_yplus

    ! Wall-cell eddy diffusivity of the thermal wall function: the total
    ! nu y+/theta+ minus the molecular nu/Pr, so the transport kernel adds
    ! it to its unchanged molecular term. EXACTLY zero on the conductive
    ! branch (theta+ = Pr y+ there), and continuous at y+_T by the
    ! definition of y+_T -- the structure of T3's nut_wall_value.
    real(C_DOUBLE) function wall_diffusivity(yplus, pr, prt, re, p, ypt) result(a)
        !$omp declare target
        real(C_DOUBLE), intent(in) :: yplus, pr, prt, re, p, ypt

        real(C_DOUBLE) :: tp

        if (yplus < ypt) then
            a = 0.0d0
        else
            tp = prt*(log(WF_E*yplus)/WF_KAPPA + p)
            a = max((yplus/tp - 1.0d0/pr)/re, 0.0d0)
        end if
    end function wall_diffusivity

    ! Cell eddy diffusivity for one scalar under wall functions: the thermal
    ! wall function at a wall cell (y+ > 0 marks them, rans_wall_yplus), the
    ! ordinary nu_t/Pr_t everywhere else.
    real(C_DOUBLE) function cell_diffusivity(nutc, yplus, pr, prt, prtModel, re, p, ypt) result(a)
        !$omp declare target
        real(C_DOUBLE), intent(in) :: nutc, yplus, pr, prt, re, p, ypt
        integer(C_INT), intent(in) :: prtModel

        if (yplus > 0.0d0) then
            a = wall_diffusivity(yplus, pr, prt, re, p, ypt)
        else
            a = eddy_diffusivity(nutc, pr, prt, prtModel, re)
        end if
    end function cell_diffusivity

    ! Face eddy diffusivity under wall functions: the 1/2 (L + R) average of
    ! the two CELL diffusivities. (Resolved mode averages nu_t first and
    ! divides by Pr_t afterwards; that is not available here, because a wall
    ! cell's diffusivity is not a function of its nu_t at all. The two agree
    ! identically for a constant Pr_t away from wall cells, and differ only
    ! in rounding -- but the wall-function branch is a separate code path,
    ! so resolved runs keep their arithmetic byte for byte.)
    real(C_DOUBLE) function wall_face_diffusivity(nutL, nutR, ypL, ypR, &
            pr, prt, prtModel, re, p, ypt) result(d)
        !$omp declare target
        real(C_DOUBLE), intent(in) :: nutL, nutR, ypL, ypR, pr, prt, re, p, ypt
        integer(C_INT), intent(in) :: prtModel

        d = 0.5d0*(cell_diffusivity(nutL, ypL, pr, prt, prtModel, re, p, ypt) &
                 + cell_diffusivity(nutR, ypR, pr, prt, prtModel, re, p, ypt))
    end function wall_face_diffusivity

    ! THE PLAIN FAMILY'S FACE DIFFUSIVITY (dirichlet / adiabatic scalars, and
    ! every scalar of a body-free run): the molecular dm plus, under a
    ! turbulence model, the eddy part -- the thermal wall function's (S5a)
    ! or nu_t,face/Pr_t on the 1/2 (L + R) face eddy viscosity. ONE
    ! definition for the transport kernel and the two statistics kernels of
    ! scalar_stats.f90, so the S4 design statement "the statistics report
    ! the flux the kernel APPLIED" holds by construction, not by copy (it
    ! was three copies until 2026-10-20). The two cells' values come in
    ! LOW -> HIGH index order and the face average is formed here in that
    ! order, which is the order the kernels used. With useNut off the face
    ! keeps dm exactly (the S1 arithmetic).
    real(C_DOUBLE) function plain_face_diffusivity(dm, useNut, wallfn, nutL, nutR, &
            ypL, ypR, pr, prt, prtModel, re, p, ypt) result(d)
        !$omp declare target
        real(C_DOUBLE), intent(in) :: dm, nutL, nutR, ypL, ypR, pr, prt, re, p, ypt
        logical, intent(in) :: useNut, wallfn
        integer(C_INT), intent(in) :: prtModel

        if (useNut .and. wallfn) then
            ! Thermal wall function: the face value is built from the two
            ! CELL diffusivities, because a wall cell's is not a function
            ! of its nu_t at all.
            d = dm + wall_face_diffusivity(nutL, nutR, ypL, ypR, pr, prt, prtModel, re, p, ypt)
        else if (useNut) then
            d = dm + eddy_diffusivity(0.5d0*(nutL + nutR), pr, prt, prtModel, re)
        else
            d = dm
        end if
    end function plain_face_diffusivity

end module scalar_wallfn
