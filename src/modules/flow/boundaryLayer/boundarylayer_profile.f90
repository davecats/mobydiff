module boundarylayer_profile
    ! Blasius similarity initial field for the boundary-layer case: fills the
    ! whole domain with the laminar flat-plate solution at the inflow
    ! displacement thickness, so the run starts on the profile (no laminar
    ! development transient upstream of the trip). Host-only. Same construction
    ! as the tutorial's make_blasius_ic.py: shoot the Blasius ODE (RK4 + secant
    ! on f'(inf)=1), then u = U_inf f'(eta), v the entrainment component, at the
    ! similarity station x + x_v with x_v = Re_theta*theta/beta^2.
    use, intrinsic :: iso_c_binding
    use :: init, only: dns_type, VAR_U, VAR_V, VAR_W
    use :: blocks, only: block_set_type
    implicit none

    private
    public :: initialise_boundarylayer_field

contains

    subroutine initialise_boundarylayer_field(blk, dns, u_inf, theta_in)
        type(block_set_type), intent(inout) :: blk
        type(dns_type), intent(in) :: dns
        real(C_DOUBLE), intent(in) :: u_inf, theta_in

        integer, parameter :: NTAB = 12000
        real(C_DOUBLE), parameter :: ETA_MAX = 12.0d0
        real(C_DOUBLE) :: fTab(0:NTAB), fpTab(0:NTAB), h, beta, disp
        real(C_DOUBLE) :: nu, re_theta, x_v, x, y, xtot, eta, f, fp, vscale
        integer :: i, j, k, b, nx, ny, nz

        h = ETA_MAX/real(NTAB, C_DOUBLE)
        call shoot_blasius(fTab, fpTab, h, NTAB, beta, disp)

        nu = 1.0d0/dns%re
        re_theta = u_inf*theta_in*dns%re          ! U_inf*theta/nu (re = 1/nu)
        x_v = re_theta*theta_in/(beta*beta)

        nx = int(blk%nb(1)); ny = int(blk%nb(2)); nz = int(blk%nb(3))
        do b = 1, int(blk%nBlocks)
            do k = 1, nz
                do j = 1, ny
                    do i = 1, nx
                        ! u at its own staggered coordinate.
                        x = blk%x(i, VAR_U, b); y = blk%y(j, VAR_U, b)
                        xtot = x + x_v
                        eta = beta*y/(beta*sqrt(nu*xtot/u_inf))    ! = y*sqrt(U/(nu*xtot))
                        call blasius_lookup(fTab, fpTab, h, NTAB, disp, eta, f, fp)
                        blk%q(i,j,k,VAR_U,b) = u_inf*fp

                        ! v at its own staggered coordinate.
                        x = blk%x(i, VAR_V, b); y = blk%y(j, VAR_V, b)
                        xtot = x + x_v
                        eta = y*sqrt(u_inf/(nu*xtot))
                        vscale = 0.5d0*sqrt(nu*u_inf/xtot)
                        call blasius_lookup(fTab, fpTab, h, NTAB, disp, eta, f, fp)
                        blk%q(i,j,k,VAR_V,b) = vscale*(eta*fp - f)

                        blk%q(i,j,k,VAR_W,b) = 0.0d0
                    end do
                end do
            end do
        end do
    end subroutine initialise_boundarylayer_field

    ! Solve f''' = -f f''/2 by RK4 + secant shooting on f'(eta_max) = 1;
    ! beta = 2 f''(0) (momentum-thickness constant), disp = displacement
    ! constant (eta_max - f(eta_max), the outer asymptote offset).
    subroutine shoot_blasius(fTab, fpTab, h, ntab, beta, disp)
        real(C_DOUBLE), intent(out) :: fTab(0:), fpTab(0:), beta, disp
        real(C_DOUBLE), intent(in) :: h
        integer, intent(in) :: ntab

        real(C_DOUBLE) :: fpp0, fpp0_prev, err, err_prev, fpp0_next
        integer :: it

        fpp0_prev = 0.33d0
        fpp0 = 0.34d0
        err_prev = integrate_blasius(fTab, fpTab, h, ntab, fpp0_prev)
        do it = 1, 50
            err = integrate_blasius(fTab, fpTab, h, ntab, fpp0)
            if (abs(err) < 1.0d-13 .or. err == err_prev) exit
            fpp0_next = fpp0 - err*(fpp0 - fpp0_prev)/(err - err_prev)
            fpp0_prev = fpp0
            err_prev = err
            fpp0 = fpp0_next
        end do
        beta = 2.0d0*fpp0
        disp = real(ntab, C_DOUBLE)*h - fTab(ntab)
    end subroutine shoot_blasius

    real(C_DOUBLE) function integrate_blasius(fTab, fpTab, h, ntab, fpp0) result(err)
        real(C_DOUBLE), intent(out) :: fTab(0:), fpTab(0:)
        real(C_DOUBLE), intent(in) :: h, fpp0
        integer, intent(in) :: ntab

        real(C_DOUBLE) :: y(3), k1(3), k2(3), k3(3), k4(3)
        integer :: i

        y = [0.0d0, 0.0d0, fpp0]
        fTab(0) = 0.0d0; fpTab(0) = 0.0d0
        do i = 1, ntab
            k1 = blasius_rhs(y)
            k2 = blasius_rhs(y + 0.5d0*h*k1)
            k3 = blasius_rhs(y + 0.5d0*h*k2)
            k4 = blasius_rhs(y + h*k3)
            y = y + h*(k1 + 2.0d0*k2 + 2.0d0*k3 + k4)/6.0d0
            fTab(i) = y(1); fpTab(i) = y(2)
        end do
        err = y(2) - 1.0d0
    end function integrate_blasius

    pure function blasius_rhs(y) result(dy)
        real(C_DOUBLE), intent(in) :: y(3)
        real(C_DOUBLE) :: dy(3)
        dy = [y(2), y(3), -0.5d0*y(1)*y(3)]
    end function blasius_rhs

    subroutine blasius_lookup(fTab, fpTab, h, ntab, disp, eta, f, fp)
        real(C_DOUBLE), intent(in) :: fTab(0:), fpTab(0:), h, disp, eta
        integer, intent(in) :: ntab
        real(C_DOUBLE), intent(out) :: f, fp

        real(C_DOUBLE) :: s
        integer :: i

        s = eta/h
        if (s >= real(ntab, C_DOUBLE)) then
            f = eta - disp        ! outer asymptote
            fp = 1.0d0
            return
        end if
        i = int(s)
        s = s - real(i, C_DOUBLE)
        f = (1.0d0 - s)*fTab(i) + s*fTab(i+1)
        fp = (1.0d0 - s)*fpTab(i) + s*fpTab(i+1)
    end subroutine blasius_lookup

end module boundarylayer_profile
