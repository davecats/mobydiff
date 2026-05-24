!--------------------------!
!                          !
!       Time-stepper       !
!          module          !
!                          !
!--------------------------! 
! 
! authors: Dr.-Ing. Davide Gatti
!          B.Sc. Ahmet Cumhur
! 
! date:    28.04.26
! 

module step
    use, intrinsic :: iso_c_binding
    use :: init, only: dns_type, grid_type, field_type, VAR_U, VAR_V, VAR_W, VAR_P
    use :: ibmm, only: ibm_type
    use :: boundary, only: boundary_type
    implicit none

    real(C_DOUBLE), parameter :: rk_alpha(3) = [64.0d0/120.0d0,  50.0d0/120.0d0,  90.0d0/120.0d0]
    real(C_DOUBLE), parameter :: rk_beta(3)  = [ 0.0d0,         -34.0d0/120.0d0, -50.0d0/120.0d0]
    real(C_DOUBLE), parameter :: rk_gamma(3) = [64.0d0/120.0d0,  16.0d0/120.0d0,  40.0d0/120.0d0]

contains

    subroutine momentum(f, dns, g, dt_alpha, dt_beta, dt_gamma, ibm, bc)
        type(field_type), intent(inout) :: f
        type(dns_type),   intent(in)    :: dns
        type(grid_type),  intent(in)    :: g
        real(C_DOUBLE),   intent(in)    :: dt_alpha, dt_beta, dt_gamma
        type(ibm_type),   intent(in)    :: ibm
        type(boundary_type),   intent(in)  :: bc

        integer :: i,j,k,ip,im,kp,km,jp,jm
        integer :: nx, ny, nz, uStartX, vStartY, wStartZ 

        real(C_DOUBLE) :: diff_ux,diff_uy,diff_uz
        real(C_DOUBLE) :: diff_vx,diff_vy,diff_vz
        real(C_DOUBLE) :: diff_wx,diff_wy,diff_wz

        real(C_DOUBLE) :: uu_p,uu_m,uv_p,uv_m,uw_p,uw_m
        real(C_DOUBLE) :: vu_p,vu_m,vv_p,vv_m,vw_p,vw_m
        real(C_DOUBLE) :: wu_p,wu_m,ww_p,ww_m,wv_p,wv_m

        real(C_DOUBLE) :: dpx,dpy,dpz,rhsu,rhsv,rhsw
        real(C_DOUBLE) :: mu_u, mu_v, mu_w
        real(C_DOUBLE) :: ire
        real(C_DOUBLE) :: forcing(1:3)

        nx = int(dns%localSize(1,2))
        ny = int(dns%localSize(2,2))
        nz = int(dns%localSize(3,2))
        ire = 1.0d0/dns%re
        forcing = dns%forcing
        uStartX = 1
        vStartY = 1
        wStartZ = 1
        if (.not. bc%isPeriodic(1) .and. dns%localSize(1,0) == 1_C_INT) uStartX = 2
        if (.not. bc%isPeriodic(2) .and. dns%localSize(2,0) == 1_C_INT) vStartY = 2
        if (.not. bc%isPeriodic(3) .and. dns%localSize(3,0) == 1_C_INT) wStartZ = 2

        ! Predictor for all staggered velocity components.
        !$omp target teams distribute parallel do collapse(3) &
        !$omp& map(to: dt_alpha, dt_beta, dt_gamma, uStartX, vStartY, wStartZ, &
        !$omp& ire, forcing(1:3), &
        !$omp& g%d1x, g%d1y, g%d1z, &
        !$omp& g%lapXm, g%lapX0, g%lapXp, g%lapYm, g%lapY0, g%lapYp, &
        !$omp& g%lapZm, g%lapZ0, g%lapZp, f%q, ibm%coef) &
        !$omp& map(tofrom: f%qs, f%oldrhs) &
        !$omp& private(i,j,k,ip,im,jp,jm,kp,km,uu_p,uu_m,uv_p,uv_m,uw_p,uw_m, &
        !$omp& vu_p,vu_m,vv_p,vv_m,vw_p,vw_m,wu_p,wu_m,ww_p,ww_m,wv_p,wv_m, &
        !$omp& diff_ux,diff_uy,diff_uz,diff_vx,diff_vy,diff_vz,diff_wx,diff_wy,diff_wz, &
        !$omp& dpx,dpy,dpz,rhsu,rhsv,rhsw,mu_u,mu_v,mu_w)
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    ip = i+1
                    im = i-1
                    jp = j+1
                    jm = j-1
                    kp = k+1
                    km = k-1

                    if (i >= uStartX) then
                        uu_p = (f%q(i,j,k,VAR_U) + f%q(ip,j,k,VAR_U))**2
                        uu_m = (f%q(im,j,k,VAR_U) + f%q(i,j,k,VAR_U))**2

                        uv_p = (f%q(i,j,k,VAR_U) + f%q(i,jp,k,VAR_U)) &
                             * (f%q(im,jp,k,VAR_V) + f%q(i,jp,k,VAR_V))
                        uv_m = (f%q(i,jm,k,VAR_U) + f%q(i,j,k,VAR_U)) &
                             * (f%q(im,j,k,VAR_V) + f%q(i,j,k,VAR_V))

                        uw_p = (f%q(i,j,k,VAR_U) + f%q(i,j,kp,VAR_U)) &
                             * (f%q(im,j,kp,VAR_W) + f%q(i,j,kp,VAR_W))
                        uw_m = (f%q(i,j,km,VAR_U) + f%q(i,j,k,VAR_U)) &
                             * (f%q(im,j,k,VAR_W) + f%q(i,j,k,VAR_W))

                        diff_ux = g%lapXm(i,VAR_U)*f%q(im,j,k,VAR_U) &
                                + g%lapX0(i,VAR_U)*f%q(i,j,k,VAR_U) &
                                + g%lapXp(i,VAR_U)*f%q(ip,j,k,VAR_U)
                        diff_uy = g%lapYm(j,VAR_U)*f%q(i,jm,k,VAR_U) &
                                + g%lapY0(j,VAR_U)*f%q(i,j,k,VAR_U) &
                                + g%lapYp(j,VAR_U)*f%q(i,jp,k,VAR_U)
                        diff_uz = g%lapZm(k,VAR_U)*f%q(i,j,km,VAR_U) &
                                + g%lapZ0(k,VAR_U)*f%q(i,j,k,VAR_U) &
                                + g%lapZp(k,VAR_U)*f%q(i,j,kp,VAR_U)

                        dpx = (f%q(i,j,k,VAR_P)-f%q(im,j,k,VAR_P))*g%d1x(i,VAR_U)

                        rhsu = ( &
                            -0.25d0*( (uu_p-uu_m)*g%d1x(i,VAR_U) &
                                     +(uv_p-uv_m)*g%d1y(j,VAR_U) &
                                     +(uw_p-uw_m)*g%d1z(k,VAR_U)) &
                            + forcing(VAR_U) &
                            + ire*(diff_ux + diff_uy + diff_uz) )

                        f%qs(i,j,k,VAR_U) = f%q(i,j,k,VAR_U) + dt_alpha*rhsu &
                            + dt_beta*f%oldrhs(i,j,k,VAR_U) - dt_gamma*dpx

                        mu_u = 1.0d0/(1.0d0 + dt_gamma*ibm%coef(i,j,k,VAR_U))
                        f%qs(i,j,k,VAR_U) = f%qs(i,j,k,VAR_U)*mu_u

                        f%oldrhs(i,j,k,VAR_U) = rhsu
                    end if

                    if (j >= vStartY) then
                        vu_p = (f%q(i,j,k,VAR_V) + f%q(ip,j,k,VAR_V)) &
                             * (f%q(ip,jm,k,VAR_U) + f%q(ip,j,k,VAR_U))
                        vu_m = (f%q(im,j,k,VAR_V) + f%q(i,j,k,VAR_V)) &
                             * (f%q(i,jm,k,VAR_U) + f%q(i,j,k,VAR_U))

                        vv_p = (f%q(i,j,k,VAR_V) + f%q(i,jp,k,VAR_V))**2
                        vv_m = (f%q(i,jm,k,VAR_V) + f%q(i,j,k,VAR_V))**2

                        vw_p = (f%q(i,j,k,VAR_V) + f%q(i,j,kp,VAR_V)) &
                             * (f%q(i,jm,kp,VAR_W) + f%q(i,j,kp,VAR_W))
                        vw_m = (f%q(i,j,km,VAR_V) + f%q(i,j,k,VAR_V)) &
                             * (f%q(i,jm,k,VAR_W) + f%q(i,j,k,VAR_W))

                        diff_vx = g%lapXm(i,VAR_V)*f%q(im,j,k,VAR_V) &
                                + g%lapX0(i,VAR_V)*f%q(i,j,k,VAR_V) &
                                + g%lapXp(i,VAR_V)*f%q(ip,j,k,VAR_V)
                        diff_vy = g%lapYm(j,VAR_V)*f%q(i,jm,k,VAR_V) &
                                + g%lapY0(j,VAR_V)*f%q(i,j,k,VAR_V) &
                                + g%lapYp(j,VAR_V)*f%q(i,jp,k,VAR_V)
                        diff_vz = g%lapZm(k,VAR_V)*f%q(i,j,km,VAR_V) &
                                + g%lapZ0(k,VAR_V)*f%q(i,j,k,VAR_V) &
                                + g%lapZp(k,VAR_V)*f%q(i,j,kp,VAR_V)

                        dpy = (f%q(i,j,k,VAR_P)-f%q(i,jm,k,VAR_P))*g%d1y(j,VAR_V)

                        rhsv = ( &
                            -0.25d0*((vu_p-vu_m)*g%d1x(i,VAR_V) &
                                    +(vv_p-vv_m)*g%d1y(j,VAR_V) &
                                    +(vw_p-vw_m)*g%d1z(k,VAR_V)) &
                            + forcing(VAR_V) &
                            + ire*(diff_vx + diff_vy + diff_vz) )

                        f%qs(i,j,k,VAR_V) = f%q(i,j,k,VAR_V) + dt_alpha*rhsv &
                            + dt_beta*f%oldrhs(i,j,k,VAR_V) - dt_gamma*dpy

                        mu_v = 1.0d0/(1.0d0 + dt_gamma*ibm%coef(i,j,k,VAR_V))
                        f%qs(i,j,k,VAR_V) = f%qs(i,j,k,VAR_V)*mu_v

                        f%oldrhs(i,j,k,VAR_V) = rhsv
                    end if

                    if (k >= wStartZ) then
                        wu_p = (f%q(i,j,k,VAR_W) + f%q(ip,j,k,VAR_W)) &
                             * (f%q(ip,j,km,VAR_U) + f%q(ip,j,k,VAR_U))
                        wu_m = (f%q(im,j,k,VAR_W) + f%q(i,j,k,VAR_W)) &
                             * (f%q(i,j,km,VAR_U) + f%q(i,j,k,VAR_U))

                        ww_p = (f%q(i,j,k,VAR_W) + f%q(i,j,kp,VAR_W))**2
                        ww_m = (f%q(i,j,km,VAR_W) + f%q(i,j,k,VAR_W))**2

                        wv_p = (f%q(i,j,k,VAR_W) + f%q(i,jp,k,VAR_W)) &
                             * (f%q(i,jp,km,VAR_V) + f%q(i,jp,k,VAR_V))
                        wv_m = (f%q(i,jm,k,VAR_W) + f%q(i,j,k,VAR_W)) &
                             * (f%q(i,j,km,VAR_V) + f%q(i,j,k,VAR_V))

                        diff_wx = g%lapXm(i,VAR_W)*f%q(im,j,k,VAR_W) &
                                + g%lapX0(i,VAR_W)*f%q(i,j,k,VAR_W) &
                                + g%lapXp(i,VAR_W)*f%q(ip,j,k,VAR_W)
                        diff_wy = g%lapYm(j,VAR_W)*f%q(i,jm,k,VAR_W) &
                                + g%lapY0(j,VAR_W)*f%q(i,j,k,VAR_W) &
                                + g%lapYp(j,VAR_W)*f%q(i,jp,k,VAR_W)
                        diff_wz = g%lapZm(k,VAR_W)*f%q(i,j,km,VAR_W) &
                                + g%lapZ0(k,VAR_W)*f%q(i,j,k,VAR_W) &
                                + g%lapZp(k,VAR_W)*f%q(i,j,kp,VAR_W)

                        dpz = (f%q(i,j,k,VAR_P)-f%q(i,j,km,VAR_P))*g%d1z(k,VAR_W)

                        rhsw = ( &
                            -0.25d0*((wu_p-wu_m)*g%d1x(i,VAR_W) &
                                     (wv_p-wv_m)*g%d1y(j,VAR_W) &
                                     (ww_p-ww_m)*g%d1z(k,VAR_W)) &
                            + forcing(VAR_W) &
                            + ire*(diff_wx + diff_wy + diff_wz) )

                        f%qs(i,j,k,VAR_W) = f%q(i,j,k,VAR_W) + dt_alpha*rhsw &
                            + dt_beta*f%oldrhs(i,j,k,VAR_W) - dt_gamma*dpz

                        mu_w = 1.0d0/(1.0d0 + dt_gamma*ibm%coef(i,j,k,VAR_W))
                        f%qs(i,j,k,VAR_W) = f%qs(i,j,k,VAR_W)*mu_w

                        f%oldrhs(i,j,k,VAR_W) = rhsw
                    end if

                end do
            end do
        end do
        !$omp end target teams distribute parallel do

        !$omp target teams distribute parallel do collapse(3) &
        !$omp& map(to: f%qs) &
        !$omp& map(tofrom: f%q) &
        !$omp& private(i,j,k)
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    f%q(i,j,k,VAR_U) = f%qs(i,j,k,VAR_U)
                    f%q(i,j,k,VAR_V) = f%qs(i,j,k,VAR_V)
                    f%q(i,j,k,VAR_W) = f%qs(i,j,k,VAR_W)
                end do
            end do
        end do
        !$omp end target teams distribute parallel do

    end subroutine momentum


    real(C_DOUBLE) function get_cfl(f, dns, g)
        type(field_type), intent(inout) :: f
        type(dns_type),   intent(in)    :: dns
        type(grid_type),  intent(in)    :: g

        integer :: i,j,k
        integer :: nx, ny, nz

        nx = int(dns%localSize(1,2))
        ny = int(dns%localSize(2,2))
        nz = int(dns%localSize(3,2))
        get_cfl = 0.0d0

        !$omp target teams distribute parallel do collapse(3) reduction(max:get_cfl) &
        !$omp& map(to: g%d1x, g%d1y, g%d1z, f%q) &
        !$omp& private(i,j,k)
        do i = 0, nx+1
            do j = 0, ny+1
                do k = 0, nz+1
                    get_cfl = max(get_cfl, abs(f%q(i,j,k,VAR_U)*g%d1x(i,VAR_U)))
                    get_cfl = max(get_cfl, abs(f%q(i,j,k,VAR_V)*g%d1y(j,VAR_V)))
                    get_cfl = max(get_cfl, abs(f%q(i,j,k,VAR_W)*g%d1z(k,VAR_W)))
                end do
            end do
        end do
        !$omp end target teams distribute parallel do
    end function get_cfl

end module step
