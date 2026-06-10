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
    use :: init, only: dns_type, grid_type, field_type, VAR_U, VAR_V, VAR_W, VAR_P, &
        CFL_COURANT, CFL_PECLET, NCFL
    use :: ibmm, only: ibm_type
    use :: boundary, only: boundary_type
    use :: comm, only: comm_type, comm_allreduce_max
    use :: les_model, only: les_type, les_is_enabled, les_profile_type, &
        les_wall_seconds, add_les_profile, LES_PROF_SGS
    implicit none

    real(C_DOUBLE), parameter :: rk_alpha(3) = [64.0d0/120.0d0,  50.0d0/120.0d0,  90.0d0/120.0d0]
    real(C_DOUBLE), parameter :: rk_beta(3)  = [ 0.0d0,         -34.0d0/120.0d0, -50.0d0/120.0d0]
    real(C_DOUBLE), parameter :: rk_gamma(3) = [64.0d0/120.0d0,  16.0d0/120.0d0,  40.0d0/120.0d0]

contains

    subroutine precompute_peclet_rate(dns, g, c)
        type(dns_type), intent(inout) :: dns
        type(grid_type), intent(in) :: g
        type(comm_type), intent(in) :: c

        integer :: i, nx, ny, nz
        real(C_DOUBLE) :: ire, local_rate(1)

        nx = int(dns%localSize(1,2))
        ny = int(dns%localSize(2,2))
        nz = int(dns%localSize(3,2))
        ire = 1.0d0/dns%re

        local_rate = 0.0d0
        do i = 0, nx+1
            local_rate(1) = max(local_rate(1), ire*g%d1x(i,VAR_P)**2)
        end do
        do i = 0, ny+1
            local_rate(1) = max(local_rate(1), ire*g%d1y(i,VAR_P)**2)
        end do
        do i = 0, nz+1
            local_rate(1) = max(local_rate(1), ire*g%d1z(i,VAR_P)**2)
        end do

        call comm_allreduce_max(c, local_rate)
        dns%peclet_rate = local_rate(1)
    end subroutine precompute_peclet_rate

    subroutine momentum(f, dns, g, dt_alpha, dt_beta, dt_gamma, ibm, bc, les, les_prof)
        type(field_type), intent(inout) :: f
        type(dns_type),   intent(in)    :: dns
        type(grid_type),  intent(in)    :: g
        real(C_DOUBLE),   intent(in)    :: dt_alpha, dt_beta, dt_gamma
        type(ibm_type),   intent(in)    :: ibm
        type(boundary_type),   intent(in)  :: bc
        type(les_type),   intent(in), optional :: les
        type(les_profile_type), intent(inout), optional :: les_prof

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
        logical :: use_les

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
        use_les = .false.
        if (present(les)) use_les = les_is_enabled(les) .and. allocated(les%nut)

        ! Predictor for all staggered velocity components.
        !$omp target teams distribute parallel do collapse(3) &
        !$omp& map(to: dt_alpha, dt_beta, dt_gamma, uStartX, vStartY, wStartZ, &
        !$omp& ire, forcing(1:3), &
        !$omp& g%d1x, g%d1y, g%d1z, &
        !$omp& g%lapXm, g%lapX0, g%lapXp, g%lapYm, g%lapY0, g%lapYp, &
        !$omp& g%lapZm, g%lapZ0, g%lapZp, f%q, ibm%mu) &
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

                        mu_u = ibm%mu(i,j,k,VAR_U)
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

                        mu_v = ibm%mu(i,j,k,VAR_V)
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
                            -0.25d0*( (wu_p-wu_m)*g%d1x(i,VAR_W) &
                                     +(wv_p-wv_m)*g%d1y(j,VAR_W) &
                                     +(ww_p-ww_m)*g%d1z(k,VAR_W)) &
                            + forcing(VAR_W) &
                            + ire*(diff_wx + diff_wy + diff_wz) )

                        f%qs(i,j,k,VAR_W) = f%q(i,j,k,VAR_W) + dt_alpha*rhsw &
                            + dt_beta*f%oldrhs(i,j,k,VAR_W) - dt_gamma*dpz

                        mu_w = ibm%mu(i,j,k,VAR_W)
                        f%qs(i,j,k,VAR_W) = f%qs(i,j,k,VAR_W)*mu_w

                        f%oldrhs(i,j,k,VAR_W) = rhsw
                    end if
                end do
            end do
        end do
        !$omp end target teams distribute parallel do

        if (present(les)) then
            if (use_les) then
                if (present(les_prof)) then
                    call add_les_momentum_correction(f, dns, g, dt_alpha, ibm, les, &
                        uStartX, vStartY, wStartZ, les_prof)
                else
                    call add_les_momentum_correction(f, dns, g, dt_alpha, ibm, les, &
                        uStartX, vStartY, wStartZ)
                end if
            end if
        end if

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


    subroutine add_les_momentum_correction(f, dns, g, dt_alpha, ibm, les, uStartX, vStartY, wStartZ, les_prof)
        type(field_type), intent(inout) :: f
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        real(C_DOUBLE), intent(in) :: dt_alpha
        type(ibm_type), intent(in) :: ibm
        type(les_type), intent(in) :: les
        integer, intent(in) :: uStartX, vStartY, wStartZ
        type(les_profile_type), intent(inout), optional :: les_prof

        integer :: i, j, k, ip, im, jp, jm, kp, km
        integer :: nx, ny, nz
        real(C_DOUBLE) :: sgs_u, sgs_v, sgs_w
        real(C_DOUBLE) :: tau_xp, tau_xm, tau_yp, tau_ym, tau_zp, tau_zm
        real(C_DOUBLE) :: nut_edge, nut0, nut1, wx, wy, mu_u, mu_v, mu_w
        real(C_DOUBLE) :: profile_start

        nx = int(dns%localSize(1,2))
        ny = int(dns%localSize(2,2))
        nz = int(dns%localSize(3,2))

        if (present(les_prof)) profile_start = les_wall_seconds()

        !$omp target teams distribute parallel do collapse(3) &
        !$omp& map(to: dt_alpha, uStartX, vStartY, wStartZ, &
        !$omp& g%d1x, g%d1y, g%d1z, f%q, ibm%mu, les%nut, &
        !$omp& les%u_from_p_x, les%v_from_p_y, les%w_from_p_z, &
        !$omp& les%inv_dx, les%inv_dy, les%inv_dz) &
        !$omp& map(tofrom: f%qs, f%oldrhs) &
        !$omp& private(i,j,k,ip,im,jp,jm,kp,km,sgs_u,sgs_v,sgs_w,tau_xp,tau_xm, &
        !$omp& tau_yp,tau_ym,tau_zp,tau_zm,nut_edge,nut0,nut1,wx,wy,mu_u,mu_v,mu_w)
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    ip = i + 1
                    im = i - 1
                    jp = j + 1
                    jm = j - 1
                    kp = k + 1
                    km = k - 1

                    if (i >= uStartX) then
                        tau_xp = 2.0d0*les%nut(i,j,k) &
                               * (f%q(ip,j,k,VAR_U) - f%q(i,j,k,VAR_U))*g%d1x(i,VAR_P)
                        tau_xm = 2.0d0*les%nut(im,j,k) &
                               * (f%q(i,j,k,VAR_U) - f%q(im,j,k,VAR_U))*g%d1x(im,VAR_P)

                        wx = les%u_from_p_x(i)
                        wy = les%v_from_p_y(jp)
                        nut0 = (1.0d0 - wx)*les%nut(im,j,k) + wx*les%nut(i,j,k)
                        nut1 = (1.0d0 - wx)*les%nut(im,jp,k) + wx*les%nut(i,jp,k)
                        nut_edge = (1.0d0 - wy)*nut0 + wy*nut1
                        tau_yp = nut_edge*( &
                            (f%q(i,jp,k,VAR_U) - f%q(i,j,k,VAR_U))*les%inv_dy(jp,VAR_U) &
                          + (f%q(i,jp,k,VAR_V) - f%q(im,jp,k,VAR_V))*les%inv_dx(i,VAR_V) )

                        wx = les%u_from_p_x(i)
                        wy = les%v_from_p_y(j)
                        nut0 = (1.0d0 - wx)*les%nut(im,jm,k) + wx*les%nut(i,jm,k)
                        nut1 = (1.0d0 - wx)*les%nut(im,j,k) + wx*les%nut(i,j,k)
                        nut_edge = (1.0d0 - wy)*nut0 + wy*nut1
                        tau_ym = nut_edge*( &
                            (f%q(i,j,k,VAR_U) - f%q(i,jm,k,VAR_U))*les%inv_dy(j,VAR_U) &
                          + (f%q(i,j,k,VAR_V) - f%q(im,j,k,VAR_V))*les%inv_dx(i,VAR_V) )

                        wx = les%u_from_p_x(i)
                        wy = les%w_from_p_z(kp)
                        nut0 = (1.0d0 - wx)*les%nut(im,j,k) + wx*les%nut(i,j,k)
                        nut1 = (1.0d0 - wx)*les%nut(im,j,kp) + wx*les%nut(i,j,kp)
                        nut_edge = (1.0d0 - wy)*nut0 + wy*nut1
                        tau_zp = nut_edge*( &
                            (f%q(i,j,kp,VAR_U) - f%q(i,j,k,VAR_U))*les%inv_dz(kp,VAR_U) &
                          + (f%q(i,j,kp,VAR_W) - f%q(im,j,kp,VAR_W))*les%inv_dx(i,VAR_W) )

                        wx = les%u_from_p_x(i)
                        wy = les%w_from_p_z(k)
                        nut0 = (1.0d0 - wx)*les%nut(im,j,km) + wx*les%nut(i,j,km)
                        nut1 = (1.0d0 - wx)*les%nut(im,j,k) + wx*les%nut(i,j,k)
                        nut_edge = (1.0d0 - wy)*nut0 + wy*nut1
                        tau_zm = nut_edge*( &
                            (f%q(i,j,k,VAR_U) - f%q(i,j,km,VAR_U))*les%inv_dz(k,VAR_U) &
                          + (f%q(i,j,k,VAR_W) - f%q(im,j,k,VAR_W))*les%inv_dx(i,VAR_W) )

                        sgs_u = (tau_xp - tau_xm)*g%d1x(i,VAR_U) &
                              + (tau_yp - tau_ym)*g%d1y(j,VAR_U) &
                              + (tau_zp - tau_zm)*g%d1z(k,VAR_U)

                        mu_u = ibm%mu(i,j,k,VAR_U)
                        f%qs(i,j,k,VAR_U) = f%qs(i,j,k,VAR_U) + dt_alpha*sgs_u*mu_u
                        f%oldrhs(i,j,k,VAR_U) = f%oldrhs(i,j,k,VAR_U) + sgs_u
                    end if

                    if (j >= vStartY) then
                        wx = les%u_from_p_x(ip)
                        wy = les%v_from_p_y(j)
                        nut0 = (1.0d0 - wx)*les%nut(i,jm,k) + wx*les%nut(ip,jm,k)
                        nut1 = (1.0d0 - wx)*les%nut(i,j,k) + wx*les%nut(ip,j,k)
                        nut_edge = (1.0d0 - wy)*nut0 + wy*nut1
                        tau_xp = nut_edge*( &
                            (f%q(ip,j,k,VAR_U) - f%q(ip,jm,k,VAR_U))*les%inv_dy(j,VAR_U) &
                          + (f%q(ip,j,k,VAR_V) - f%q(i,j,k,VAR_V))*les%inv_dx(ip,VAR_V) )

                        wx = les%u_from_p_x(i)
                        wy = les%v_from_p_y(j)
                        nut0 = (1.0d0 - wx)*les%nut(im,jm,k) + wx*les%nut(i,jm,k)
                        nut1 = (1.0d0 - wx)*les%nut(im,j,k) + wx*les%nut(i,j,k)
                        nut_edge = (1.0d0 - wy)*nut0 + wy*nut1
                        tau_xm = nut_edge*( &
                            (f%q(i,j,k,VAR_U) - f%q(i,jm,k,VAR_U))*les%inv_dy(j,VAR_U) &
                          + (f%q(i,j,k,VAR_V) - f%q(im,j,k,VAR_V))*les%inv_dx(i,VAR_V) )

                        tau_yp = 2.0d0*les%nut(i,j,k) &
                               * (f%q(i,jp,k,VAR_V) - f%q(i,j,k,VAR_V))*g%d1y(j,VAR_P)
                        tau_ym = 2.0d0*les%nut(i,jm,k) &
                               * (f%q(i,j,k,VAR_V) - f%q(i,jm,k,VAR_V))*g%d1y(jm,VAR_P)

                        wx = les%v_from_p_y(j)
                        wy = les%w_from_p_z(kp)
                        nut0 = (1.0d0 - wx)*les%nut(i,jm,k) + wx*les%nut(i,j,k)
                        nut1 = (1.0d0 - wx)*les%nut(i,jm,kp) + wx*les%nut(i,j,kp)
                        nut_edge = (1.0d0 - wy)*nut0 + wy*nut1
                        tau_zp = nut_edge*( &
                            (f%q(i,j,kp,VAR_V) - f%q(i,j,k,VAR_V))*les%inv_dz(kp,VAR_V) &
                          + (f%q(i,j,kp,VAR_W) - f%q(i,jm,kp,VAR_W))*les%inv_dy(j,VAR_W) )

                        wx = les%v_from_p_y(j)
                        wy = les%w_from_p_z(k)
                        nut0 = (1.0d0 - wx)*les%nut(i,jm,km) + wx*les%nut(i,j,km)
                        nut1 = (1.0d0 - wx)*les%nut(i,jm,k) + wx*les%nut(i,j,k)
                        nut_edge = (1.0d0 - wy)*nut0 + wy*nut1
                        tau_zm = nut_edge*( &
                            (f%q(i,j,k,VAR_V) - f%q(i,j,km,VAR_V))*les%inv_dz(k,VAR_V) &
                          + (f%q(i,j,k,VAR_W) - f%q(i,jm,k,VAR_W))*les%inv_dy(j,VAR_W) )

                        sgs_v = (tau_xp - tau_xm)*g%d1x(i,VAR_V) &
                              + (tau_yp - tau_ym)*g%d1y(j,VAR_V) &
                              + (tau_zp - tau_zm)*g%d1z(k,VAR_V)

                        mu_v = ibm%mu(i,j,k,VAR_V)
                        f%qs(i,j,k,VAR_V) = f%qs(i,j,k,VAR_V) + dt_alpha*sgs_v*mu_v
                        f%oldrhs(i,j,k,VAR_V) = f%oldrhs(i,j,k,VAR_V) + sgs_v
                    end if

                    if (k >= wStartZ) then
                        wx = les%u_from_p_x(ip)
                        wy = les%w_from_p_z(k)
                        nut0 = (1.0d0 - wx)*les%nut(i,j,km) + wx*les%nut(ip,j,km)
                        nut1 = (1.0d0 - wx)*les%nut(i,j,k) + wx*les%nut(ip,j,k)
                        nut_edge = (1.0d0 - wy)*nut0 + wy*nut1
                        tau_xp = nut_edge*( &
                            (f%q(ip,j,k,VAR_U) - f%q(ip,j,km,VAR_U))*les%inv_dz(k,VAR_U) &
                          + (f%q(ip,j,k,VAR_W) - f%q(i,j,k,VAR_W))*les%inv_dx(ip,VAR_W) )

                        wx = les%u_from_p_x(i)
                        wy = les%w_from_p_z(k)
                        nut0 = (1.0d0 - wx)*les%nut(im,j,km) + wx*les%nut(i,j,km)
                        nut1 = (1.0d0 - wx)*les%nut(im,j,k) + wx*les%nut(i,j,k)
                        nut_edge = (1.0d0 - wy)*nut0 + wy*nut1
                        tau_xm = nut_edge*( &
                            (f%q(i,j,k,VAR_U) - f%q(i,j,km,VAR_U))*les%inv_dz(k,VAR_U) &
                          + (f%q(i,j,k,VAR_W) - f%q(im,j,k,VAR_W))*les%inv_dx(i,VAR_W) )

                        wx = les%v_from_p_y(jp)
                        wy = les%w_from_p_z(k)
                        nut0 = (1.0d0 - wx)*les%nut(i,j,km) + wx*les%nut(i,jp,km)
                        nut1 = (1.0d0 - wx)*les%nut(i,j,k) + wx*les%nut(i,jp,k)
                        nut_edge = (1.0d0 - wy)*nut0 + wy*nut1
                        tau_yp = nut_edge*( &
                            (f%q(i,jp,k,VAR_V) - f%q(i,jp,km,VAR_V))*les%inv_dz(k,VAR_V) &
                          + (f%q(i,jp,k,VAR_W) - f%q(i,j,k,VAR_W))*les%inv_dy(jp,VAR_W) )

                        wx = les%v_from_p_y(j)
                        wy = les%w_from_p_z(k)
                        nut0 = (1.0d0 - wx)*les%nut(i,jm,km) + wx*les%nut(i,j,km)
                        nut1 = (1.0d0 - wx)*les%nut(i,jm,k) + wx*les%nut(i,j,k)
                        nut_edge = (1.0d0 - wy)*nut0 + wy*nut1
                        tau_ym = nut_edge*( &
                            (f%q(i,j,k,VAR_V) - f%q(i,j,km,VAR_V))*les%inv_dz(k,VAR_V) &
                          + (f%q(i,j,k,VAR_W) - f%q(i,jm,k,VAR_W))*les%inv_dy(j,VAR_W) )

                        tau_zp = 2.0d0*les%nut(i,j,k) &
                               * (f%q(i,j,kp,VAR_W) - f%q(i,j,k,VAR_W))*g%d1z(k,VAR_P)
                        tau_zm = 2.0d0*les%nut(i,j,km) &
                               * (f%q(i,j,k,VAR_W) - f%q(i,j,km,VAR_W))*g%d1z(km,VAR_P)

                        sgs_w = (tau_xp - tau_xm)*g%d1x(i,VAR_W) &
                              + (tau_yp - tau_ym)*g%d1y(j,VAR_W) &
                              + (tau_zp - tau_zm)*g%d1z(k,VAR_W)

                        mu_w = ibm%mu(i,j,k,VAR_W)
                        f%qs(i,j,k,VAR_W) = f%qs(i,j,k,VAR_W) + dt_alpha*sgs_w*mu_w
                        f%oldrhs(i,j,k,VAR_W) = f%oldrhs(i,j,k,VAR_W) + sgs_w
                    end if
                end do
            end do
        end do
        !$omp end target teams distribute parallel do

        if (present(les_prof)) call add_les_profile(les_prof, LES_PROF_SGS, &
            les_wall_seconds() - profile_start)
    end subroutine add_les_momentum_correction

    subroutine get_timestep_rates(f, dns, g, rates, les)
        type(field_type), intent(inout) :: f
        type(dns_type),   intent(in)    :: dns
        type(grid_type),  intent(in)    :: g
        real(C_DOUBLE), intent(out) :: rates(1:NCFL)
        type(les_type), intent(in), optional :: les

        integer :: i,j,k
        integer :: nx, ny, nz
        real(C_DOUBLE) :: cfl_rate, peclet_rate, ire, nu_eff
        logical :: use_les

        nx = int(dns%localSize(1,2))
        ny = int(dns%localSize(2,2))
        nz = int(dns%localSize(3,2))
        cfl_rate = 0.0d0
        use_les = .false.
        if (present(les)) use_les = les_is_enabled(les) .and. allocated(les%nut)

        !$omp target teams distribute parallel do collapse(3) reduction(max:cfl_rate) &
        !$omp& map(to: g%d1x, g%d1y, g%d1z, f%q) &
        !$omp& private(i,j,k)
        do k = 0, nz+1
            do j = 0, ny+1
                do i = 0, nx+1
                    cfl_rate = max(cfl_rate, abs(f%q(i,j,k,VAR_U)*g%d1x(i,VAR_U)))
                    cfl_rate = max(cfl_rate, abs(f%q(i,j,k,VAR_V)*g%d1y(j,VAR_V)))
                    cfl_rate = max(cfl_rate, abs(f%q(i,j,k,VAR_W)*g%d1z(k,VAR_W)))
                end do
            end do
        end do
        !$omp end target teams distribute parallel do

        rates(CFL_COURANT) = cfl_rate
        rates(CFL_PECLET) = dns%peclet_rate
        if (.not. use_les) return

        peclet_rate = dns%peclet_rate
        ire = 1.0d0/dns%re

        !$omp target teams distribute parallel do collapse(3) reduction(max:peclet_rate) &
        !$omp& map(to: g%d1x, g%d1y, g%d1z, les%nut) &
        !$omp& private(i,j,k,nu_eff)
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    nu_eff = ire + max(0.0d0, les%nut(i,j,k))
                    peclet_rate = max(peclet_rate, nu_eff*g%d1x(i,VAR_P)**2)
                    peclet_rate = max(peclet_rate, nu_eff*g%d1y(j,VAR_P)**2)
                    peclet_rate = max(peclet_rate, nu_eff*g%d1z(k,VAR_P)**2)
                end do
            end do
        end do
        !$omp end target teams distribute parallel do

        rates(CFL_PECLET) = peclet_rate
    end subroutine get_timestep_rates

    logical function run_should_continue(dns, loop_steps)
        type(dns_type), intent(in) :: dns
        integer(C_INT), intent(in) :: loop_steps

        real(C_DOUBLE) :: time_tol

        run_should_continue = .true.
        if (dns%nsteps > 0_C_INT .and. loop_steps >= dns%nsteps) run_should_continue = .false.
        if (dns%t_final > 0.0d0) then
            time_tol = max(1.0d-12, 100.0d0*epsilon(1.0d0)*max(1.0d0, abs(dns%t_final)))
            if (dns%t_current >= dns%t_final - time_tol) run_should_continue = .false.
        end if
    end function run_should_continue

    subroutine trim_dt_for_final_time(dns)
        type(dns_type), intent(inout) :: dns

        real(C_DOUBLE) :: remaining

        if (dns%t_final <= 0.0d0) return

        remaining = dns%t_final - dns%t_current
        dns%dt = min(dns%dt, max(0.0d0, remaining))
    end subroutine trim_dt_for_final_time

    subroutine update_timestep_limits(f, dns, g, c, les)
        type(field_type), intent(inout) :: f
        type(dns_type), intent(inout) :: dns
        type(grid_type), intent(in) :: g
        type(comm_type), intent(in) :: c
        type(les_type), intent(in), optional :: les

        real(C_DOUBLE) :: rates(1:NCFL), next_dt
        logical :: have_limit

        if (dns%cflmax <= 0.0d0 .and. dns%pecletmax <= 0.0d0) return

        if (present(les)) then
            call get_timestep_rates(f, dns, g, rates, les)
        else
            call get_timestep_rates(f, dns, g, rates)
        end if
        call comm_allreduce_max(c, rates)

        next_dt = dns%dt
        have_limit = .false.
        if (dns%cflmax > 0.0d0 .and. rates(CFL_COURANT) > 0.0d0) then
            if (.not. have_limit) next_dt = dns%dtmax
            next_dt = min(next_dt, dns%cflmax/rates(CFL_COURANT))
            have_limit = .true.
        end if
        if (dns%pecletmax > 0.0d0 .and. rates(CFL_PECLET) > 0.0d0) then
            if (.not. have_limit) next_dt = dns%dtmax
            next_dt = min(next_dt, dns%pecletmax/rates(CFL_PECLET))
            have_limit = .true.
        end if

        if (have_limit) dns%dt = next_dt
        dns%cfl(CFL_COURANT) = dns%dt*rates(CFL_COURANT)
        dns%cfl(CFL_PECLET) = dns%dt*rates(CFL_PECLET)
    end subroutine update_timestep_limits

end module step
