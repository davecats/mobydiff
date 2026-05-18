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
    use :: init, only: grid_type, field_type
    use :: ibmm, only: ibm_type
    use :: boundary, only: boundary_type
    implicit none

    real(C_DOUBLE), parameter :: rk_alpha(3) = [64.0d0/120.0d0,  50.0d0/120.0d0,  90.0d0/120.0d0]
    real(C_DOUBLE), parameter :: rk_beta(3)  = [ 0.0d0,         -34.0d0/120.0d0, -50.0d0/120.0d0]
    real(C_DOUBLE), parameter :: rk_gamma(3) = [64.0d0/120.0d0,  16.0d0/120.0d0,  40.0d0/120.0d0]

contains

    subroutine momentum(f, g, dt_alpha, dt_beta, dt_gamma, ibm, bc)
        type(field_type), intent(inout) :: f
        type(grid_type),  intent(in)    :: g
        real(C_DOUBLE),   intent(in)    :: dt_alpha, dt_beta, dt_gamma
        type(ibm_type),   intent(in)    :: ibm
        type(boundary_type),   intent(in)  :: bc

        integer :: i,j,k,ip,im,kp,km,jp,jm
        integer :: nx, ny, nz, nx0, ny0, nz0 

        real(C_DOUBLE) :: diff_ux,diff_uy,diff_uz
        real(C_DOUBLE) :: diff_vx,diff_vy,diff_vz
        real(C_DOUBLE) :: diff_wx,diff_wy,diff_wz

        real(C_DOUBLE) :: uu_p,uu_m,uv_p,uv_m,uw_p,uw_m
        real(C_DOUBLE) :: vu_p,vu_m,vv_p,vv_m,vw_p,vw_m
        real(C_DOUBLE) :: wu_p,wu_m,ww_p,ww_m,wv_p,wv_m

        real(C_DOUBLE) :: dpx,dpy,dpz,rhsu,rhsv,rhsw
        real(C_DOUBLE) :: mu_u, mu_v, mu_w
        real(C_DOUBLE) :: dx, dy, dz, dx2, dy2, dz2, re
        real(C_DOUBLE) :: idx, idy, idz, idx2, idy2, idz2, ire
        real(C_DOUBLE) :: forcing_x, forcing_y, forcing_z

        nx = g%nx
        ny = g%ny
        nz = g%nz
        dx = g%dx
        dy = g%dy
        dz = g%dz
        dx2 = dx*dx
        dy2 = dy*dy
        dz2 = dz*dz
        re = g%re
        idx = 1.0d0/dx
        idy = 1.0d0/dy
        idz = 1.0d0/dz
        idx2 = 1.0d0/dx2
        idy2 = 1.0d0/dy2
        idz2 = 1.0d0/dz2
        ire = 1.0d0/re
        forcing_x = g%forcing_x
        forcing_y = g%forcing_y
        forcing_z = g%forcing_z
        nx0 = 2; if (bc%isPeriodic(1)) nx0 = nx0 - 1
        ny0 = 2; if (bc%isPeriodic(2)) ny0 = ny0 - 1
        nz0 = 2; if (bc%isPeriodic(3)) nz0 = nz0 - 1

        ! Predictor for all staggered velocity components.
        !$omp target teams distribute parallel do collapse(3) &
        !$omp& map(to: dt_alpha, dt_beta, dt_gamma, nx0, ny0, nz0, &
        !$omp& idx, idy, idz, idx2, idy2, idz2, ire, forcing_x, forcing_y, forcing_z, &
        !$omp& f%un(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%vn(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%wn(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%pn(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& ibm%coef_u(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& ibm%coef_v(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& ibm%coef_w(0:nx+1,0:ny+1,0:nz+1)) &
        !$omp& map(tofrom: f%us(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%vs(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%ws(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%oldrhsu(1:nx,1:ny,1:nz), &
        !$omp& f%oldrhsv(1:nx,1:ny,1:nz), &
        !$omp& f%oldrhsw(1:nx,1:ny,1:nz)) &
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

                    if (i >= nx0) then
                    uu_p = 0.25d0*(f%un(ip,j,k)+f%un(i,j,k))**2
                    uu_m = 0.25d0*(f%un(i,j,k)+f%un(im,j,k))**2

                    uv_p = 0.25d0*(f%un(i,jp,k)+f%un(i,j,k))*(f%vn(i,jp,k)+f%vn(im,jp,k))
                    uv_m = 0.25d0*(f%un(i,j,k)+f%un(i,jm,k))*(f%vn(i,j,k)+f%vn(im,j,k))

                    uw_p = 0.25d0*(f%un(i,j,k)+f%un(i,j,kp))*(f%wn(i,j,kp)+f%wn(im,j,kp))
                    uw_m = 0.25d0*(f%un(i,j,k)+f%un(i,j,km))*(f%wn(i,j,k)+f%wn(im,j,k))

                    diff_ux = (f%un(im,j,k)-2.0d0*f%un(i,j,k)+f%un(ip,j,k))*idx2
                    diff_uy = (f%un(i,jm,k)-2.0d0*f%un(i,j,k)+f%un(i,jp,k))*idy2
                    diff_uz = (f%un(i,j,km)-2.0d0*f%un(i,j,k)+f%un(i,j,kp))*idz2

                    dpx = (f%pn(i,j,k)-f%pn(im,j,k))*idx

                    rhsu = ( &
                        -(uu_p-uu_m)*idx &
                        -(uv_p-uv_m)*idy &
                        -(uw_p-uw_m)*idz &
                        + forcing_x &
                        + ire*(diff_ux + diff_uy + diff_uz) )

                    f%us(i,j,k) = f%un(i,j,k) + dt_alpha*rhsu &
                        + dt_beta*f%oldrhsu(i,j,k) - dt_gamma*dpx

                    mu_u = 1.0d0/(1.0d0 + dt_gamma*ibm%coef_u(i,j,k))
                    f%us(i,j,k) = f%us(i,j,k)*mu_u

                    f%oldrhsu(i,j,k) = rhsu
                    end if

                    if (j >= ny0) then
                    vu_p = 0.25d0*(f%vn(i,j,k)+f%vn(ip,j,k))*(f%un(ip,j,k)+f%un(ip,jm,k))
                    vu_m = 0.25d0*(f%vn(i,j,k)+f%vn(im,j,k))*(f%un(i,j,k)+f%un(i,jm,k))

                    vv_p = 0.25d0*(f%vn(i,j,k)+f%vn(i,jp,k))**2
                    vv_m = 0.25d0*(f%vn(i,j,k)+f%vn(i,jm,k))**2

                    vw_p = 0.25d0*(f%vn(i,j,kp)+f%vn(i,j,k))*(f%wn(i,j,kp)+f%wn(i,jm,kp))
                    vw_m = 0.25d0*(f%vn(i,j,km)+f%vn(i,j,k))*(f%wn(i,j,k)+f%wn(i,jm,k))

                    diff_vx = (f%vn(im,j,k)-2.0d0*f%vn(i,j,k)+f%vn(ip,j,k))*idx2
                    diff_vy = (f%vn(i,jm,k)-2.0d0*f%vn(i,j,k)+f%vn(i,jp,k))*idy2
                    diff_vz = (f%vn(i,j,km)-2.0d0*f%vn(i,j,k)+f%vn(i,j,kp))*idz2

                    dpy = (f%pn(i,j,k)-f%pn(i,jm,k))*idy

                    rhsv = ( &
                        -(vu_p-vu_m)*idx &
                        -(vv_p-vv_m)*idy &
                        -(vw_p-vw_m)*idz &
                        + forcing_y &
                        + ire*(diff_vx + diff_vy + diff_vz) )

                    f%vs(i,j,k) = f%vn(i,j,k) + dt_alpha*rhsv &
                        + dt_beta*f%oldrhsv(i,j,k) - dt_gamma*dpy

                    mu_v = 1.0d0/(1.0d0 + dt_gamma*ibm%coef_v(i,j,k))
                    f%vs(i,j,k) = f%vs(i,j,k)*mu_v

                    f%oldrhsv(i,j,k) = rhsv
                    end if

                    if (k >= nz0) then
                    wu_p = 0.25d0*(f%wn(i,j,k)+f%wn(ip,j,k))*(f%un(ip,j,k)+f%un(ip,j,km))
                    wu_m = 0.25d0*(f%wn(i,j,k)+f%wn(im,j,k))*(f%un(i,j,k)+f%un(i,j,km))

                    ww_p = 0.25d0*(f%wn(i,j,k)+f%wn(i,j,kp))**2
                    ww_m = 0.25d0*(f%wn(i,j,k)+f%wn(i,j,km))**2

                    wv_p = 0.25d0*(f%wn(i,j,k)+f%wn(i,jp,k))*(f%vn(i,jp,k)+f%vn(i,jp,km))
                    wv_m = 0.25d0*(f%wn(i,j,k)+f%wn(i,jm,k))*(f%vn(i,j,k)+f%vn(i,j,km))

                    diff_wx = (f%wn(im,j,k)-2.0d0*f%wn(i,j,k)+f%wn(ip,j,k))*idx2
                    diff_wy = (f%wn(i,jm,k)-2.0d0*f%wn(i,j,k)+f%wn(i,jp,k))*idy2
                    diff_wz = (f%wn(i,j,km)-2.0d0*f%wn(i,j,k)+f%wn(i,j,kp))*idz2

                    dpz = (f%pn(i,j,k)-f%pn(i,j,km))*idz

                    rhsw = ( &
                        -(wu_p-wu_m)*idx &
                        -(wv_p-wv_m)*idy &
                        -(ww_p-ww_m)*idz &
                        + forcing_z &
                        + ire*(diff_wx + diff_wy + diff_wz) )

                    f%ws(i,j,k) = f%wn(i,j,k) + dt_alpha*rhsw &
                        + dt_beta*f%oldrhsw(i,j,k) - dt_gamma*dpz

                    mu_w = 1.0d0/(1.0d0 + dt_gamma*ibm%coef_w(i,j,k))
                    f%ws(i,j,k) = f%ws(i,j,k)*mu_w

                    f%oldrhsw(i,j,k) = rhsw
                    end if

                end do
            end do
        end do
        !$omp end target teams distribute parallel do

    end subroutine momentum


    real(C_DOUBLE) function get_cfl(f,g)
        type(field_type), intent(inout) :: f
        type(grid_type),  intent(in)    :: g

        integer :: i,j,k
        integer :: nx, ny, nz
        real(C_DOUBLE) :: dx, dy, dz

        nx = g%nx
        ny = g%ny
        nz = g%nz
        dx = g%dx
        dy = g%dy
        dz = g%dz
        get_cfl = 0.0d0

        !$omp target teams distribute parallel do collapse(3) reduction(max:get_cfl) &
        !$omp& map(to: f%un(0:nx+1,0:ny+1,0:nz+1)) private(i,j,k)
        do i = 0, nx+1
            do j = 0, ny+1
                do k = 0, nz+1
                    get_cfl = max(get_cfl, abs(f%un(i,j,k)/dx))
                end do
            end do
        end do
        !$omp end target teams distribute parallel do

        !$omp target teams distribute parallel do collapse(3) reduction(max:get_cfl) &
        !$omp& map(to: f%vn(0:nx+1,0:ny+1,0:nz+1)) private(i,j,k)
        do i = 0, nx+1
            do j = 0, ny+1
                do k = 0, nz+1
                    get_cfl = max(get_cfl, abs(f%vn(i,j,k)/dy))
                end do
            end do
        end do
        !$omp end target teams distribute parallel do

        !$omp target teams distribute parallel do collapse(3) reduction(max:get_cfl) &
        !$omp& map(to: f%wn(0:nx+1,0:ny+1,0:nz+1)) private(i,j,k)
        do i = 0, nx+1
            do j = 0, ny+1
                do k = 0, nz+1
                    get_cfl = max(get_cfl, abs(f%wn(i,j,k)/dz))
                end do
            end do
        end do
        !$omp end target teams distribute parallel do
    end function get_cfl

end module step
