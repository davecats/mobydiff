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
    use :: poisson_workspace, only: poisson_fft_workspace
    implicit none

    real(C_DOUBLE), parameter :: rk_alpha(3) = [64.0d0/120.0d0,  50.0d0/120.0d0,  90.0d0/120.0d0]
    real(C_DOUBLE), parameter :: rk_beta(3)  = [ 0.0d0,         -34.0d0/120.0d0, -50.0d0/120.0d0]
    real(C_DOUBLE), parameter :: rk_gamma(3) = [64.0d0/120.0d0,  16.0d0/120.0d0,  40.0d0/120.0d0]

contains

    subroutine momentum(f, g, dt_alpha, dt_beta, dt_gamma)
        type(field_type), intent(inout) :: f
        type(grid_type),  intent(in)    :: g
        real(C_DOUBLE),   intent(in)    :: dt_alpha, dt_beta, dt_gamma

        integer :: i,j,k,ip,im,kp,km,jp,jm
        integer :: nx, ny, nz

        real(C_DOUBLE) :: diff_ux,diff_uy,diff_uz
        real(C_DOUBLE) :: diff_vx,diff_vy,diff_vz
        real(C_DOUBLE) :: diff_wx,diff_wy,diff_wz

        real(C_DOUBLE) :: uu_p,uu_m,uv_p,uv_m,uw_p,uw_m
        real(C_DOUBLE) :: vu_p,vu_m,vv_p,vv_m,vw_p,vw_m
        real(C_DOUBLE) :: wu_p,wu_m,ww_p,ww_m,wv_p,wv_m

        real(C_DOUBLE) :: dpx,dpy,dpz,rhs
        real(C_DOUBLE) :: dx, dy, dz, dx2, dy2, dz2, re

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

        ! x-momentum
        !$omp target teams distribute parallel do collapse(3) &
        !$omp& map(to: dt_alpha, dt_beta, dt_gamma, f%un(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%vn(0:nx+1,1:ny+1,0:nz+1), &
        !$omp& f%wn(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%pn(0:nx+1,1:ny,0:nz+1)) &
        !$omp& map(tofrom: f%us(0:nx+1,0:ny+1,0:nz+1), f%oldrhsu(1:nx,1:ny,1:nz)) &
        !$omp& private(i,j,k,ip,im,jp,jm,kp,km,uu_p,uu_m,uv_p,uv_m,uw_p,uw_m, &
        !$omp& diff_ux,diff_uy,diff_uz,dpx,rhs)
        do i = 1, nx
            do j = 1, ny
                do k = 1, nz

                    ip = i+1
                    im = i-1
                    jp = j+1
                    jm = j-1
                    kp = k+1
                    km = k-1
                     
                    uu_p = 0.25d0*(f%un(ip,j,k)+f%un(i,j,k))**2
                    uu_m = 0.25d0*(f%un(i,j,k)+f%un(im,j,k))**2

                    uv_p = 0.25d0*(f%un(i,jp,k)+f%un(i,j,k))*(f%vn(i,jp,k)+f%vn(im,jp,k))
                    uv_m = 0.25d0*(f%un(i,j,k)+f%un(i,jm,k))*(f%vn(i,j,k)+f%vn(im,j,k))

                    uw_p = 0.25d0*(f%un(i,j,k)+f%un(i,j,kp))*(f%wn(i,j,kp)+f%wn(im,j,kp))
                    uw_m = 0.25d0*(f%un(i,j,k)+f%un(i,j,km))*(f%wn(i,j,k)+f%wn(im,j,k))

                    diff_ux = (f%un(im,j,k)-2.0d0*f%un(i,j,k)+f%un(ip,j,k))/dx2
                    diff_uy = (f%un(i,jm,k)-2.0d0*f%un(i,j,k)+f%un(i,jp,k))/dy2
                    diff_uz = (f%un(i,j,km)-2.0d0*f%un(i,j,k)+f%un(i,j,kp))/dz2

                    dpx = (f%pn(i,j,k)-f%pn(im,j,k))/dx 

                    rhs = ( &
                        -(uu_p-uu_m)/dx &
                        -(uv_p-uv_m)/dy &
                        -(uw_p-uw_m)/dz &
                        + 1 &
                        + (1.0d0/re)*(diff_ux + diff_uy + diff_uz) )

                    f%us(i,j,k) = f%un(i,j,k) + dt_alpha*rhs &
                        + dt_beta*f%oldrhsu(i,j,k) - dt_gamma*dpx

                    f%oldrhsu(i,j,k) = rhs

                end do
            end do
        end do
        !$omp end target teams distribute parallel do

        ! y-momentum
        !$omp target teams distribute parallel do collapse(3) &
        !$omp& map(to: dt_alpha, dt_beta, dt_gamma, f%un(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%vn(0:nx+1,1:ny+1,0:nz+1), &
        !$omp& f%wn(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%pn(0:nx+1,1:ny,0:nz+1)) &
        !$omp& map(tofrom: f%vs(0:nx+1,1:ny+1,0:nz+1), f%oldrhsv(1:nx,2:ny,1:nz)) &
        !$omp& private(i,j,k,ip,im,jp,jm,kp,km,vu_p,vu_m,vv_p,vv_m,vw_p,vw_m, &
        !$omp& diff_vx,diff_vy,diff_vz,dpy,rhs)
        do i = 1, nx
            do j = 2, ny
                do k = 1, nz

                    ip = i+1
                    im = i-1
                    jp = j+1
                    jm = j-1
                    kp = k+1
                    km = k-1

                    vu_p = 0.25d0*(f%vn(i,j,k)+f%vn(ip,j,k))*(f%un(ip,j,k)+f%un(ip,jm,k))
                    vu_m = 0.25d0*(f%vn(i,j,k)+f%vn(im,j,k))*(f%un(i,j,k)+f%un(i,jm,k))

                    vv_p = 0.25d0*(f%vn(i,j,k)+f%vn(i,jp,k))**2
                    vv_m = 0.25d0*(f%vn(i,j,k)+f%vn(i,jm,k))**2

                    vw_p = 0.25d0*(f%vn(i,j,kp)+f%vn(i,j,k))*(f%wn(i,j,kp)+f%wn(i,jm,kp))
                    vw_m = 0.25d0*(f%vn(i,j,km)+f%vn(i,j,k))*(f%wn(i,j,k)+f%wn(i,jm,k))

                    diff_vx = (f%vn(im,j,k)-2.0d0*f%vn(i,j,k)+f%vn(ip,j,k))/dx2
                    diff_vy = (f%vn(i,jm,k)-2.0d0*f%vn(i,j,k)+f%vn(i,jp,k))/dy2
                    diff_vz = (f%vn(i,j,km)-2.0d0*f%vn(i,j,k)+f%vn(i,j,kp))/dz2

                    dpy = (f%pn(i,j,k)-f%pn(i,jm,k))/dy

                    rhs = ( &
                        -(vu_p-vu_m)/dx &
                        -(vv_p-vv_m)/dy &
                        -(vw_p-vw_m)/dz &
                        + (1.0d0/re)*(diff_vx + diff_vy + diff_vz) )

                    f%vs(i,j,k) = f%vn(i,j,k) + dt_alpha*rhs &
                        + dt_beta*f%oldrhsv(i,j,k) - dt_gamma*dpy

                    f%oldrhsv(i,j,k) = rhs

                end do
            end do
        end do
        !$omp end target teams distribute parallel do

        ! z-momentum
        !$omp target teams distribute parallel do collapse(3) &
        !$omp& map(to: dt_alpha, dt_beta, dt_gamma, f%un(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%vn(0:nx+1,1:ny+1,0:nz+1), &
        !$omp& f%wn(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%pn(0:nx+1,1:ny,0:nz+1)) &
        !$omp& map(tofrom: f%ws(0:nx+1,0:ny+1,0:nz+1), f%oldrhsw(1:nx,1:ny,1:nz)) &
        !$omp& private(i,j,k,ip,im,jp,jm,kp,km,wu_p,wu_m,ww_p,ww_m,wv_p,wv_m, &
        !$omp& diff_wx,diff_wy,diff_wz,dpz,rhs)
        do i = 1, nx
            do j = 1, ny
                do k = 1, nz

                    ip = i+1
                    im = i-1
                    jp = j+1
                    jm = j-1
                    kp = k+1
                    km = k-1

                    wu_p = 0.25d0*(f%wn(i,j,k)+f%wn(ip,j,k))*(f%un(ip,j,k)+f%un(ip,j,km))
                    wu_m = 0.25d0*(f%wn(i,j,k)+f%wn(im,j,k))*(f%un(i,j,k)+f%un(i,j,km))

                    ww_p = 0.25d0*(f%wn(i,j,k)+f%wn(i,j,kp))**2
                    ww_m = 0.25d0*(f%wn(i,j,k)+f%wn(i,j,km))**2

                    wv_p = 0.25d0*(f%wn(i,j,k)+f%wn(i,jp,k))*(f%vn(i,jp,k)+f%vn(i,jp,km))
                    wv_m = 0.25d0*(f%wn(i,j,k)+f%wn(i,jm,k))*(f%vn(i,j,k)+f%vn(i,j,km))

                    diff_wx = (f%wn(im,j,k)-2.0d0*f%wn(i,j,k)+f%wn(ip,j,k))/dx2
                    diff_wy = (f%wn(i,jm,k)-2.0d0*f%wn(i,j,k)+f%wn(i,jp,k))/dy2
                    diff_wz = (f%wn(i,j,km)-2.0d0*f%wn(i,j,k)+f%wn(i,j,kp))/dz2

                    dpz = (f%pn(i,j,k)-f%pn(i,j,km))/dz

                    rhs = ( &
                        -(wu_p-wu_m)/dx &
                        -(wv_p-wv_m)/dy &
                        -(ww_p-ww_m)/dz &
                        + (1.0d0/re)*(diff_wx + diff_wy + diff_wz) )

                    f%ws(i,j,k) = f%wn(i,j,k) + dt_alpha*rhs &
                        + dt_beta*f%oldrhsw(i,j,k) - dt_gamma*dpz

                    f%oldrhsw(i,j,k) = rhs

                end do
            end do
        end do
        !$omp end target teams distribute parallel do

    end subroutine momentum


    subroutine corrector(f, g, dt_gamma)
        type(field_type), intent(inout) :: f
        type(grid_type),  intent(in)    :: g
        real(C_DOUBLE),   intent(in)    :: dt_gamma

        integer :: i,j,k,im,km
        integer :: nx, ny, nz
        real(C_DOUBLE) :: dx, dy, dz

        nx = g%nx
        ny = g%ny
        nz = g%nz
        dx = g%dx
        dy = g%dy
        dz = g%dz

        !$omp target teams distribute parallel do collapse(3) &
        !$omp& map(to: dt_gamma, f%us(0:nx+1,0:ny+1,0:nz+1), f%ws(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%pc(0:nx+1,1:ny,0:nz+1)) &
        !$omp& map(tofrom: f%un(0:nx+1,0:ny+1,0:nz+1), f%wn(0:nx+1,0:ny+1,0:nz+1)) &
        !$omp& private(i,j,k)
        do i = 1, nx
            do j = 1, ny
                do k = 1, nz
                    f%un(i,j,k) = f%us(i,j,k) - dt_gamma*(f%pc(i,j,k)-f%pc(i-1,j,k))/dx
                    f%wn(i,j,k) = f%ws(i,j,k) - dt_gamma*(f%pc(i,j,k)-f%pc(i,j,k-1))/dz
                end do
            end do
        end do
        !$omp end target teams distribute parallel do

        !$omp target teams distribute parallel do collapse(3) &
        !$omp& map(to: dt_gamma, f%vs(0:nx+1,1:ny+1,0:nz+1), &
        !$omp& f%pc(0:nx+1,1:ny,0:nz+1)) &
        !$omp& map(tofrom: f%vn(0:nx+1,1:ny+1,0:nz+1)) &
        !$omp& private(i,j,k)
        do i = 1, nx
            do j = 2, ny
                do k = 1, nz
                    f%vn(i,j,k) = f%vs(i,j,k) - dt_gamma*(f%pc(i,j,k)-f%pc(i,j-1,k))/dy
                end do
            end do
        end do
        !$omp end target teams distribute parallel do

        !$omp target teams distribute parallel do collapse(3) &
        !$omp& map(to: f%pc(0:nx+1,1:ny,0:nz+1)) &
        !$omp& map(tofrom: f%pn(0:nx+1,1:ny,0:nz+1)) &
        !$omp& private(i,j,k)
        do i = 0, nx+1
            do j = 1, ny
                do k = 0, nz+1
                    f%pn(i,j,k) = f%pn(i,j,k) + f%pc(i,j,k)
                end do
            end do
        end do
        !$omp end target teams distribute parallel do

    end subroutine corrector


    subroutine divU(f, g, ws, dt_gamma)
        type(field_type), intent(inout) :: f
        type(grid_type),  intent(in)    :: g
        type(poisson_fft_workspace), intent(inout) :: ws
        real(C_DOUBLE),   intent(in)    :: dt_gamma

        integer :: i,j,k
        integer :: nx, ny, nz
        real(C_DOUBLE) :: dx, dy, dz

        nx = g%nx
        ny = g%ny
        nz = g%nz
        dx = g%dx
        dy = g%dy
        dz = g%dz

#ifdef USE_CUFFT
        !$omp target teams distribute parallel do collapse(3) &
        !$omp& map(to: dt_gamma, f%us(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%vs(0:nx+1,1:ny+1,0:nz+1), &
        !$omp& f%ws(0:nx+1,0:ny+1,0:nz+1)) &
        !$omp& map(tofrom: ws%plane(1:nx,1:nz,1:ny)) &
        !$omp& private(i,j,k)
        do i = 1, nx
            do j = 1, ny
                do k = 1, nz

                    ws%plane(i,k,j) = ( &
                        (f%us(i+1,j,k)-f%us(i,j,k))/dx &
                      + (f%vs(i,j+1,k)-f%vs(i,j,k))/dy &
                      + (f%ws(i,j,k+1)-f%ws(i,j,k))/dz ) / dt_gamma

                end do
            end do
        end do
        !$omp end target teams distribute parallel do
#else
        do i = 1, nx
            do j = 1, ny
                do k = 1, nz
                    ws%rhs(i,j,k) = ( &
                        (f%us(i+1,j,k)-f%us(i,j,k))/dx &
                      + (f%vs(i,j+1,k)-f%vs(i,j,k))/dy &
                      + (f%ws(i,j,k+1)-f%ws(i,j,k))/dz ) / dt_gamma
                end do
            end do
        end do
#endif

    end subroutine divU
    
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
        !$omp& map(to: f%vn(0:nx+1,1:ny+1,0:nz+1)) private(i,j,k)
        do i = 0, nx+1
            do j = 1, ny+1
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
