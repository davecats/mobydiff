module boundary
    use, intrinsic :: iso_c_binding
    use :: init, only: grid_type, field_type
    implicit none

contains

    subroutine apply_bc(f, g)
        type(field_type), intent(inout) :: f
        type(grid_type),  intent(in)    :: g
        integer :: i, j, k
        integer :: nx, ny, nz

        nx = g%nx
        ny = g%ny
        nz = g%nz

        ! No-slip walls in y.
#ifdef USE_REDBLACK
        !$omp target teams distribute parallel do collapse(2) &
        !$omp& map(tofrom: f%un(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%vn(0:nx+1,1:ny+1,0:nz+1), &
        !$omp& f%wn(0:nx+1,0:ny+1,0:nz+1)) private(i,k)
#else
        !$omp target teams distribute parallel do collapse(2) &
        !$omp& map(tofrom: f%un(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%vn(0:nx+1,1:ny+1,0:nz+1), &
        !$omp& f%wn(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%us(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%vs(0:nx+1,1:ny+1,0:nz+1), &
        !$omp& f%ws(0:nx+1,0:ny+1,0:nz+1)) private(i,k)
#endif
        do i = 0, nx+1
            do k = 0, nz+1
                f%un(i,ny+1,k) = -f%un(i,ny,k)
                f%vn(i,ny+1,k) = 0.0d0
                f%wn(i,ny+1,k) = -f%wn(i,ny,k)

                f%un(i,0,k) = -f%un(i,1,k)
                f%vn(i,1,k) = 0.0d0
                f%wn(i,0,k) = -f%wn(i,1,k)
#ifndef USE_REDBLACK
                f%us(i,ny+1,k) = -f%us(i,ny,k)
                f%vs(i,ny+1,k) = 0.0d0
                f%ws(i,ny+1,k) = -f%ws(i,ny,k)

                f%us(i,0,k) = -f%us(i,1,k)
                f%vs(i,1,k) = 0.0d0
                f%ws(i,0,k) = -f%ws(i,1,k)
#endif                
            end do
        end do
        !$omp end target teams distribute parallel do

        ! Periodicity in x for velocity and pressure.
#ifdef USE_REDBLACK
        !$omp target teams distribute parallel do collapse(2) &
        !$omp& map(tofrom: f%un(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%vn(0:nx+1,1:ny+1,0:nz+1), &
        !$omp& f%wn(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%pn(0:nx+1,1:ny,0:nz+1)) private(j,k)
#else
        !$omp target teams distribute parallel do collapse(2) &
        !$omp& map(tofrom: f%un(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%vn(0:nx+1,1:ny+1,0:nz+1), &
        !$omp& f%wn(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%us(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%vs(0:nx+1,1:ny+1,0:nz+1), &
        !$omp& f%ws(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%pn(0:nx+1,1:ny,0:nz+1)) private(j,k)
#endif
        do j = 0, ny+1
            do k = 0, nz+1
                f%un(0,j,k) = f%un(nx,j,k)
                f%wn(0,j,k) = f%wn(nx,j,k)

                f%un(nx+1,j,k) = f%un(1,j,k)
                f%wn(nx+1,j,k) = f%wn(1,j,k)

#ifndef USE_REDBLACK
                f%us(0,j,k) = f%us(nx,j,k)
                f%ws(0,j,k) = f%ws(nx,j,k)

                f%us(nx+1,j,k) = f%us(1,j,k)
                f%ws(nx+1,j,k) = f%ws(1,j,k)
#endif

                if (j >= 1) then
                    f%vn(0,j,k) = f%vn(nx,j,k)
                    f%vn(nx+1,j,k) = f%vn(1,j,k)
#ifndef USE_REDBLACK
                    f%vs(0,j,k) = f%vs(nx,j,k)
                    f%vs(nx+1,j,k) = f%vs(1,j,k)
#endif
                end if

                if (j >= 1 .and. j <= ny) then
                    f%pn(0,j,k) = f%pn(nx,j,k)
                    f%pn(nx+1,j,k) = f%pn(1,j,k)
                end if
            end do
        end do
        !$omp end target teams distribute parallel do

        ! Periodicity in z for velocity and pressure.
#ifdef USE_REDBLACK
        !$omp target teams distribute parallel do collapse(2) &
        !$omp& map(tofrom: f%un(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%vn(0:nx+1,1:ny+1,0:nz+1), &
        !$omp& f%wn(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%pn(0:nx+1,1:ny,0:nz+1)) private(i,j)
#else
        !$omp target teams distribute parallel do collapse(2) &
        !$omp& map(tofrom: f%un(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%vn(0:nx+1,1:ny+1,0:nz+1), &
        !$omp& f%us(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%vs(0:nx+1,1:ny+1,0:nz+1), &
        !$omp& f%ws(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%wn(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%pn(0:nx+1,1:ny,0:nz+1)) private(i,j)
#endif
        do i = 0, nx+1
            do j = 0, ny+1
                f%un(i,j,0) = f%un(i,j,nz)
                f%wn(i,j,0) = f%wn(i,j,nz)

                f%un(i,j,nz+1) = f%un(i,j,1)
                f%wn(i,j,nz+1) = f%wn(i,j,1)

#ifndef USE_REDBLACK
                f%us(i,j,0) = f%us(i,j,nz)
                f%ws(i,j,0) = f%ws(i,j,nz)
                
                f%us(i,j,nz+1) = f%us(i,j,1)
                f%ws(i,j,nz+1) = f%ws(i,j,1)
#endif

                if (j >= 1) then
                    f%vn(i,j,0) = f%vn(i,j,nz)
                    f%vn(i,j,nz+1) = f%vn(i,j,1)

#ifndef USE_REDBLACK
                    f%vs(i,j,0) = f%vs(i,j,nz)
                    f%vs(i,j,nz+1) = f%vs(i,j,1)
#endif                    
                end if

                if (j >= 1 .and. j <= ny) then
                    f%pn(i,j,0) = f%pn(i,j,nz)
                    f%pn(i,j,nz+1) = f%pn(i,j,1)
                end if
            end do
        end do
        !$omp end target teams distribute parallel do
    end subroutine apply_bc

#ifdef USE_REDBLACK
    subroutine apply_bc_redblack(f, g)
        type(field_type), intent(inout) :: f
        type(grid_type),  intent(in)    :: g

        integer :: n, m, i, j, k
        integer :: nx, ny, nz
        integer :: n_y, n_x, n_z, offset

        nx = g%nx
        ny = g%ny
        nz = g%nz
        n_y = nx*nz
        n_x = ny*nz
        n_z = nx*ny

        !$omp target teams distribute parallel do &
        !$omp& map(to: nx, ny, nz, n_y, n_x, n_z) &
        !$omp& map(tofrom: f%un(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%vn(0:nx+1,1:ny+1,0:nz+1), &
        !$omp& f%wn(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%pn(0:nx+1,1:ny,0:nz+1)) private(n,m,i,j,k,offset)
        do n = 1, 2*n_y + 2*n_x + 2*n_z
            if (n <= n_y) then
                m = n - 1
                i = modulo(m, nx) + 1
                k = m/nx + 1
                f%un(i,0,k) = -f%un(i,1,k)
                f%vn(i,1,k) = 0.0d0
                f%wn(i,0,k) = -f%wn(i,1,k)
            else if (n <= 2*n_y) then
                m = n - n_y - 1
                i = modulo(m, nx) + 1
                k = m/nx + 1
                f%un(i,ny+1,k) = -f%un(i,ny,k)
                f%vn(i,ny+1,k) = 0.0d0
                f%wn(i,ny+1,k) = -f%wn(i,ny,k)
            else if (n <= 2*n_y + n_x) then
                m = n - 2*n_y - 1
                j = modulo(m, ny) + 1
                k = m/ny + 1
                f%un(0,j,k) = f%un(nx,j,k)
                f%vn(0,j,k) = f%vn(nx,j,k)
                f%wn(0,j,k) = f%wn(nx,j,k)
                f%pn(0,j,k) = f%pn(nx,j,k)
            else if (n <= 2*n_y + 2*n_x) then
                m = n - 2*n_y - n_x - 1
                j = modulo(m, ny) + 1
                k = m/ny + 1
                f%un(nx+1,j,k) = f%un(1,j,k)
                f%vn(nx+1,j,k) = f%vn(1,j,k)
                f%wn(nx+1,j,k) = f%wn(1,j,k)
                f%pn(nx+1,j,k) = f%pn(1,j,k)
            else if (n <= 2*n_y + 2*n_x + n_z) then
                offset = 2*n_y + 2*n_x
                m = n - offset - 1
                i = modulo(m, nx) + 1
                j = m/nx + 1
                f%un(i,j,0) = f%un(i,j,nz)
                f%vn(i,j,0) = f%vn(i,j,nz)
                f%wn(i,j,0) = f%wn(i,j,nz)
                f%pn(i,j,0) = f%pn(i,j,nz)
            else
                offset = 2*n_y + 2*n_x + n_z
                m = n - offset - 1
                i = modulo(m, nx) + 1
                j = m/nx + 1
                f%un(i,j,nz+1) = f%un(i,j,1)
                f%vn(i,j,nz+1) = f%vn(i,j,1)
                f%wn(i,j,nz+1) = f%wn(i,j,1)
                f%pn(i,j,nz+1) = f%pn(i,j,1)
            end if
        end do
        !$omp end target teams distribute parallel do
    end subroutine apply_bc_redblack
#endif

end module boundary
