module boundary
    use, intrinsic :: iso_c_binding
    use :: init, only: grid_type, field_type
    implicit none

    type :: boundary_type
        logical(C_BOOL) :: isPeriodic(1:3)    
        integer(C_INT) :: bcType(1:3,0:1,0:3)     ! first index is direction(x,y,z), second index is boundary (min,max), third index is variable (p,u,v,w), 0 for Dirichlet, 1 for Neumann
        real(C_DOUBLE) :: bcValue(1:3,0:1,0:3)
    end type boundary_type

contains

    subroutine init_bc(bc)
        type(boundary_type), intent(inout) :: bc

        bc%isPeriodic(1:3) = .True.
        bc%bcType(:,:,:) = 0;  
        bc%bcValue(:,:,:) = 0;
    end subroutine init_bc


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
        !$omp& f%vn(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%wn(0:nx+1,0:ny+1,0:nz+1)) private(i,k)
#else
        !$omp target teams distribute parallel do collapse(2) &
        !$omp& map(tofrom: f%un(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%vn(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%wn(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%us(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%vs(0:nx+1,0:ny+1,0:nz+1), &
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
        !$omp& f%vn(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%wn(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%pn(0:nx+1,0:ny+1,0:nz+1)) private(j,k)
#else
        !$omp target teams distribute parallel do collapse(2) &
        !$omp& map(tofrom: f%un(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%vn(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%wn(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%us(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%vs(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%ws(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%pn(0:nx+1,0:ny+1,0:nz+1)) private(j,k)
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
        !$omp& f%vn(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%wn(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%pn(0:nx+1,0:ny+1,0:nz+1)) private(i,j)
#else
        !$omp target teams distribute parallel do collapse(2) &
        !$omp& map(tofrom: f%un(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%vn(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%us(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%vs(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%ws(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%wn(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%pn(0:nx+1,0:ny+1,0:nz+1)) private(i,j)
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
    subroutine apply_bc_redblack(f, g, bc)
        type(field_type), intent(inout)  :: f
        type(grid_type),  intent(in)     :: g
        type(boundary_type),  intent(in) :: bc

        integer :: n, m, i, j, k
        integer :: nx, ny, nz
        integer :: n_y, n_x, n_z, offset
        real(C_DOUBLE) :: dx, dy, dz
        logical(C_BOOL) :: isPeriodic(1:3)
        integer(C_INT) :: bcType(1:3,0:1,0:3)
        real(C_DOUBLE) :: bcValue(1:3,0:1,0:3)

        nx = g%nx
        ny = g%ny
        nz = g%nz
        dx = g%dx
        dy = g%dy
        dz = g%dz
        isPeriodic = bc%isPeriodic
        bcType = bc%bcType
        bcValue = bc%bcValue
        n_y = nx*nz
        n_x = ny*nz
        n_z = nx*ny

        !$omp target teams distribute parallel do &
        !$omp& map(to: nx, ny, nz, n_y, n_x, n_z, dx, dy, dz, &
        !$omp& isPeriodic(1:3), bcType(1:3,0:1,0:3), bcValue(1:3,0:1,0:3)) &
        !$omp& map(tofrom: f%un(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%vn(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%wn(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%pn(0:nx+1,0:ny+1,0:nz+1)) private(n,m,i,j,k,offset)
        do n = 1, 2*n_y + 2*n_x + 2*n_z
            if (n <= n_y) then                  ! Bottom boundary         
                m = n - 1
                i = modulo(m, nx) + 1
                k = m/nx + 1
                if (isPeriodic(2)) then
                    f%un(i,0,k) = f%un(i,ny,k)
                    f%vn(i,0,k) = f%vn(i,ny,k)
                    f%wn(i,0,k) = f%wn(i,ny,k)
                    f%pn(i,0,k) = f%pn(i,ny,k)
                else
                    if ( bcType(2,0,1)==0 ) then
                        f%un(i,0,k) = 2*bcValue(2,0,1)-f%un(i,1,k)
                    else
                        f%un(i,0,k) = f%un(i,1,k) - bcValue(2,0,1)*dy
                    end if
                    if ( bcType(2,0,2)==0 ) then
                        f%vn(i,1,k) = bcValue(2,0,2)
                    else
                        f%vn(i,1,k) = f%vn(i,2,k) - bcValue(2,0,2)*dy  ! first-order approx.
                    end if
                    if ( bcType(2,0,3)==0 ) then
                        f%wn(i,0,k) = 2*bcValue(2,0,3)-f%wn(i,1,k)
                    else
                        f%wn(i,0,k) = f%wn(i,1,k) - bcValue(2,0,3)*dy
                    end if
                end if
            else if (n <= 2*n_y) then           ! Top boundary 
                m = n - n_y - 1
                i = modulo(m, nx) + 1
                k = m/nx + 1
                if (isPeriodic(2)) then
                    f%un(i,ny+1,k) = f%un(i,1,k)
                    f%vn(i,ny+1,k) = f%vn(i,1,k)
                    f%wn(i,ny+1,k) = f%wn(i,1,k)
                    f%pn(i,ny+1,k) = f%pn(i,1,k)
                else
                    if ( bcType(2,1,1)==0 ) then
                        f%un(i,ny+1,k) = 2*bcValue(2,1,1)-f%un(i,ny,k)
                    else
                        f%un(i,ny+1,k) = f%un(i,ny,k) + bcValue(2,1,1)*dy
                    end if
                    if ( bcType(2,1,2)==0 ) then
                        f%vn(i,ny+1,k) = bcValue(2,1,2)
                    else
                        f%vn(i,ny+1,k) = f%vn(i,ny,k) + bcValue(2,1,2)*dy  ! first-order approx.
                    end if
                    if ( bcType(2,1,3)==0 ) then
                        f%wn(i,ny+1,k) = 2*bcValue(2,1,3)-f%wn(i,ny,k)
                    else
                        f%wn(i,ny+1,k) = f%wn(i,ny,k) + bcValue(2,1,3)*dy
                    end if
                end if
            else if (n <= 2*n_y + n_x) then         ! Left boundary 
                m = n - 2*n_y - 1
                j = modulo(m, ny) + 1
                k = m/ny + 1
                if (isPeriodic(1)) then
                    f%un(0,j,k) = f%un(nx,j,k)
                    f%vn(0,j,k) = f%vn(nx,j,k)
                    f%wn(0,j,k) = f%wn(nx,j,k)
                    f%pn(0,j,k) = f%pn(nx,j,k)
                else
                    if ( bcType(1,0,1)==0 ) then
                        f%un(1,j,k) = bcValue(1,0,1)
                    else
                        f%un(1,j,k) = f%un(2,j,k) - bcValue(1,0,1)*dx ! first-order approx.
                    end if
                    if ( bcType(1,0,2)==0 ) then
                        f%vn(0,j,k) = 2*bcValue(1,0,2) - f%vn(1,j,k)
                    else
                        f%vn(0,j,k) = f%vn(1,j,k) - bcValue(1,0,2)*dx
                    end if
                    if ( bcType(1,0,3)==0 ) then
                        f%wn(0,j,k) = 2*bcValue(1,0,3)-f%wn(1,j,k)
                    else
                        f%wn(0,j,k) = f%wn(1,j,k) - bcValue(1,0,3)*dx
                    end if
                end if
            else if (n <= 2*n_y + 2*n_x) then     ! Right boundary 
                m = n - 2*n_y - n_x - 1
                j = modulo(m, ny) + 1
                k = m/ny + 1
                if (isPeriodic(1)) then
                    f%un(nx+1,j,k) = f%un(1,j,k)
                    f%vn(nx+1,j,k) = f%vn(1,j,k)
                    f%wn(nx+1,j,k) = f%wn(1,j,k)
                    f%pn(nx+1,j,k) = f%pn(1,j,k)
                else
                    if ( bcType(1,1,1)==0 ) then
                        f%un(nx+1,j,k) = bcValue(1,1,1)
                    else
                        f%un(nx+1,j,k) = f%un(nx,j,k) + bcValue(1,1,1)*dx ! first-order approx.
                    end if
                    if ( bcType(1,1,2)==0 ) then
                        f%vn(nx+1,j,k) = 2*bcValue(1,1,2) - f%vn(nx,j,k)
                    else
                        f%vn(nx+1,j,k) = f%vn(nx,j,k) + bcValue(1,1,2)*dx
                    end if
                    if ( bcType(1,1,3)==0 ) then
                        f%wn(nx+1,j,k) = 2*bcValue(1,1,3)-f%wn(nx,j,k)
                    else
                        f%wn(nx+1,j,k) = f%wn(nx,j,k) + bcValue(1,1,3)*dx
                    end if
                end if
            else if (n <= 2*n_y + 2*n_x + n_z) then      ! Front boundary 
                offset = 2*n_y + 2*n_x
                m = n - offset - 1
                i = modulo(m, nx) + 1
                j = m/nx + 1
                if (isPeriodic(3)) then
                    f%un(i,j,0) = f%un(i,j,nz)
                    f%vn(i,j,0) = f%vn(i,j,nz)
                    f%wn(i,j,0) = f%wn(i,j,nz)
                    f%pn(i,j,0) = f%pn(i,j,nz)
                else
                    if ( bcType(3,0,1)==0 ) then
                        f%un(i,j,0) = 2*bcValue(3,0,1) - f%un(i,j,1)
                    else
                        f%un(i,j,0) = f%un(i,j,1) - bcValue(3,0,1)*dz
                    end if
                    if ( bcType(3,0,2)==0 ) then
                        f%vn(i,j,0) = 2*bcValue(3,0,2) - f%vn(i,j,1)
                    else
                        f%vn(i,j,0) = f%vn(i,j,1) - bcValue(3,0,2)*dz
                    end if
                    if ( bcType(3,0,3)==0 ) then
                        f%wn(i,j,1) = bcValue(3,0,3)
                    else
                        f%wn(i,j,1) = f%wn(i,j,2) - bcValue(3,0,3)*dz    ! first-order approx.
                    end if
                end if
            else                                ! back boundary
                offset = 2*n_y + 2*n_x + n_z
                m = n - offset - 1
                i = modulo(m, nx) + 1
                j = m/nx + 1
                if (isPeriodic(3)) then
                    f%un(i,j,nz+1) = f%un(i,j,1)
                    f%vn(i,j,nz+1) = f%vn(i,j,1)
                    f%wn(i,j,nz+1) = f%wn(i,j,1)
                    f%pn(i,j,nz+1) = f%pn(i,j,1)
                else
                    if ( bcType(3,1,1)==0 ) then
                        f%un(i,j,nz+1) = 2*bcValue(3,1,1) - f%un(i,j,nz)
                    else
                        f%un(i,j,nz+1) = f%un(i,j,nz) + bcValue(3,1,1)*dz
                    end if
                    if ( bcType(3,1,2)==0 ) then
                        f%vn(i,j,nz+1) = 2*bcValue(3,1,2) - f%vn(i,j,nz)
                    else
                        f%vn(i,j,nz+1) = f%vn(i,j,nz) + bcValue(3,1,2)*dz
                    end if
                    if ( bcType(3,1,3)==0 ) then
                        f%wn(i,j,nz+1) = bcValue(3,1,3)
                    else
                        f%wn(i,j,nz+1) = f%wn(i,j,nz) + bcValue(3,1,3)*dz    ! first-order approx.
                    end if
                end if
            end if
        end do
        !$omp end target teams distribute parallel do
    end subroutine apply_bc_redblack
#endif

end module boundary
