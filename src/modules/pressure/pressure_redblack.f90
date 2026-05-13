module pressure_redblack
    use, intrinsic :: iso_c_binding
    use :: init, only: grid_type, field_type
    use :: pressure_workspace, only: pressure_solver_type
    use :: ibmm, only: ibm_type

    implicit none

contains

    subroutine init_redblack_solver(ps, g)
        type(pressure_solver_type), intent(inout) :: ps
        type(grid_type), intent(in) :: g

        if (mod(g%nx,2) /= 0 .or. mod(g%nz,2) /= 0) then
            error stop "red-black pressure solver requires even nx and nz for periodic coloring"
        end if

    end subroutine init_redblack_solver

    subroutine destroy_redblack_solver(ps)
        type(pressure_solver_type), intent(inout) :: ps

    end subroutine destroy_redblack_solver

    subroutine pressure_projection_redblack(ps, f, g, dt_gamma, ibm)
        type(pressure_solver_type), intent(in) :: ps
        type(field_type), intent(inout) :: f
        type(grid_type), intent(in) :: g
        real(C_DOUBLE), intent(in) :: dt_gamma
        type(ibm_type), intent(in) :: ibm

        real(C_DOUBLE) :: phi,denom,dx,dy,dz,idx,idy,idz,idt,sor
        integer(C_INT) :: i,ip,iLo,ni,j,k,kp,iIter,parity,nx,ny,nz,nActive

        nx = g%nx
        ny = g%ny
        nz = g%nz
        dx = g%dx
        dy = g%dy
        dz = g%dz
        nActive = nx/2

        idx = 1.0_C_DOUBLE/dx
        idy = 1.0_C_DOUBLE/dy
        idz = 1.0_C_DOUBLE/dz
        idt = 1.0_C_DOUBLE/dt_gamma
        sor = ps%sor

#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do collapse(3) &
        !$omp& map(to: f%us(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%vs(0:nx+1,1:ny+1,0:nz+1), &
        !$omp& f%ws(0:nx+1,0:ny+1,0:nz+1)) &
        !$omp& map(tofrom: f%un(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%vn(0:nx+1,1:ny+1,0:nz+1), &
        !$omp& f%wn(0:nx+1,0:ny+1,0:nz+1)) &
        !$omp& private(i,j,k)
#endif
        do k = 0, nz+1
            do j = 0, ny+1
                do i = 0, nx+1
                    f%un(i,j,k) = f%us(i,j,k)
                    f%wn(i,j,k) = f%ws(i,j,k)
                    if (j == 1 .or. j == ny+1) then
                        f%vn(i,j,k) = 0.0d0
                    else if (j >= 2) then
                        f%vn(i,j,k) = f%vs(i,j,k)
                    end if
                end do
            end do
        end do
#ifdef USE_OPENMP_OFFLOAD
        !$omp end target teams distribute parallel do
#endif
        
        DO iIter=1,ps%nIter
            DO parity=0,1

#ifdef USE_OPENMP_OFFLOAD
                !$omp target teams distribute parallel do collapse(3) &
                !$omp& map(to: parity, sor, idx, idy, idz, idt, dt_gamma, &
                !$omp& ibm%coef_u(0:nx+1,0:ny+1,0:nz+1), &
                !$omp& ibm%coef_v(0:nx+1,1:ny+1,0:nz+1), &
                !$omp& ibm%coef_w(0:nx+1,0:ny+1,0:nz+1)) &
                !$omp& map(tofrom: f%un(0:nx+1,0:ny+1,0:nz+1), &
                !$omp& f%vn(0:nx+1,1:ny+1,0:nz+1), &
                !$omp& f%wn(0:nx+1,0:ny+1,0:nz+1), &
                !$omp& f%pn(0:nx+1,1:ny,0:nz+1)) &
                !$omp& private(i,ip,iLo,ni,j,k,kp,phi,denom)
#endif
                DO j=1,ny
                    DO k=1,nz
                        DO ni=1,nActive
                            kp = k + 1
                            IF (kp > nz) kp = 1
                            iLo = 1 + modulo(j+k+parity,2)
                            i = iLo + 2*(ni-1)
                            ip = i + 1
                            IF (ip > nx) ip = 1
                            denom = 2.0_C_DOUBLE*idx**2 + 2.0_C_DOUBLE*idz**2 + 2.0_C_DOUBLE*idy**2

                            phi = ( &
                                   (f%un(ip,j,k)-f%un(i,j,k))*idx &
                                 + (f%vn(i,j+1,k)-f%vn(i,j,k))*idy &
                                 + (f%wn(i,j,kp)-f%wn(i,j,k))*idz ) * (-sor/denom)

                            f%pn(i,j,k) = f%pn(i,j,k) + phi*idt

                            f%un(i,j,k) = f%un(i,j,k) - phi*idx
                            f%un(ip,j,k) = f%un(ip,j,k) + phi*idx

                            f%wn(i,j,k) = f%wn(i,j,k) - phi*idz
                            f%wn(i,j,kp) = f%wn(i,j,kp) + phi*idz

                            IF (j>1) THEN
                                f%vn(i,j,k) = f%vn(i,j,k) - phi*idy
                            END IF
                            IF (j<ny) THEN
                                f%vn(i,j+1,k) = f%vn(i,j+1,k) + phi*idy
                            END IF
                        END DO
                    END DO
                END DO
#ifdef USE_OPENMP_OFFLOAD
                !$omp end target teams distribute parallel do
#endif
            END DO
        END DO 

#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do collapse(3) &
        !$omp& map(to: dt_gamma, &
        !$omp& ibm%coef_u(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& ibm%coef_v(0:nx+1,1:ny+1,0:nz+1), &
        !$omp& ibm%coef_w(0:nx+1,0:ny+1,0:nz+1)) &
        !$omp& map(tofrom: f%un(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%vn(0:nx+1,1:ny+1,0:nz+1), &
        !$omp& f%wn(0:nx+1,0:ny+1,0:nz+1)) &
        !$omp& private(i,j,k)
#endif
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    f%un(i,j,k) = f%un(i,j,k) / (1.0d0 + dt_gamma*ibm%coef_u(i,j,k))
                    f%wn(i,j,k) = f%wn(i,j,k) / (1.0d0 + dt_gamma*ibm%coef_w(i,j,k))
                    if (j >= 2) then
                        f%vn(i,j,k) = f%vn(i,j,k) / (1.0d0 + dt_gamma*ibm%coef_v(i,j,k))
                    end if
                end do
            end do
        end do
#ifdef USE_OPENMP_OFFLOAD
        !$omp end target teams distribute parallel do
#endif

    end subroutine pressure_projection_redblack

end module pressure_redblack
