module pressure_redblack
    use, intrinsic :: iso_c_binding
    use :: init, only: grid_type, field_type
    use :: pressure_workspace, only: pressure_solver_type
    use :: ibmm, only: ibm_type
    use :: boundary, only: boundary_type

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

    subroutine pressure_projection_redblack(ps, f, g, dt_gamma, ibm, bc)
        type(pressure_solver_type), intent(in) :: ps
        type(field_type), intent(inout) :: f
        type(grid_type), intent(in) :: g
        real(C_DOUBLE), intent(in) :: dt_gamma
        type(ibm_type), intent(in) :: ibm
        type(boundary_type), intent(in) :: bc

        real(C_DOUBLE) :: phi,denom,dx,dy,dz,idx,idy,idz,idx2,idy2,idz2,idt,sor
        real(C_DOUBLE) :: mu_u_i,mu_u_ip,mu_v_j,mu_v_jp,mu_w_k,mu_w_kp
        integer(C_INT) :: i,ip,iLo,ni,j,jp,k,kp,iIter,parity,nx,ny,nz,nActive
        logical(C_BOOL) :: isPeriodic(1:3)
        integer(C_INT) :: bcType(1:3,0:1,0:3)
        real(C_DOUBLE) :: bcValue(1:3,0:1,0:3)

        nx = g%nx
        ny = g%ny
        nz = g%nz
        dx = g%dx
        dy = g%dy
        dz = g%dz
        nActive = nx/2

        isPeriodic = bc%isPeriodic
        bcType = bc%bcType
        bcValue = bc%bcValue

        idx = 1.0_C_DOUBLE/dx
        idy = 1.0_C_DOUBLE/dy
        idz = 1.0_C_DOUBLE/dz
        idx2 = idx*idx
        idy2 = idy*idy
        idz2 = idz*idz
        idt = 1.0_C_DOUBLE/dt_gamma
        sor = ps%sor

#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do collapse(3) &
        !$omp& map(to: f%us(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%vs(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%ws(0:nx+1,0:ny+1,0:nz+1)) &
        !$omp& map(tofrom: f%un(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%vn(0:nx+1,0:ny+1,0:nz+1), &
        !$omp& f%wn(0:nx+1,0:ny+1,0:nz+1)) &
        !$omp& private(i,j,k)
#endif
        do k = 0, nz+1
            do j = 0, ny+1
                do i = 0, nx+1
                        f%un(i,j,k) = f%us(i,j,k)
                        f%vn(i,j,k) = f%vs(i,j,k)
                        f%wn(i,j,k) = f%ws(i,j,k)
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
                !$omp& map(to: parity, sor, dx, dy, dz, idx, idy, idz, idx2, idy2, idz2, idt, dt_gamma, &
                !$omp& isPeriodic(1:3), bcType(1:3,0:1,0:3), bcValue(1:3,0:1,0:3), &
                !$omp& ibm%coef_u(0:nx+1,0:ny+1,0:nz+1), &
                !$omp& ibm%coef_v(0:nx+1,0:ny+1,0:nz+1), &
                !$omp& ibm%coef_w(0:nx+1,0:ny+1,0:nz+1)) &
                !$omp& map(tofrom: f%un(0:nx+1,0:ny+1,0:nz+1), &
                !$omp& f%vn(0:nx+1,0:ny+1,0:nz+1), &
                !$omp& f%wn(0:nx+1,0:ny+1,0:nz+1), &
                !$omp& f%pn(0:nx+1,0:ny+1,0:nz+1)) &
                !$omp& private(i,ip,iLo,ni,j,jp,k,kp,phi,denom, &
                !$omp& mu_u_i,mu_u_ip,mu_v_j,mu_v_jp,mu_w_k,mu_w_kp)
#endif
            DO k=1,nz
                DO j=1,ny
                        DO ni=1,nActive

                            ! Compute indices (account for periodicity)
                            jp = j + 1
                            IF (isPeriodic(2) .and. jp > ny) jp = 1
                            kp = k + 1
                            IF (isPeriodic(3) .and. kp > nz) kp = 1
                            iLo = 1 + modulo(j+k+parity,2)
                            i = iLo + 2*(ni-1)
                            ip = i + 1
                            IF (isPeriodic(1) .and. ip > nx) ip = 1

                            mu_u_i  = 1.0d0/(1.0d0+dt_gamma*ibm%coef_u(i,j,k))
                            mu_u_ip = 1.0d0/(1.0d0+dt_gamma*ibm%coef_u(ip,j,k))
                            mu_v_j  = 1.0d0/(1.0d0+dt_gamma*ibm%coef_v(i,j,k))
                            mu_v_jp = 1.0d0/(1.0d0+dt_gamma*ibm%coef_v(i,jp,k))
                            mu_w_k  = 1.0d0/(1.0d0+dt_gamma*ibm%coef_w(i,j,k))
                            mu_w_kp = 1.0d0/(1.0d0+dt_gamma*ibm%coef_w(i,j,kp))

                            denom=(mu_u_i+mu_u_ip)*idx2 & 
                                 +(mu_w_k+mu_w_kp)*idz2 &
                                 +(mu_v_j+mu_v_jp)*idy2

                            ! Left and right boundary
                            if (.not. isPeriodic(1) ) then
                            if (i<2) then
                                if ( bcType(1,0,1)==0 ) then
                                    f%un(i,j,k) = bcValue(1,0,1)
                                else
                                    f%un(i,j,k) = f%un(ip,j,k) - bcValue(1,0,1)*dx ! first-order approx.
                                end if
                            end if
                            if (i>=nx) then
                                if ( bcType(1,1,1)==0 ) then
                                    f%un(ip,j,k) = bcValue(1,1,1)
                                else
                                    f%un(ip,j,k) = f%un(nx,j,k) + bcValue(1,1,1)*dx ! first-order approx.
                                end if
                            end if
                            end if

                            ! Bottom and top boundary
                            if (.not. isPeriodic(2)) then
                            if (j<2) then
                                if ( bcType(2,0,2)==0 ) then
                                    f%vn(i,j,k) = bcValue(2,0,2)
                                else
                                    f%vn(i,j,k) = f%vn(i,jp,k) - bcValue(2,0,2)*dy  ! first-order approx.
                                end if
                            end if
                            if (j>=ny) then
                                if ( bcType(2,1,2)==0 ) then
                                    f%vn(i,jp,k) = bcValue(2,1,2)
                                else
                                    f%vn(i,jp,k) = f%vn(i,j,k) + bcValue(2,1,2)*dy  ! first-order approx.
                                end if
                            end if
                            end if

                            ! Front and back boundary
                            if (.not. isPeriodic(3)) then
                            if (k<2) then
                                if ( bcType(3,0,3)==0 ) then
                                    f%wn(i,j,k) = bcValue(3,0,3)
                                else
                                    f%wn(i,j,k) = f%wn(i,j,kp) - bcValue(3,0,3)*dz    ! first-order approx.
                                end if
                            end if
                            if (k>=nz) then
                                if ( bcType(3,1,3)==0 ) then
                                    f%wn(i,j,kp) = bcValue(3,1,3)
                                else
                                    f%wn(i,j,kp) = f%wn(i,j,k) + bcValue(3,1,3)*dz    ! first-order approx.
                                end if
                            end if
                            end if
                            
                            ! Pressure point interation (account for pressure BC)
                            phi = ( &
                                   (f%un(ip,j,k)-f%un(i,j,k))*idx &
                                 + (f%vn(i,jp,k)-f%vn(i,j,k))*idy &
                                 + (f%wn(i,j,kp)-f%wn(i,j,k))*idz ) * (-sor/denom)

                            f%pn(i,j,k) = f%pn(i,j,k) + phi*idt
                            
                            f%un(i,j,k) = f%un(i,j,k) - phi*idx*mu_u_i
                            IF ( (.not. isPeriodic(1)) .and. i==1 ) THEN
                                IF ( bcType(1,0,0)==0 ) THEN
                                    f%un(i,j,k) = f%un(i,j,k) + bcValue(1,0,0)*idx*mu_u_i
                                ELSE
                                    f%un(i,j,k) = f%un(i,j,k) + (phi*idx-bcValue(1,0,0))*mu_u_i
                                END IF
                            END IF

                            f%un(ip,j,k) = f%un(ip,j,k) + phi*idx*mu_u_ip
                            IF ( (.not. isPeriodic(1)) .and. i==nx ) THEN
                                IF ( bcType(1,1,0)==0 ) THEN
                                    f%un(ip,j,k) = f%un(ip,j,k) - bcValue(1,1,0)*idx*mu_u_ip
                                ELSE
                                    f%un(ip,j,k) = f%un(ip,j,k) + (bcValue(1,1,0)-phi*idx)*mu_u_ip
                                END IF
                            END IF

                            f%vn(i,j,k) = f%vn(i,j,k) - phi*idy*mu_v_j
                            IF ( (.not. isPeriodic(2)) .and. j==1 ) THEN
                                IF ( bcType(2,0,0)==0 ) THEN
                                    f%vn(i,j,k) = f%vn(i,j,k) + bcValue(2,0,0)*idy*mu_v_j
                                ELSE
                                    f%vn(i,j,k) = f%vn(i,j,k) + (phi*idy-bcValue(2,0,0))*mu_v_j
                                END IF
                            END IF

                            f%vn(i,jp,k) = f%vn(i,jp,k) + phi*idy*mu_v_jp
                            IF ( (.not. isPeriodic(2)) .and. j==ny ) THEN
                                IF ( bcType(2,1,0)==0 ) THEN
                                    f%vn(i,jp,k) = f%vn(i,jp,k) - bcValue(2,1,0)*idy*mu_v_jp
                                ELSE
                                    f%vn(i,jp,k) = f%vn(i,jp,k) +(bcValue(2,1,0)- phi*idy)*mu_v_jp
                                END IF
                            END IF

                            f%wn(i,j,k) = f%wn(i,j,k) - phi*idz*mu_w_k
                            IF ( (.not. isPeriodic(3)) .and. k==1 ) THEN
                                IF ( bcType(3,0,0)==0 ) THEN
                                    f%wn(i,j,k) = f%wn(i,j,k) + bcValue(3,0,0)*idz*mu_w_k
                                ELSE
                                    f%wn(i,j,k) = f%wn(i,j,k) + (phi*idz-bcValue(3,0,0))*mu_w_k
                                END IF
                            END IF
                            f%wn(i,j,kp) = f%wn(i,j,kp) + phi*idz*mu_w_kp
                            IF ( (.not. isPeriodic(3)) .and. k==nz ) THEN
                                IF ( bcType(3,1,0)==0 ) THEN
                                    f%wn(i,j,kp) = f%wn(i,j,kp) - bcValue(3,1,0)*idz
                                ELSE
                                    f%wn(i,j,kp) = f%wn(i,j,kp) +(bcValue(3,1,0)- phi*idz)*mu_w_kp
                                END IF
                            END IF

                        END DO
                    END DO
                END DO
#ifdef USE_OPENMP_OFFLOAD
                !$omp end target teams distribute parallel do
#endif
            END DO
        END DO 



    end subroutine pressure_projection_redblack

end module pressure_redblack
