!--------------------------!
!                          !
!    Immersed Boundary     !
!         Module           !
!                          !
!--------------------------! 
! 
! authors: Dr.-Ing. Davide Gatti
!          B.Sc. Ahmet Cumhur
! 
! date:    28.04.26
! 


module ibmm
    use, intrinsic :: iso_c_binding
    use :: init, only: grid_type
    implicit none

    real(C_DOUBLE), parameter :: SOLID = 1.0d30
    real(C_DOUBLE), parameter :: DEFAULT_TOL = 1.0d-10
    integer(C_INT), parameter :: MAX_ITER = 200

    !========================
    ! IBM TYPE
    !========================
    type :: ibm_type
        integer :: n_wave_x, n_wave_z
        real(C_DOUBLE) :: amp_x, phase_x
        real(C_DOUBLE) :: amp_z, phase_z

        real(C_DOUBLE), allocatable :: coef_u(:,:,:), coef_v(:,:,:), coef_w(:,:,:)

    end type ibm_type

contains

!========================
! INITIALIZE IBM
!========================
    subroutine init_ibm(ibm, g)
        type(ibm_type), intent(inout) :: ibm
        type(grid_type), intent(in)   :: g

        ibm%n_wave_x = 1
        ibm%n_wave_z = 1
        ibm%amp_x = 5*g%dy
        ibm%amp_z = 5*g%dy
        ibm%phase_x = 0.0d0
        ibm%phase_z = 0.0d0

        allocate(ibm%coef_u(0:g%nx+1,0:g%ny+1,0:g%nz+1))
        allocate(ibm%coef_v(0:g%nx+1,0:g%ny+1,0:g%nz+1))
        allocate(ibm%coef_w(0:g%nx+1,0:g%ny+1,0:g%nz+1))
    end subroutine init_ibm

    subroutine enter_ibm_data(ibm, g)
        type(ibm_type), intent(inout) :: ibm
        type(grid_type), intent(in)   :: g

#ifdef USE_OPENMP_OFFLOAD
        !$omp target enter data map(to: ibm)
        !$omp target enter data map(to: &
        !$omp& ibm%coef_u(0:g%nx+1,0:g%ny+1,0:g%nz+1), &
        !$omp& ibm%coef_v(0:g%nx+1,0:g%ny+1,0:g%nz+1), &
        !$omp& ibm%coef_w(0:g%nx+1,0:g%ny+1,0:g%nz+1))
#endif
    end subroutine enter_ibm_data

    subroutine exit_ibm_data(ibm, g)
        type(ibm_type), intent(inout) :: ibm
        type(grid_type), intent(in)   :: g

#ifdef USE_OPENMP_OFFLOAD
        !$omp target exit data map(delete: &
        !$omp& ibm%coef_u(0:g%nx+1,0:g%ny+1,0:g%nz+1), &
        !$omp& ibm%coef_v(0:g%nx+1,0:g%ny+1,0:g%nz+1), &
        !$omp& ibm%coef_w(0:g%nx+1,0:g%ny+1,0:g%nz+1))
        !$omp target exit data map(delete: ibm)
#endif
    end subroutine exit_ibm_data


    logical function isInBody(xIN, ibm, g)
        implicit none

        real(C_DOUBLE), intent(in) :: xIN(1:3)
        type(ibm_type), intent(in) :: ibm
        type(grid_type), intent(in) :: g

        real(C_DOUBLE), parameter :: pi = 3.141592653589793d0
        real(C_DOUBLE) :: y_body

        y_body = ibm%amp_x * 0.5d0 * &
                 (1.0d0 + sin(2.0d0*pi*real(ibm%n_wave_x,C_DOUBLE)*xIN(1)/g%lx + ibm%phase_x)) + 2*g%dy 

        isInBody = (xIN(2) < y_body)

    end function isInBody

    subroutine bisection(xAin,xB,ibm,g) 
        real(C_DOUBLE), intent(in) :: xAin(1:3)
        real(C_DOUBLE), intent(inout):: xB(1:3)
        type(ibm_type), intent(in) :: ibm
        type(grid_type), intent(in) :: g
        real(C_DOUBLE) :: xA(1:3),xM(1:3)
        logical :: la, lb, lm

        integer(C_INT) :: it

        xA = xAin

        DO it=1,MAX_ITER

            xm = 0.5*(xa+xb)

            if (NORM2(xb-xa) < DEFAULT_TOL) then
                exit
            end if

            la = isInBody(xa,ibm,g)
            lm = isInBody(xm,ibm,g)
            IF (lm .eqv. la) THEN
                xa = xm
            ELSE
                xb = xm
            END IF
        END DO
        xb = xm
    end subroutine  bisection


    subroutine set_ibm_coeff(g, ibm, coeff, dix, diy, diz)
        implicit none

        type(grid_type), intent(in) :: g
        type(ibm_type), intent(in) :: ibm

        real(C_DOUBLE), intent(inout) :: coeff(:,:,:)

        integer, intent(in) :: dix, diy, diz
        integer :: ix, iy, iz
        real(C_DOUBLE) :: xA(1:3)
#ifdef USE_IBM_SECONDORDER
        integer(C_INT) :: neighbours(1:3,1:6), iN
        real(C_DOUBLE) :: xB(1:3), d0, d

        neighbours(1:3,1) = (/-1, 0, 0 /)
        neighbours(1:3,2) = (/ 1, 0, 0 /)
        neighbours(1:3,3) = (/ 0,-1, 0 /)
        neighbours(1:3,4) = (/ 0, 1, 0 /)
        neighbours(1:3,5) = (/ 0, 0,-1 /)
        neighbours(1:3,6) = (/ 0, 0, 1 /)
#endif

        coeff = 0.0d0

        do iz = 1, size(coeff,3)
            do iy = 1, size(coeff,2)
                do ix = 1, size(coeff,1)
                    xA(1) = (real(ix,C_DOUBLE) - real(dix,C_DOUBLE)*0.5d0 - 1.5d0)*g%dx
                    xA(2) = (real(iy,C_DOUBLE) - real(diy,C_DOUBLE)*0.5d0 - 1.5d0)*g%dy
                    xA(3) = (real(iz,C_DOUBLE) - real(diz,C_DOUBLE)*0.5d0 - 1.5d0)*g%dz
                    if (isInBody(xA, ibm, g)) then
                        coeff(ix,iy,iz) = SOLID
#ifdef USE_IBM_SECONDORDER
                    else
                        do iN = 1,6
                            xB(1) = (real(ix+neighbours(1,iN),C_DOUBLE) - real(dix,C_DOUBLE)*0.5d0 - 1.5d0)*g%dx
                            xB(2) = (real(iy+neighbours(2,iN),C_DOUBLE) - real(diy,C_DOUBLE)*0.5d0 - 1.5d0)*g%dy
                            xB(3) = (real(iz+neighbours(3,iN),C_DOUBLE) - real(diz,C_DOUBLE)*0.5d0 - 1.5d0)*g%dz
                            !coordMod will be put here
                            d0 = norm2(xB-xA)
                            if (isInBody(xB, ibm, g)) then
                                call bisection(xA,xB,ibm,g)
                                d = norm2(xB-xA)
                                coeff(ix,iy,iz) =  coeff(ix,iy,iz) + ((d0-d)/d)/d0**2  ! adjust for noneq. grid
                            end if
                        end do
#endif
                    end if

                end do
            end do
        end do
        coeff = coeff/g%re
    end subroutine set_ibm_coeff
    
end module ibmm
