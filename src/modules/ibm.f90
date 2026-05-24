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
    use :: init, only: dns_type, grid_type, VAR_U, VAR_V, VAR_W
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

        real(C_DOUBLE), allocatable :: coef(:,:,:,:)

    end type ibm_type

contains

!========================
! INITIALIZE IBM
!========================
    subroutine init_ibm(ibm, dns, g)
        type(ibm_type), intent(inout) :: ibm
        type(dns_type), intent(in)    :: dns
        type(grid_type), intent(in)   :: g
        integer :: nx, ny, nz

        nx = int(dns%localSize(1,2))
        ny = int(dns%localSize(2,2))
        nz = int(dns%localSize(3,2))

        ibm%n_wave_x = 1
        ibm%n_wave_z = 1
        ibm%amp_x = 5*g%havg(2)
        ibm%amp_z = 5*g%havg(2)
        ibm%phase_x = 0.0d0
        ibm%phase_z = 0.0d0

        allocate(ibm%coef(0:nx+1,0:ny+1,0:nz+1,VAR_U:VAR_W))
        ibm%coef = 0.0d0
    end subroutine init_ibm

    subroutine enter_ibm_data(ibm, dns)
        type(ibm_type), intent(inout) :: ibm
        type(dns_type), intent(in)    :: dns

#ifdef USE_OPENMP_OFFLOAD
        !$omp target enter data map(to: ibm)
        !$omp target enter data map(to: ibm%coef)
#endif
    end subroutine enter_ibm_data

    subroutine exit_ibm_data(ibm, dns)
        type(ibm_type), intent(inout) :: ibm
        type(dns_type), intent(in)    :: dns

#ifdef USE_OPENMP_OFFLOAD
        !$omp target exit data map(delete: ibm%coef)
        !$omp target exit data map(delete: ibm)
#endif
    end subroutine exit_ibm_data


    logical function isInBody(xIN, ibm, dns, g)
        implicit none

        real(C_DOUBLE), intent(in) :: xIN(1:3)
        type(ibm_type), intent(in) :: ibm
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g

        real(C_DOUBLE), parameter :: pi = 3.141592653589793d0
        real(C_DOUBLE) :: y_body

        y_body = ibm%amp_x * 0.5d0 * &
                 (1.0d0 + sin(2.0d0*pi*real(ibm%n_wave_x,C_DOUBLE)*xIN(1)/dns%leng(1) + ibm%phase_x)) + 2*g%havg(2)

        isInBody = (xIN(2) < y_body)

    end function isInBody

    subroutine bisection(xAin,xB,ibm,dns,g)
        real(C_DOUBLE), intent(in) :: xAin(1:3)
        real(C_DOUBLE), intent(inout):: xB(1:3)
        type(ibm_type), intent(in) :: ibm
        type(dns_type), intent(in) :: dns
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

            la = isInBody(xa,ibm,dns,g)
            lm = isInBody(xm,ibm,dns,g)
            IF (lm .eqv. la) THEN
                xa = xm
            ELSE
                xb = xm
            END IF
        END DO
        xb = xm
    end subroutine  bisection


    subroutine set_ibm_coeff(dns, g, ibm, var)
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        type(ibm_type), intent(inout) :: ibm
        integer(C_INT), intent(in) :: var

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

        if (var < VAR_U .or. var > VAR_W) error stop "invalid IBM coefficient variable"

        ibm%coef(:,:,:,var) = 0.0d0
        if (.not. dns%ibm_enabled) return

        do iz = lbound(ibm%coef,3), ubound(ibm%coef,3)
            do iy = lbound(ibm%coef,2), ubound(ibm%coef,2)
                do ix = lbound(ibm%coef,1), ubound(ibm%coef,1)
                    xA(1) = g%x(ix,var)
                    xA(2) = g%y(iy,var)
                    xA(3) = g%z(iz,var)
                    if (isInBody(xA, ibm, dns, g)) then
                        ibm%coef(ix,iy,iz,var) = SOLID
#ifdef USE_IBM_SECONDORDER
                    else
                        do iN = 1,6
                            xB(1) = g%x(ix+neighbours(1,iN),var)
                            xB(2) = g%y(iy+neighbours(2,iN),var)
                            xB(3) = g%z(iz+neighbours(3,iN),var)
                            d0 = norm2(xB-xA)
                            if (isInBody(xB, ibm, dns, g)) then
                                call bisection(xA,xB,ibm,dns,g)
                                d = norm2(xB-xA)
                                ibm%coef(ix,iy,iz,var) = ibm%coef(ix,iy,iz,var) + ((d0-d)/d)/d0**2
                            end if
                        end do
#endif
                    end if
                end do
            end do
        end do
        ibm%coef(:,:,:,var) = ibm%coef(:,:,:,var)/dns%re
    end subroutine set_ibm_coeff
    
end module ibmm
