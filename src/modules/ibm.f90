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

!$omp declare target(isInBody, distance3, bisection, add_neighbor_coeff)

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
                 (1.0d0 + sin(2.0d0*pi*real(ibm%n_wave_x,C_DOUBLE)*xIN(1)/dns%leng(1) + ibm%phase_x)) + &
                 2.0d0*g%havg(2)
        isInBody = (xIN(2) < y_body)
    end function isInBody

    real(C_DOUBLE) function distance3(xA, xB) result(d)
        real(C_DOUBLE), intent(in) :: xA(1:3), xB(1:3)

        d = sqrt((xB(1)-xA(1))**2 + (xB(2)-xA(2))**2 + (xB(3)-xA(3))**2)
    end function distance3

    subroutine bisection(xAin,xB,ibm,dns,g)
        real(C_DOUBLE), intent(in) :: xAin(1:3)
        real(C_DOUBLE), intent(inout):: xB(1:3)
        type(ibm_type), intent(in) :: ibm
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g

        real(C_DOUBLE) :: xA(1:3), xM(1:3)
        logical :: la, lm
        integer(C_INT) :: it

        xA = xAin
        xM = xA
        la = isInBody(xA, ibm, dns, g)

        do it = 1, MAX_ITER
            xM = 0.5d0*(xA + xB)
            if (distance3(xA, xB) < DEFAULT_TOL) exit

            lm = isInBody(xM, ibm, dns, g)
            if (lm .eqv. la) then
                xA = xM
            else
                xB = xM
            end if
        end do
        xB = xM
    end subroutine bisection

    subroutine add_neighbor_coeff(coeff, xA, xB, ibm, dns, g)
        real(C_DOUBLE), intent(inout) :: coeff
        real(C_DOUBLE), intent(in) :: xA(1:3)
        real(C_DOUBLE), intent(inout) :: xB(1:3)
        type(ibm_type), intent(in) :: ibm
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g

        real(C_DOUBLE) :: d0, d

        if (isInBody(xB, ibm, dns, g)) then
            d0 = distance3(xA, xB)
            call bisection(xA, xB, ibm, dns, g)
            d = distance3(xA, xB)
            coeff = coeff + ((d0-d)/d)/d0**2
        end if
    end subroutine add_neighbor_coeff

    subroutine set_ibm_coeff(dns, g, ibm, var)
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        type(ibm_type), intent(inout) :: ibm
        integer(C_INT), intent(in) :: var

        integer :: ix, iy, iz
        integer :: ilo, ihi, jlo, jhi, klo, khi
        real(C_DOUBLE) :: xA(1:3), xB(1:3), coeff
        real(C_DOUBLE) :: re_inv, solid_coef
        logical(C_BOOL) :: enabled

        if (var < VAR_U .or. var > VAR_W) error stop "invalid IBM coefficient variable"

        ilo = lbound(ibm%coef,1)
        ihi = ubound(ibm%coef,1)
        jlo = lbound(ibm%coef,2)
        jhi = ubound(ibm%coef,2)
        klo = lbound(ibm%coef,3)
        khi = ubound(ibm%coef,3)

        enabled = dns%ibm_enabled
        re_inv = 1.0d0/dns%re
        solid_coef = SOLID*re_inv
#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do collapse(3) &
        !$omp& map(to: var, enabled, re_inv, solid_coef, ilo, ihi, jlo, jhi, klo, khi, dns, ibm, g, &
        !$omp& g%x, g%y, g%z) &
        !$omp& map(tofrom: ibm%coef) &
        !$omp& private(ix,iy,iz,xA,xB,coeff)
#endif
        do iz = klo, khi
            do iy = jlo, jhi
                do ix = ilo, ihi
                    ibm%coef(ix,iy,iz,var) = 0.0d0
                    if (.not. enabled) cycle

                    xA(1) = g%x(ix,var)
                    xA(2) = g%y(iy,var)
                    xA(3) = g%z(iz,var)
                    if (isInBody(xA, ibm, dns, g)) then
                        ibm%coef(ix,iy,iz,var) = solid_coef
#ifdef USE_IBM_SECONDORDER
                    else
                        coeff = 0.0d0

                        xB(1) = g%x(ix-1,var); xB(2) = xA(2);          xB(3) = xA(3)
                        call add_neighbor_coeff(coeff, xA, xB, ibm, dns, g)
                        xB(1) = g%x(ix+1,var); xB(2) = xA(2);          xB(3) = xA(3)
                        call add_neighbor_coeff(coeff, xA, xB, ibm, dns, g)
                        xB(1) = xA(1);          xB(2) = g%y(iy-1,var); xB(3) = xA(3)
                        call add_neighbor_coeff(coeff, xA, xB, ibm, dns, g)
                        xB(1) = xA(1);          xB(2) = g%y(iy+1,var); xB(3) = xA(3)
                        call add_neighbor_coeff(coeff, xA, xB, ibm, dns, g)
                        xB(1) = xA(1);          xB(2) = xA(2);          xB(3) = g%z(iz-1,var)
                        call add_neighbor_coeff(coeff, xA, xB, ibm, dns, g)
                        xB(1) = xA(1);          xB(2) = xA(2);          xB(3) = g%z(iz+1,var)
                        call add_neighbor_coeff(coeff, xA, xB, ibm, dns, g)

                        ibm%coef(ix,iy,iz,var) = coeff*re_inv
#endif
                    end if
                end do
            end do
        end do
#ifdef USE_OPENMP_OFFLOAD
        !$omp end target teams distribute parallel do
#endif
    end subroutine set_ibm_coeff
    
end module ibmm
