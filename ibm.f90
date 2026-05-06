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
    allocate(ibm%coef_v(0:g%nx+1,1:g%ny+1,0:g%nz+1))
    allocate(ibm%coef_w(0:g%nx+1,0:g%ny+1,0:g%nz+1))


end subroutine init_ibm

    subroutine enter_ibm_data(ibm, g)
        type(ibm_type), intent(inout) :: ibm
        type(grid_type), intent(in)   :: g

#ifdef USE_OPENMP_OFFLOAD
        !$omp target enter data map(to: ibm)
        !$omp target enter data map(to: &
        !$omp& ibm%coef_u(0:g%nx+1,0:g%ny+1,0:g%nz+1), &
        !$omp& ibm%coef_v(0:g%nx+1,1:g%ny+1,0:g%nz+1), &
        !$omp& ibm%coef_w(0:g%nx+1,0:g%ny+1,0:g%nz+1))
#endif
    end subroutine enter_ibm_data

    subroutine exit_ibm_data(ibm, g)
        type(ibm_type), intent(inout) :: ibm
        type(grid_type), intent(in)   :: g

#ifdef USE_OPENMP_OFFLOAD
        !$omp target exit data map(delete: &
        !$omp& ibm%coef_u(0:g%nx+1,0:g%ny+1,0:g%nz+1), &
        !$omp& ibm%coef_v(0:g%nx+1,1:g%ny+1,0:g%nz+1), &
        !$omp& ibm%coef_w(0:g%nx+1,0:g%ny+1,0:g%nz+1))
        !$omp target exit data map(delete: ibm)
#endif
    end subroutine exit_ibm_data


    logical function isInBody(x, y, z, ibm, g)
        implicit none

        real(C_DOUBLE), intent(in) :: x, y, z
        type(ibm_type), intent(in) :: ibm
        type(grid_type), intent(in) :: g

        real(C_DOUBLE), parameter :: pi = 3.141592653589793d0
        real(C_DOUBLE) :: y_body

        y_body = ibm%amp_x * 0.5d0 * &
                 (1.0d0 + sin(2.0d0*pi*real(ibm%n_wave_x,C_DOUBLE)*x/g%lx + ibm%phase_x)) + 2*g%dy !* &
!                 ibm%amp_z * 0.5d0 * &
                !(1.0d0 + sin(2.0d0*pi*real(ibm%n_wave_z,C_DOUBLE)*z/g%lz + ibm%phase_z))

        isInBody = (y < y_body)

    end function isInBody


    subroutine set_ibm_coeff(g, ibm, coeff, dix, diy, diz)
        implicit none

        type(grid_type), intent(in) :: g
        type(ibm_type), intent(in) :: ibm

        real(C_DOUBLE), intent(inout) :: coeff(:,:,:)

        integer, intent(in) :: dix, diy, diz
        integer :: ix, iy, iz
        real(C_DOUBLE) :: x, y, z

        coeff = 0.0d0

        do iz = 1, size(coeff,3)
            do iy = 1, size(coeff,2)
                do ix = 1, size(coeff,1)

                    x = (real(ix,C_DOUBLE) - real(dix,C_DOUBLE)*0.5d0        )*g%dx
                    y = (real(iy,C_DOUBLE) - real(diy,C_DOUBLE)*0.5d0 + 0.5d0)*g%dy
                    z = (real(iz,C_DOUBLE) - real(diz,C_DOUBLE)*0.5d0        )*g%dz

                    if (isInBody(x, y, z, ibm, g)) then
                        coeff(ix,iy,iz) = SOLID
                    end if

                end do
            end do
        end do

    end subroutine set_ibm_coeff
    
    subroutine apply_ibm(field, coeff, g)
        implicit none

        real(C_DOUBLE), intent(inout) :: field(:,:,:)
        real(C_DOUBLE), intent(in)    :: coeff(:,:,:)
        type(grid_type), intent(in) :: g

        integer :: ix, iy, iz
        real(C_DOUBLE) :: dt

        dt = g%dt

        !$omp target teams distribute parallel do collapse(3) &
        !$omp& map(to: coeff(1:size(coeff,1),1:size(coeff,2),1:size(coeff,3))) &
        !$omp& map(tofrom: field(1:size(field,1),1:size(field,2),1:size(field,3))) &
        !$omp& private(ix,iy,iz)
        do iz = 1, size(field,3)
       	   do iy = 1, size(field,2)
	      do ix = 1, size(field,1)

	        field(ix,iy,iz) = field(ix,iy,iz) / &
	                          (1.0d0 + dt*coeff(ix,iy,iz))

	      end do
	   end do
        end do
        !$omp end target teams distribute parallel do

end subroutine apply_ibm
end module ibmm
