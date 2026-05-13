module pressure_solver
    use, intrinsic :: iso_c_binding
    use :: init, only: grid_type, field_type
    use :: pressure_workspace, only: pressure_solver_type
    use :: ibmm, only: ibm_type
#ifdef USE_REDBLACK
    use :: pressure_redblack, only: &
        init_redblack_solver, &
        pressure_projection_redblack, &
        destroy_redblack_solver
#else
    use :: pressure_fft, only: init_tridiag_coefficients, deallocate_workspace_arrays
#ifdef USE_CUFFT
    use :: pressure_backend_cufft, only: &
        allocate_backend_workspace => allocate_cufft_workspace, &
        create_backend_plans => create_cufft_plans, &
        solve_pressure_fft_backend => solve_pressure_cufft, &
        destroy_backend_workspace => destroy_cufft_workspace
#else
    use :: pressure_backend_fftw, only: &
        allocate_backend_workspace => allocate_fftw_workspace, &
        create_backend_plans => create_fftw_plans, &
        solve_pressure_fft_backend => solve_pressure_fftw, &
        destroy_backend_workspace => destroy_fftw_workspace
#endif
#endif
    implicit none

    private
    public :: pressure_solver_type, init_pressure_solver, pressure_projection, destroy_pressure_solver

contains

subroutine init_pressure_solver(ps, g)
    type(pressure_solver_type), intent(inout) :: ps
    type(grid_type), intent(in) :: g

#ifdef USE_REDBLACK
    call init_redblack_solver(ps, g)
#else
    call allocate_backend_workspace(ps, g)
    call create_backend_plans(ps, g)
    call init_tridiag_coefficients(ps, g)
#endif
end subroutine init_pressure_solver

subroutine pressure_projection(ps, f, g, dt_gamma, ibm)
    type(pressure_solver_type), intent(inout) :: ps
    type(field_type), intent(inout) :: f
    type(grid_type), intent(in) :: g
    real(C_DOUBLE), intent(in) :: dt_gamma
    type(ibm_type), intent(in) :: ibm

#ifdef USE_REDBLACK
    call pressure_projection_redblack(ps, f, g, dt_gamma, ibm)
#else
    call build_fft_rhs(ps, f, g, dt_gamma)
    call solve_pressure_fft_backend(g, ps)
    call apply_fft_pressure_correction(ps, f, g, dt_gamma, ibm)
#endif
end subroutine pressure_projection

subroutine destroy_pressure_solver(ps)
    type(pressure_solver_type), intent(inout) :: ps

#ifdef USE_REDBLACK
    call destroy_redblack_solver(ps)
#else
    call destroy_backend_workspace(ps)
    call deallocate_workspace_arrays(ps)
#endif
end subroutine destroy_pressure_solver

#ifndef USE_REDBLACK
subroutine build_fft_rhs(ps, f, g, dt_gamma)
    type(pressure_solver_type), intent(inout) :: ps
    type(field_type), intent(inout) :: f
    type(grid_type), intent(in) :: g
    real(C_DOUBLE), intent(in) :: dt_gamma

    integer :: i,j,k
    integer :: nx, ny, nz
    real(C_DOUBLE) :: dx, dy, dz

    nx = g%nx
    ny = g%ny
    nz = g%nz
    dx = g%dx
    dy = g%dy
    dz = g%dz

#ifdef USE_OPENMP_OFFLOAD
    !$omp target teams distribute parallel do collapse(3) &
    !$omp& map(to: dt_gamma, f%us(0:nx+1,0:ny+1,0:nz+1), &
    !$omp& f%vs(0:nx+1,1:ny+1,0:nz+1), &
    !$omp& f%ws(0:nx+1,0:ny+1,0:nz+1)) &
    !$omp& map(tofrom: ps%rhs(1:nx,1:nz,1:ny)) &
    !$omp& private(i,j,k)
#endif
    do i = 1, nx
        do j = 1, ny
            do k = 1, nz
                ps%rhs(i,k,j) = ( &
                    (f%us(i+1,j,k)-f%us(i,j,k))/dx &
                  + (f%vs(i,j+1,k)-f%vs(i,j,k))/dy &
                  + (f%ws(i,j,k+1)-f%ws(i,j,k))/dz ) / dt_gamma
            end do
        end do
    end do
#ifdef USE_OPENMP_OFFLOAD
    !$omp end target teams distribute parallel do
#endif
end subroutine build_fft_rhs

subroutine apply_fft_pressure_correction(ps, f, g, dt_gamma, ibm)
    type(pressure_solver_type), intent(inout) :: ps
    type(field_type), intent(inout) :: f
    type(grid_type), intent(in) :: g
    real(C_DOUBLE), intent(in) :: dt_gamma
    type(ibm_type), intent(in) :: ibm

    integer :: i,j,k,im,km
    integer :: nx, ny, nz
    real(C_DOUBLE) :: dx, dy, dz

    nx = g%nx
    ny = g%ny
    nz = g%nz
    dx = g%dx
    dy = g%dy
    dz = g%dz

#ifdef USE_OPENMP_OFFLOAD
    !$omp target teams distribute parallel do collapse(3) &
    !$omp& map(to: dt_gamma, ps%rhs(1:nx,1:nz,1:ny), &
    !$omp& f%us(0:nx+1,0:ny+1,0:nz+1), &
    !$omp& f%vs(0:nx+1,1:ny+1,0:nz+1), &
    !$omp& f%ws(0:nx+1,0:ny+1,0:nz+1), &
    !$omp& ibm%coef_u(0:nx+1,0:ny+1,0:nz+1), &
    !$omp& ibm%coef_v(0:nx+1,1:ny+1,0:nz+1), &
    !$omp& ibm%coef_w(0:nx+1,0:ny+1,0:nz+1)) &
    !$omp& map(tofrom: f%un(0:nx+1,0:ny+1,0:nz+1), &
    !$omp& f%vn(0:nx+1,1:ny+1,0:nz+1), &
    !$omp& f%wn(0:nx+1,0:ny+1,0:nz+1), &
    !$omp& f%pn(0:nx+1,1:ny,0:nz+1)) &
    !$omp& private(i,j,k,im,km)
#endif
    do i = 1, nx
        do j = 1, ny
            do k = 1, nz
                im = i - 1
                km = k - 1
                if (im < 1) im = nx
                if (km < 1) km = nz

                f%un(i,j,k) = (f%us(i,j,k) - dt_gamma*(ps%rhs(i,k,j)-ps%rhs(im,k,j))/dx) &
                    / (1.0d0 + dt_gamma*ibm%coef_u(i,j,k))
                f%wn(i,j,k) = (f%ws(i,j,k) - dt_gamma*(ps%rhs(i,k,j)-ps%rhs(i,km,j))/dz) &
                    / (1.0d0 + dt_gamma*ibm%coef_w(i,j,k))
                f%pn(i,j,k) = f%pn(i,j,k) + ps%rhs(i,k,j)
                if (j >= 2) then
                    f%vn(i,j,k) = (f%vs(i,j,k) - dt_gamma*(ps%rhs(i,k,j)-ps%rhs(i,k,j-1))/dy) &
                        / (1.0d0 + dt_gamma*ibm%coef_v(i,j,k))
                end if
            end do
        end do
    end do
#ifdef USE_OPENMP_OFFLOAD
    !$omp end target teams distribute parallel do
#endif
end subroutine apply_fft_pressure_correction
#endif

end module pressure_solver
