module poisson_fftw_backend
    use, intrinsic :: iso_c_binding
    use :: init, only: grid_type, field_type
    use :: poisson_workspace, only: poisson_fft_workspace
    use :: poisson_common, only: allocate_common_workspace, solve_tridiagonal_y
    implicit none

    include "fftw3.f03"

contains

subroutine allocate_fftw_workspace(ws, g)
    type(poisson_fft_workspace), intent(inout) :: ws
    type(grid_type), intent(in) :: g
    integer :: nxh

    nxh = g%nx/2 + 1

    call allocate_common_workspace(ws, nxh, g%ny, g%nz)
    allocate(ws%rhs(g%nx, g%ny, g%nz))
    allocate(ws%plane(g%nx, g%nz))
    allocate(ws%plane_hat(nxh, g%nz))
end subroutine allocate_fftw_workspace

subroutine create_fftw_plans(ws, g)
    type(poisson_fft_workspace), intent(inout) :: ws
    type(grid_type), intent(in) :: g

    ws%plan_fwd = fftw_plan_dft_r2c_2d( &
        g%nz, g%nx, ws%plane, ws%plane_hat, FFTW_ESTIMATE)

    ws%plan_bwd = fftw_plan_dft_c2r_2d( &
        g%nz, g%nx, ws%plane_hat, ws%plane, FFTW_ESTIMATE)
end subroutine create_fftw_plans

subroutine poisson_fftw(g, f, ws)
    type(grid_type), intent(in) :: g
    type(field_type), intent(inout) :: f
    type(poisson_fft_workspace), intent(inout) :: ws

    integer :: nxh
    real(C_DOUBLE) :: dyi2, scale

    nxh = g%nx/2 + 1
    dyi2 = 1.0d0/g%dy**2
    scale = 1.0d0/real(g%nx*g%nz, C_DOUBLE)

    call execute_fftw_forward_y_planes(g, ws)
    call solve_tridiagonal_y(ws, nxh, g%ny, g%nz, dyi2)
    call execute_fftw_backward_y_planes(g, f, ws, scale)
end subroutine poisson_fftw

subroutine execute_fftw_forward_y_planes(g, ws)
    type(grid_type), intent(in) :: g
    type(poisson_fft_workspace), intent(inout) :: ws
    integer :: j

    do j = 1, g%ny
        ws%plane = ws%rhs(:,j,:)
        call fftw_execute_dft_r2c(ws%plan_fwd, ws%plane, ws%plane_hat)
        ws%p_hat(:,j,:) = ws%plane_hat
    end do
end subroutine execute_fftw_forward_y_planes

subroutine execute_fftw_backward_y_planes(g, f, ws, scale)
    type(grid_type), intent(in) :: g
    type(field_type), intent(inout) :: f
    type(poisson_fft_workspace), intent(inout) :: ws
    real(C_DOUBLE), intent(in) :: scale
    integer :: j

    do j = 1, g%ny
        ws%plane_hat = ws%p_hat(:,j,:)
        call fftw_execute_dft_c2r(ws%plan_bwd, ws%plane_hat, ws%plane)
        f%pc(1:g%nx,j,1:g%nz) = ws%plane * scale
    end do
end subroutine execute_fftw_backward_y_planes

subroutine destroy_fftw_workspace(ws)
    type(poisson_fft_workspace), intent(inout) :: ws

    if (c_associated(ws%plan_fwd)) call fftw_destroy_plan(ws%plan_fwd)
    if (c_associated(ws%plan_bwd)) call fftw_destroy_plan(ws%plan_bwd)

    ws%plan_fwd = C_NULL_PTR
    ws%plan_bwd = C_NULL_PTR
end subroutine destroy_fftw_workspace

end module poisson_fftw_backend
