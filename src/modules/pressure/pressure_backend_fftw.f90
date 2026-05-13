module pressure_backend_fftw
    use, intrinsic :: iso_c_binding
    use :: init, only: grid_type
    use :: pressure_workspace, only: pressure_solver_type
    use :: pressure_fft, only: allocate_common_workspace, solve_tridiagonal_y
    implicit none

    include "fftw3.f03"

contains

subroutine allocate_fftw_workspace(ws, g)
    type(pressure_solver_type), intent(inout) :: ws
    type(grid_type), intent(in) :: g
    integer :: nxh

    nxh = g%nx/2 + 1

    call allocate_common_workspace(ws, g%nx, nxh, g%ny, g%nz)
end subroutine allocate_fftw_workspace

subroutine create_fftw_plans(ws, g)
    type(pressure_solver_type), intent(inout) :: ws
    type(grid_type), intent(in) :: g

    integer(C_INT) :: rank, istride, ostride, idist, odist, batch
    integer(C_INT), target :: fft_dims(2), inembed(2), onembed(2)
    integer :: nxh

    nxh = g%nx/2 + 1

    fft_dims = [int(g%nz, C_INT), int(g%nx, C_INT)]
    inembed = fft_dims
    onembed = [int(g%nz, C_INT), int(nxh, C_INT)]
    rank = 2_C_INT
    istride = 1_C_INT
    ostride = 1_C_INT
    idist = int(g%nx*g%nz, C_INT)
    odist = int(nxh*g%nz, C_INT)
    batch = int(g%ny, C_INT)

    ws%plan_fwd = fftw_plan_many_dft_r2c(rank, fft_dims, batch, &
        ws%rhs, inembed, istride, idist, ws%plane_hat, onembed, &
        ostride, odist, FFTW_ESTIMATE)

    ws%plan_bwd = fftw_plan_many_dft_c2r(rank, fft_dims, batch, &
        ws%plane_hat, onembed, ostride, odist, ws%rhs, inembed, &
        istride, idist, FFTW_ESTIMATE)
end subroutine create_fftw_plans

subroutine solve_pressure_fftw(g, ws)
    type(grid_type), intent(in) :: g
    type(pressure_solver_type), intent(inout) :: ws

    integer :: nxh
    real(C_DOUBLE) :: dyi2, scale

    nxh = g%nx/2 + 1
    dyi2 = 1.0d0/g%dy**2
    scale = 1.0d0/real(g%nx*g%nz, C_DOUBLE)

    call execute_fftw_forward(ws)
    call solve_tridiagonal_y(ws, nxh, g%ny, g%nz, dyi2)
    call execute_fftw_backward(ws)
    call scale_pressure_correction_fftw(g, ws, scale)
end subroutine solve_pressure_fftw

subroutine execute_fftw_forward(ws)
    type(pressure_solver_type), intent(inout) :: ws

    call fftw_execute_dft_r2c(ws%plan_fwd, ws%rhs, ws%plane_hat)
end subroutine execute_fftw_forward

subroutine execute_fftw_backward(ws)
    type(pressure_solver_type), intent(inout) :: ws

    call fftw_execute_dft_c2r(ws%plan_bwd, ws%plane_hat, ws%rhs)
end subroutine execute_fftw_backward

subroutine scale_pressure_correction_fftw(g, ws, scale)
    type(grid_type), intent(in) :: g
    type(pressure_solver_type), intent(inout) :: ws
    real(C_DOUBLE), intent(in) :: scale
    integer :: i, j, k

    do i = 1, g%nx
        do j = 1, g%ny
            do k = 1, g%nz
                ws%rhs(i,k,j) = ws%rhs(i,k,j) * scale
            end do
        end do
    end do
end subroutine scale_pressure_correction_fftw

subroutine destroy_fftw_workspace(ws)
    type(pressure_solver_type), intent(inout) :: ws

    if (c_associated(ws%plan_fwd)) call fftw_destroy_plan(ws%plan_fwd)
    if (c_associated(ws%plan_bwd)) call fftw_destroy_plan(ws%plan_bwd)

    ws%plan_fwd = C_NULL_PTR
    ws%plan_bwd = C_NULL_PTR
end subroutine destroy_fftw_workspace

end module pressure_backend_fftw
