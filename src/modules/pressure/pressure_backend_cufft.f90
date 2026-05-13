module pressure_backend_cufft
    use, intrinsic :: iso_c_binding
    use :: init, only: grid_type
    use :: pressure_workspace, only: pressure_solver_type
    use :: pressure_fft, only: allocate_common_workspace, solve_tridiagonal_y
    use :: cufft_bindings
    implicit none

contains

subroutine allocate_cufft_workspace(ws, g)
    type(pressure_solver_type), intent(inout) :: ws
    type(grid_type), intent(in) :: g
    integer :: nxh

    nxh = g%nx/2 + 1

    call allocate_common_workspace(ws, g%nx, nxh, g%ny, g%nz)

    !$omp target enter data map(to: ws)
    !$omp target enter data map(alloc: &
    !$omp& ws%cp_hat(1:nxh,1:g%ny,1:g%nz), &
    !$omp& ws%dp_hat(1:nxh,1:g%ny,1:g%nz), &
    !$omp& ws%den_inv_hat(1:nxh,1:g%ny,1:g%nz), &
    !$omp& ws%rhs(1:g%nx,1:g%nz,1:g%ny), &
    !$omp& ws%plane_hat(1:nxh,1:g%nz,1:g%ny))
end subroutine allocate_cufft_workspace

subroutine create_cufft_plans(ws, g)
    type(pressure_solver_type), intent(inout) :: ws
    type(grid_type), intent(in) :: g

    integer(C_INT), target :: fft_dims(2), inembed(2), onembed(2)
    integer(C_INT) :: rank, istride, ostride, idist, odist, batch
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

    call check_cufft(cufftPlanMany(ws%plan_fwd, rank, c_loc(fft_dims(1)), &
        c_loc(inembed(1)), istride, idist, c_loc(onembed(1)), &
        ostride, odist, CUFFT_D2Z, batch), "cufftPlanMany D2Z")

    call check_cufft(cufftPlanMany(ws%plan_bwd, rank, c_loc(fft_dims(1)), &
        c_loc(onembed(1)), ostride, odist, c_loc(inembed(1)), &
        istride, idist, CUFFT_Z2D, batch), "cufftPlanMany Z2D")
end subroutine create_cufft_plans

subroutine solve_pressure_cufft(g, ws)
    type(grid_type), intent(in) :: g
    type(pressure_solver_type), intent(inout) :: ws

    integer :: nx, ny, nz, nxh
    real(C_DOUBLE) :: dyi2, scale

    nx = g%nx
    ny = g%ny
    nz = g%nz
    nxh = nx/2 + 1
    dyi2 = 1.0d0/g%dy**2
    scale = 1.0d0/real(nx*nz, C_DOUBLE)

    call execute_cufft_forward(ws)
    call solve_tridiagonal_y(ws, nxh, ny, nz, dyi2)
    call execute_cufft_backward(ws)
    call scale_pressure_correction_cufft(ws, nx, ny, nz, scale)
end subroutine solve_pressure_cufft

subroutine execute_cufft_forward(ws)
    type(pressure_solver_type), intent(inout) :: ws
    integer(C_INT) :: ierr
    type(C_PTR) :: rhs_ptr, plane_hat_ptr

    !$omp target data use_device_addr(ws%rhs, ws%plane_hat)
    rhs_ptr = c_loc(ws%rhs(1,1,1))
    plane_hat_ptr = c_loc(ws%plane_hat(1,1,1))
    ierr = cufftExecD2Z(ws%plan_fwd, rhs_ptr, plane_hat_ptr)
    !$omp end target data
    call check_cufft(ierr, "batched cufftExecD2Z")
    call check_cuda(cudaStreamSynchronize(C_NULL_PTR), "cudaStreamSynchronize after batched cufftExecD2Z")
end subroutine execute_cufft_forward

subroutine execute_cufft_backward(ws)
    type(pressure_solver_type), intent(inout) :: ws
    integer(C_INT) :: ierr
    type(C_PTR) :: rhs_ptr, plane_hat_ptr

    !$omp target data use_device_addr(ws%plane_hat, ws%rhs)
    plane_hat_ptr = c_loc(ws%plane_hat(1,1,1))
    rhs_ptr = c_loc(ws%rhs(1,1,1))
    ierr = cufftExecZ2D(ws%plan_bwd, plane_hat_ptr, rhs_ptr)
    !$omp end target data
    call check_cufft(ierr, "batched cufftExecZ2D")
    call check_cuda(cudaStreamSynchronize(C_NULL_PTR), "cudaStreamSynchronize after batched cufftExecZ2D")
end subroutine execute_cufft_backward

subroutine scale_pressure_correction_cufft(ws, nx, ny, nz, scale)
    type(pressure_solver_type), intent(inout) :: ws
    integer, intent(in) :: nx, ny, nz
    real(C_DOUBLE), intent(in) :: scale
    integer :: i, j, ikz

    !$omp target teams distribute parallel do collapse(3) &
    !$omp& map(tofrom: ws%rhs(1:nx,1:nz,1:ny)) private(i,j,ikz)
    do j = 1, ny
        do ikz = 1, nz
            do i = 1, nx
                ws%rhs(i,ikz,j) = ws%rhs(i,ikz,j) * scale
            end do
        end do
    end do
    !$omp end target teams distribute parallel do
end subroutine scale_pressure_correction_cufft

subroutine destroy_cufft_workspace(ws)
    type(pressure_solver_type), intent(inout) :: ws

    if (ws%plan_fwd /= 0) call check_cufft(cufftDestroy(ws%plan_fwd), "cufftDestroy forward")
    if (ws%plan_bwd /= 0) call check_cufft(cufftDestroy(ws%plan_bwd), "cufftDestroy backward")

    if (allocated(ws%cp_hat)) then
        !$omp target exit data map(delete: &
        !$omp& ws%cp_hat, ws%dp_hat, ws%den_inv_hat, &
        !$omp& ws%rhs, ws%plane_hat)
    end if
    !$omp target exit data map(delete: ws)

    ws%plan_fwd = 0
    ws%plan_bwd = 0
end subroutine destroy_cufft_workspace

end module pressure_backend_cufft
