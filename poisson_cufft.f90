module poisson_cufft_backend
    use, intrinsic :: iso_c_binding
    use :: init, only: grid_type, field_type
    use :: poisson_workspace, only: poisson_fft_workspace
    use :: poisson_common, only: allocate_common_workspace, solve_tridiagonal_y
    use :: cufft_bindings
    implicit none

contains

subroutine allocate_cufft_workspace(ws, g)
    type(poisson_fft_workspace), intent(inout) :: ws
    type(grid_type), intent(in) :: g
    integer :: nxh

    nxh = g%nx/2 + 1

    call allocate_common_workspace(ws, nxh, g%ny, g%nz)
    allocate(ws%plane(g%nx, g%nz, g%ny))
    allocate(ws%plane_hat(nxh, g%nz, g%ny))

    !$omp target enter data map(to: ws)
    !$omp target enter data map(alloc: &
    !$omp& ws%cp_hat(1:nxh,1:g%ny,1:g%nz), &
    !$omp& ws%dp_hat(1:nxh,1:g%ny,1:g%nz), &
    !$omp& ws%den_inv_hat(1:nxh,1:g%ny,1:g%nz), &
    !$omp& ws%plane(1:g%nx,1:g%nz,1:g%ny), &
    !$omp& ws%plane_hat(1:nxh,1:g%nz,1:g%ny))
end subroutine allocate_cufft_workspace

subroutine create_cufft_plans(ws, g)
    type(poisson_fft_workspace), intent(inout) :: ws
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

subroutine poisson_cufft(g, f, ws)
    type(grid_type), intent(in) :: g
    type(field_type), intent(inout) :: f
    type(poisson_fft_workspace), intent(inout) :: ws

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
    call unpack_pressure_from_cufft(f, ws, nx, ny, nz, scale)
end subroutine poisson_cufft

subroutine execute_cufft_forward(ws)
    type(poisson_fft_workspace), intent(inout) :: ws
    integer(C_INT) :: ierr
    type(C_PTR) :: plane_ptr, plane_hat_ptr

    !$omp target data use_device_addr(ws%plane, ws%plane_hat)
    plane_ptr = c_loc(ws%plane(1,1,1))
    plane_hat_ptr = c_loc(ws%plane_hat(1,1,1))
    ierr = cufftExecD2Z(ws%plan_fwd, plane_ptr, plane_hat_ptr)
    !$omp end target data
    call check_cufft(ierr, "batched cufftExecD2Z")
    call check_cuda(cudaStreamSynchronize(C_NULL_PTR), "cudaStreamSynchronize after batched cufftExecD2Z")
end subroutine execute_cufft_forward

subroutine execute_cufft_backward(ws)
    type(poisson_fft_workspace), intent(inout) :: ws
    integer(C_INT) :: ierr
    type(C_PTR) :: plane_ptr, plane_hat_ptr

    !$omp target data use_device_addr(ws%plane_hat, ws%plane)
    plane_hat_ptr = c_loc(ws%plane_hat(1,1,1))
    plane_ptr = c_loc(ws%plane(1,1,1))
    ierr = cufftExecZ2D(ws%plan_bwd, plane_hat_ptr, plane_ptr)
    !$omp end target data
    call check_cufft(ierr, "batched cufftExecZ2D")
    call check_cuda(cudaStreamSynchronize(C_NULL_PTR), "cudaStreamSynchronize after batched cufftExecZ2D")
end subroutine execute_cufft_backward

subroutine unpack_pressure_from_cufft(f, ws, nx, ny, nz, scale)
    type(field_type), intent(inout) :: f
    type(poisson_fft_workspace), intent(inout) :: ws
    integer, intent(in) :: nx, ny, nz
    real(C_DOUBLE), intent(in) :: scale
    integer :: i, j, ikz

    !$omp target teams distribute parallel do collapse(3) &
    !$omp& map(to: ws%plane(1:nx,1:nz,1:ny)) &
    !$omp& map(tofrom: f%pc(0:nx+1,1:ny,0:nz+1)) private(i,j,ikz)
    do j = 1, ny
        do ikz = 1, nz
            do i = 1, nx
                f%pc(i,j,ikz) = ws%plane(i,ikz,j) * scale
            end do
        end do
    end do
    !$omp end target teams distribute parallel do
end subroutine unpack_pressure_from_cufft

subroutine destroy_cufft_workspace(ws)
    type(poisson_fft_workspace), intent(inout) :: ws

    if (ws%plan_fwd /= 0) call check_cufft(cufftDestroy(ws%plan_fwd), "cufftDestroy forward")
    if (ws%plan_bwd /= 0) call check_cufft(cufftDestroy(ws%plan_bwd), "cufftDestroy backward")

    if (allocated(ws%cp_hat)) then
        !$omp target exit data map(delete: &
        !$omp& ws%cp_hat, ws%dp_hat, ws%den_inv_hat, &
        !$omp& ws%plane, ws%plane_hat)
    end if
    !$omp target exit data map(delete: ws)

    ws%plan_fwd = 0
    ws%plan_bwd = 0
end subroutine destroy_cufft_workspace

end module poisson_cufft_backend
