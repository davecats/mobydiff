module poisson_common
    use, intrinsic :: iso_c_binding
    use :: init, only: grid_type
    use :: poisson_workspace, only: poisson_fft_workspace
    implicit none

contains

subroutine allocate_common_workspace(ws, nxh, ny, nz)
    type(poisson_fft_workspace), intent(inout) :: ws
    integer, intent(in) :: nxh, ny, nz

    ! GPU cuFFT solves directly in plane_hat layout, so p_hat is CPU-only.
#ifndef USE_CUFFT
    allocate(ws%p_hat(nxh, ny, nz))
#endif
    allocate(ws%cp_hat(nxh, ny, nz))
    allocate(ws%dp_hat(nxh, ny, nz))
    allocate(ws%den_inv_hat(nxh, ny, nz))
end subroutine allocate_common_workspace

subroutine init_tridiag_coefficients(ws, g)
    type(poisson_fft_workspace), intent(inout) :: ws
    type(grid_type), intent(in) :: g

    integer :: ikx, ikz, jj, kx, kz
    integer :: nx, ny, nz, nxh
    real(C_DOUBLE) :: dx, dy, dz, dyi2
    real(C_DOUBLE) :: kx_ph, kz_ph, k_tot
    real(C_DOUBLE), parameter :: pi = 3.141592653589793d0
    complex(C_DOUBLE_COMPLEX) :: a_val, b_val, c_val, denom
    complex(C_DOUBLE_COMPLEX) :: zero_c, one_c
    logical :: zero_mode

    nx = g%nx
    ny = g%ny
    nz = g%nz
    nxh = nx/2 + 1
    dx = g%dx
    dy = g%dy
    dz = g%dz
    dyi2 = 1.0d0/dy**2
    zero_c = cmplx(0.0d0, 0.0d0, kind=C_DOUBLE_COMPLEX)
    one_c = cmplx(1.0d0, 0.0d0, kind=C_DOUBLE_COMPLEX)

    ! Precompute the mode-dependent Thomas coefficients for all (kx,kz) systems.
#ifdef USE_CUFFT
    !$omp target teams distribute parallel do collapse(2) &
    !$omp& map(tofrom: ws%cp_hat(1:nxh,1:ny,1:nz), ws%den_inv_hat(1:nxh,1:ny,1:nz)) &
    !$omp& private(ikx,ikz,jj,kx,kz,kx_ph,kz_ph,k_tot,a_val,b_val,c_val,denom,zero_mode)
#endif
    do ikx = 1, nxh
        do ikz = 1, nz
            kx = ikx - 1
            kx_ph = (4.0d0/dx**2) * sin(pi*real(kx,C_DOUBLE)/real(nx,C_DOUBLE))**2

            if (ikz <= nz/2 + 1) then
                kz = ikz - 1
            else
                kz = ikz - nz - 1
            end if

            kz_ph = (4.0d0/dz**2) * sin(pi*real(kz,C_DOUBLE)/real(nz,C_DOUBLE))**2
            k_tot = kx_ph + kz_ph
            zero_mode = (kx == 0 .and. kz == 0)

            ! Fix the null pressure mode by forcing the first RHS entry to zero later.
            if (zero_mode) then
                ws%den_inv_hat(ikx,1,ikz) = one_c
                ws%cp_hat(ikx,1,ikz) = zero_c
            else
                b_val = cmplx(-dyi2 - k_tot, 0.0d0, kind=C_DOUBLE_COMPLEX)
                c_val = cmplx( dyi2,         0.0d0, kind=C_DOUBLE_COMPLEX)
                ws%den_inv_hat(ikx,1,ikz) = one_c/b_val
                ws%cp_hat(ikx,1,ikz) = c_val*ws%den_inv_hat(ikx,1,ikz)
            end if

            a_val = cmplx(dyi2, 0.0d0, kind=C_DOUBLE_COMPLEX)
            do jj = 2, ny
                if (jj < ny) then
                    b_val = cmplx(-2.0d0*dyi2 - k_tot, 0.0d0, kind=C_DOUBLE_COMPLEX)
                    c_val = cmplx( dyi2,             0.0d0, kind=C_DOUBLE_COMPLEX)
                else
                    b_val = cmplx(-dyi2 - k_tot, 0.0d0, kind=C_DOUBLE_COMPLEX)
                    c_val = zero_c
                end if

                denom = b_val - a_val*ws%cp_hat(ikx,jj-1,ikz)
                ws%den_inv_hat(ikx,jj,ikz) = one_c/denom

                if (jj < ny) then
                    ws%cp_hat(ikx,jj,ikz) = c_val*ws%den_inv_hat(ikx,jj,ikz)
                else
                    ws%cp_hat(ikx,jj,ikz) = zero_c
                end if
            end do
        end do
    end do
#ifdef USE_CUFFT
    !$omp end target teams distribute parallel do
#endif
end subroutine init_tridiag_coefficients

subroutine solve_tridiagonal_y(ws, nxh, ny, nz, dyi2)
    type(poisson_fft_workspace), intent(inout) :: ws
    integer, intent(in) :: nxh, ny, nz
    real(C_DOUBLE), intent(in) :: dyi2

#ifdef USE_CUFFT
    call solve_tridiagonal_y_cufft(ws, nxh, ny, nz, dyi2)
#else
    call solve_tridiagonal_y_fftw(ws, nxh, ny, nz, dyi2)
#endif
end subroutine solve_tridiagonal_y

#ifdef USE_CUFFT
subroutine solve_tridiagonal_y_cufft(ws, nxh, ny, nz, dyi2)
    type(poisson_fft_workspace), intent(inout) :: ws
    integer, intent(in) :: nxh, ny, nz
    real(C_DOUBLE), intent(in) :: dyi2

    integer :: ikx, ikz, jj
    complex(C_DOUBLE_COMPLEX) :: a_val, rhs_val, zero_c

    zero_c = cmplx(0.0d0, 0.0d0, kind=C_DOUBLE_COMPLEX)

    ! Each Fourier mode is an independent tridiagonal solve in y.
    !$omp target teams distribute parallel do collapse(2) &
    !$omp& map(to: ws%cp_hat(1:nxh,1:ny,1:nz), ws%den_inv_hat(1:nxh,1:ny,1:nz)) &
    !$omp& map(tofrom: ws%plane_hat(1:nxh,1:nz,1:ny), ws%dp_hat(1:nxh,1:ny,1:nz)) &
    !$omp& private(ikx,ikz,jj,a_val,rhs_val)
    do ikx = 1, nxh
        do ikz = 1, nz
            a_val = cmplx(dyi2, 0.0d0, kind=C_DOUBLE_COMPLEX)

            ! The zero mode is singular; choose zero mean pressure correction.
            if (ikx == 1 .and. ikz == 1) then
                ws%dp_hat(ikx,1,ikz) = zero_c
            else
                rhs_val = ws%plane_hat(ikx,ikz,1)
                ws%dp_hat(ikx,1,ikz) = rhs_val*ws%den_inv_hat(ikx,1,ikz)
            end if

            do jj = 2, ny
                rhs_val = ws%plane_hat(ikx,ikz,jj)
                ws%dp_hat(ikx,jj,ikz) = (rhs_val &
                    - a_val*ws%dp_hat(ikx,jj-1,ikz))*ws%den_inv_hat(ikx,jj,ikz)
            end do

            ! Store the back substitution in the layout expected by the active FFT backend.
            ws%plane_hat(ikx,ikz,ny) = ws%dp_hat(ikx,ny,ikz)
            do jj = ny-1, 1, -1
                ws%plane_hat(ikx,ikz,jj) = ws%dp_hat(ikx,jj,ikz) &
                    - ws%cp_hat(ikx,jj,ikz)*ws%plane_hat(ikx,ikz,jj+1)
            end do
        end do
    end do
    !$omp end target teams distribute parallel do
end subroutine solve_tridiagonal_y_cufft
#else
subroutine solve_tridiagonal_y_fftw(ws, nxh, ny, nz, dyi2)
    type(poisson_fft_workspace), intent(inout) :: ws
    integer, intent(in) :: nxh, ny, nz
    real(C_DOUBLE), intent(in) :: dyi2

    integer :: ikx, ikz, jj
    complex(C_DOUBLE_COMPLEX) :: a_val, rhs_val, zero_c

    zero_c = cmplx(0.0d0, 0.0d0, kind=C_DOUBLE_COMPLEX)

    ! Each Fourier mode is an independent tridiagonal solve in y.
    do ikx = 1, nxh
        do ikz = 1, nz
            a_val = cmplx(dyi2, 0.0d0, kind=C_DOUBLE_COMPLEX)

            ! The zero mode is singular; choose zero mean pressure correction.
            if (ikx == 1 .and. ikz == 1) then
                ws%dp_hat(ikx,1,ikz) = zero_c
            else
                rhs_val = ws%p_hat(ikx,1,ikz)
                ws%dp_hat(ikx,1,ikz) = rhs_val*ws%den_inv_hat(ikx,1,ikz)
            end if

            do jj = 2, ny
                rhs_val = ws%p_hat(ikx,jj,ikz)
                ws%dp_hat(ikx,jj,ikz) = (rhs_val &
                    - a_val*ws%dp_hat(ikx,jj-1,ikz))*ws%den_inv_hat(ikx,jj,ikz)
            end do

            ! Store the back substitution in the layout expected by the active FFT backend.
            ws%p_hat(ikx,ny,ikz) = ws%dp_hat(ikx,ny,ikz)
            do jj = ny-1, 1, -1
                ws%p_hat(ikx,jj,ikz) = ws%dp_hat(ikx,jj,ikz) &
                    - ws%cp_hat(ikx,jj,ikz)*ws%p_hat(ikx,jj+1,ikz)
            end do
        end do
    end do
end subroutine solve_tridiagonal_y_fftw
#endif

subroutine deallocate_workspace_arrays(ws)
    type(poisson_fft_workspace), intent(inout) :: ws

    if (allocated(ws%den_inv_hat)) deallocate(ws%den_inv_hat)
    if (allocated(ws%dp_hat)) deallocate(ws%dp_hat)
    if (allocated(ws%cp_hat)) deallocate(ws%cp_hat)
#ifndef USE_CUFFT
    if (allocated(ws%rhs)) deallocate(ws%rhs)
#endif
    if (allocated(ws%p_hat)) deallocate(ws%p_hat)
    if (allocated(ws%plane_hat)) deallocate(ws%plane_hat)
    if (allocated(ws%plane)) deallocate(ws%plane)
end subroutine deallocate_workspace_arrays

end module poisson_common
