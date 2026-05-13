module pressure_fft
    use, intrinsic :: iso_c_binding
    use :: init, only: grid_type
    use :: pressure_workspace, only: pressure_solver_type
    implicit none

contains

subroutine allocate_common_workspace(ps, nx, nxh, ny, nz)
    type(pressure_solver_type), intent(inout) :: ps
    integer, intent(in) :: nx, nxh, ny, nz

    allocate(ps%cp_hat(nxh, ny, nz))
    allocate(ps%dp_hat(nxh, ny, nz))
    allocate(ps%den_inv_hat(nxh, ny, nz))
    allocate(ps%rhs(nx, nz, ny))
    allocate(ps%plane_hat(nxh, nz, ny))
end subroutine allocate_common_workspace

subroutine init_tridiag_coefficients(ps, g)
    type(pressure_solver_type), intent(inout) :: ps
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
    !$omp& map(tofrom: ps%cp_hat(1:nxh,1:ny,1:nz), ps%den_inv_hat(1:nxh,1:ny,1:nz)) &
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
                ps%den_inv_hat(ikx,1,ikz) = one_c
                ps%cp_hat(ikx,1,ikz) = zero_c
            else
                b_val = cmplx(-dyi2 - k_tot, 0.0d0, kind=C_DOUBLE_COMPLEX)
                c_val = cmplx( dyi2,         0.0d0, kind=C_DOUBLE_COMPLEX)
                ps%den_inv_hat(ikx,1,ikz) = one_c/b_val
                ps%cp_hat(ikx,1,ikz) = c_val*ps%den_inv_hat(ikx,1,ikz)
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

                denom = b_val - a_val*ps%cp_hat(ikx,jj-1,ikz)
                ps%den_inv_hat(ikx,jj,ikz) = one_c/denom

                if (jj < ny) then
                    ps%cp_hat(ikx,jj,ikz) = c_val*ps%den_inv_hat(ikx,jj,ikz)
                else
                    ps%cp_hat(ikx,jj,ikz) = zero_c
                end if
            end do
        end do
    end do
#ifdef USE_CUFFT
    !$omp end target teams distribute parallel do
#endif
end subroutine init_tridiag_coefficients

subroutine solve_tridiagonal_y(ps, nxh, ny, nz, dyi2)
    type(pressure_solver_type), intent(inout) :: ps
    integer, intent(in) :: nxh, ny, nz
    real(C_DOUBLE), intent(in) :: dyi2

    integer :: ikx, ikz, jj
    complex(C_DOUBLE_COMPLEX) :: a_val, rhs_val, zero_c

    zero_c = cmplx(0.0d0, 0.0d0, kind=C_DOUBLE_COMPLEX)

    ! Each Fourier mode is an independent tridiagonal solve in y.
#ifdef USE_OPENMP_OFFLOAD
    !$omp target teams distribute parallel do collapse(2) &
    !$omp& map(to: ps%cp_hat(1:nxh,1:ny,1:nz), ps%den_inv_hat(1:nxh,1:ny,1:nz)) &
    !$omp& map(tofrom: ps%plane_hat(1:nxh,1:nz,1:ny), ps%dp_hat(1:nxh,1:ny,1:nz)) &
    !$omp& private(ikx,ikz,jj,a_val,rhs_val)
#endif
    do ikx = 1, nxh
        do ikz = 1, nz
            a_val = cmplx(dyi2, 0.0d0, kind=C_DOUBLE_COMPLEX)

            ! The zero mode is singular; choose zero mean pressure correction.
            if (ikx == 1 .and. ikz == 1) then
                ps%dp_hat(ikx,1,ikz) = zero_c
            else
                rhs_val = ps%plane_hat(ikx,ikz,1)
                ps%dp_hat(ikx,1,ikz) = rhs_val*ps%den_inv_hat(ikx,1,ikz)
            end if

            do jj = 2, ny
                rhs_val = ps%plane_hat(ikx,ikz,jj)
                ps%dp_hat(ikx,jj,ikz) = (rhs_val &
                    - a_val*ps%dp_hat(ikx,jj-1,ikz))*ps%den_inv_hat(ikx,jj,ikz)
            end do

            ! Store the back substitution in the layout expected by the active FFT backend.
            ps%plane_hat(ikx,ikz,ny) = ps%dp_hat(ikx,ny,ikz)
            do jj = ny-1, 1, -1
                ps%plane_hat(ikx,ikz,jj) = ps%dp_hat(ikx,jj,ikz) &
                    - ps%cp_hat(ikx,jj,ikz)*ps%plane_hat(ikx,ikz,jj+1)
            end do
        end do
    end do
#ifdef USE_OPENMP_OFFLOAD
    !$omp end target teams distribute parallel do
#endif
end subroutine solve_tridiagonal_y

subroutine deallocate_workspace_arrays(ps)
    type(pressure_solver_type), intent(inout) :: ps

    if (allocated(ps%den_inv_hat)) deallocate(ps%den_inv_hat)
    if (allocated(ps%dp_hat)) deallocate(ps%dp_hat)
    if (allocated(ps%cp_hat)) deallocate(ps%cp_hat)
    if (allocated(ps%rhs)) deallocate(ps%rhs)
    if (allocated(ps%plane_hat)) deallocate(ps%plane_hat)
end subroutine deallocate_workspace_arrays

end module pressure_fft
