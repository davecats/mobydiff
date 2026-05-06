module poisson_solver
    use :: init, only: grid_type, field_type
    use :: poisson_workspace, only: poisson_fft_workspace
    use :: poisson_common, only: init_tridiag_coefficients, deallocate_workspace_arrays
#ifdef USE_CUFFT
    use :: poisson_cufft_backend, only: &
        allocate_backend_workspace => allocate_cufft_workspace, &
        create_backend_plans => create_cufft_plans, &
        poisson_backend => poisson_cufft, &
        destroy_backend_workspace => destroy_cufft_workspace
#else
    use :: poisson_fftw_backend, only: &
        allocate_backend_workspace => allocate_fftw_workspace, &
        create_backend_plans => create_fftw_plans, &
        poisson_backend => poisson_fftw, &
        destroy_backend_workspace => destroy_fftw_workspace
#endif
    implicit none

contains

subroutine init_poisson_fft_workspace(ws, g)
    type(poisson_fft_workspace), intent(inout) :: ws
    type(grid_type), intent(in) :: g

    call allocate_backend_workspace(ws, g)
    call create_backend_plans(ws, g)
    call init_tridiag_coefficients(ws, g)
end subroutine init_poisson_fft_workspace

subroutine poisson(g, f, ws)
    type(grid_type), intent(in) :: g
    type(field_type), intent(inout) :: f
    type(poisson_fft_workspace), intent(inout) :: ws

    call poisson_backend(g, f, ws)
end subroutine poisson

subroutine destroy_poisson_fft_workspace(ws)
    type(poisson_fft_workspace), intent(inout) :: ws

    call destroy_backend_workspace(ws)
    call deallocate_workspace_arrays(ws)
end subroutine destroy_poisson_fft_workspace

end module poisson_solver
