module gpu_runtime
    use, intrinsic :: iso_c_binding
    use :: init, only: dns_type, grid_type
#ifdef USE_OPENMP_OFFLOAD
    use omp_lib
#endif
    implicit none

contains

    subroutine init_openmp_offload(has_terminal)
        logical, intent(in), optional :: has_terminal
        logical :: terminal

        terminal = .true.
        if (present(has_terminal)) terminal = has_terminal
#ifdef USE_OPENMP_OFFLOAD
        if (terminal) then
            if (omp_get_num_devices() <= 0) then
                print *, "WARNING: OpenMP offload enabled, but no target device was reported."
                print *, "         The OpenMP runtime may execute target regions on the host."
            else
                print *, "OpenMP target devices available:", omp_get_num_devices()
            end if
        end if
#endif
    end subroutine init_openmp_offload

    subroutine enter_grid_data(g, dns)
        type(grid_type), intent(inout) :: g
        type(dns_type), intent(in)     :: dns

#ifdef USE_OPENMP_OFFLOAD
        !$omp target enter data map(to: g)
        !$omp target enter data map(to: g%xNode, g%yNode, g%zNode)
#endif
    end subroutine enter_grid_data

    subroutine exit_grid_data(g, dns)
        type(grid_type), intent(inout) :: g
        type(dns_type), intent(in)     :: dns

#ifdef USE_OPENMP_OFFLOAD
        !$omp target exit data map(delete: g%xNode, g%yNode, g%zNode)
        !$omp target exit data map(delete: g)
#endif
    end subroutine exit_grid_data

end module gpu_runtime
