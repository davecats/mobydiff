module gpu_runtime
    use, intrinsic :: iso_c_binding
    use :: init, only: dns_type, field_type, grid_type
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
        !$omp target enter data map(to: &
        !$omp& g%xNode, g%yNode, g%zNode, &
        !$omp& g%x, g%y, g%z, g%d1x, g%d1y, g%d1z, &
        !$omp& g%lapXm, g%lapX0, g%lapXp, &
        !$omp& g%lapYm, g%lapY0, g%lapYp, &
        !$omp& g%lapZm, g%lapZ0, g%lapZp)
#endif
    end subroutine enter_grid_data

    subroutine exit_grid_data(g, dns)
        type(grid_type), intent(inout) :: g
        type(dns_type), intent(in)     :: dns

#ifdef USE_OPENMP_OFFLOAD
        !$omp target exit data map(delete: &
        !$omp& g%xNode, g%yNode, g%zNode, &
        !$omp& g%x, g%y, g%z, g%d1x, g%d1y, g%d1z, &
        !$omp& g%lapXm, g%lapX0, g%lapXp, &
        !$omp& g%lapYm, g%lapY0, g%lapYp, &
        !$omp& g%lapZm, g%lapZ0, g%lapZp)
        !$omp target exit data map(delete: g)
#endif
    end subroutine exit_grid_data

    subroutine enter_field_data(f, dns)
        type(field_type), intent(inout) :: f
        type(dns_type),   intent(in)    :: dns

#ifdef USE_OPENMP_OFFLOAD
        !$omp target enter data map(to: f)
        !$omp target enter data map(to: &
        !$omp& f%q, f%qs, f%oldrhs)
#endif
    end subroutine enter_field_data

    subroutine exit_field_data(f, dns)
        type(field_type), intent(inout) :: f
        type(dns_type),   intent(in)    :: dns

#ifdef USE_OPENMP_OFFLOAD
        !$omp target exit data map(delete: &
        !$omp& f%q, f%qs, f%oldrhs)
        !$omp target exit data map(delete: f)
#endif
    end subroutine exit_field_data

end module gpu_runtime
