module gpu_runtime
    use, intrinsic :: iso_c_binding
    use :: init, only: grid_type, field_type
#ifdef USE_OPENMP_OFFLOAD
    use omp_lib
#endif
    implicit none

contains

    subroutine init_openmp_offload()
#ifdef USE_OPENMP_OFFLOAD
        if (omp_get_num_devices() <= 0) then
            print *, "WARNING: OpenMP offload enabled, but no target device was reported."
            print *, "         The OpenMP runtime may execute target regions on the host."
        else
            print *, "OpenMP target devices available:", omp_get_num_devices()
        end if
#endif
    end subroutine init_openmp_offload

    subroutine enter_field_data(f, g)
        type(field_type), intent(inout) :: f
        type(grid_type),  intent(in)    :: g

#ifdef USE_OPENMP_OFFLOAD
        !$omp target enter data map(to: f)
        !$omp target enter data map(to: &
        !$omp& f%un(0:g%nx+1,0:g%ny+1,0:g%nz+1), f%us(0:g%nx+1,0:g%ny+1,0:g%nz+1), f%oldrhsu(1:g%nx,1:g%ny,1:g%nz), &
        !$omp& f%vn(0:g%nx+1,1:g%ny+1,0:g%nz+1), f%vs(0:g%nx+1,1:g%ny+1,0:g%nz+1), f%oldrhsv(1:g%nx,2:g%ny,1:g%nz),&
        !$omp& f%wn(0:g%nx+1,0:g%ny+1,0:g%nz+1), f%ws(0:g%nx+1,0:g%ny+1,0:g%nz+1), f%oldrhsw(1:g%nx,1:g%ny,1:g%nz),&
        !$omp& f%pn(0:g%nx+1,1:g%ny,0:g%nz+1))
#endif
    end subroutine enter_field_data

    subroutine exit_field_data(f, g)
        type(field_type), intent(inout) :: f
        type(grid_type),  intent(in)    :: g

#ifdef USE_OPENMP_OFFLOAD
        !$omp target exit data map(delete: &
        !$omp& f%un(0:g%nx+1,0:g%ny+1,0:g%nz+1), f%us(0:g%nx+1,0:g%ny+1,0:g%nz+1), f%oldrhsu(1:g%nx,1:g%ny,1:g%nz),&
        !$omp& f%vn(0:g%nx+1,1:g%ny+1,0:g%nz+1), f%vs(0:g%nx+1,1:g%ny+1,0:g%nz+1), f%oldrhsv(1:g%nx,2:g%ny,1:g%nz),&
        !$omp& f%wn(0:g%nx+1,0:g%ny+1,0:g%nz+1), f%ws(0:g%nx+1,0:g%ny+1,0:g%nz+1), f%oldrhsw(1:g%nx,1:g%ny,1:g%nz),&
        !$omp& f%pn(0:g%nx+1,1:g%ny,0:g%nz+1))
        !$omp target exit data map(delete: f)
#endif
    end subroutine exit_field_data

end module gpu_runtime
