module poisson_workspace
    use, intrinsic :: iso_c_binding
    implicit none

    type :: poisson_fft_workspace
        complex(C_DOUBLE_COMPLEX), allocatable :: p_hat(:,:,:)
        complex(C_DOUBLE_COMPLEX), allocatable :: cp_hat(:,:,:)
        complex(C_DOUBLE_COMPLEX), allocatable :: dp_hat(:,:,:)
        complex(C_DOUBLE_COMPLEX), allocatable :: den_inv_hat(:,:,:)

#ifdef USE_CUFFT
        real(C_DOUBLE), allocatable :: plane(:,:,:)
        complex(C_DOUBLE_COMPLEX), allocatable :: plane_hat(:,:,:)

        integer(C_INT) :: plan_fwd = 0
        integer(C_INT) :: plan_bwd = 0
#else
        real(C_DOUBLE), allocatable :: rhs(:,:,:)
        real(C_DOUBLE), allocatable :: plane(:,:)
        complex(C_DOUBLE_COMPLEX), allocatable :: plane_hat(:,:)

        type(C_PTR) :: plan_fwd = C_NULL_PTR
        type(C_PTR) :: plan_bwd = C_NULL_PTR
#endif
    end type poisson_fft_workspace
end module poisson_workspace
