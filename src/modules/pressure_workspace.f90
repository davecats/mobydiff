module pressure_workspace
    use, intrinsic :: iso_c_binding
    implicit none

    type :: pressure_solver_type
#ifdef USE_REDBLACK
        integer(C_INT) :: nIter=3
        real(C_DOUBLE) :: sor=1.5
#else
        complex(C_DOUBLE_COMPLEX), allocatable :: cp_hat(:,:,:)
        complex(C_DOUBLE_COMPLEX), allocatable :: dp_hat(:,:,:)
        complex(C_DOUBLE_COMPLEX), allocatable :: den_inv_hat(:,:,:)
        real(C_DOUBLE), allocatable :: rhs(:,:,:)
        complex(C_DOUBLE_COMPLEX), allocatable :: plane_hat(:,:,:)

#ifdef USE_CUFFT
        integer(C_INT) :: plan_fwd = 0
        integer(C_INT) :: plan_bwd = 0
#else
        type(C_PTR) :: plan_fwd = C_NULL_PTR
        type(C_PTR) :: plan_bwd = C_NULL_PTR
#endif
#endif
    end type pressure_solver_type
end module pressure_workspace
