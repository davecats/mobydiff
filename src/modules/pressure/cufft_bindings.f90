module cufft_bindings
    use, intrinsic :: iso_c_binding
    implicit none

    integer(C_INT), parameter :: CUFFT_SUCCESS = 0
    integer(C_INT), parameter :: CUFFT_D2Z = 106
    integer(C_INT), parameter :: CUFFT_Z2D = 108

    interface
        function cufftPlanMany(plan, rank, n, inembed, istride, idist, &
                               onembed, ostride, odist, fft_type, batch) &
                               bind(C, name="cufftPlanMany") result(ierr)
            import :: C_INT, C_PTR
            integer(C_INT) :: plan
            integer(C_INT), value :: rank
            type(C_PTR), value :: n
            type(C_PTR), value :: inembed
            integer(C_INT), value :: istride, idist
            type(C_PTR), value :: onembed
            integer(C_INT), value :: ostride, odist
            integer(C_INT), value :: fft_type, batch
            integer(C_INT) :: ierr
        end function cufftPlanMany

        function cufftDestroy(plan) bind(C, name="cufftDestroy") result(ierr)
            import :: C_INT
            integer(C_INT), value :: plan
            integer(C_INT) :: ierr
        end function cufftDestroy

        function cufftExecD2Z(plan, idata, odata) bind(C, name="cufftExecD2Z") result(ierr)
            import :: C_INT, C_PTR
            integer(C_INT), value :: plan
            type(C_PTR), value :: idata
            type(C_PTR), value :: odata
            integer(C_INT) :: ierr
        end function cufftExecD2Z

        function cufftExecZ2D(plan, idata, odata) bind(C, name="cufftExecZ2D") result(ierr)
            import :: C_INT, C_PTR
            integer(C_INT), value :: plan
            type(C_PTR), value :: idata
            type(C_PTR), value :: odata
            integer(C_INT) :: ierr
        end function cufftExecZ2D

        function cudaStreamSynchronize(stream) bind(C, name="cudaStreamSynchronize") result(ierr)
            import :: C_INT, C_PTR
            type(C_PTR), value :: stream
            integer(C_INT) :: ierr
        end function cudaStreamSynchronize
    end interface

contains

subroutine check_cufft(ierr, where)
    integer(C_INT), intent(in) :: ierr
    character(len=*), intent(in) :: where

    if (ierr /= CUFFT_SUCCESS) then
        write(*,'(A,1X,A,1X,I0)') "cuFFT error in", trim(where), ierr
        error stop
	end if
end subroutine check_cufft

subroutine check_cuda(ierr, where)
    integer(C_INT), intent(in) :: ierr
    character(len=*), intent(in) :: where

    if (ierr /= 0_C_INT) then
        write(*,'(A,1X,A,1X,I0)') "CUDA error in", trim(where), ierr
        error stop
    end if
end subroutine check_cuda

end module cufft_bindings
