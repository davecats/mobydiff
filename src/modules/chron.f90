module chron
    use, intrinsic :: iso_c_binding, only: C_DOUBLE, C_INT
    use, intrinsic :: iso_fortran_env, only: int64
    implicit none

    private

    type, public :: chron_type
        integer(int64) :: start_count = 0_int64
        integer(int64) :: end_count = 0_int64
        integer(int64) :: count_rate = 0_int64
        integer(C_INT) :: nsteps = 0_C_INT
        real(C_DOUBLE) :: elapsed_seconds = 0.0d0
        real(C_DOUBLE) :: seconds_per_step = 0.0d0
        logical :: running = .false.
    end type chron_type

    public :: init_chron, start_chron, stop_chron, write_chron

contains

    subroutine init_chron(timer)
        type(chron_type), intent(inout) :: timer

        timer%start_count = 0_int64
        timer%end_count = 0_int64
        timer%nsteps = 0_C_INT
        timer%elapsed_seconds = 0.0d0
        timer%seconds_per_step = 0.0d0
        timer%running = .false.
        call system_clock(count_rate=timer%count_rate)
    end subroutine init_chron

    subroutine start_chron(timer)
        type(chron_type), intent(inout) :: timer

        call init_chron(timer)
        call system_clock(count=timer%start_count)
        timer%running = .true.
    end subroutine start_chron

    subroutine stop_chron(timer, nsteps)
        type(chron_type), intent(inout) :: timer
        integer(C_INT), intent(in) :: nsteps

        if (.not. timer%running) error stop "chron timer stopped before it was started"

        call system_clock(count=timer%end_count)
        timer%running = .false.
        timer%nsteps = nsteps
        timer%elapsed_seconds = real(timer%end_count - timer%start_count, C_DOUBLE) / &
            real(timer%count_rate, C_DOUBLE)

        if (nsteps > 0_C_INT) then
            timer%seconds_per_step = timer%elapsed_seconds / real(nsteps, C_DOUBLE)
        else
            timer%seconds_per_step = 0.0d0
        end if
    end subroutine stop_chron

    subroutine write_chron(timer)
        type(chron_type), intent(in) :: timer

        write(*,'(A,1X,I0,1X,A,1X,ES16.8,1X,A,1X,ES16.8)') &
            "timing: nsteps", timer%nsteps, "loop_seconds", timer%elapsed_seconds, &
            "seconds_per_step", timer%seconds_per_step
    end subroutine write_chron

end module chron
