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

    ! Generic multi-bucket phase profiler: accumulate wall time and call counts
    ! into named categories and print one line each plus a total. The category
    ! labels and the output-line tag are supplied by the caller, so any
    ! subsystem (LES today, the projection profiling to come) can reuse it.
    integer, parameter, public :: PROFILER_MAX_CATS = 8

    type, public :: profiler_type
        character(len=24) :: tag = "timing"                 ! output line prefix
        integer :: ncats = 0
        character(len=24) :: labels(PROFILER_MAX_CATS) = ""
        real(C_DOUBLE) :: seconds(PROFILER_MAX_CATS) = 0.0d0
        integer(C_INT) :: calls(PROFILER_MAX_CATS) = 0_C_INT
    end type profiler_type

    public :: init_chron, start_chron, stop_chron, write_chron
    public :: wall_seconds, init_profiler, reset_profiler, profiler_add, write_profiler
    public :: profiler_total

contains

    ! Wall-clock reading in seconds, shared by the stopwatch and the profiler.
    real(C_DOUBLE) function wall_seconds() result(seconds)
        integer(int64) :: count, rate

        call system_clock(count=count, count_rate=rate)
        seconds = real(count, C_DOUBLE)/real(rate, C_DOUBLE)
    end function wall_seconds

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

    ! Set the output tag and the per-category labels; zero the accumulators.
    subroutine init_profiler(profile, tag, labels)
        type(profiler_type), intent(out) :: profile
        character(len=*), intent(in) :: tag
        character(len=*), intent(in) :: labels(:)

        if (size(labels) > PROFILER_MAX_CATS) error stop "profiler: too many categories"
        profile%tag = tag
        profile%ncats = size(labels)
        profile%labels(1:profile%ncats) = labels
        profile%seconds = 0.0d0
        profile%calls = 0_C_INT
    end subroutine init_profiler

    ! Zero the accumulators, keeping the tag and labels.
    subroutine reset_profiler(profile)
        type(profiler_type), intent(inout) :: profile

        profile%seconds = 0.0d0
        profile%calls = 0_C_INT
    end subroutine reset_profiler

    subroutine profiler_add(profile, category, elapsed_seconds)
        type(profiler_type), intent(inout) :: profile
        integer, intent(in) :: category
        real(C_DOUBLE), intent(in) :: elapsed_seconds

        if (category < 1 .or. category > profile%ncats) return
        profile%seconds(category) = profile%seconds(category) + elapsed_seconds
        profile%calls(category) = profile%calls(category) + 1_C_INT
    end subroutine profiler_add

    ! Wall time accumulated over all categories (the caller decides whether that
    ! is a partition of something or a subset of it).
    real(C_DOUBLE) function profiler_total(profile) result(seconds)
        type(profiler_type), intent(in) :: profile

        seconds = sum(profile%seconds(1:profile%ncats))
    end function profiler_total

    subroutine write_profiler(profile, nsteps)
        type(profiler_type), intent(in) :: profile
        integer(C_INT), intent(in) :: nsteps

        integer :: i

        do i = 1, profile%ncats
            call write_profiler_line(profile%tag, profile%labels(i), &
                profile%seconds(i), profile%calls(i), nsteps)
        end do
        call write_profiler_line(profile%tag, "total_measured", &
            sum(profile%seconds(1:profile%ncats)), sum(profile%calls(1:profile%ncats)), nsteps)
    end subroutine write_profiler

    subroutine write_profiler_line(tag, label, seconds, calls, nsteps)
        character(len=*), intent(in) :: tag, label
        real(C_DOUBLE), intent(in) :: seconds
        integer(C_INT), intent(in) :: calls, nsteps

        real(C_DOUBLE) :: seconds_per_step, seconds_per_call

        seconds_per_step = 0.0d0
        seconds_per_call = 0.0d0
        if (nsteps > 0_C_INT) seconds_per_step = seconds/real(nsteps, C_DOUBLE)
        if (calls > 0_C_INT) seconds_per_call = seconds/real(calls, C_DOUBLE)

        write(*,'(A,1X,A,1X,A,1X,I0,1X,A,1X,I0,1X,A,1X,ES16.8,1X,A,1X,ES16.8,1X,A,1X,ES16.8)') &
            trim(tag)//":", trim(label), "calls", calls, "nsteps", nsteps, &
            "seconds", seconds, "seconds_per_step", seconds_per_step, &
            "seconds_per_call", seconds_per_call
    end subroutine write_profiler_line

end module chron
