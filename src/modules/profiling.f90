! Step-phase timing for the solver loop ([output] profile = true).
!
! WHY MODULE-LEVEL STATE: the timers have to sit inside comm.f90 and
! pressure_solver.f90, whose instrumented routines are called from a dozen
! places each. Threading a profiler argument through all of them would cost
! more readability than a diagnostic is worth, so the three profilers live here
! as module state and every hook is
!
!     t0 = prof_tic()
!     call <the work>
!     call prof_toc(step_prof, PROF_MOMENTUM, t0)
!
! Disabled (the default) prof_tic/prof_toc are one branch on a saved logical and
! the clock is never read, so production runs are unaffected. Timing only reads
! clocks: results are bit-exact either way.
!
! THE THREE PROFILERS ARE NESTED, NOT DISJOINT:
!
!   step_timing   one bucket per phase of an RK substage; its buckets ARE
!                 disjoint and their sum should cover the chron loop time
!                 (write_step_profilers prints that coverage -- if it is not
!                 close to 1 the instrumentation is lying and nothing derived
!                 from it can be trusted).
!   proj_timing   the inside of step_timing's `projection` bucket.
!   exch_timing   the inside of every halo exchange, wherever called from
!                 (step_timing's `vel_exchange` and `turbulence`, proj_timing's
!                 `phi_exchange` and `vel_exchange`).
!
! At ONE rank there are no peers, so pack/mpi_post/mpi_wait/unpack are exactly
! zero and the whole exchange shows up as local_copy -- which is the intended
! reading: local_copy is device-local halo traffic, the MPI buckets are network.
module profiling
    use, intrinsic :: iso_c_binding, only: C_DOUBLE, C_INT
    use :: chron, only: profiler_type, init_profiler, profiler_add, write_profiler, &
        profiler_total, wall_seconds

    implicit none
    private

    ! step_timing: the phases of one RK substage, plus the per-step tail.
    integer, parameter, public :: PROF_MOMENTUM = 1     ! fused predictor kernel
    integer, parameter, public :: PROF_IBM_MU = 2       ! update_ibm_mu (pointwise, halo-carrying)
    integer, parameter, public :: PROF_BODYFORCE = 3
    integer, parameter, public :: PROF_TURBULENCE = 4   ! nut producers + their scalar exchange
    integer, parameter, public :: PROF_APPLY_BC = 5     ! post-predictor apply_bc
    integer, parameter, public :: PROF_VEL_EXCHANGE = 6 ! post-predictor exchange (syncface)
    integer, parameter, public :: PROF_PROJECTION = 7
    integer, parameter, public :: PROF_IO_STATS = 8     ! dt limits, snapshots, case after_step

    ! proj_timing: the inside of PROF_PROJECTION.
    integer, parameter, public :: PROF_SWEEP = 1        ! compute-phi / red-black sweep (+cheb combine)
    integer, parameter, public :: PROF_APPLY = 2        ! jacobi_apply
    integer, parameter, public :: PROF_PHI_EXCHANGE = 3 ! scalar phi halos (+ outlet ghosts)
    integer, parameter, public :: PROF_PROJ_VEL_EXCHANGE = 4
    integer, parameter, public :: PROF_PROJ_BC = 5
    integer, parameter, public :: PROF_PROJ_SETUP = 6   ! buffer allocation, face flags

    ! exch_timing: the inside of every halo exchange, vector and scalar alike.
    integer, parameter, public :: PROF_PACK = 1
    integer, parameter, public :: PROF_MPI_POST = 2     ! Irecv/Isend (incl. use_device_addr)
    integer, parameter, public :: PROF_MPI_WAIT = 3
    integer, parameter, public :: PROF_UNPACK = 4
    integer, parameter, public :: PROF_LOCAL_COPY = 5   ! same-rank block-pair copies

    type(profiler_type), save, public :: step_prof, proj_prof, exch_prof
    logical, save :: profEnabled = .false.

    public :: init_step_profilers, write_step_profilers, prof_tic, prof_toc, profiling_enabled

contains

    subroutine init_step_profilers(enabled)
        logical, intent(in) :: enabled

        profEnabled = enabled
        call init_profiler(step_prof, "step_timing", &
            [character(len=24) :: "momentum", "ibm_mu", "bodyforce", "turbulence", &
             "apply_bc", "vel_exchange", "projection", "io_stats"])
        call init_profiler(proj_prof, "proj_timing", &
            [character(len=24) :: "sweep", "apply", "phi_exchange", "vel_exchange", &
             "apply_bc", "setup"])
        call init_profiler(exch_prof, "exch_timing", &
            [character(len=24) :: "pack", "mpi_post", "mpi_wait", "unpack", "local_copy"])
    end subroutine init_step_profilers

    logical function profiling_enabled()
        profiling_enabled = profEnabled
    end function profiling_enabled

    ! Clock reading that costs nothing when profiling is off.
    real(C_DOUBLE) function prof_tic() result(t)
        t = 0.0d0
        if (profEnabled) t = wall_seconds()
    end function prof_tic

    subroutine prof_toc(profile, category, t0)
        type(profiler_type), intent(inout) :: profile
        integer, intent(in) :: category
        real(C_DOUBLE), intent(in) :: t0

        if (.not. profEnabled) return
        call profiler_add(profile, category, wall_seconds() - t0)
    end subroutine prof_toc

    ! Print the three profilers plus the coverage of step_timing against the
    ! loop timer. Coverage is THE validity check on the whole measurement: the
    ! target regions carry no `nowait`, so the host blocks on each kernel and
    ! wall-clock brackets are meaningful -- but only if they add up.
    subroutine write_step_profilers(nsteps, loop_seconds)
        integer(C_INT), intent(in) :: nsteps
        real(C_DOUBLE), intent(in) :: loop_seconds

        real(C_DOUBLE) :: covered

        if (.not. profEnabled) return

        call write_profiler(step_prof, nsteps)
        covered = profiler_total(step_prof)
        if (loop_seconds > 0.0d0) then
            write(*,'(A,1X,ES16.8,1X,A,1X,ES16.8,1X,A,1X,F8.4)') &
                "step_timing: coverage measured", covered, "loop_seconds", loop_seconds, &
                "fraction", covered/loop_seconds
        end if
        ! proj_timing sits inside step_timing's projection bucket and exch_timing
        ! inside every exchange -- their totals are subsets, not additions.
        call write_profiler(proj_prof, nsteps)
        call write_profiler(exch_prof, nsteps)
    end subroutine write_step_profilers

end module profiling
