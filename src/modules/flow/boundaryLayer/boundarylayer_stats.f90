module boundarylayer_stats
    ! Statistics for the spatially developing boundary layer. Unlike the
    ! channel (profiles of y, averaged in x, z, time), the boundary layer is
    ! INHOMOGENEOUS in x, so profiles are functions of (x, y) and averaging is
    ! only over the spanwise direction (z) and time. Accumulators live on the
    ! global base (x, y) plane (single level: the boundary-layer case does not
    ! refine), flattened y-fastest (row = (gx-1)*ny + gy), allreduce-summed and
    ! written by rank 0 -- the same recipe as channel_stats, one dimension up.
    !
    ! It also owns the runtime health line requested for this case: iteration,
    ! time, L2 and global (net) divergence residual, Linf velocity, CFL,
    ! Peclet, and wall-clock per step.
    use, intrinsic :: iso_c_binding
    use, intrinsic :: iso_fortran_env, only: int64
    use :: init, only: dns_type, grid_type, VAR_U, VAR_V, VAR_W, VAR_P, &
        CFL_COURANT, CFL_PECLET
    use :: blocks, only: block_set_type
    use :: comm, only: comm_type, comm_allreduce_sum, comm_allreduce_max
    use :: io, only: to_c_string
    implicit none

    private

    integer, parameter :: STAT_U = 1
    integer, parameter :: STAT_V = 2
    integer, parameter :: STAT_W = 3
    integer, parameter :: STAT_UU = 4
    integer, parameter :: STAT_VV = 5
    integer, parameter :: STAT_WW = 6
    integer, parameter :: STAT_UV = 7
    integer, parameter :: STAT_UW = 8
    integer, parameter :: STAT_VW = 9
    integer, parameter :: STAT_P = 10
    integer, parameter :: BL_NSTAT = 10

    type, public :: bl_stats_type
        integer :: sample_interval = -1
        integer :: write_interval = -1
        integer :: runtime_interval = -1
        character(len=256) :: file = "bl_stats.h5"
        character(len=256) :: runtime_file = "runtimedata.txt"
        logical :: runtime_header_written = .false.
        integer(int64) :: clock_start = 0_int64
        integer(int64) :: clock_rate = 0_int64
        integer(C_INT) :: clock_step_start = 0_C_INT
        integer(C_INT) :: last_write_step = -1_C_INT
        integer :: nx = 0
        integer :: ny = 0
        real(C_DOUBLE), allocatable :: sum(:)       ! (BL_NSTAT*nx*ny)
        real(C_DOUBLE), allocatable :: count(:)     ! (nx*ny)
        real(C_DOUBLE), allocatable :: profile(:)   ! (BL_NSTAT*nx*ny)
        real(C_DOUBLE), allocatable :: xcoord(:)    ! (nx)
        real(C_DOUBLE), allocatable :: ycoord(:)    ! (ny)
    contains
        procedure :: setup => bl_stats_setup
        procedure :: after_step => bl_stats_after_step
        procedure :: write_hdf5 => bl_stats_write_hdf5
        procedure :: finalize => bl_stats_finalize
    end type bl_stats_type

    interface
        function fdm_h5_write_bl_stats(file_name, nx, ny, nstat, step, t_current, re, &
                xcoord, ycoord, profile, raw_sum, count) &
                bind(C, name="fdm_h5_write_bl_stats") result(ierr)
            import :: C_CHAR, C_INT, C_DOUBLE
            character(kind=C_CHAR), intent(in) :: file_name(*)
            integer(C_INT), value :: nx, ny, nstat, step
            real(C_DOUBLE), value :: t_current, re
            real(C_DOUBLE), intent(in) :: xcoord(*), ycoord(*), profile(*), raw_sum(*), count(*)
            integer(C_INT) :: ierr
        end function fdm_h5_write_bl_stats

        function fdm_h5_read_bl_stats(file_name, nx, ny, nstat, step, t_current, raw_sum, count) &
                bind(C, name="fdm_h5_read_bl_stats") result(ierr)
            import :: C_CHAR, C_INT, C_DOUBLE
            character(kind=C_CHAR), intent(in) :: file_name(*)
            integer(C_INT), value :: nx, ny, nstat
            integer(C_INT), intent(inout) :: step
            real(C_DOUBLE), intent(inout) :: t_current
            real(C_DOUBLE), intent(inout) :: raw_sum(*), count(*)
            integer(C_INT) :: ierr
        end function fdm_h5_read_bl_stats
    end interface

contains

    subroutine bl_stats_setup(this, blk, dns, g, c)
        class(bl_stats_type), intent(inout) :: this
        type(block_set_type), intent(in) :: blk
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        type(comm_type), intent(in) :: c

        integer :: n, npt

        ! The runtime line and the profile sampling are independently gated.
        if (this%sample_interval <= 0 .and. this%write_interval <= 0 .and. &
            this%runtime_interval <= 0) return

        call system_clock(count=this%clock_start, count_rate=this%clock_rate)
        this%clock_step_start = dns%step_current
        this%runtime_header_written = .false.
        this%last_write_step = -1_C_INT

        if (blk%nLevels > 1_C_INT .and. c%has_terminal) &
            print *, "warning: boundaryLayer stats assume a single level; refinement ignored"

        this%nx = int(dns%globalSize(1))
        this%ny = int(dns%globalSize(2))
        npt = this%nx*this%ny

        if (this%sample_interval > 0 .or. this%write_interval > 0) then
            allocate(this%sum(BL_NSTAT*npt), this%count(npt), this%profile(BL_NSTAT*npt))
            this%sum = 0.0d0
            this%count = 0.0d0
            this%profile = 0.0d0
            allocate(this%xcoord(this%nx), this%ycoord(this%ny))
            do n = 1, this%nx
                this%xcoord(n) = 0.5d0*(g%xNode(n-1) + g%xNode(n))
            end do
            do n = 1, this%ny
                this%ycoord(n) = 0.5d0*(g%yNode(n-1) + g%yNode(n))
            end do
            if (len_trim(dns%restart_file) > 0) call read_bl_stats_restart(this, dns, c)
        end if
    end subroutine bl_stats_setup

    subroutine bl_stats_after_step(this, blk, dns, g, c)
        class(bl_stats_type), intent(inout) :: this
        type(block_set_type), intent(inout) :: blk
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        type(comm_type), intent(in) :: c

        real(C_DOUBLE), allocatable :: sample_sum(:), sample_count(:)

        ! Runtime health line: the requested diagnostics.
        if (interval_is_due(dns%step_current, this%runtime_interval)) &
            call write_runtime_line(this, blk, dns, c)

        if (.not. allocated(this%sum)) return

        if (interval_is_due(dns%step_current, this%sample_interval)) then
            allocate(sample_sum(size(this%sum)), sample_count(size(this%count)))
            call collect_bl_sample(this, blk, sample_sum, sample_count)
            this%sum = this%sum + sample_sum
            this%count = this%count + sample_count
        end if

        if (interval_is_due(dns%step_current, this%write_interval)) call this%write_hdf5(dns, c)
    end subroutine bl_stats_after_step

    logical function interval_is_due(step, interval) result(is_due)
        integer(C_INT), intent(in) :: step
        integer, intent(in) :: interval
        is_due = .false.
        if (interval <= 0) return
        is_due = modulo(int(step), interval) == 0
    end function interval_is_due

    ! Per (x,y) cell: sum over the spanwise direction (z) with the cell's
    ! z-width weight; time averaging accumulates across samples.
    subroutine collect_bl_sample(this, blk, sample_sum, sample_count)
        class(bl_stats_type), intent(in) :: this
        type(block_set_type), intent(inout) :: blk
        real(C_DOUBLE), intent(inout) :: sample_sum(:), sample_count(:)

        integer :: i, j, k, b, s, nx, ny, nz, nBlocks, gx, gy, row, base, ny_g
        real(C_DOUBLE) :: dz, p, velocity(3), sample(BL_NSTAT)

        nx = int(blk%nb(1)); ny = int(blk%nb(2)); nz = int(blk%nb(3))
        nBlocks = int(blk%nBlocks)
        ny_g = this%ny
        sample_sum = 0.0d0
        sample_count = 0.0d0

#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do collapse(4) &
        !$omp& map(to: ny_g, blk%origin, blk%q, blk%z) &
        !$omp& map(tofrom: sample_sum, sample_count) &
        !$omp& private(i,j,k,b,s,gx,gy,row,base,dz,p,velocity,sample)
#endif
        do b = 1, nBlocks
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    gx = int(blk%origin(1,b)) + i
                    gy = int(blk%origin(2,b)) + j
                    row = (gx - 1)*ny_g + gy
                    base = BL_NSTAT*(row - 1)
                    dz = blk%z(k+1,VAR_W,b) - blk%z(k,VAR_W,b)

                    velocity(1) = 0.5d0*(blk%q(i,j,k,VAR_U,b) + blk%q(i+1,j,k,VAR_U,b))
                    velocity(2) = 0.5d0*(blk%q(i,j,k,VAR_V,b) + blk%q(i,j+1,k,VAR_V,b))
                    velocity(3) = 0.5d0*(blk%q(i,j,k,VAR_W,b) + blk%q(i,j,k+1,VAR_W,b))
                    p = blk%q(i,j,k,VAR_P,b)
                    sample(STAT_U) = velocity(1)
                    sample(STAT_V) = velocity(2)
                    sample(STAT_W) = velocity(3)
                    sample(STAT_UU) = velocity(1)*velocity(1)
                    sample(STAT_VV) = velocity(2)*velocity(2)
                    sample(STAT_WW) = velocity(3)*velocity(3)
                    sample(STAT_UV) = velocity(1)*velocity(2)
                    sample(STAT_UW) = velocity(1)*velocity(3)
                    sample(STAT_VW) = velocity(2)*velocity(3)
                    sample(STAT_P) = p

                    !$omp atomic update
                    sample_count(row) = sample_count(row) + dz
                    do s = 1, BL_NSTAT
                        !$omp atomic update
                        sample_sum(base + s) = sample_sum(base + s) + dz*sample(s)
                    end do
                end do
            end do
        end do
        end do
#ifdef USE_OPENMP_OFFLOAD
        !$omp end target teams distribute parallel do
#endif
    end subroutine collect_bl_sample

    subroutine bl_stats_write_hdf5(this, dns, c)
        class(bl_stats_type), intent(inout) :: this
        type(dns_type), intent(in) :: dns
        type(comm_type), intent(in) :: c

        character(kind=C_CHAR,len=:), allocatable :: c_file_name
        integer(C_INT) :: ierr
        real(C_DOUBLE), allocatable :: raw_sum(:), reduced_count(:)
        integer :: n, s, base

        if (.not. allocated(this%sum)) return
        if (this%last_write_step == dns%step_current) return

        allocate(raw_sum(size(this%sum)), reduced_count(size(this%count)))
        raw_sum = this%sum
        reduced_count = this%count
        call comm_allreduce_sum(c, raw_sum)
        call comm_allreduce_sum(c, reduced_count)
        if (sum(reduced_count) <= 0.0d0) return

        ! Time+span mean per (x,y) cell (raw_sum / count).
        this%profile = raw_sum
        do n = 1, size(reduced_count)
            if (reduced_count(n) <= 0.0d0) cycle
            base = BL_NSTAT*(n - 1)
            do s = 1, BL_NSTAT
                this%profile(base+s) = this%profile(base+s)/reduced_count(n)
            end do
        end do

        if (c%has_terminal) then
            c_file_name = to_c_string(this%file)
            ierr = fdm_h5_write_bl_stats(c_file_name, int(this%nx, C_INT), int(this%ny, C_INT), &
                int(BL_NSTAT, C_INT), dns%step_current, dns%t_current, dns%re, &
                this%xcoord, this%ycoord, this%profile, raw_sum, reduced_count)
            if (ierr /= 0_C_INT) then
                print *, "error: could not write boundaryLayer statistics file: ", trim(this%file)
                error stop
            end if
        end if
        this%last_write_step = dns%step_current
    end subroutine bl_stats_write_hdf5

    subroutine bl_stats_finalize(this, dns, c)
        class(bl_stats_type), intent(inout) :: this
        type(dns_type), intent(in) :: dns
        type(comm_type), intent(in) :: c

        call this%write_hdf5(dns, c)
    end subroutine bl_stats_finalize

    ! Runtime health line. The divergence residual (volume-weighted L2 and the
    ! net global integral = mass-conservation check), the Linf velocity, plus
    ! CFL/Peclet (from the last update_timestep_limits) and wall-clock/step.
    subroutine write_runtime_line(this, blk, dns, c)
        class(bl_stats_type), intent(inout) :: this
        type(block_set_type), intent(inout) :: blk
        type(dns_type), intent(in) :: dns
        type(comm_type), intent(in) :: c

        real(C_DOUBLE) :: sums(3), maxes(1)
        real(C_DOUBLE) :: l2_div, global_div, linf_vel, cfl, peclet, seconds_per_step
        integer(int64) :: clock_now
        integer(C_INT) :: nsteps
        integer :: unit, stat
        character(len=*), parameter :: header = &
            "iteration time L2_div global_div Linf_vel CFL Peclet wall_seconds_per_step"

        call divergence_reductions(blk, sums, maxes)
        ! sums = [sum(div^2*vol), sum(vol), sum(div*vol)]; maxes = [Linf vel]
        call comm_allreduce_sum(c, sums)
        call comm_allreduce_max(c, maxes)

        l2_div = 0.0d0
        if (sums(2) > 0.0d0) l2_div = sqrt(max(0.0d0, sums(1)/sums(2)))
        global_div = sums(3)
        linf_vel = maxes(1)
        cfl = dns%cfl(CFL_COURANT)
        peclet = dns%cfl(CFL_PECLET)

        seconds_per_step = 0.0d0
        if (this%clock_rate > 0_int64) then
            nsteps = max(1_C_INT, dns%step_current - this%clock_step_start)
            call system_clock(count=clock_now)
            seconds_per_step = real(clock_now - this%clock_start, C_DOUBLE) / &
                (real(this%clock_rate, C_DOUBLE)*real(nsteps, C_DOUBLE))
        end if

        if (.not. c%has_terminal) return

        if (.not. this%runtime_header_written) then
            write(*,'(A)') header
            open(newunit=unit, file=trim(this%runtime_file), status="replace", action="write", iostat=stat)
            if (stat == 0) then
                write(unit,'(A)') header
                close(unit)
            else
                print *, "warning: could not open runtime data file: ", trim(this%runtime_file)
            end if
            this%runtime_header_written = .true.
        end if

        write(*,'(I10,7(1X,ES16.8))') int(dns%step_current), dns%t_current, &
            l2_div, global_div, linf_vel, cfl, peclet, seconds_per_step
        open(newunit=unit, file=trim(this%runtime_file), status="old", position="append", &
            action="write", iostat=stat)
        if (stat == 0) then
            write(unit,'(I10,7(1X,ES16.8))') int(dns%step_current), dns%t_current, &
                l2_div, global_div, linf_vel, cfl, peclet, seconds_per_step
            close(unit)
        end if
    end subroutine write_runtime_line

    ! Per-rank reductions of the discrete divergence and velocity magnitude
    ! over interior cells (the same divergence stencil as the projection).
    subroutine divergence_reductions(blk, sums, maxes)
        type(block_set_type), intent(inout) :: blk
        real(C_DOUBLE), intent(out) :: sums(3), maxes(1)

        integer :: i, j, k, b, nx, ny, nz, nBlocks
        real(C_DOUBLE) :: div, vol, uc, vc, wc, vmag
        real(C_DOUBLE) :: s_div2, s_vol, s_div, m_vel

        nx = int(blk%nb(1)); ny = int(blk%nb(2)); nz = int(blk%nb(3))
        nBlocks = int(blk%nBlocks)
        s_div2 = 0.0d0; s_vol = 0.0d0; s_div = 0.0d0; m_vel = 0.0d0

#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do collapse(4) &
        !$omp& reduction(+: s_div2, s_vol, s_div) reduction(max: m_vel) &
        !$omp& map(to: blk%q, blk%d1x, blk%d1y, blk%d1z) &
        !$omp& private(i,j,k,b,div,vol,uc,vc,wc,vmag)
#endif
        do b = 1, nBlocks
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    div = (blk%q(i+1,j,k,VAR_U,b)-blk%q(i,j,k,VAR_U,b))*blk%d1x(i,VAR_P,b) &
                        + (blk%q(i,j+1,k,VAR_V,b)-blk%q(i,j,k,VAR_V,b))*blk%d1y(j,VAR_P,b) &
                        + (blk%q(i,j,k+1,VAR_W,b)-blk%q(i,j,k,VAR_W,b))*blk%d1z(k,VAR_P,b)
                    vol = 1.0d0/(blk%d1x(i,VAR_P,b)*blk%d1y(j,VAR_P,b)*blk%d1z(k,VAR_P,b))
                    uc = 0.5d0*(blk%q(i,j,k,VAR_U,b) + blk%q(i+1,j,k,VAR_U,b))
                    vc = 0.5d0*(blk%q(i,j,k,VAR_V,b) + blk%q(i,j+1,k,VAR_V,b))
                    wc = 0.5d0*(blk%q(i,j,k,VAR_W,b) + blk%q(i,j,k+1,VAR_W,b))
                    vmag = sqrt(uc*uc + vc*vc + wc*wc)
                    s_div2 = s_div2 + div*div*vol
                    s_vol  = s_vol + vol
                    s_div  = s_div + div*vol
                    m_vel  = max(m_vel, vmag)
                end do
            end do
        end do
        end do
#ifdef USE_OPENMP_OFFLOAD
        !$omp end target teams distribute parallel do
#endif
        sums(1) = s_div2; sums(2) = s_vol; sums(3) = s_div
        maxes(1) = m_vel
    end subroutine divergence_reductions

    subroutine read_bl_stats_restart(this, dns, c)
        class(bl_stats_type), intent(inout) :: this
        type(dns_type), intent(in) :: dns
        type(comm_type), intent(in) :: c

        character(kind=C_CHAR,len=:), allocatable :: c_file_name
        integer(C_INT) :: ierr, restart_step
        real(C_DOUBLE) :: restart_time
        logical :: exists

        if (c%world_rank /= 0) return
        inquire(file=trim(this%file), exist=exists)
        if (.not. exists) return

        restart_step = 0_C_INT
        restart_time = 0.0d0
        c_file_name = to_c_string(this%file)
        ierr = fdm_h5_read_bl_stats(c_file_name, int(this%nx, C_INT), int(this%ny, C_INT), &
            int(BL_NSTAT, C_INT), restart_step, restart_time, this%sum, this%count)
        if (ierr /= 0_C_INT) then
            if (c%has_terminal) print *, "warning: could not read boundaryLayer statistics; starting fresh: ", &
                trim(this%file)
            this%sum = 0.0d0
            this%count = 0.0d0
        else if (c%has_terminal) then
            print *, "continuing boundaryLayer statistics from: ", trim(this%file)
        end if
    end subroutine read_bl_stats_restart

end module boundarylayer_stats
