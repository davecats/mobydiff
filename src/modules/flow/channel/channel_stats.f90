module channel_stats
    use, intrinsic :: iso_c_binding
    use, intrinsic :: iso_fortran_env, only: int64
    use :: init, only: dns_type, grid_type, VAR_U, VAR_V, VAR_W, VAR_P
    use :: blocks, only: block_set_type
    use :: comm, only: comm_type, comm_allreduce_sum
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
    integer, parameter :: STAT_UP = 10
    integer, parameter :: STAT_VP = 11
    integer, parameter :: STAT_WP = 12
    integer, parameter :: STAT_K = 13
    integer, parameter :: STAT_EPSILON = 14
    integer, parameter :: CHANNEL_NSTAT = 14

    type, public :: channel_stats_type
        integer :: sample_interval = -1
        integer :: write_interval = -1
        character(len=256) :: file = "channel_stats.h5"
        character(len=256) :: runtime_file = "runtimedata.txt"
        logical :: runtime_header_written = .false.
        integer(int64) :: clock_start = 0_int64
        integer(int64) :: clock_rate = 0_int64
        integer(C_INT) :: clock_step_start = 0_C_INT
        integer(C_INT) :: last_write_step = -1_C_INT
        real(C_DOUBLE), allocatable :: sum(:)
        real(C_DOUBLE), allocatable :: count(:)
        real(C_DOUBLE), allocatable :: profile(:)
        real(C_DOUBLE), allocatable :: coord(:)
    contains
        procedure :: setup => channel_stats_setup
        procedure :: after_step => channel_stats_after_step
        procedure :: write_hdf5 => channel_stats_write_hdf5
        procedure :: finalize => channel_stats_finalize
    end type channel_stats_type

    interface
        function fdm_h5_write_channel_stats(file_name, nwall, nstat, step, t_current, wall_dir, re, &
                forcing, coord, profile, raw_sum, count) &
                bind(C, name="fdm_h5_write_channel_stats") result(ierr)
            import :: C_CHAR, C_INT, C_DOUBLE
            character(kind=C_CHAR), intent(in) :: file_name(*)
            integer(C_INT), value :: nwall, nstat, step, wall_dir
            real(C_DOUBLE), value :: t_current, re
            real(C_DOUBLE), intent(in) :: forcing(*), coord(*), profile(*), raw_sum(*), count(*)
            integer(C_INT) :: ierr
        end function fdm_h5_write_channel_stats

        function fdm_h5_read_channel_stats(file_name, nwall, nstat, step, t_current, raw_sum, count) &
                bind(C, name="fdm_h5_read_channel_stats") result(ierr)
            import :: C_CHAR, C_INT, C_DOUBLE
            character(kind=C_CHAR), intent(in) :: file_name(*)
            integer(C_INT), value :: nwall, nstat
            integer(C_INT), intent(inout) :: step
            real(C_DOUBLE), intent(inout) :: t_current
            real(C_DOUBLE), intent(inout) :: raw_sum(*), count(*)
            integer(C_INT) :: ierr
        end function fdm_h5_read_channel_stats
    end interface

contains

    subroutine channel_stats_setup(this, dns, g, c)
        class(channel_stats_type), intent(inout) :: this
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        type(comm_type), intent(in) :: c

        integer :: nwall, n

        if (this%sample_interval <= 0 .and. this%write_interval <= 0) return

        call system_clock(count=this%clock_start, count_rate=this%clock_rate)
        this%clock_step_start = dns%step_current
        this%runtime_header_written = .false.
        this%last_write_step = -1_C_INT

        nwall = int(dns%globalSize(2))
        allocate(this%sum(CHANNEL_NSTAT*nwall))
        allocate(this%count(nwall))
        allocate(this%profile(CHANNEL_NSTAT*nwall))
        allocate(this%coord(nwall))

        this%sum = 0.0d0
        this%count = 0.0d0
        this%profile = 0.0d0

        do n = 1, nwall
            this%coord(n) = 0.5d0*(g%yNode(n-1) + g%yNode(n))
        end do

        if (len_trim(dns%restart_file) > 0) then
            call read_channel_stats_restart(this, dns, c)
        end if

    end subroutine channel_stats_setup

    subroutine channel_stats_after_step(this, blk, dns, g, c)
        class(channel_stats_type), intent(inout) :: this
        type(block_set_type), intent(inout) :: blk
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        type(comm_type), intent(in) :: c

        logical :: sample_stats, write_stats
        real(C_DOUBLE), allocatable :: sample_sum(:), sample_count(:)

        if (.not. allocated(this%sum)) return

        sample_stats = interval_is_due(dns%step_current, this%sample_interval)
        write_stats = interval_is_due(dns%step_current, this%write_interval)

        if (sample_stats) then
            allocate(sample_sum(size(this%sum)))
            allocate(sample_count(size(this%count)))
            call collect_channel_sample(blk, dns, g, sample_sum, sample_count)
            call add_channel_sample(this, sample_sum, sample_count)
            call write_runtime_sample(this, dns, g, c, sample_sum, sample_count)
        end if

        if (write_stats) call this%write_hdf5(dns, c)
    end subroutine channel_stats_after_step

    logical function interval_is_due(step, interval) result(is_due)
        integer(C_INT), intent(in) :: step
        integer, intent(in) :: interval

        is_due = .false.
        if (interval <= 0) return
        is_due = modulo(int(step), interval) == 0
    end function interval_is_due

    ! Samples are gathered over this rank's box, i.e. block 1 in Phase 0.
    subroutine collect_channel_sample(blk, dns, g, sample_sum, sample_count)
        type(block_set_type), intent(inout) :: blk
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        real(C_DOUBLE), intent(inout) :: sample_sum(:), sample_count(:)

        integer :: i, j, k, s, nx, ny, nz, global_i, global_k, global_wall_idx, base
        real(C_DOUBLE) :: cell_area
        real(C_DOUBLE) :: p, kin, eps, velocity(3), sample(CHANNEL_NSTAT)

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))
        sample_sum = 0.0d0
        sample_count = 0.0d0

#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do collapse(3) &
        !$omp& map(to: dns, g, blk%q, blk%d1x, blk%d1y, blk%d1z) map(tofrom: sample_sum, sample_count) &
        !$omp& private(i,j,k,s,global_i,global_k,global_wall_idx,base,cell_area,p,kin,eps,velocity,sample)
#endif
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    global_i = int(dns%localSize(1,0)) + i - 1
                    global_k = int(dns%localSize(3,0)) + k - 1
                    global_wall_idx = int(dns%localSize(2,0)) + j - 1
                    base = CHANNEL_NSTAT*(global_wall_idx - 1)
                    cell_area = (g%xNode(global_i) - g%xNode(global_i - 1)) * &
                                (g%zNode(global_k) - g%zNode(global_k - 1))

                    call centered_velocity(blk, i, j, k, velocity)
                    p = blk%q(i,j,k,VAR_P,1)
                    kin = 0.5d0*sum(velocity*velocity)
                    eps = channel_dissipation(blk, dns, i, j, k)
                    sample = channel_sample(velocity, p, kin, eps)

                    !$omp atomic update
                    sample_count(global_wall_idx) = sample_count(global_wall_idx) + cell_area
                    do s = 1, CHANNEL_NSTAT
                        !$omp atomic update
                        sample_sum(base + s) = sample_sum(base + s) + cell_area*sample(s)
                    end do
                end do
            end do
        end do
#ifdef USE_OPENMP_OFFLOAD
        !$omp end target teams distribute parallel do
#endif
    end subroutine collect_channel_sample

    subroutine add_channel_sample(this, sample_sum, sample_count)
        class(channel_stats_type), intent(inout) :: this
        real(C_DOUBLE), intent(in) :: sample_sum(:), sample_count(:)

        this%sum = this%sum + sample_sum
        this%count = this%count + sample_count
    end subroutine add_channel_sample

    subroutine write_runtime_sample(this, dns, g, c, sample_sum, sample_count)
        class(channel_stats_type), intent(inout) :: this
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        type(comm_type), intent(in) :: c
        real(C_DOUBLE), intent(in) :: sample_sum(:), sample_count(:)

        real(C_DOUBLE), allocatable :: reduced_sum(:), reduced_count(:)
        real(C_DOUBLE) :: flow_values(2), volume_values(2)
        real(C_DOUBLE) :: wall_seconds_per_step

        allocate(reduced_sum(size(sample_sum)))
        allocate(reduced_count(size(sample_count)))
        reduced_sum = sample_sum
        reduced_count = sample_count
        call comm_allreduce_sum(c, reduced_sum)
        call comm_allreduce_sum(c, reduced_count)

        call channel_runtime_values(this, dns, g, reduced_sum, reduced_count, flow_values, volume_values)
        wall_seconds_per_step = channel_wall_seconds_per_step(this, dns)
        if (c%has_terminal) then
            call write_runtime_output(this, dns, flow_values, volume_values, wall_seconds_per_step)
        end if
    end subroutine write_runtime_sample

    subroutine channel_stats_write_hdf5(this, dns, c)
        class(channel_stats_type), intent(inout) :: this
        type(dns_type), intent(in) :: dns
        type(comm_type), intent(in) :: c

        character(kind=C_CHAR,len=:), allocatable :: c_file_name
        integer(C_INT) :: ierr
        integer :: nwall
        real(C_DOUBLE), allocatable :: raw_sum(:), reduced_count(:)

        if (.not. allocated(this%sum)) return
        if (this%last_write_step == dns%step_current) return

        nwall = size(this%count)
        allocate(raw_sum(size(this%sum)))
        allocate(reduced_count(size(this%count)))
        raw_sum = this%sum
        reduced_count = this%count
        call comm_allreduce_sum(c, raw_sum)
        call comm_allreduce_sum(c, reduced_count)
        if (sum(reduced_count) <= 0.0d0) return

        call build_channel_profile(this, raw_sum, reduced_count, dns)

        if (c%has_terminal) then
            c_file_name = to_c_string(this%file)
            ierr = fdm_h5_write_channel_stats(c_file_name, int(nwall, C_INT), int(CHANNEL_NSTAT, C_INT), &
                dns%step_current, dns%t_current, int(2, C_INT), dns%re, dns%forcing, &
                this%coord, this%profile, raw_sum, reduced_count)
            if (ierr /= 0_C_INT) then
                print *, "error: could not write channel statistics file: ", trim(this%file)
                error stop
            end if
        end if

        this%last_write_step = dns%step_current
    end subroutine channel_stats_write_hdf5

    subroutine channel_stats_finalize(this, dns, c)
        class(channel_stats_type), intent(inout) :: this
        type(dns_type), intent(in) :: dns
        type(comm_type), intent(in) :: c

        call this%write_hdf5(dns, c)
    end subroutine channel_stats_finalize

    subroutine centered_velocity(blk, i, j, k, velocity)
        type(block_set_type), intent(in) :: blk
        integer, intent(in) :: i, j, k
        real(C_DOUBLE), intent(out) :: velocity(3)

        velocity(1) = 0.5d0*(blk%q(i,j,k,VAR_U,1) + blk%q(i+1,j,k,VAR_U,1))
        velocity(2) = 0.5d0*(blk%q(i,j,k,VAR_V,1) + blk%q(i,j+1,k,VAR_V,1))
        velocity(3) = 0.5d0*(blk%q(i,j,k,VAR_W,1) + blk%q(i,j,k+1,VAR_W,1))
    end subroutine centered_velocity

    function channel_sample(velocity, p, kin, eps) result(sample)
        real(C_DOUBLE), intent(in) :: velocity(3), p, kin, eps
        real(C_DOUBLE) :: sample(CHANNEL_NSTAT)

        sample(STAT_U) = velocity(1)
        sample(STAT_V) = velocity(2)
        sample(STAT_W) = velocity(3)
        sample(STAT_UU) = velocity(1)*velocity(1)
        sample(STAT_VV) = velocity(2)*velocity(2)
        sample(STAT_WW) = velocity(3)*velocity(3)
        sample(STAT_UV) = velocity(1)*velocity(2)
        sample(STAT_UW) = velocity(1)*velocity(3)
        sample(STAT_VW) = velocity(2)*velocity(3)
        sample(STAT_UP) = velocity(1)*p
        sample(STAT_VP) = velocity(2)*p
        sample(STAT_WP) = velocity(3)*p
        sample(STAT_K) = kin
        sample(STAT_EPSILON) = eps
    end function channel_sample

    real(C_DOUBLE) function channel_dissipation(blk, dns, i, j, k) result(eps)
        type(block_set_type), intent(in) :: blk
        type(dns_type), intent(in) :: dns
        integer, intent(in) :: i, j, k

        real(C_DOUBLE) :: dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
        real(C_DOUBLE) :: s11, s22, s33, s12, s13, s23

        dudx = (blk%q(i+1,j,k,VAR_U,1) - blk%q(i,j,k,VAR_U,1))*blk%d1x(i,VAR_P,1)
        dvdy = (blk%q(i,j+1,k,VAR_V,1) - blk%q(i,j,k,VAR_V,1))*blk%d1y(j,VAR_P,1)
        dwdz = (blk%q(i,j,k+1,VAR_W,1) - blk%q(i,j,k,VAR_W,1))*blk%d1z(k,VAR_P,1)

        dudy = 0.25d0*((blk%q(i,j+1,k,VAR_U,1) + blk%q(i+1,j+1,k,VAR_U,1)) - &
                       (blk%q(i,j-1,k,VAR_U,1) + blk%q(i+1,j-1,k,VAR_U,1)))*blk%d1y(j,VAR_U,1)
        dudz = 0.25d0*((blk%q(i,j,k+1,VAR_U,1) + blk%q(i+1,j,k+1,VAR_U,1)) - &
                       (blk%q(i,j,k-1,VAR_U,1) + blk%q(i+1,j,k-1,VAR_U,1)))*blk%d1z(k,VAR_U,1)
        dvdx = 0.25d0*((blk%q(i+1,j,k,VAR_V,1) + blk%q(i+1,j+1,k,VAR_V,1)) - &
                       (blk%q(i-1,j,k,VAR_V,1) + blk%q(i-1,j+1,k,VAR_V,1)))*blk%d1x(i,VAR_V,1)
        dvdz = 0.25d0*((blk%q(i,j,k+1,VAR_V,1) + blk%q(i,j+1,k+1,VAR_V,1)) - &
                       (blk%q(i,j,k-1,VAR_V,1) + blk%q(i,j+1,k-1,VAR_V,1)))*blk%d1z(k,VAR_V,1)
        dwdx = 0.25d0*((blk%q(i+1,j,k,VAR_W,1) + blk%q(i+1,j,k+1,VAR_W,1)) - &
                       (blk%q(i-1,j,k,VAR_W,1) + blk%q(i-1,j,k+1,VAR_W,1)))*blk%d1x(i,VAR_W,1)
        dwdy = 0.25d0*((blk%q(i,j+1,k,VAR_W,1) + blk%q(i,j+1,k+1,VAR_W,1)) - &
                       (blk%q(i,j-1,k,VAR_W,1) + blk%q(i,j-1,k+1,VAR_W,1)))*blk%d1y(j,VAR_W,1)

        s11 = dudx
        s22 = dvdy
        s33 = dwdz
        s12 = 0.5d0*(dudy + dvdx)
        s13 = 0.5d0*(dudz + dwdx)
        s23 = 0.5d0*(dvdz + dwdy)
        eps = 2.0d0/dns%re*(s11*s11 + s22*s22 + s33*s33 + &
            2.0d0*(s12*s12 + s13*s13 + s23*s23))
    end function channel_dissipation

    subroutine build_channel_profile(this, raw_sum, reduced_count, dns)
        class(channel_stats_type), intent(inout) :: this
        real(C_DOUBLE), intent(in) :: raw_sum(:), reduced_count(:)
        type(dns_type), intent(in) :: dns

        integer :: n, s, base

        this%profile = raw_sum
        do n = 1, size(reduced_count)
            if (reduced_count(n) <= 0.0d0) cycle
            base = CHANNEL_NSTAT*(n - 1)
            do s = 1, CHANNEL_NSTAT
                this%profile(base+s) = this%profile(base+s)/reduced_count(n)
            end do
            this%profile(base+STAT_K) = 0.5d0*(this%profile(base+STAT_UU) + &
                this%profile(base+STAT_VV) + this%profile(base+STAT_WW) - &
                this%profile(base+STAT_U)**2 - this%profile(base+STAT_V)**2 - &
                this%profile(base+STAT_W)**2)
        end do
        call subtract_mean_profile_dissipation(this, dns, reduced_count)
    end subroutine build_channel_profile

    subroutine subtract_mean_profile_dissipation(this, dns, reduced_count)
        class(channel_stats_type), intent(inout) :: this
        type(dns_type), intent(in) :: dns
        real(C_DOUBLE), intent(in) :: reduced_count(:)

        integer :: n, nwall, base, mean_base
        real(C_DOUBLE), allocatable :: mean_velocity(:)
        real(C_DOUBLE) :: mean_eps

        nwall = size(reduced_count)
        if (nwall < 2) return

        allocate(mean_velocity(3*nwall))
        mean_velocity = 0.0d0
        do n = 1, nwall
            base = CHANNEL_NSTAT*(n - 1)
            mean_base = 3*(n - 1)
            mean_velocity(mean_base+1) = this%profile(base+STAT_U)
            mean_velocity(mean_base+2) = this%profile(base+STAT_V)
            mean_velocity(mean_base+3) = this%profile(base+STAT_W)
        end do

        do n = 1, nwall
            if (reduced_count(n) <= 0.0d0) cycle
            base = CHANNEL_NSTAT*(n - 1)
            mean_eps = mean_profile_dissipation_at(mean_velocity, this%coord, reduced_count, n, dns%re)
            this%profile(base+STAT_EPSILON) = max(0.0d0, this%profile(base+STAT_EPSILON) - mean_eps)
        end do
    end subroutine subtract_mean_profile_dissipation

    real(C_DOUBLE) function channel_wall_seconds_per_step(this, dns) result(seconds_per_step)
        class(channel_stats_type), intent(in) :: this
        type(dns_type), intent(in) :: dns

        integer(int64) :: clock_now
        integer(C_INT) :: nsteps

        seconds_per_step = 0.0d0
        if (this%clock_rate <= 0_int64) return

        nsteps = max(1_C_INT, dns%step_current - this%clock_step_start)
        call system_clock(count=clock_now)
        seconds_per_step = real(clock_now - this%clock_start, C_DOUBLE) / &
            (real(this%clock_rate, C_DOUBLE)*real(nsteps, C_DOUBLE))
    end function channel_wall_seconds_per_step

    subroutine write_runtime_output(this, dns, flow_values, volume_values, wall_seconds_per_step)
        class(channel_stats_type), intent(inout) :: this
        type(dns_type), intent(in) :: dns
        real(C_DOUBLE), intent(in) :: flow_values(2), volume_values(2), wall_seconds_per_step

        integer :: unit, stat
        character(len=*), parameter :: header = &
            "iteration time meanpx meanpz flowratex flowratez volume_k volume_epsilon wall_seconds_per_step"

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

        write(*,'(I10,8(1X,ES16.8))') int(dns%step_current), dns%t_current, &
            dns%forcing(1), dns%forcing(3), flow_values(1), flow_values(2), &
            volume_values(1), volume_values(2), wall_seconds_per_step

        open(newunit=unit, file=trim(this%runtime_file), status="old", position="append", action="write", &
            iostat=stat)
        if (stat == 0) then
            write(unit,'(I10,8(1X,ES16.8))') int(dns%step_current), dns%t_current, &
                dns%forcing(1), dns%forcing(3), flow_values(1), flow_values(2), &
                volume_values(1), volume_values(2), wall_seconds_per_step
            close(unit)
        else
            print *, "warning: could not append runtime data file: ", trim(this%runtime_file)
        end if
    end subroutine write_runtime_output

    subroutine channel_runtime_values(this, dns, g, sample_sum, sample_count, flow_values, volume_values)
        class(channel_stats_type), intent(in) :: this
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        real(C_DOUBLE), intent(in) :: sample_sum(:), sample_count(:)
        real(C_DOUBLE), intent(out) :: flow_values(2)
        real(C_DOUBLE), intent(out) :: volume_values(2)

        integer :: nwall, base, mean_base, n
        real(C_DOUBLE) :: dy, total_volume, mean_eps, k_turb, eps_turb
        real(C_DOUBLE) :: flow_integral(2), volume_integral(2)
        real(C_DOUBLE), allocatable :: mean_velocity(:)

        nwall = size(sample_count)
        allocate(mean_velocity(3*nwall))
        mean_velocity = 0.0d0

        flow_values = 0.0d0
        volume_values = 0.0d0
        flow_integral = 0.0d0
        volume_integral = 0.0d0
        total_volume = 0.0d0

        do n = 1, nwall
            if (sample_count(n) <= 0.0d0) cycle
            base = CHANNEL_NSTAT*(n - 1)
            mean_base = 3*(n - 1)
            mean_velocity(mean_base+1) = sample_sum(base+STAT_U)/sample_count(n)
            mean_velocity(mean_base+2) = sample_sum(base+STAT_V)/sample_count(n)
            mean_velocity(mean_base+3) = sample_sum(base+STAT_W)/sample_count(n)
        end do

        do n = 1, nwall
            if (sample_count(n) <= 0.0d0) cycle
            base = CHANNEL_NSTAT*(n - 1)
            mean_base = 3*(n - 1)
            dy = g%yNode(n) - g%yNode(n - 1)
            total_volume = total_volume + dy*sample_count(n)

            flow_integral(1) = flow_integral(1) + dy*sample_sum(base+STAT_U)
            flow_integral(2) = flow_integral(2) + dy*sample_sum(base+STAT_W)

            k_turb = max(0.0d0, sample_sum(base+STAT_K)/sample_count(n) - &
                0.5d0*(mean_velocity(mean_base+1)**2 + mean_velocity(mean_base+2)**2 + &
                mean_velocity(mean_base+3)**2))
            mean_eps = mean_profile_dissipation_at(mean_velocity, this%coord, sample_count, n, dns%re)
            eps_turb = max(0.0d0, sample_sum(base+STAT_EPSILON)/sample_count(n) - mean_eps)

            volume_integral(1) = volume_integral(1) + dy*sample_count(n)*k_turb
            volume_integral(2) = volume_integral(2) + dy*sample_count(n)*eps_turb
        end do

        if (total_volume > 0.0d0) then
            flow_values = flow_integral/total_volume
            volume_values = volume_integral/total_volume
        end if
    end subroutine channel_runtime_values

    real(C_DOUBLE) function mean_profile_dissipation_at(mean_velocity, coord, count, n, re) result(mean_eps)
        real(C_DOUBLE), intent(in) :: mean_velocity(:), coord(:), count(:)
        integer, intent(in) :: n
        real(C_DOUBLE), intent(in) :: re

        integer :: nwall, n0, n1, dir, base0, base1
        real(C_DOUBLE) :: dy, grad(3)

        mean_eps = 0.0d0
        nwall = size(count)
        if (nwall < 2 .or. re <= 0.0d0) return

        if (n <= 1) then
            n0 = 1
            n1 = 2
        else if (n >= nwall) then
            n0 = nwall - 1
            n1 = nwall
        else
            n0 = n - 1
            n1 = n + 1
        end if
        if (count(n0) <= 0.0d0 .or. count(n1) <= 0.0d0) return

        dy = coord(n1) - coord(n0)
        if (abs(dy) <= tiny(1.0d0)) return

        base0 = 3*(n0 - 1)
        base1 = 3*(n1 - 1)
        do dir = 1, 3
            grad(dir) = (mean_velocity(base1+dir) - mean_velocity(base0+dir))/dy
        end do

        mean_eps = (grad(1)**2 + grad(3)**2)/re + 2.0d0*grad(2)**2/re
    end function mean_profile_dissipation_at

    subroutine read_channel_stats_restart(this, dns, c)
        class(channel_stats_type), intent(inout) :: this
        type(dns_type), intent(in) :: dns
        type(comm_type), intent(in) :: c

        character(kind=C_CHAR,len=:), allocatable :: c_file_name
        integer(C_INT) :: ierr, restart_step
        real(C_DOUBLE) :: restart_time
        logical :: exists

        inquire(file=trim(this%file), exist=exists)
        if (.not. exists) return
        if (c%world_rank /= 0) return

        restart_step = 0_C_INT
        restart_time = 0.0d0
        c_file_name = to_c_string(this%file)
        ierr = fdm_h5_read_channel_stats(c_file_name, int(size(this%count), C_INT), &
            int(CHANNEL_NSTAT, C_INT), restart_step, restart_time, this%sum, this%count)
        if (ierr /= 0_C_INT) then
            if (c%has_terminal) print *, "error: could not read channel statistics file: ", trim(this%file)
            error stop
        else if (c%has_terminal) then
            print *, "continuing channel statistics from: ", trim(this%file)
        end if
    end subroutine read_channel_stats_restart

end module channel_stats
