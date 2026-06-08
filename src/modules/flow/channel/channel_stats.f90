module channel_stats
    use, intrinsic :: iso_c_binding
    use, intrinsic :: iso_fortran_env, only: int64
    use :: init, only: dns_type, grid_type, field_type, VAR_U, VAR_V, VAR_W, VAR_P, &
        CFL_COURANT, CFL_PECLET
    use :: comm, only: comm_type, comm_allreduce_sum
    use :: io, only: to_c_string
    use :: channel_profile, only: channel_plane_coord
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
        integer :: interval = 100
        integer :: runtime_interval = -1
        character(len=256) :: file = "channel_stats.h5"
        character(len=256) :: runtime_file = "runtimedata.txt"
        logical :: on_device = .false.
        logical :: runtime_header_written = .false.
        integer(int64) :: clock_start = 0_int64
        integer(int64) :: clock_rate = 0_int64
        integer(C_INT) :: clock_step_start = 0_C_INT
        real(C_DOUBLE), allocatable :: sum(:)
        real(C_DOUBLE), allocatable :: count(:)
        real(C_DOUBLE), allocatable :: profile(:)
        real(C_DOUBLE), allocatable :: coord(:)
    contains
        procedure :: setup => channel_stats_setup
        procedure :: accumulate => channel_stats_accumulate
        procedure :: write => channel_stats_write
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

    subroutine channel_stats_setup(this, dns, g, wall_dir, c)
        class(channel_stats_type), intent(inout) :: this
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        integer(C_INT), intent(in) :: wall_dir
        type(comm_type), intent(in) :: c

        integer :: nwall, n

        if (this%interval <= 0 .and. effective_runtime_interval(this) <= 0) return

        call system_clock(count=this%clock_start, count_rate=this%clock_rate)
        this%clock_step_start = dns%step_current
        this%runtime_header_written = .false.

        nwall = int(dns%globalSize(wall_dir))
        allocate(this%sum(CHANNEL_NSTAT*nwall))
        allocate(this%count(nwall))
        allocate(this%profile(CHANNEL_NSTAT*nwall))
        allocate(this%coord(nwall))

        this%sum = 0.0d0
        this%count = 0.0d0
        this%profile = 0.0d0

        do n = 1, nwall
            this%coord(n) = channel_plane_coord(g, wall_dir, n)
        end do

        if (this%interval > 0) call read_channel_stats_restart(this, dns, c)

#ifdef USE_OPENMP_OFFLOAD
        !$omp target enter data map(to: this%sum, this%count)
        this%on_device = .true.
#endif
    end subroutine channel_stats_setup

    subroutine channel_stats_accumulate(this, f, dns, g, wall_dir)
        class(channel_stats_type), intent(inout) :: this
        type(field_type), intent(inout) :: f
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        integer(C_INT), intent(in) :: wall_dir

        integer :: i, j, k, s, nx, ny, nz, wall_idx, global_wall_idx, base
        real(C_DOUBLE) :: p, kin, eps, velocity(3), sample(CHANNEL_NSTAT)

        if (.not. allocated(this%sum)) return

        nx = int(dns%localSize(1,2))
        ny = int(dns%localSize(2,2))
        nz = int(dns%localSize(3,2))

#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do collapse(3) &
        !$omp& map(to: this, dns, g, f%q) map(tofrom: this%sum, this%count) &
        !$omp& private(i,j,k,s,wall_idx,global_wall_idx,base,p,kin,eps,velocity,sample)
#endif
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    select case (wall_dir)
                    case (1)
                        wall_idx = i
                    case (2)
                        wall_idx = j
                    case default
                        wall_idx = k
                    end select
                    global_wall_idx = int(dns%localSize(wall_dir,0)) + wall_idx - 1
                    base = CHANNEL_NSTAT*(global_wall_idx - 1)

                    call centered_velocity(f, i, j, k, velocity)
                    p = f%q(i,j,k,VAR_P)
                    kin = 0.5d0*sum(velocity*velocity)
                    eps = channel_dissipation(f, dns, g, i, j, k)
                    sample = channel_sample(velocity, p, kin, eps)

                    !$omp atomic update
                    this%count(global_wall_idx) = this%count(global_wall_idx) + 1.0d0
                    do s = 1, CHANNEL_NSTAT
                        !$omp atomic update
                        this%sum(base + s) = this%sum(base + s) + sample(s)
                    end do
                end do
            end do
        end do
#ifdef USE_OPENMP_OFFLOAD
        !$omp end target teams distribute parallel do
#endif
    end subroutine channel_stats_accumulate

    subroutine channel_stats_write(this, f, dns, g, c, wall_dir, write_hdf5, write_runtime)
        class(channel_stats_type), intent(inout) :: this
        type(field_type), intent(inout) :: f
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        type(comm_type), intent(in) :: c
        integer(C_INT), intent(in) :: wall_dir
        logical, intent(in) :: write_hdf5, write_runtime

        character(kind=C_CHAR,len=:), allocatable :: c_file_name
        integer(C_INT) :: ierr
        integer :: nwall
        real(C_DOUBLE), allocatable :: raw_sum(:), reduced_count(:)
        real(C_DOUBLE) :: flow_values(2), volume_values(2)
        real(C_DOUBLE) :: wall_seconds_per_step

        if (.not. (write_hdf5 .or. write_runtime)) return
        if (.not. allocated(this%sum)) return

#ifdef USE_OPENMP_OFFLOAD
        !$omp target update from(this%sum, this%count)
#endif

        nwall = size(this%count)
        allocate(raw_sum(size(this%sum)))
        allocate(reduced_count(size(this%count)))
        raw_sum = this%sum
        reduced_count = this%count
        call comm_allreduce_sum(c, raw_sum)
        call comm_allreduce_sum(c, reduced_count)

        call build_channel_profile(this, raw_sum, reduced_count)

        if (write_runtime) then
            call channel_current_diagnostics(f, dns, g, c, flow_values)
            call channel_volume_averages(this, reduced_count, volume_values)
            wall_seconds_per_step = channel_wall_seconds_per_step(this, dns)
            if (c%has_terminal) then
                call write_runtime_output(this, dns, flow_values, volume_values, wall_seconds_per_step)
            end if
        end if

        if (write_hdf5 .and. c%has_terminal) then
            c_file_name = to_c_string(this%file)
            ierr = fdm_h5_write_channel_stats(c_file_name, int(nwall, C_INT), int(CHANNEL_NSTAT, C_INT), &
                dns%step_current, dns%t_current, wall_dir, dns%re, dns%forcing, &
                this%coord, this%profile, raw_sum, reduced_count)
            if (ierr /= 0_C_INT) then
                print *, "error: could not write channel statistics file: ", trim(this%file)
                error stop
            end if
        end if
    end subroutine channel_stats_write

    subroutine channel_stats_finalize(this)
        class(channel_stats_type), intent(inout) :: this

#ifdef USE_OPENMP_OFFLOAD
        if (this%on_device) then
            !$omp target exit data map(delete: this%sum, this%count)
        end if
#endif
        this%on_device = .false.
    end subroutine channel_stats_finalize

    subroutine centered_velocity(f, i, j, k, velocity)
        type(field_type), intent(in) :: f
        integer, intent(in) :: i, j, k
        real(C_DOUBLE), intent(out) :: velocity(3)

        velocity(1) = 0.5d0*(f%q(i,j,k,VAR_U) + f%q(i+1,j,k,VAR_U))
        velocity(2) = 0.5d0*(f%q(i,j,k,VAR_V) + f%q(i,j+1,k,VAR_V))
        velocity(3) = 0.5d0*(f%q(i,j,k,VAR_W) + f%q(i,j,k+1,VAR_W))
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

    real(C_DOUBLE) function channel_dissipation(f, dns, g, i, j, k) result(eps)
        type(field_type), intent(in) :: f
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        integer, intent(in) :: i, j, k

        real(C_DOUBLE) :: dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
        real(C_DOUBLE) :: s11, s22, s33, s12, s13, s23

        dudx = (f%q(i+1,j,k,VAR_U) - f%q(i,j,k,VAR_U))*g%d1x(i,VAR_P)
        dvdy = (f%q(i,j+1,k,VAR_V) - f%q(i,j,k,VAR_V))*g%d1y(j,VAR_P)
        dwdz = (f%q(i,j,k+1,VAR_W) - f%q(i,j,k,VAR_W))*g%d1z(k,VAR_P)

        dudy = 0.25d0*((f%q(i,j+1,k,VAR_U) + f%q(i+1,j+1,k,VAR_U)) - &
                       (f%q(i,j-1,k,VAR_U) + f%q(i+1,j-1,k,VAR_U)))*g%d1y(j,VAR_U)
        dudz = 0.25d0*((f%q(i,j,k+1,VAR_U) + f%q(i+1,j,k+1,VAR_U)) - &
                       (f%q(i,j,k-1,VAR_U) + f%q(i+1,j,k-1,VAR_U)))*g%d1z(k,VAR_U)
        dvdx = 0.25d0*((f%q(i+1,j,k,VAR_V) + f%q(i+1,j+1,k,VAR_V)) - &
                       (f%q(i-1,j,k,VAR_V) + f%q(i-1,j+1,k,VAR_V)))*g%d1x(i,VAR_V)
        dvdz = 0.25d0*((f%q(i,j,k+1,VAR_V) + f%q(i,j+1,k+1,VAR_V)) - &
                       (f%q(i,j,k-1,VAR_V) + f%q(i,j+1,k-1,VAR_V)))*g%d1z(k,VAR_V)
        dwdx = 0.25d0*((f%q(i+1,j,k,VAR_W) + f%q(i+1,j,k+1,VAR_W)) - &
                       (f%q(i-1,j,k,VAR_W) + f%q(i-1,j,k+1,VAR_W)))*g%d1x(i,VAR_W)
        dwdy = 0.25d0*((f%q(i,j+1,k,VAR_W) + f%q(i,j+1,k+1,VAR_W)) - &
                       (f%q(i,j-1,k,VAR_W) + f%q(i,j-1,k+1,VAR_W)))*g%d1y(j,VAR_W)

        s11 = dudx
        s22 = dvdy
        s33 = dwdz
        s12 = 0.5d0*(dudy + dvdx)
        s13 = 0.5d0*(dudz + dwdx)
        s23 = 0.5d0*(dvdz + dwdy)
        eps = 2.0d0/dns%re*(s11*s11 + s22*s22 + s33*s33 + &
            2.0d0*(s12*s12 + s13*s13 + s23*s23))
    end function channel_dissipation

    subroutine build_channel_profile(this, raw_sum, reduced_count)
        class(channel_stats_type), intent(inout) :: this
        real(C_DOUBLE), intent(in) :: raw_sum(:), reduced_count(:)

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
    end subroutine build_channel_profile

    subroutine channel_volume_averages(this, reduced_count, volume_values)
        class(channel_stats_type), intent(in) :: this
        real(C_DOUBLE), intent(in) :: reduced_count(:)
        real(C_DOUBLE), intent(out) :: volume_values(2)

        integer :: n, base
        real(C_DOUBLE) :: count_sum

        volume_values = 0.0d0
        count_sum = 0.0d0
        do n = 1, size(reduced_count)
            base = CHANNEL_NSTAT*(n - 1)
            volume_values(1) = volume_values(1) + this%profile(base+STAT_K)*reduced_count(n)
            volume_values(2) = volume_values(2) + this%profile(base+STAT_EPSILON)*reduced_count(n)
            count_sum = count_sum + reduced_count(n)
        end do
        if (count_sum > 0.0d0) volume_values = volume_values/count_sum
    end subroutine channel_volume_averages

    integer function effective_runtime_interval(this) result(interval)
        class(channel_stats_type), intent(in) :: this

        interval = this%runtime_interval
        if (interval < 0) interval = this%interval
    end function effective_runtime_interval

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

    subroutine channel_current_diagnostics(f, dns, g, c, flow_values)
        type(field_type), intent(inout) :: f
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        type(comm_type), intent(in) :: c
        real(C_DOUBLE), intent(out) :: flow_values(2)

        integer :: i, j, k, nx, ny, nz
        real(C_DOUBLE) :: velocity(3), diagnostics(3)

        nx = int(dns%localSize(1,2))
        ny = int(dns%localSize(2,2))
        nz = int(dns%localSize(3,2))
        diagnostics = 0.0d0

#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do collapse(3) reduction(+:diagnostics) &
        !$omp& map(to: dns, f%q) private(i,j,k,velocity)
#endif
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    call centered_velocity(f, i, j, k, velocity)
                    diagnostics(1) = diagnostics(1) + velocity(1)
                    diagnostics(2) = diagnostics(2) + velocity(3)
                    diagnostics(3) = diagnostics(3) + 1.0d0
                end do
            end do
        end do
#ifdef USE_OPENMP_OFFLOAD
        !$omp end target teams distribute parallel do
#endif

        call comm_allreduce_sum(c, diagnostics)
        flow_values = 0.0d0
        if (diagnostics(3) > 0.0d0) then
            flow_values(1) = diagnostics(1)/diagnostics(3)
            flow_values(2) = diagnostics(2)/diagnostics(3)
        end if
    end subroutine channel_current_diagnostics

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
            if (c%has_terminal) print *, "warning: could not read channel statistics file: ", trim(this%file)
            this%sum = 0.0d0
            this%count = 0.0d0
        else if (c%has_terminal) then
            print *, "continuing channel statistics from: ", trim(this%file)
        end if
    end subroutine read_channel_stats_restart

end module channel_stats
