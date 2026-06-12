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
        ! Wall-normal tables are per refinement level, concatenated: level
        ! l owns rows lvlOff(l)+1 .. lvlOff(l+1) (ny*2^l rows), each block
        ! accumulating into its own level's rows. Level 0 is written to
        ! `file`, level l to the same name with an _l<l> suffix, so
        ! single-level runs are unchanged.
        integer :: nLevels = 1
        integer, allocatable :: lvlOff(:)
        real(C_DOUBLE), allocatable :: sum(:)
        real(C_DOUBLE), allocatable :: count(:)
        real(C_DOUBLE), allocatable :: profile(:)
        real(C_DOUBLE), allocatable :: coord(:)
        real(C_DOUBLE), allocatable :: width(:)
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

    subroutine channel_stats_setup(this, blk, dns, g, c)
        class(channel_stats_type), intent(inout) :: this
        type(block_set_type), intent(in) :: blk
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        type(comm_type), intent(in) :: c

        integer :: total, l, n, nwall, off

        if (this%sample_interval <= 0 .and. this%write_interval <= 0) return

        call system_clock(count=this%clock_start, count_rate=this%clock_rate)
        this%clock_step_start = dns%step_current
        this%runtime_header_written = .false.
        this%last_write_step = -1_C_INT

        this%nLevels = int(blk%nLevels)
        allocate(this%lvlOff(0:this%nLevels))
        this%lvlOff(0) = 0
        do l = 1, this%nLevels
            this%lvlOff(l) = this%lvlOff(l-1) + int(dns%globalSize(2))*2**(l-1)
        end do
        total = this%lvlOff(this%nLevels)

        allocate(this%sum(CHANNEL_NSTAT*total))
        allocate(this%count(total))
        allocate(this%profile(CHANNEL_NSTAT*total))
        allocate(this%coord(total))
        allocate(this%width(total))

        this%sum = 0.0d0
        this%count = 0.0d0
        this%profile = 0.0d0

        ! Per-level y lines come from the block set; without [blocks] nb
        ! (rank-box mode) only level 0 exists and the global line serves.
        do l = 0, this%nLevels - 1
            nwall = int(dns%globalSize(2))*2**l
            off = this%lvlOff(l)
            do n = 1, nwall
                if (allocated(blk%lineY)) then
                    this%coord(off+n) = 0.5d0*(blk%lineY(n-1, l+1) + blk%lineY(n, l+1))
                    this%width(off+n) = blk%lineY(n, l+1) - blk%lineY(n-1, l+1)
                else
                    this%coord(off+n) = 0.5d0*(g%yNode(n-1) + g%yNode(n))
                    this%width(off+n) = g%yNode(n) - g%yNode(n-1)
                end if
            end do
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
            call collect_channel_sample(this, blk, dns, sample_sum, sample_count)
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

    subroutine collect_channel_sample(this, blk, dns, sample_sum, sample_count)
        class(channel_stats_type), intent(in) :: this
        type(block_set_type), intent(inout) :: blk
        type(dns_type), intent(in) :: dns
        real(C_DOUBLE), intent(inout) :: sample_sum(:), sample_count(:)

        integer :: i, j, k, b, s, nx, ny, nz, nBlocks, global_wall_idx, base, row
        integer :: lvlOff(0:this%nLevels)
        real(C_DOUBLE) :: cell_area
        real(C_DOUBLE) :: p, kin, eps, velocity(3), sample(CHANNEL_NSTAT)

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))
        nBlocks = int(blk%nBlocks)
        lvlOff = this%lvlOff
        sample_sum = 0.0d0
        sample_count = 0.0d0

        ! Each block accumulates into its own level's rows; cell areas come
        ! from the block's level-l coordinate slices (bitwise the global
        ! line values at level 0).
#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do collapse(4) &
        !$omp& map(to: dns, lvlOff, blk%origin, blk%level, blk%q, blk%x, blk%z, &
        !$omp& blk%d1x, blk%d1y, blk%d1z) &
        !$omp& map(tofrom: sample_sum, sample_count) &
        !$omp& private(i,j,k,b,s,global_wall_idx,base,row,cell_area,p,kin,eps,velocity,sample)
#endif
        do b = 1, nBlocks
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    global_wall_idx = int(blk%origin(2,b)) + j
                    row = lvlOff(int(blk%level(b))) + global_wall_idx
                    base = CHANNEL_NSTAT*(row - 1)
                    cell_area = (blk%x(i+1,VAR_U,b) - blk%x(i,VAR_U,b)) * &
                                (blk%z(k+1,VAR_W,b) - blk%z(k,VAR_W,b))

                    call centered_velocity(blk, b, i, j, k, velocity)
                    p = blk%q(i,j,k,VAR_P,b)
                    kin = 0.5d0*sum(velocity*velocity)
                    eps = channel_dissipation(blk, dns, b, i, j, k)
                    sample = channel_sample(velocity, p, kin, eps)

                    !$omp atomic update
                    sample_count(row) = sample_count(row) + cell_area
                    do s = 1, CHANNEL_NSTAT
                        !$omp atomic update
                        sample_sum(base + s) = sample_sum(base + s) + cell_area*sample(s)
                    end do
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

        integer :: l, off, base

        if (.not. allocated(this%sum)) return
        if (this%last_write_step == dns%step_current) return

        allocate(raw_sum(size(this%sum)))
        allocate(reduced_count(size(this%count)))
        raw_sum = this%sum
        reduced_count = this%count
        call comm_allreduce_sum(c, raw_sum)
        call comm_allreduce_sum(c, reduced_count)
        if (sum(reduced_count) <= 0.0d0) return

        call build_channel_profile(this, raw_sum, reduced_count, dns)

        ! One file per level (level 0 keeps the configured name), each in
        ! the single-level layout the existing tooling reads.
        if (c%has_terminal) then
            do l = 0, this%nLevels - 1
                off = this%lvlOff(l)
                nwall = this%lvlOff(l+1) - off
                base = CHANNEL_NSTAT*off
                if (sum(reduced_count(off+1:off+nwall)) <= 0.0d0) cycle
                c_file_name = to_c_string(level_file_name(this%file, l))
                ierr = fdm_h5_write_channel_stats(c_file_name, int(nwall, C_INT), int(CHANNEL_NSTAT, C_INT), &
                    dns%step_current, dns%t_current, int(2, C_INT), dns%re, dns%forcing, &
                    this%coord(off+1:off+nwall), this%profile(base+1:base+CHANNEL_NSTAT*nwall), &
                    raw_sum(base+1:base+CHANNEL_NSTAT*nwall), reduced_count(off+1:off+nwall))
                if (ierr /= 0_C_INT) then
                    print *, "error: could not write channel statistics file: ", &
                        trim(level_file_name(this%file, l))
                    error stop
                end if
            end do
        end if

        this%last_write_step = dns%step_current
    end subroutine channel_stats_write_hdf5

    subroutine channel_stats_finalize(this, dns, c)
        class(channel_stats_type), intent(inout) :: this
        type(dns_type), intent(in) :: dns
        type(comm_type), intent(in) :: c

        call this%write_hdf5(dns, c)
    end subroutine channel_stats_finalize

    subroutine centered_velocity(blk, b, i, j, k, velocity)
        type(block_set_type), intent(in) :: blk
        integer, intent(in) :: b, i, j, k
        real(C_DOUBLE), intent(out) :: velocity(3)

        velocity(1) = 0.5d0*(blk%q(i,j,k,VAR_U,b) + blk%q(i+1,j,k,VAR_U,b))
        velocity(2) = 0.5d0*(blk%q(i,j,k,VAR_V,b) + blk%q(i,j+1,k,VAR_V,b))
        velocity(3) = 0.5d0*(blk%q(i,j,k,VAR_W,b) + blk%q(i,j,k+1,VAR_W,b))
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

    real(C_DOUBLE) function channel_dissipation(blk, dns, b, i, j, k) result(eps)
        type(block_set_type), intent(in) :: blk
        type(dns_type), intent(in) :: dns
        integer, intent(in) :: b, i, j, k

        real(C_DOUBLE) :: dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
        real(C_DOUBLE) :: s11, s22, s33, s12, s13, s23

        dudx = (blk%q(i+1,j,k,VAR_U,b) - blk%q(i,j,k,VAR_U,b))*blk%d1x(i,VAR_P,b)
        dvdy = (blk%q(i,j+1,k,VAR_V,b) - blk%q(i,j,k,VAR_V,b))*blk%d1y(j,VAR_P,b)
        dwdz = (blk%q(i,j,k+1,VAR_W,b) - blk%q(i,j,k,VAR_W,b))*blk%d1z(k,VAR_P,b)

        dudy = 0.25d0*((blk%q(i,j+1,k,VAR_U,b) + blk%q(i+1,j+1,k,VAR_U,b)) - &
                       (blk%q(i,j-1,k,VAR_U,b) + blk%q(i+1,j-1,k,VAR_U,b)))*blk%d1y(j,VAR_U,b)
        dudz = 0.25d0*((blk%q(i,j,k+1,VAR_U,b) + blk%q(i+1,j,k+1,VAR_U,b)) - &
                       (blk%q(i,j,k-1,VAR_U,b) + blk%q(i+1,j,k-1,VAR_U,b)))*blk%d1z(k,VAR_U,b)
        dvdx = 0.25d0*((blk%q(i+1,j,k,VAR_V,b) + blk%q(i+1,j+1,k,VAR_V,b)) - &
                       (blk%q(i-1,j,k,VAR_V,b) + blk%q(i-1,j+1,k,VAR_V,b)))*blk%d1x(i,VAR_V,b)
        dvdz = 0.25d0*((blk%q(i,j,k+1,VAR_V,b) + blk%q(i,j+1,k+1,VAR_V,b)) - &
                       (blk%q(i,j,k-1,VAR_V,b) + blk%q(i,j+1,k-1,VAR_V,b)))*blk%d1z(k,VAR_V,b)
        dwdx = 0.25d0*((blk%q(i+1,j,k,VAR_W,b) + blk%q(i+1,j,k+1,VAR_W,b)) - &
                       (blk%q(i-1,j,k,VAR_W,b) + blk%q(i-1,j,k+1,VAR_W,b)))*blk%d1x(i,VAR_W,b)
        dwdy = 0.25d0*((blk%q(i,j+1,k,VAR_W,b) + blk%q(i,j+1,k+1,VAR_W,b)) - &
                       (blk%q(i,j-1,k,VAR_W,b) + blk%q(i,j-1,k+1,VAR_W,b)))*blk%d1y(j,VAR_W,b)

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

        integer :: n, nwall, base, mean_base, l, off
        real(C_DOUBLE), allocatable :: mean_velocity(:)
        real(C_DOUBLE) :: mean_eps

        ! Per level: the wall-normal finite differences below must not mix
        ! rows of different levels.
        do l = 0, this%nLevels - 1
            off = this%lvlOff(l)
            nwall = this%lvlOff(l+1) - off
            if (nwall < 2) cycle

            if (allocated(mean_velocity)) deallocate(mean_velocity)
            allocate(mean_velocity(3*nwall))
            mean_velocity = 0.0d0
            do n = 1, nwall
                base = CHANNEL_NSTAT*(off + n - 1)
                mean_base = 3*(n - 1)
                mean_velocity(mean_base+1) = this%profile(base+STAT_U)
                mean_velocity(mean_base+2) = this%profile(base+STAT_V)
                mean_velocity(mean_base+3) = this%profile(base+STAT_W)
            end do

            do n = 1, nwall
                if (reduced_count(off+n) <= 0.0d0) cycle
                base = CHANNEL_NSTAT*(off + n - 1)
                mean_eps = mean_profile_dissipation_at(mean_velocity, &
                    this%coord(off+1:off+nwall), reduced_count(off+1:off+nwall), n, dns%re)
                this%profile(base+STAT_EPSILON) = max(0.0d0, this%profile(base+STAT_EPSILON) - mean_eps)
            end do
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

        integer :: nwall, base, mean_base, n, l, off
        real(C_DOUBLE) :: dy, total_volume, mean_eps, k_turb, eps_turb
        real(C_DOUBLE) :: flow_integral(2), volume_integral(2)
        real(C_DOUBLE), allocatable :: mean_velocity(:)

        flow_values = 0.0d0
        volume_values = 0.0d0
        flow_integral = 0.0d0
        volume_integral = 0.0d0
        total_volume = 0.0d0

        ! Levels tile the volume (each cell is counted at its block's
        ! level), so the volume integrals just accumulate over all levels.
        do l = 0, this%nLevels - 1
            off = this%lvlOff(l)
            nwall = this%lvlOff(l+1) - off

            if (allocated(mean_velocity)) deallocate(mean_velocity)
            allocate(mean_velocity(3*nwall))
            mean_velocity = 0.0d0

            do n = 1, nwall
                if (sample_count(off+n) <= 0.0d0) cycle
                base = CHANNEL_NSTAT*(off + n - 1)
                mean_base = 3*(n - 1)
                mean_velocity(mean_base+1) = sample_sum(base+STAT_U)/sample_count(off+n)
                mean_velocity(mean_base+2) = sample_sum(base+STAT_V)/sample_count(off+n)
                mean_velocity(mean_base+3) = sample_sum(base+STAT_W)/sample_count(off+n)
            end do

            do n = 1, nwall
                if (sample_count(off+n) <= 0.0d0) cycle
                base = CHANNEL_NSTAT*(off + n - 1)
                mean_base = 3*(n - 1)
                dy = this%width(off+n)
                total_volume = total_volume + dy*sample_count(off+n)

                flow_integral(1) = flow_integral(1) + dy*sample_sum(base+STAT_U)
                flow_integral(2) = flow_integral(2) + dy*sample_sum(base+STAT_W)

                k_turb = max(0.0d0, sample_sum(base+STAT_K)/sample_count(off+n) - &
                    0.5d0*(mean_velocity(mean_base+1)**2 + mean_velocity(mean_base+2)**2 + &
                    mean_velocity(mean_base+3)**2))
                mean_eps = mean_profile_dissipation_at(mean_velocity, &
                    this%coord(off+1:off+nwall), sample_count(off+1:off+nwall), n, dns%re)
                eps_turb = max(0.0d0, sample_sum(base+STAT_EPSILON)/sample_count(off+n) - mean_eps)

                volume_integral(1) = volume_integral(1) + dy*sample_count(off+n)*k_turb
                volume_integral(2) = volume_integral(2) + dy*sample_count(off+n)*eps_turb
            end do
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

        integer :: l, off, nwall, base

        if (c%world_rank /= 0) return

        do l = 0, this%nLevels - 1
            off = this%lvlOff(l)
            nwall = this%lvlOff(l+1) - off
            base = CHANNEL_NSTAT*off
            inquire(file=trim(level_file_name(this%file, l)), exist=exists)
            if (.not. exists) cycle

            restart_step = 0_C_INT
            restart_time = 0.0d0
            c_file_name = to_c_string(level_file_name(this%file, l))
            ierr = fdm_h5_read_channel_stats(c_file_name, int(nwall, C_INT), &
                int(CHANNEL_NSTAT, C_INT), restart_step, restart_time, &
                this%sum(base+1:base+CHANNEL_NSTAT*nwall), this%count(off+1:off+nwall))
            if (ierr /= 0_C_INT) then
                if (c%has_terminal) print *, "error: could not read channel statistics file: ", &
                    trim(level_file_name(this%file, l))
                error stop
            else if (c%has_terminal) then
                print *, "continuing channel statistics from: ", trim(level_file_name(this%file, l))
            end if
        end do
    end subroutine read_channel_stats_restart

    ! channel_stats.h5 -> channel_stats.h5 (level 0), channel_stats_l1.h5, ...
    function level_file_name(file, l) result(name)
        character(len=*), intent(in) :: file
        integer, intent(in) :: l
        character(len=300) :: name

        integer :: dot
        character(len=8) :: tag

        if (l == 0) then
            name = file
            return
        end if
        write(tag, '("_l",I0)') l
        dot = index(file, ".", back=.true.)
        if (dot > 0) then
            name = file(1:dot-1)//trim(tag)//file(dot:)
        else
            name = trim(file)//trim(tag)
        end if
    end function level_file_name

end module channel_stats
