!--------------------------------!
!                                !
!   Passive-scalar statistics    !
!                                !
!--------------------------------!
!
! Increment S4 of docs/next_session_scalar.md: the in-solver form of the
! statistics validation/scalar/check_scalar_turb.py computes in post
! (mean profile, rms, turbulent flux, wall flux -> Nusselt), plus the
! immersed body's heat release.
!
! WHY THIS IS A SOLVER-LEVEL FACILITY and not a case component like
! channel_stats/bl_stats: the case `after_step` interface carries neither the
! scalar configuration (Pr, Pr_t) nor turb%nut, and the same statistics have
! to serve the channel, the boundary layer and body cases (the heated
! cylinder runs the generic case). So the accumulator lives in moby_solve.f90
! and is driven by [scalar] keys. Everything here is gated on
! `stats_sample_interval`/`stats_write_interval`/`heat_interval` (all off by
! default) AND on scalars being configured, so a run without them calls no
! kernel at all -- the by-construction bit-exactness argument S0-S3 use.
!
! TWO layouts, one kernel ([scalar] stats_layout):
!   profile  wall-normal rows, x-z averaged, per refinement level -- the
!            channel form (channel_stats' lvlOff tables verbatim);
!   plane    rows on the global (x,y) plane, z averaged -- the boundary-layer
!            form (bl_stats' flattening verbatim), single level.
! The files are written with the channel_stats / bl_stats HDF5 writers, so
! there is no new C code and the existing readers work: `nstat` is
! SCALAR_NSTAT*nScalar and stat s of scalar is column
! SCALAR_NSTAT*(is-1) + s. tools/scalar_stats.py does the reading.
!
! The face fluxes are the TRANSPORT KERNEL's own expressions (same face
! diffusivity 1/(Re Pr) + nut/Pr_t(face), same FACE_CLOSED / adiabatic-body
! masks), which is what makes the wall row an exact discrete wall flux and
! therefore an exact theta_tau / Nusselt number.
!
! BODY HEAT (`heat_interval`): a Dirichlet body's heat release CANNOT be
! measured as int coef_p (s_body - s) dV -- a solid cell holds the body value
! to the LAST BIT, so the product is 1e28 x 0 = 0 and ~37 % of the heat is
! invisible (the S3 FINDING, docs/next_session_scalar.md). What is written is
! the cancellation-free pair check_scalar_ibm.py `surface` uses: the flux
! across every staircase face separating a solid cell from a fluid one, plus
! the penalization delivered into the GRADED fluid cells.
module scalar_stats
    use, intrinsic :: iso_c_binding
    use :: init, only: dns_type, grid_type, VAR_U, VAR_V, VAR_W, VAR_P, VAR_S0
    use :: blocks, only: block_set_type, FACE_CLOSED
    use :: comm, only: comm_type, comm_allreduce_sum
    use :: io, only: to_c_string
    use :: turbulence, only: turb_type, turbulence_is_enabled
    use :: ibmm, only: ibm_type
    use :: scalar, only: scalar_type, scalars_enabled, eddy_diffusivity, &
        wall_face_diffusivity, &
        SC_IBM_ADIABATIC, SC_LAYOUT_PROFILE, SC_LAYOUT_PLANE
    implicit none

    private

    ! Per-scalar statistics, in the order they are stored in a row.
    integer, parameter :: SSTAT_S   = 1     ! <s>
    integer, parameter :: SSTAT_SS  = 2     ! <s^2>            (rms = sqrt(SS - S^2))
    integer, parameter :: SSTAT_US  = 3     ! <u_c s>          (streamwise flux)
    integer, parameter :: SSTAT_CLO = 4     ! <v s>|face       , cell's LOW  y face
    integer, parameter :: SSTAT_JLO = 5     ! <v s - D ds/dy>| , cell's LOW  y face
    integer, parameter :: SSTAT_CHI = 6     ! the same two on the HIGH y face; a
    integer, parameter :: SSTAT_JHI = 7     ! wall row then carries the wall flux
    integer, parameter :: SCALAR_NSTAT = 7

    ! The solid-coefficient test the transport kernel uses (SOLID/Re = 1e30/Re).
    real(C_DOUBLE), parameter :: SOLID_FACE_THRESHOLD = 1.0d20

    type, public :: scalar_stats_type
        integer :: layout = SC_LAYOUT_PROFILE
        integer :: sample_interval = -1
        integer :: write_interval = -1
        integer :: heat_interval = -1
        character(len=256) :: file = "scalar_stats.h5"
        character(len=256) :: heat_file = "scalar_heat.txt"
        logical :: heat_header_written = .false.
        integer(C_INT) :: last_write_step = -1_C_INT
        integer :: nScal = 0
        integer :: nStat = 0                       ! SCALAR_NSTAT*nScal
        ! profile layout: wall-normal rows per level, concatenated (level l
        ! owns rows lvlOff(l)+1 .. lvlOff(l+1)), the channel_stats layout.
        integer :: nLevels = 1
        integer, allocatable :: lvlOff(:)
        ! plane layout: the global (x,y) plane flattened y-fastest.
        integer :: nx = 0, ny = 0
        real(C_DOUBLE), allocatable :: sum(:), count(:), profile(:)
        real(C_DOUBLE), allocatable :: coord(:), xcoord(:)
    end type scalar_stats_type

    public :: scalar_stats_setup, scalar_stats_after_step, scalar_stats_finalize
    public :: SCALAR_NSTAT

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

    !--------------------------------------------------------------------
    ! Setup
    !--------------------------------------------------------------------

    subroutine scalar_stats_setup(this, sc, blk, dns, g, c)
        type(scalar_stats_type), intent(inout) :: this
        type(scalar_type), intent(in) :: sc
        type(block_set_type), intent(in) :: blk
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        type(comm_type), intent(in) :: c

        integer :: total, l, n, nwall, off

        ! The intervals come from [scalar]; nothing here runs without them.
        this%layout = int(sc%statsLayout)
        this%sample_interval = int(sc%statsSample)
        this%write_interval = int(sc%statsWrite)
        this%heat_interval = int(sc%heatInterval)
        this%file = sc%statsFile
        this%heat_file = sc%heatFile
        this%heat_header_written = .false.
        this%last_write_step = -1_C_INT

        if (.not. scalars_enabled(sc)) then
            this%sample_interval = -1
            this%write_interval = -1
            this%heat_interval = -1
            return
        end if
        if (this%heat_interval > 0 .and. .not. dns%ibm_enabled) then
            if (c%has_terminal) print *, "warning: [scalar] heat_interval needs an immersed body; ignored"
            this%heat_interval = -1
        end if
        if (this%sample_interval <= 0 .and. this%write_interval <= 0) return

        this%nScal = int(sc%n)
        this%nStat = SCALAR_NSTAT*this%nScal

        if (this%layout == SC_LAYOUT_PLANE) then
            if (blk%nLevels > 1_C_INT .and. c%has_terminal) &
                print *, "warning: [scalar] stats_layout = plane assumes a single level; refinement ignored"
            this%nx = int(dns%globalSize(1))
            this%ny = int(dns%globalSize(2))
            total = this%nx*this%ny
            allocate(this%xcoord(this%nx), this%coord(this%ny))
            do n = 1, this%nx
                this%xcoord(n) = 0.5d0*(g%xNode(n-1) + g%xNode(n))
            end do
            do n = 1, this%ny
                this%coord(n) = 0.5d0*(g%yNode(n-1) + g%yNode(n))
            end do
        else
            ! Wall-normal rows, one table per level (the channel_stats layout).
            ! The row count per level is the PER-DIRECTION scaling
            ! 2**(l*refMask(2)) (blocks.f90): under [blocks] refine_dims = xz
            ! the y line is the same at every level, and sizing the tables as
            ! if it doubled would index blk%lineY out of bounds.
            this%nLevels = int(blk%nLevels)
            allocate(this%lvlOff(0:this%nLevels))
            this%lvlOff(0) = 0
            do l = 1, this%nLevels
                this%lvlOff(l) = this%lvlOff(l-1) + level_rows(dns, blk, l-1)
            end do
            total = this%lvlOff(this%nLevels)
            allocate(this%coord(total))
            do l = 0, this%nLevels - 1
                nwall = level_rows(dns, blk, l)
                off = this%lvlOff(l)
                do n = 1, nwall
                    if (allocated(blk%lineY)) then
                        this%coord(off+n) = 0.5d0*(blk%lineY(n-1, l+1) + blk%lineY(n, l+1))
                    else
                        this%coord(off+n) = 0.5d0*(g%yNode(n-1) + g%yNode(n))
                    end if
                end do
            end do
        end if

        allocate(this%sum(this%nStat*total), this%count(total), this%profile(this%nStat*total))
        this%sum = 0.0d0
        this%count = 0.0d0
        this%profile = 0.0d0

        if (len_trim(dns%restart_file) > 0) call read_stats_restart(this, c)
    end subroutine scalar_stats_setup

    !--------------------------------------------------------------------
    ! Per-step entry point
    !--------------------------------------------------------------------

    subroutine scalar_stats_after_step(this, sc, blk, dns, turb, ibm, c)
        type(scalar_stats_type), intent(inout) :: this
        type(scalar_type), intent(in) :: sc
        type(block_set_type), intent(inout) :: blk
        type(dns_type), intent(in) :: dns
        type(turb_type), intent(in) :: turb
        type(ibm_type), intent(in) :: ibm
        type(comm_type), intent(in) :: c

        real(C_DOUBLE), allocatable :: sample_sum(:), sample_count(:)

        if (this%heat_interval > 0) then
            if (interval_is_due(dns%step_current, this%heat_interval)) &
                call sample_body_heat(this, sc, blk, dns, turb, ibm, c)
        end if

        if (.not. allocated(this%sum)) return

        if (interval_is_due(dns%step_current, this%sample_interval)) then
            allocate(sample_sum(size(this%sum)), sample_count(size(this%count)))
            call collect_scalar_sample(this, sc, blk, dns, turb, ibm, sample_sum, sample_count)
            this%sum = this%sum + sample_sum
            this%count = this%count + sample_count
        end if

        if (interval_is_due(dns%step_current, this%write_interval)) &
            call write_stats(this, dns, c)
    end subroutine scalar_stats_after_step

    subroutine scalar_stats_finalize(this, dns, c)
        type(scalar_stats_type), intent(inout) :: this
        type(dns_type), intent(in) :: dns
        type(comm_type), intent(in) :: c

        call write_stats(this, dns, c)
    end subroutine scalar_stats_finalize

    ! Wall-normal rows of refinement level l.
    integer function level_rows(dns, blk, l) result(n)
        type(dns_type), intent(in) :: dns
        type(block_set_type), intent(in) :: blk
        integer, intent(in) :: l

        n = int(dns%globalSize(2))*2**(l*int(blk%refMask(2)))
    end function level_rows

    logical function interval_is_due(step, interval) result(is_due)
        integer(C_INT), intent(in) :: step
        integer, intent(in) :: interval

        is_due = .false.
        if (interval <= 0) return
        is_due = modulo(int(step), interval) == 0
    end function interval_is_due

    !--------------------------------------------------------------------
    ! Sampling
    !--------------------------------------------------------------------

    ! One sample of every scalar into the row accumulators. The nut argument
    ! follows the transport kernel: with no turbulence model turb%nut does not
    ! exist, so the 1-cell dummy is mapped and the eddy term switched off.
    subroutine collect_scalar_sample(this, sc, blk, dns, turb, ibm, sample_sum, sample_count)
        type(scalar_stats_type), intent(in) :: this
        type(scalar_type), intent(in) :: sc
        type(block_set_type), intent(inout) :: blk
        type(dns_type), intent(in) :: dns
        type(turb_type), intent(in) :: turb
        type(ibm_type), intent(in) :: ibm
        real(C_DOUBLE), intent(inout) :: sample_sum(:), sample_count(:)

        if (turbulence_is_enabled(turb) .and. allocated(turb%nut)) then
            call scalar_sample_kernel(this, sc, blk, dns, turb%nut, .true., ibm%coef, &
                sample_sum, sample_count)
        else
            call scalar_sample_kernel(this, sc, blk, dns, sc%nutNone, .false., ibm%coef, &
                sample_sum, sample_count)
        end if
    end subroutine collect_scalar_sample

    subroutine scalar_sample_kernel(this, sc, blk, dns, nut, useNut, coef, sample_sum, sample_count)
        type(scalar_stats_type), intent(in) :: this
        type(scalar_type), intent(in) :: sc
        type(block_set_type), intent(inout) :: blk
        type(dns_type), intent(in) :: dns
        real(C_DOUBLE), intent(in) :: nut(0:,0:,0:,1:)
        logical, intent(in) :: useNut
        real(C_DOUBLE), intent(in) :: coef(0:,0:,0:,1:,1:)
        real(C_DOUBLE), intent(inout) :: sample_sum(:), sample_count(:)

        integer :: i, j, k, b, is, s, nx, ny, nz, nBlocks, nScal, var, row, base, gx, gy
        integer :: nStat, plane, ny_g
        integer :: lvlOff(0:this%nLevels)
        real(C_DOUBLE) :: weight, ire, re, dm, uc, s0, sS, sN, vs, vn, dys, dyn
        real(C_DOUBLE) :: nts, ntn, clo, jlo, chi, jhi
        logical :: useIbm, adiab, cls, cln, sols, soln, ms, mn, wallfn

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))
        nBlocks = int(blk%nBlocks)
        nScal = this%nScal
        nStat = this%nStat
        ny_g = this%ny
        plane = merge(1, 0, this%layout == SC_LAYOUT_PLANE)
        lvlOff = 0
        if (allocated(this%lvlOff)) lvlOff = this%lvlOff
        re = dns%re
        ire = 1.0d0/re
        useIbm = logical(dns%ibm_enabled)
        ! S5a: under wall functions the wall row's flux is the THERMAL WALL
        ! FUNCTION's, not nu_t/Pr_t's -- the statistics must report the flux
        ! the transport kernel actually applied, or theta_tau is not the
        ! solver's theta_tau. Same branch, same helper as scalar.f90.
        wallfn = dns%rans_wall_treatment == 1_C_INT
        sample_sum = 0.0d0
        sample_count = 0.0d0

        !$omp target teams distribute parallel do collapse(4) &
        !$omp& map(to: ire, re, useNut, useIbm, nScal, nStat, nx, ny, nz, plane, ny_g, lvlOff, &
        !$omp& blk%q, blk%x, blk%z, blk%origin, blk%level, blk%physLow, blk%physHigh, nut, coef, &
        !$omp& sc%pr, sc%prt, sc%prtModel, sc%invDy, sc%ibmMode, &
        !$omp& sc%wfP, sc%wfYpt, sc%wfYplus, wallfn) &
        !$omp& map(tofrom: sample_sum, sample_count) &
        !$omp& private(i,j,k,b,is,s,var,row,base,gx,gy,weight,dm,uc,s0,sS,sN,vs,vn, &
        !$omp& dys,dyn,nts,ntn,clo,jlo,chi,jhi,adiab,cls,cln,sols,soln,ms,mn)
        do b = 1, nBlocks
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    if (plane == 1) then
                        gx = int(blk%origin(1,b)) + i
                        gy = int(blk%origin(2,b)) + j
                        row = (gx - 1)*ny_g + gy
                        weight = blk%z(k+1,VAR_W,b) - blk%z(k,VAR_W,b)
                    else
                        row = lvlOff(int(blk%level(b))) + int(blk%origin(2,b)) + j
                        weight = (blk%x(i+1,VAR_U,b) - blk%x(i,VAR_U,b))* &
                                 (blk%z(k+1,VAR_W,b) - blk%z(k,VAR_W,b))
                    end if
                    base = nStat*(row - 1)

                    uc = 0.5d0*(blk%q(i,j,k,VAR_U,b) + blk%q(i+1,j,k,VAR_U,b))
                    vs = blk%q(i,j,  k,VAR_V,b)
                    vn = blk%q(i,j+1,k,VAR_V,b)

                    ! The transport kernel's y-face masks, verbatim.
                    cls = j == 1  .and. blk%physLow(2,b)  == FACE_CLOSED
                    cln = j == ny .and. blk%physHigh(2,b) == FACE_CLOSED
                    sols = .false.
                    soln = .false.
                    if (useIbm) then
                        sols = abs(coef(i,j,  k,VAR_V,b)) > SOLID_FACE_THRESHOLD
                        soln = abs(coef(i,j+1,k,VAR_V,b)) > SOLID_FACE_THRESHOLD
                    end if
                    nts = 0.0d0
                    ntn = 0.0d0
                    if (useNut) then
                        nts = 0.5d0*(nut(i,j-1,k,b) + nut(i,j,k,b))
                        ntn = 0.5d0*(nut(i,j,k,b) + nut(i,j+1,k,b))
                    end if

                    !$omp atomic update
                    sample_count(row) = sample_count(row) + weight

                    do is = 1, nScal
                        var = VAR_S0 + is
                        s0 = blk%q(i,j,  k,var,b)
                        sS = blk%q(i,j-1,k,var,b)
                        sN = blk%q(i,j+1,k,var,b)

                        adiab = useIbm .and. sc%ibmMode(is) == SC_IBM_ADIABATIC
                        ms = cls .or. (adiab .and. sols)
                        mn = cln .or. (adiab .and. soln)

                        dm = ire/sc%pr(is)
                        dys = dm
                        dyn = dm
                        if (useNut .and. wallfn) then
                            dys = dm + wall_face_diffusivity(nut(i,j-1,k,b), nut(i,j,k,b), &
                                sc%wfYplus(i,j-1,k,b), sc%wfYplus(i,j,k,b), &
                                sc%pr(is), sc%prt(is), sc%prtModel(is), re, &
                                sc%wfP(is), sc%wfYpt(is))
                            dyn = dm + wall_face_diffusivity(nut(i,j,k,b), nut(i,j+1,k,b), &
                                sc%wfYplus(i,j,k,b), sc%wfYplus(i,j+1,k,b), &
                                sc%pr(is), sc%prt(is), sc%prtModel(is), re, &
                                sc%wfP(is), sc%wfYpt(is))
                        else if (useNut) then
                            dys = dm + eddy_diffusivity(nts, sc%pr(is), sc%prt(is), &
                                sc%prtModel(is), re)
                            dyn = dm + eddy_diffusivity(ntn, sc%pr(is), sc%prt(is), &
                                sc%prtModel(is), re)
                        end if

                        ! Fluxes in the +y sense: convection minus D ds/dy.
                        clo = merge(0.0d0, vs*0.5d0*(sS + s0), ms)
                        jlo = clo - merge(0.0d0, dys*(s0 - sS)*sc%invDy(j,b), ms)
                        chi = merge(0.0d0, vn*0.5d0*(s0 + sN), mn)
                        jhi = chi - merge(0.0d0, dyn*(sN - s0)*sc%invDy(j+1,b), mn)

                        s = base + SCALAR_NSTAT*(is - 1)
                        !$omp atomic update
                        sample_sum(s+SSTAT_S) = sample_sum(s+SSTAT_S) + weight*s0
                        !$omp atomic update
                        sample_sum(s+SSTAT_SS) = sample_sum(s+SSTAT_SS) + weight*s0*s0
                        !$omp atomic update
                        sample_sum(s+SSTAT_US) = sample_sum(s+SSTAT_US) + weight*uc*s0
                        !$omp atomic update
                        sample_sum(s+SSTAT_CLO) = sample_sum(s+SSTAT_CLO) + weight*clo
                        !$omp atomic update
                        sample_sum(s+SSTAT_JLO) = sample_sum(s+SSTAT_JLO) + weight*jlo
                        !$omp atomic update
                        sample_sum(s+SSTAT_CHI) = sample_sum(s+SSTAT_CHI) + weight*chi
                        !$omp atomic update
                        sample_sum(s+SSTAT_JHI) = sample_sum(s+SSTAT_JHI) + weight*jhi
                    end do
                end do
            end do
        end do
        end do
        !$omp end target teams distribute parallel do
    end subroutine scalar_sample_kernel

    !--------------------------------------------------------------------
    ! Output
    !--------------------------------------------------------------------

    subroutine write_stats(this, dns, c)
        type(scalar_stats_type), intent(inout) :: this
        type(dns_type), intent(in) :: dns
        type(comm_type), intent(in) :: c

        character(kind=C_CHAR,len=:), allocatable :: c_file_name
        real(C_DOUBLE), allocatable :: raw_sum(:), reduced_count(:)
        integer(C_INT) :: ierr
        integer :: n, s, l, off, nwall, base

        if (.not. allocated(this%sum)) return
        if (this%last_write_step == dns%step_current) return

        allocate(raw_sum(size(this%sum)), reduced_count(size(this%count)))
        raw_sum = this%sum
        reduced_count = this%count
        call comm_allreduce_sum(c, raw_sum)
        call comm_allreduce_sum(c, reduced_count)
        if (sum(reduced_count) <= 0.0d0) return

        this%profile = raw_sum
        do n = 1, size(reduced_count)
            if (reduced_count(n) <= 0.0d0) cycle
            base = this%nStat*(n - 1)
            do s = 1, this%nStat
                this%profile(base+s) = this%profile(base+s)/reduced_count(n)
            end do
        end do

        if (.not. c%has_terminal) then
            this%last_write_step = dns%step_current
            return
        end if

        if (this%layout == SC_LAYOUT_PLANE) then
            c_file_name = to_c_string(trim(this%file))
            ierr = fdm_h5_write_bl_stats(c_file_name, int(this%nx, C_INT), int(this%ny, C_INT), &
                int(this%nStat, C_INT), dns%step_current, dns%t_current, dns%re, &
                this%xcoord, this%coord, this%profile, raw_sum, reduced_count)
            if (ierr /= 0_C_INT) then
                print *, "error: could not write scalar statistics file: ", trim(this%file)
                error stop
            end if
        else
            ! One file per level, level 0 keeping the configured name.
            do l = 0, this%nLevels - 1
                off = this%lvlOff(l)
                nwall = this%lvlOff(l+1) - off
                base = this%nStat*off
                if (sum(reduced_count(off+1:off+nwall)) <= 0.0d0) cycle
                c_file_name = to_c_string(level_file_name(this%file, l))
                ierr = fdm_h5_write_channel_stats(c_file_name, int(nwall, C_INT), &
                    int(this%nStat, C_INT), dns%step_current, dns%t_current, int(2, C_INT), &
                    dns%re, dns%forcing, this%coord(off+1:off+nwall), &
                    this%profile(base+1:base+this%nStat*nwall), &
                    raw_sum(base+1:base+this%nStat*nwall), reduced_count(off+1:off+nwall))
                if (ierr /= 0_C_INT) then
                    print *, "error: could not write scalar statistics file: ", &
                        trim(level_file_name(this%file, l))
                    error stop
                end if
            end do
        end if

        this%last_write_step = dns%step_current
    end subroutine write_stats

    ! Continue the accumulators from a previous run's files (channel_stats'
    ! recipe: rank 0 reads, the sums are already rank-reduced on disk, so the
    ! other ranks stay at zero and the next allreduce is still correct).
    subroutine read_stats_restart(this, c)
        type(scalar_stats_type), intent(inout) :: this
        type(comm_type), intent(in) :: c

        character(kind=C_CHAR,len=:), allocatable :: c_file_name
        integer(C_INT) :: ierr, restart_step
        real(C_DOUBLE) :: restart_time
        logical :: exists
        integer :: l, off, nwall, base

        if (c%world_rank /= 0) return

        if (this%layout == SC_LAYOUT_PLANE) then
            inquire(file=trim(this%file), exist=exists)
            if (.not. exists) return
            restart_step = 0_C_INT
            restart_time = 0.0d0
            c_file_name = to_c_string(trim(this%file))
            ierr = fdm_h5_read_bl_stats(c_file_name, int(this%nx, C_INT), int(this%ny, C_INT), &
                int(this%nStat, C_INT), restart_step, restart_time, this%sum, this%count)
            call report_restart(ierr, trim(this%file), c)
            return
        end if

        do l = 0, this%nLevels - 1
            off = this%lvlOff(l)
            nwall = this%lvlOff(l+1) - off
            base = this%nStat*off
            inquire(file=trim(level_file_name(this%file, l)), exist=exists)
            if (.not. exists) cycle
            restart_step = 0_C_INT
            restart_time = 0.0d0
            c_file_name = to_c_string(level_file_name(this%file, l))
            ierr = fdm_h5_read_channel_stats(c_file_name, int(nwall, C_INT), &
                int(this%nStat, C_INT), restart_step, restart_time, &
                this%sum(base+1:base+this%nStat*nwall), this%count(off+1:off+nwall))
            call report_restart(ierr, level_file_name(this%file, l), c)
        end do
    end subroutine read_stats_restart

    subroutine report_restart(ierr, name, c)
        integer(C_INT), intent(in) :: ierr
        character(len=*), intent(in) :: name
        type(comm_type), intent(in) :: c

        if (ierr /= 0_C_INT) then
            if (c%has_terminal) print *, "error: could not read scalar statistics file: ", trim(name)
            error stop
        else if (c%has_terminal) then
            print *, "continuing scalar statistics from: ", trim(name)
        end if
    end subroutine report_restart

    ! scalar_stats.h5 -> scalar_stats.h5 (level 0), scalar_stats_l1.h5, ...
    ! (channel_stats.f90 builds its per-level names the same way; the twelve
    ! lines are repeated rather than shared because this module deliberately
    ! does not depend on the flow-case layer.)
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

    !--------------------------------------------------------------------
    ! Immersed-body heat release (the cancellation-free form)
    !--------------------------------------------------------------------

    subroutine sample_body_heat(this, sc, blk, dns, turb, ibm, c)
        type(scalar_stats_type), intent(inout) :: this
        type(scalar_type), intent(in) :: sc
        type(block_set_type), intent(inout) :: blk
        type(dns_type), intent(in) :: dns
        type(turb_type), intent(in) :: turb
        type(ibm_type), intent(in) :: ibm
        type(comm_type), intent(in) :: c

        real(C_DOUBLE), allocatable :: heat(:)
        integer :: is, unit, stat, nScal
        character(len=32) :: tag

        nScal = int(sc%n)
        allocate(heat(2*nScal))
        if (turbulence_is_enabled(turb) .and. allocated(turb%nut)) then
            call body_heat_kernel(sc, blk, dns, turb%nut, .true., ibm%coef, heat)
        else
            call body_heat_kernel(sc, blk, dns, sc%nutNone, .false., ibm%coef, heat)
        end if
        call comm_allreduce_sum(c, heat)

        if (.not. c%has_terminal) return

        if (.not. this%heat_header_written) then
            open(newunit=unit, file=trim(this%heat_file), status="replace", action="write", iostat=stat)
            if (stat /= 0) then
                print *, "warning: could not open scalar heat file: ", trim(this%heat_file)
                return
            end if
            write(unit,'(A)', advance="no") "# step time"
            do is = 1, nScal
                write(tag,'(A)') trim(sc%name(is))
                write(unit,'(1X,A)', advance="no") &
                    trim(tag)//"_staircase "//trim(tag)//"_graded "//trim(tag)//"_total"
            end do
            write(unit,'(A)') ""
            close(unit)
            this%heat_header_written = .true.
        end if

        open(newunit=unit, file=trim(this%heat_file), status="old", position="append", &
            action="write", iostat=stat)
        if (stat /= 0) then
            print *, "warning: could not append scalar heat file: ", trim(this%heat_file)
            return
        end if
        ! Full double precision: this file is compared with an independent
        ! Python transcription of the same integral (validation/scalar), and
        ! ES16.8 would put a 1e-9 floor under that comparison.
        write(unit,'(I10,1X,ES24.16)', advance="no") int(dns%step_current), dns%t_current
        do is = 1, nScal
            write(unit,'(3(1X,ES24.16))', advance="no") heat(2*is-1), heat(2*is), &
                heat(2*is-1) + heat(2*is)
        end do
        write(unit,'(A)') ""
        close(unit)
    end subroutine sample_body_heat

    ! Per scalar: heat(2*is-1) = the flux across every staircase face
    ! separating a solid cell from a fluid one, heat(2*is) = the penalization
    ! delivered into the NON-SOLID cells (a solid cell's own share is the
    ! staircase term -- see the guard below, and the S3 FINDING for why the
    ! split exists at all).
    !
    ! Each interior face is visited ONCE, as the LOW face of the cell that
    ! owns it, so blocks and ranks never double count. A FACE_CLOSED face
    ! borders a block removed inside the body: both sides are then solid (that
    ! is the removal criterion), so skipping it loses no flux.
    subroutine body_heat_kernel(sc, blk, dns, nut, useNut, coef, heat)
        type(scalar_type), intent(in) :: sc
        type(block_set_type), intent(inout) :: blk
        type(dns_type), intent(in) :: dns
        real(C_DOUBLE), intent(in) :: nut(0:,0:,0:,1:)
        logical, intent(in) :: useNut
        real(C_DOUBLE), intent(in) :: coef(0:,0:,0:,1:,1:)
        real(C_DOUBLE), intent(inout) :: heat(:)

        integer :: i, j, k, b, is, nx, ny, nz, nBlocks, nScal, var
        real(C_DOUBLE) :: ire, re, dm, dx, dy, dz, s0, sW, sS, sB, flux, sgn
        real(C_DOUBLE) :: ntw, nts, ntb, dxw, dys, dzb, cp
        logical :: solc, solw, sols, solb, clw, cls, clb, wallfn

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))
        nBlocks = int(blk%nBlocks)
        nScal = int(sc%n)
        re = dns%re
        ire = 1.0d0/re
        wallfn = dns%rans_wall_treatment == 1_C_INT      ! S5a, as above
        heat = 0.0d0

        !$omp target teams distribute parallel do collapse(4) &
        !$omp& map(to: ire, re, useNut, nScal, nx, ny, nz, &
        !$omp& blk%q, blk%x, blk%y, blk%z, blk%physLow, nut, coef, &
        !$omp& sc%pr, sc%prt, sc%prtModel, sc%invDx, sc%invDy, sc%invDz, &
        !$omp& sc%ibmValue, sc%ibmMode, sc%wfP, sc%wfYpt, sc%wfYplus, wallfn) &
        !$omp& map(tofrom: heat) &
        !$omp& private(i,j,k,b,is,var,dm,dx,dy,dz,s0,sW,sS,sB,flux,sgn,cp, &
        !$omp& ntw,nts,ntb,dxw,dys,dzb,solc,solw,sols,solb,clw,cls,clb)
        do b = 1, nBlocks
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    dx = blk%x(i+1,VAR_U,b) - blk%x(i,VAR_U,b)
                    dy = blk%y(j+1,VAR_V,b) - blk%y(j,VAR_V,b)
                    dz = blk%z(k+1,VAR_W,b) - blk%z(k,VAR_W,b)

                    cp = coef(i,j,k,VAR_P,b)
                    solc = abs(cp) > SOLID_FACE_THRESHOLD
                    solw = abs(coef(i-1,j,k,VAR_P,b)) > SOLID_FACE_THRESHOLD
                    sols = abs(coef(i,j-1,k,VAR_P,b)) > SOLID_FACE_THRESHOLD
                    solb = abs(coef(i,j,k-1,VAR_P,b)) > SOLID_FACE_THRESHOLD
                    clw = i == 1 .and. blk%physLow(1,b) == FACE_CLOSED
                    cls = j == 1 .and. blk%physLow(2,b) == FACE_CLOSED
                    clb = k == 1 .and. blk%physLow(3,b) == FACE_CLOSED

                    ntw = 0.0d0
                    nts = 0.0d0
                    ntb = 0.0d0
                    if (useNut) then
                        ntw = 0.5d0*(nut(i-1,j,k,b) + nut(i,j,k,b))
                        nts = 0.5d0*(nut(i,j-1,k,b) + nut(i,j,k,b))
                        ntb = 0.5d0*(nut(i,j,k-1,b) + nut(i,j,k,b))
                    end if

                    do is = 1, nScal
                        ! An adiabatic body exchanges nothing with the fluid
                        ! BY CONSTRUCTION: no penalization, and every solid
                        ! face masked in the transport kernel. Reporting the
                        ! two terms for it would be reporting a flux that was
                        ! never applied.
                        if (sc%ibmMode(is) == SC_IBM_ADIABATIC) cycle
                        var = VAR_S0 + is
                        s0 = blk%q(i,j,k,var,b)
                        dm = ire/sc%pr(is)
                        dxw = dm
                        dys = dm
                        dzb = dm
                        if (useNut .and. wallfn) then
                            dxw = dm + wall_face_diffusivity(nut(i-1,j,k,b), nut(i,j,k,b), &
                                sc%wfYplus(i-1,j,k,b), sc%wfYplus(i,j,k,b), &
                                sc%pr(is), sc%prt(is), sc%prtModel(is), re, sc%wfP(is), sc%wfYpt(is))
                            dys = dm + wall_face_diffusivity(nut(i,j-1,k,b), nut(i,j,k,b), &
                                sc%wfYplus(i,j-1,k,b), sc%wfYplus(i,j,k,b), &
                                sc%pr(is), sc%prt(is), sc%prtModel(is), re, sc%wfP(is), sc%wfYpt(is))
                            dzb = dm + wall_face_diffusivity(nut(i,j,k-1,b), nut(i,j,k,b), &
                                sc%wfYplus(i,j,k-1,b), sc%wfYplus(i,j,k,b), &
                                sc%pr(is), sc%prt(is), sc%prtModel(is), re, sc%wfP(is), sc%wfYpt(is))
                        else if (useNut) then
                            dxw = dm + eddy_diffusivity(ntw, sc%pr(is), sc%prt(is), sc%prtModel(is), re)
                            dys = dm + eddy_diffusivity(nts, sc%pr(is), sc%prt(is), sc%prtModel(is), re)
                            dzb = dm + eddy_diffusivity(ntb, sc%pr(is), sc%prt(is), sc%prtModel(is), re)
                        end if

                        ! The penalization delivered into this cell. SOLID
                        ! cells are EXCLUDED explicitly, not left to cancel:
                        ! their heat is the staircase term below, and counting
                        ! both double counts it. The exclusion used to be
                        ! implicit -- a solid cell holds the body value, so
                        ! cp*(ibmValue - s0) reads 1e28*0 = 0 (the S3 FINDING)
                        ! -- but that cancellation is an ARTEFACT OF THE VALUE:
                        ! it is bitwise exact only when ibmValue is large
                        ! enough to swallow the O(1e-29) penalization residual
                        ! under its own ulp. With ibmValue = 0 the residual
                        ! survives, cp*(0 - 1.6e-29) is O(0.1) per cell, and
                        ! the reported total came out 128 % high (measured
                        ! 2026-08-05 on validation/scalar/ibmwf180.ini, whose
                        ! walls are at 0). The guard makes the diagnostic
                        ! independent of ibmValue, which it must be: with
                        ! ibmValue = 1 it changes nothing (those terms are
                        ! already exactly 0), so every S3/S4 gate is unmoved.
                        if (.not. solc) then
                            !$omp atomic update
                            heat(2*is) = heat(2*is) + cp*(sc%ibmValue(is) - s0)*dx*dy*dz/sc%pr(is)
                        end if

                        ! The three LOW faces; a face counts only when exactly
                        ! one side is solid, positive INTO the fluid.
                        if ((solc .neqv. solw) .and. .not. clw) then
                            sW = blk%q(i-1,j,k,var,b)
                            flux = blk%q(i,j,k,VAR_U,b)*0.5d0*(sW + s0) &
                                 - dxw*(s0 - sW)*sc%invDx(i,b)
                            sgn = merge(1.0d0, -1.0d0, solw)
                            !$omp atomic update
                            heat(2*is-1) = heat(2*is-1) + sgn*flux*dy*dz
                        end if
                        if ((solc .neqv. sols) .and. .not. cls) then
                            sS = blk%q(i,j-1,k,var,b)
                            flux = blk%q(i,j,k,VAR_V,b)*0.5d0*(sS + s0) &
                                 - dys*(s0 - sS)*sc%invDy(j,b)
                            sgn = merge(1.0d0, -1.0d0, sols)
                            !$omp atomic update
                            heat(2*is-1) = heat(2*is-1) + sgn*flux*dx*dz
                        end if
                        if ((solc .neqv. solb) .and. .not. clb) then
                            sB = blk%q(i,j,k-1,var,b)
                            flux = blk%q(i,j,k,VAR_W,b)*0.5d0*(sB + s0) &
                                 - dzb*(s0 - sB)*sc%invDz(k,b)
                            sgn = merge(1.0d0, -1.0d0, solb)
                            !$omp atomic update
                            heat(2*is-1) = heat(2*is-1) + sgn*flux*dx*dy
                        end if
                    end do
                end do
            end do
        end do
        end do
        !$omp end target teams distribute parallel do
    end subroutine body_heat_kernel

end module scalar_stats
