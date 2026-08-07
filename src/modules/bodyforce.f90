!--------------------------!
!                          !
!     Volumetric body      !
!      force module        !
!                          !
!--------------------------!
!
! Optional spatially-varying volumetric force f(x) added to the momentum
! predictor ON TOP of the constant [flow] forcing_* term. Disabled by
! default -> the solver is byte-identical to a build without this module
! (the correction kernel in step.f90 is simply never called).
!
! The force owns its own flat contiguous array, mapped to the device once
! (per the CLAUDE.md GPU-data convention). Each component lives at ITS
! staggered face location, mirroring blk%q / dns%forcing. The array holds
! only interior cells (1:nb): the momentum predictor reads the force at the
! predicted faces, so no halos are needed.
!
! Three sources ([force] type):
!   profile  a built-in named analytic form filled once at init at each
!            component's staggered coordinate (no expression parser --
!            "constant" and "sine"; anything else uses custom).
!   file     fx/fy/fz read from an HDF5 field laid out like a velocity
!            field (un/vn/wn), via the io block-table read path.
!   custom   the user fills bf%f inside the RK loop by editing the clearly
!            marked hook update_bodyforce below (time-dependent forcing,
!            controllers, actuators).
!   trip     the Schlatter & Orlu (2012) random wall-normal trip forcing
!            that triggers laminar->turbulent transition in a boundary
!            layer (fill_trip below; refreshed each substage from the loop).

module bodyforce
    use, intrinsic :: iso_c_binding
    use :: init, only: dns_type, grid_type, VAR_U, VAR_V, VAR_W, NVAR, NVEL
    use :: blocks, only: block_set_type
    use :: io, only: read_force_file
    implicit none

    private
    public :: bodyforce_type
    public :: init_bodyforce, destroy_bodyforce
    public :: enter_bodyforce_data, exit_bodyforce_data
    public :: bodyforce_is_enabled, bodyforce_zero
    public :: bodyforce_update_to_device, bodyforce_update_from_device
    public :: update_bodyforce

    ! Force source ([force] type).
    integer(C_INT), parameter, public :: SRC_NONE    = 0_C_INT
    integer(C_INT), parameter, public :: SRC_PROFILE = 1_C_INT
    integer(C_INT), parameter, public :: SRC_FILE    = 2_C_INT
    integer(C_INT), parameter, public :: SRC_CUSTOM  = 3_C_INT
    integer(C_INT), parameter, public :: SRC_TRIP    = 4_C_INT

    ! Named analytic profiles ([force] profile). Add new forms here and in
    ! fill_profile below; keep them cheap closed-form expressions.
    integer(C_INT), parameter, public :: PROF_CONSTANT = 1_C_INT
    integer(C_INT), parameter, public :: PROF_SINE     = 2_C_INT

    type :: bodyforce_type
        logical(C_BOOL) :: enabled = .false.
        integer(C_INT)  :: source  = SRC_NONE
        integer(C_INT)  :: profile = PROF_CONSTANT
        ! Per-component amplitude and (profile) wavenumber; prof_dir selects
        ! the coordinate the profile varies along (1=x, 2=y, 3=z).
        real(C_DOUBLE)  :: amp(1:3)      = 0.0d0
        real(C_DOUBLE)  :: wavenumber(1:3) = 0.0d0
        integer(C_INT)  :: prof_dir = 1_C_INT
        ! The force at each component's staggered face, interior cells only.
        real(C_DOUBLE), allocatable :: f(:,:,:,:,:)  ! (1:nb,1:nb,1:nb,NVEL,nBlocks)

        ! Trip forcing (SRC_TRIP) state. The force acts on the wall-normal
        ! (v) momentum:
        !   f_v(x,y,z,t) = amp * exp(-((x-x0)/lx)^2 - (y/ly)^2) * g(z,t),
        !   g(z,t) = (1-b(t)) g_k(z) + b(t) g_{k+1}(z),  b = 3p^2 - 2p^3,
        !   p = t/ts - k,  k = floor(t/ts),
        ! g_k(z) a unit-rms random spanwise function (nmodes Fourier modes,
        ! period Lz), regenerated each ts; b(t) is the C^1 smooth step that
        ! makes the forcing continuous in time (Schlatter & Orlu 2012).
        real(C_DOUBLE) :: trip_x0 = 0.0d0, trip_lx = 4.0d0, trip_ly = 1.0d0
        real(C_DOUBLE) :: trip_amp = 0.0d0, trip_ts = 4.0d0, trip_lz = 0.0d0
        integer(C_INT) :: trip_nmodes = 16_C_INT, trip_seed = 1_C_INT
        integer(C_INT) :: trip_kindex = -1_C_INT
        ! Spanwise Fourier coefficients of g_k (…ak/…bk) and g_{k+1} (…akp1/…bkp1).
        real(C_DOUBLE), allocatable :: trip_ak(:), trip_bk(:)
        real(C_DOUBLE), allocatable :: trip_akp1(:), trip_bkp1(:)
    end type bodyforce_type

contains

    logical function bodyforce_is_enabled(bf)
        type(bodyforce_type), intent(in) :: bf

        bodyforce_is_enabled = bf%enabled .and. bf%source /= SRC_NONE
    end function bodyforce_is_enabled

    ! Interpret the parsed [force] config (strings on dns) and, when enabled,
    ! allocate and fill f. profile/file are filled once here; custom starts
    ! zeroed and is (re)filled by update_bodyforce each substage.
    subroutine init_bodyforce(bf, dns, blk, g, c_has_terminal)
        type(bodyforce_type), intent(inout) :: bf
        type(dns_type), intent(in) :: dns
        type(block_set_type), intent(in) :: blk
        type(grid_type), intent(in) :: g
        logical, intent(in) :: c_has_terminal

        integer :: nx, ny, nz

        call destroy_bodyforce(bf)

        bf%enabled = dns%force_enabled
        if (.not. bf%enabled) return

        bf%source     = force_source_id(dns%force_type)
        bf%profile    = force_profile_id(dns%force_profile)
        bf%amp        = dns%force_amp
        bf%wavenumber = dns%force_wavenumber
        bf%prof_dir   = dns%force_dir

        if (bf%source == SRC_NONE) error stop "[force] type must be profile, file, custom or trip"
        if (bf%prof_dir < 1_C_INT .or. bf%prof_dir > 3_C_INT) &
            error stop "[force] dir must be 1, 2 or 3"

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))
        allocate(bf%f(nx, ny, nz, NVEL, blk%nBlocks))
        bf%f = 0.0d0

        select case (bf%source)
        case (SRC_PROFILE)
            call fill_profile(bf, blk, nx, ny, nz)
        case (SRC_FILE)
            call read_force_file(bf%f, blk, dns, dns%force_file, c_has_terminal)
        case (SRC_CUSTOM)
            ! Left zeroed; the user fills it via update_bodyforce in the loop.
        case (SRC_TRIP)
            call init_trip(bf, dns, c_has_terminal)
        end select
    end subroutine init_bodyforce

    ! Copy the parsed [force] trip_* parameters and prime the random spanwise
    ! states g_0, g_1. Host-only; the field itself is filled each substage by
    ! update_bodyforce -> fill_trip.
    subroutine init_trip(bf, dns, c_has_terminal)
        type(bodyforce_type), intent(inout) :: bf
        type(dns_type), intent(in) :: dns
        logical, intent(in) :: c_has_terminal

        integer :: n, sz
        integer, allocatable :: seed(:)

        bf%trip_x0     = dns%trip_x0
        bf%trip_lx     = dns%trip_lx
        bf%trip_ly     = dns%trip_ly
        bf%trip_amp    = dns%trip_amp
        bf%trip_ts     = dns%trip_ts
        bf%trip_nmodes = dns%trip_nmodes
        bf%trip_seed   = dns%trip_seed
        bf%trip_lz     = dns%leng(3)

        if (bf%trip_ts <= 0.0d0) error stop "[force] trip_ts must be > 0"
        if (bf%trip_lx <= 0.0d0 .or. bf%trip_ly <= 0.0d0) &
            error stop "[force] trip_lx and trip_ly must be > 0"
        if (bf%trip_nmodes < 1_C_INT) error stop "[force] trip_nmodes must be >= 1"

        n = int(bf%trip_nmodes)
        allocate(bf%trip_ak(n), bf%trip_bk(n), bf%trip_akp1(n), bf%trip_bkp1(n))

        ! Deterministic seed so a run is reproducible and rank-independent
        ! (the trip is a global spanwise function evaluated identically on
        ! every rank).
        call random_seed(size=sz)
        allocate(seed(sz))
        seed = bf%trip_seed + 37*[(n, n=0, sz-1)]
        call random_seed(put=seed)
        deallocate(seed)

        ! Prime the walk: g_k for k = 0 and g_{k+1} for k = 1.
        call trip_gen_coeffs(bf%trip_ak,   bf%trip_bk)
        call trip_gen_coeffs(bf%trip_akp1, bf%trip_bkp1)
        bf%trip_kindex = 0_C_INT

        if (c_has_terminal) print '(a,es10.3,a,es10.3,a,es10.3,a,i0)', &
            " trip forcing: amp=", bf%trip_amp, " x0=", bf%trip_x0, &
            " ts=", bf%trip_ts, " nmodes=", bf%trip_nmodes
    end subroutine init_trip

    ! Draw a unit-rms random spanwise function: nmodes Fourier coefficients
    ! uniform in [-1,1], then normalized so var(g) = sum(a^2+b^2)/2 = 1.
    subroutine trip_gen_coeffs(a, b)
        real(C_DOUBLE), intent(out) :: a(:), b(:)
        real(C_DOUBLE) :: r, s
        integer :: n

        call random_number(a); a = 2.0d0*a - 1.0d0
        call random_number(b); b = 2.0d0*b - 1.0d0
        s = 0.0d0
        do n = 1, size(a)
            s = s + a(n)*a(n) + b(n)*b(n)
        end do
        s = sqrt(0.5d0*s)
        if (s > 0.0d0) then
            r = 1.0d0/s
            a = a*r
            b = b*r
        end if
    end subroutine trip_gen_coeffs

    subroutine destroy_bodyforce(bf)
        type(bodyforce_type), intent(inout) :: bf

        if (allocated(bf%f)) deallocate(bf%f)
        if (allocated(bf%trip_ak)) deallocate(bf%trip_ak, bf%trip_bk, &
            bf%trip_akp1, bf%trip_bkp1)
        bf%trip_kindex = -1_C_INT
    end subroutine destroy_bodyforce

    subroutine enter_bodyforce_data(bf)
        type(bodyforce_type), intent(inout) :: bf

        if (.not. allocated(bf%f)) return

        !$omp target enter data map(to: bf)
        !$omp target enter data map(to: bf%f)
    end subroutine enter_bodyforce_data

    subroutine exit_bodyforce_data(bf)
        type(bodyforce_type), intent(inout) :: bf

        if (.not. allocated(bf%f)) return

        !$omp target exit data map(delete: bf%f)
        !$omp target exit data map(delete: bf)
    end subroutine exit_bodyforce_data

    ! Device-side zero of f (for a custom fill that adds increments).
    subroutine bodyforce_zero(bf)
        type(bodyforce_type), intent(inout) :: bf

        integer :: i, j, k, v, b, nx, ny, nz, nBlocks

        if (.not. allocated(bf%f)) return
        nx = size(bf%f, 1); ny = size(bf%f, 2); nz = size(bf%f, 3)
        nBlocks = size(bf%f, 5)

        !$omp target teams distribute parallel do collapse(5) &
        !$omp& map(to: bf%f) private(i,j,k,v,b)
        do b = 1, nBlocks
        do v = 1, int(NVEL)
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    bf%f(i,j,k,v,b) = 0.0d0
                end do
            end do
        end do
        end do
        end do
        !$omp end target teams distribute parallel do
    end subroutine bodyforce_zero

    ! Push a host-side fill of f to the device (call after writing bf%f on the
    ! host in update_bodyforce, before momentum uses it).
    subroutine bodyforce_update_to_device(bf)
        type(bodyforce_type), intent(inout) :: bf

        if (.not. allocated(bf%f)) return
#ifdef USE_OPENMP_OFFLOAD
        !$omp target update to(bf%f)
#endif
    end subroutine bodyforce_update_to_device

    ! Pull the device copy of f back to the host (e.g. to inspect it).
    subroutine bodyforce_update_from_device(bf)
        type(bodyforce_type), intent(inout) :: bf

        if (.not. allocated(bf%f)) return
#ifdef USE_OPENMP_OFFLOAD
        !$omp target update from(bf%f)
#endif
    end subroutine bodyforce_update_from_device

    !------------------------------------------------------------------!
    !  USER HOOK -- edit this routine for [force] type = custom.        !
    !                                                                    !
    !  Called once per RK substage BEFORE the momentum predictor, with  !
    !  the current time t. Fill bf%f(i,j,k,VAR_U/V/W,b) as a function of !
    !  the staggered coordinates blk%x/y/z(:,VAR_*,b), the current field !
    !  blk%q (for controllers), and t. Each component lives at its own  !
    !  staggered face; interior indices run 1:blk%nb.                    !
    !                                                                    !
    !  Host-side fill: write bf%f on the host, then call                 !
    !  bodyforce_update_to_device(bf) (done for you at the end here).    !
    !  A bulk-velocity controller would comm_allreduce_sum blk%q over    !
    !  the domain, then write a uniform bf%f.                            !
    !                                                                    !
    !  BUT on an offload build blk%q's HOST copy is STALE inside the     !
    !  time loop -- the solver's kernels write the DEVICE copy and       !
    !  nothing pulls it back (moby_solve.f90's `target update           !
    !  from(blk%q)` runs at init only; io.f90 pulls it per snapshot).    !
    !  A controller reading blk%q here would silently integrate the      !
    !  INITIAL field on GPU and the current one on CPU. Put an explicit  !
    !  `!$omp target update from(blk%q)` (guarded by                     !
    !  USE_OPENMP_OFFLOAD) ahead of the read, or do the reduction in a   !
    !  target region.                                                    !
    !------------------------------------------------------------------!
    subroutine update_bodyforce(bf, blk, dns, g, t)
        type(bodyforce_type), intent(inout) :: bf
        type(block_set_type), intent(in) :: blk
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        real(C_DOUBLE), intent(in) :: t

        ! The trip source is refreshed each substage (time-dependent);
        ! profile/file are already filled at init.
        if (.not. bodyforce_is_enabled(bf)) return
        if (bf%source == SRC_TRIP) then
            ! fill_trip fills bf%f directly on the device (no host copy).
            call fill_trip(bf, blk, t)
            return
        end if
        if (bf%source /= SRC_CUSTOM) return

        !==== USER CODE HERE (fill bf%f on the host) ==================
        ! Example (uncomment and adapt): a streamwise force that ramps in
        ! time and varies sinusoidally across the channel height.
        ! integer :: i, j, k, b
        ! do b = 1, int(blk%nBlocks)
        !     do k = 1, int(blk%nb(3))
        !         do j = 1, int(blk%nb(2))
        !             do i = 1, int(blk%nb(1))
        !                 bf%f(i,j,k,VAR_U,b) = tanh(t) &
        !                     * sin(blk%y(j,VAR_U,b))
        !             end do
        !         end do
        !     end do
        ! end do
        !=============================================================

        ! Nothing filled by default -> a no-op until the user edits the block
        ! above. Ship the (possibly updated) host array to the device.
        call bodyforce_update_to_device(bf)
    end subroutine update_bodyforce

    ! Schlatter & Orlu trip force on the wall-normal (v) component at time t:
    !   f_v = amp * exp(-((x-x0)/lx)^2 - (y/ly)^2) * g(z,t).
    ! The random walk (regenerating the nmodes spanwise coefficients each ts)
    ! is advanced on the HOST -- cheap and rare (once per ts ~ hundreds of
    ! steps) -- and the FIELD is filled by a device kernel (fill_trip_kernel):
    ! only the small coefficient arrays cross to the device, not the whole
    ! force field every substage. No host-side full-field loop, no H2D copy of
    ! f. The arithmetic matches the old host fill exactly (CPU bit-exact).
    subroutine fill_trip(bf, blk, t)
        type(bodyforce_type), intent(inout) :: bf
        type(block_set_type), intent(in) :: blk
        real(C_DOUBLE), intent(in) :: t

        integer(C_INT) :: kidx
        real(C_DOUBLE) :: p, bstep

        ! Advance the walk until g_k / g_{k+1} bracket [k*ts, (k+1)*ts] ∋ t.
        kidx = int(floor(t/bf%trip_ts), C_INT)
        do while (bf%trip_kindex < kidx)
            bf%trip_ak = bf%trip_akp1
            bf%trip_bk = bf%trip_bkp1
            call trip_gen_coeffs(bf%trip_akp1, bf%trip_bkp1)
            bf%trip_kindex = bf%trip_kindex + 1_C_INT
        end do

        p = t/bf%trip_ts - real(kidx, C_DOUBLE)
        bstep = p*p*(3.0d0 - 2.0d0*p)        ! 3p^2 - 2p^3, C^1 smooth step

        call fill_trip_kernel(bf, blk, bstep, bf%trip_ak, bf%trip_bk, &
            bf%trip_akp1, bf%trip_bkp1)
    end subroutine fill_trip

    ! Device kernel: fill bf%f's v-component from the (small) coefficient
    ! arrays, evaluating the Gaussian envelope and the spanwise Fourier sum
    ! per cell. bf%f and blk%{x,y,z} are already device-resident; only the
    ! coefficients + scalars cross. On the CPU build this is a plain loop.
    subroutine fill_trip_kernel(bf, blk, bstep, ak, bk, akp1, bkp1)
        type(bodyforce_type), intent(inout) :: bf
        type(block_set_type), intent(in) :: blk
        real(C_DOUBLE), intent(in) :: bstep, ak(:), bk(:), akp1(:), bkp1(:)

        integer :: i, j, k, b, m, nx, ny, nz, nBlocks, nm
        real(C_DOUBLE) :: x, y, z, ex, env, arg, w, gk, gkp1
        real(C_DOUBLE) :: amp, x0, lx, ly

        nx = int(blk%nb(1)); ny = int(blk%nb(2)); nz = int(blk%nb(3))
        nBlocks = int(blk%nBlocks)
        nm = int(bf%trip_nmodes)
        amp = bf%trip_amp; x0 = bf%trip_x0; lx = bf%trip_lx; ly = bf%trip_ly
        w = 8.0d0*atan(1.0d0)/bf%trip_lz         ! 2*pi/Lz

#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do collapse(4) &
        !$omp& map(to: bstep, amp, x0, lx, ly, w, nm, nx, ny, nz, &
        !$omp& ak(1:nm), bk(1:nm), akp1(1:nm), bkp1(1:nm), blk%x, blk%y, blk%z) &
        !$omp& map(to: bf%f) &
        !$omp& private(i,j,k,b,m,x,y,z,ex,env,arg,gk,gkp1)
#endif
        do b = 1, nBlocks
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    bf%f(i,j,k,VAR_U,b) = 0.0d0
                    bf%f(i,j,k,VAR_W,b) = 0.0d0
                    ! v (VAR_V) lives at blk%{x,y,z}(:,VAR_V,b).
                    x = blk%x(i, VAR_V, b)
                    y = blk%y(j, VAR_V, b)
                    z = blk%z(k, VAR_V, b)
                    ex = -((x - x0)/lx)**2 - (y/ly)**2
                    if (ex < -50.0d0) then
                        bf%f(i,j,k,VAR_V,b) = 0.0d0
                    else
                        env = exp(ex)
                        gk = 0.0d0; gkp1 = 0.0d0
                        do m = 1, nm
                            arg = w*real(m, C_DOUBLE)*z
                            gk   = gk   + ak(m)  *cos(arg) + bk(m)  *sin(arg)
                            gkp1 = gkp1 + akp1(m)*cos(arg) + bkp1(m)*sin(arg)
                        end do
                        bf%f(i,j,k,VAR_V,b) = amp*env* &
                            ((1.0d0 - bstep)*gk + bstep*gkp1)
                    end if
                end do
            end do
        end do
        end do
#ifdef USE_OPENMP_OFFLOAD
        !$omp end target teams distribute parallel do
#endif
    end subroutine fill_trip_kernel

    ! Fill f from a named analytic profile at each component's staggered
    ! coordinate. Host-side; the caller maps f to the device afterwards.
    subroutine fill_profile(bf, blk, nx, ny, nz)
        type(bodyforce_type), intent(inout) :: bf
        type(block_set_type), intent(in) :: blk
        integer, intent(in) :: nx, ny, nz

        integer :: i, j, k, v, b, d
        real(C_DOUBLE) :: coord(1:3), kd

        do b = 1, int(blk%nBlocks)
            do v = 1, int(NVEL)
                do k = 1, nz
                    do j = 1, ny
                        do i = 1, nx
                            ! Component v sits at blk%{x,y,z}(:,v,b).
                            coord(1) = blk%x(i, v, b)
                            coord(2) = blk%y(j, v, b)
                            coord(3) = blk%z(k, v, b)
                            select case (bf%profile)
                            case (PROF_CONSTANT)
                                bf%f(i,j,k,v,b) = bf%amp(v)
                            case (PROF_SINE)
                                d = int(bf%prof_dir)
                                kd = bf%wavenumber(d)
                                bf%f(i,j,k,v,b) = bf%amp(v)*sin(kd*coord(d))
                            end select
                        end do
                    end do
                end do
            end do
        end do
    end subroutine fill_profile

    integer(C_INT) function force_source_id(name) result(id)
        character(len=*), intent(in) :: name

        select case (trim(adjustl(name)))
        case ("profile", "")
            id = SRC_PROFILE
        case ("file")
            id = SRC_FILE
        case ("custom")
            id = SRC_CUSTOM
        case ("trip")
            id = SRC_TRIP
        case default
            id = SRC_NONE
        end select
    end function force_source_id

    integer(C_INT) function force_profile_id(name) result(id)
        character(len=*), intent(in) :: name

        select case (trim(adjustl(name)))
        case ("constant", "")
            id = PROF_CONSTANT
        case ("sine", "sin")
            id = PROF_SINE
        case default
            id = PROF_CONSTANT
        end select
    end function force_profile_id

end module bodyforce
