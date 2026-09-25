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

        ! Trip forcing (SRC_TRIP) state -- an EXACT port of the CaNS/SIMSON trip
        ! (cans_js/src/trip.f90; the AMPHIBIOUS "CaNS-exact trip" port). The force
        ! acts on the wall-normal (v) momentum:
        !   f_v(x,y,z,t) = exp(-((x-x0)/lx)^2 - (y/ly)^2) * g(z,t),
        !   g(z,t) = (amp_s/nmodes) g_s(z)
        !          + (amp_t/nmodes) [ (1-b(t)) g_k(z) + b(t) g_{k+1}(z) ],
        !   g(z) = sqrt(2) cos(phi_0) + sum_{m=1}^{nmodes/2} 2 cos(2 pi m z/Lz + phi_m),
        !   b = p^2 (3-2p),  p = t/ts - k,  k = floor(t/ts).
        ! The spanwise signal is a CaNS-flat Fourier series (mode 0 + nmodes/2
        ! non-zero harmonics, each with fixed amplitude and a random phase), scaled
        ! by amp/nmodes and NOT unit-rms normalized -- so its spanwise rms is
        ! sqrt(nmodes) and the effective forcing amplitude is amp/sqrt(nmodes).
        ! (This is why CaNS amp=0.18854 with nmodes=16 gives an effective 0.047.)
        ! g_k is redrawn every ts; b(t) is the C^1 smooth step (Schlatter & Orlu 2012).
        real(C_DOUBLE) :: trip_x0 = 0.0d0, trip_lx = 4.0d0, trip_ly = 1.0d0
        real(C_DOUBLE) :: trip_amp = 0.0d0, trip_amp_s = 0.0d0, trip_ts = 4.0d0, trip_lz = 0.0d0
        integer(C_INT) :: trip_nmodes = 16_C_INT, trip_seed = 1_C_INT
        integer(C_INT) :: trip_nharm = 0_C_INT        ! nmodes/2 non-zero harmonics
        integer(C_INT) :: trip_kindex = -1_C_INT
        ! Random spanwise PHASES: index 1 = CaNS mode 0, 2..nharm+1 = harmonics.
        ! _s = steady realisation, _old/_new = the two temporal realisations.
        real(C_DOUBLE), allocatable :: trip_phase_s(:), trip_phase_old(:), trip_phase_new(:)
        ! The blocks whose cells can carry a non-zero trip force, listed once
        ! at init. The envelope is a small ellipse -- the kernel's own
        ! ex < -50 test cuts it at |x-x0| < lx*sqrt(50), |y| < ly*sqrt(50),
        ! which on a boundary-layer domain is under 1 % of the volume -- so
        ! almost every block is permanently zero. Those blocks are zeroed once
        ! at allocation and never written again, which is what keeps the
        ! per-substage refresh off them: it used to store three doubles per
        ! cell over the WHOLE domain every substage to write mostly zeros
        ! (4.4-5.5 % of the step at 16 ranks, measured).
        integer(C_INT), allocatable :: trip_blocks(:)
        integer(C_INT) :: n_trip_blocks = 0_C_INT
        ! The spanwise signals g_s(z), g_k(z) and g_{k+1}(z), tabulated per
        ! (k, listed block). They are the nharm-term Fourier sums the field
        ! kernel used to evaluate PER CELL -- 3*(nharm+1) transcendentals each
        ! -- although they depend on z alone. Worse, they only change when the
        ! random walk advances, once per trip_ts (hundreds of steps), so they
        ! are refreshed then and not per substage. The per-cell blend
        ! (1-b)g_k + b g_{k+1} still happens in the field kernel, because b(t)
        ! does change every substage.
        real(C_DOUBLE), allocatable :: trip_gs(:,:), trip_g(:,:), trip_gp1(:,:)  ! (nz, n_trip_blocks)
        logical(C_BOOL) :: trip_span_valid = .false.
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
            call init_trip(bf, dns, blk, c_has_terminal)
        end select
    end subroutine init_bodyforce

    ! Copy the parsed [force] trip_* parameters and prime the random spanwise
    ! states g_0, g_1. Host-only; the field itself is filled each substage by
    ! update_bodyforce -> fill_trip.
    subroutine init_trip(bf, dns, blk, c_has_terminal)
        type(bodyforce_type), intent(inout) :: bf
        type(dns_type), intent(in) :: dns
        type(block_set_type), intent(in) :: blk
        logical, intent(in) :: c_has_terminal

        integer :: n, sz
        integer, allocatable :: seed(:)

        bf%trip_x0     = dns%trip_x0
        bf%trip_lx     = dns%trip_lx
        bf%trip_ly     = dns%trip_ly
        bf%trip_amp    = dns%trip_amp
        bf%trip_amp_s  = dns%trip_amp_s
        bf%trip_ts     = dns%trip_ts
        bf%trip_nmodes = dns%trip_nmodes
        bf%trip_seed   = dns%trip_seed
        bf%trip_lz     = dns%leng(3)

        if (bf%trip_ts <= 0.0d0) error stop "[force] trip_ts must be > 0"
        if (bf%trip_lx <= 0.0d0 .or. bf%trip_ly <= 0.0d0) &
            error stop "[force] trip_lx and trip_ly must be > 0"
        ! CaNS uses nmodes/2 non-zero harmonics, so nmodes must be even and >= 2.
        if (bf%trip_nmodes < 2_C_INT .or. modulo(bf%trip_nmodes, 2_C_INT) /= 0_C_INT) &
            error stop "[force] trip_nmodes must be even and >= 2 (CaNS uses nmodes/2 harmonics)"
        bf%trip_nharm = bf%trip_nmodes/2_C_INT

        n = int(bf%trip_nharm) + 1        ! phase index 1 = mode 0, 2..nharm+1 = harmonics
        allocate(bf%trip_phase_s(n), bf%trip_phase_old(n), bf%trip_phase_new(n))

        ! Deterministic seed so a run is reproducible and rank-independent
        ! (the trip is a global spanwise function evaluated identically on
        ! every rank).
        call random_seed(size=sz)
        allocate(seed(sz))
        seed = bf%trip_seed + 37*[(n, n=0, sz-1)]
        call random_seed(put=seed)
        deallocate(seed)

        ! Reproduce the CaNS init RNG consumption order: steady realisation first
        ! (drawn even when amp_s = 0), then the temporal states at k = 0 and k = 1.
        call trip_gen_phases(bf%trip_phase_s)
        call trip_gen_phases(bf%trip_phase_old)
        call trip_gen_phases(bf%trip_phase_new)
        bf%trip_kindex = 0_C_INT

        call select_trip_blocks(bf, blk)

        if (c_has_terminal) print '(a,es10.3,a,es10.3,a,es10.3,a,i0,a,i0,a,i0,a,i0)', &
            " trip forcing (CaNS-exact): amp=", bf%trip_amp, " x0=", bf%trip_x0, &
            " ts=", bf%trip_ts, " nmodes=", bf%trip_nmodes, " harmonics=", bf%trip_nharm, &
            " active blocks=", bf%n_trip_blocks, "/", blk%nBlocks
    end subroutine init_trip

    ! List the blocks the trip envelope can reach. ex = -((x-x0)/lx)^2
    ! - (y/ly)^2 is separable, so its MAXIMUM over a block is minus the sum of
    ! the squared distances from x0 and 0 to the block's own x and y ranges; if
    ! that maximum is already below the kernel's -50 cutoff, every cell of the
    ! block is exactly zero. The test is therefore the kernel's own test
    ! evaluated at the block's closest point -- skipping a block writes exactly
    ! the zeros it already holds from bf%f = 0.0d0 at allocation, so this is
    ! bit-exact, not an approximation.
    subroutine select_trip_blocks(bf, blk)
        type(bodyforce_type), intent(inout) :: bf
        type(block_set_type), intent(in) :: blk

        integer(C_INT) :: b, n, nx, ny
        real(C_DOUBLE) :: xlo, xhi, ylo, yhi, dx, dy

        nx = blk%nb(1); ny = blk%nb(2)
        allocate(bf%trip_blocks(blk%nBlocks))
        n = 0_C_INT
        do b = 1_C_INT, blk%nBlocks
            xlo = minval(blk%x(1:nx,VAR_V,b)); xhi = maxval(blk%x(1:nx,VAR_V,b))
            ylo = minval(blk%y(1:ny,VAR_V,b)); yhi = maxval(blk%y(1:ny,VAR_V,b))
            dx = max(0.0d0, xlo - bf%trip_x0, bf%trip_x0 - xhi)
            dy = max(0.0d0, ylo, -yhi)
            if (-((dx/bf%trip_lx)**2 + (dy/bf%trip_ly)**2) >= -50.0d0) then
                n = n + 1_C_INT
                bf%trip_blocks(n) = b
            end if
        end do
        bf%n_trip_blocks = n
        if (n > 0_C_INT) allocate(bf%trip_gs(blk%nb(3), n), &
            bf%trip_g(blk%nb(3), n), bf%trip_gp1(blk%nb(3), n))
        bf%trip_span_valid = .false.
    end subroutine select_trip_blocks

    ! Draw the CaNS random spanwise phases: nharm+1 values (mode 0 + harmonics)
    ! uniform in [0, 2*pi]. The CaNS signal amplitudes are fixed (sqrt(2) for mode
    ! zero, 2 per harmonic); only the phase is random, so there is no unit-rms
    ! normalization (unlike the pre-port implementation).
    subroutine trip_gen_phases(phase)
        real(C_DOUBLE), intent(out) :: phase(:)
        call random_number(phase)
        phase = 8.0d0*atan(1.0d0)*phase        ! [0, 2*pi]
    end subroutine trip_gen_phases

    subroutine destroy_bodyforce(bf)
        type(bodyforce_type), intent(inout) :: bf

        if (allocated(bf%f)) deallocate(bf%f)
        if (allocated(bf%trip_phase_s)) deallocate(bf%trip_phase_s, &
            bf%trip_phase_old, bf%trip_phase_new)
        if (allocated(bf%trip_blocks)) deallocate(bf%trip_blocks)
        if (allocated(bf%trip_gs)) deallocate(bf%trip_gs, bf%trip_g, bf%trip_gp1)
        bf%n_trip_blocks = 0_C_INT
        bf%trip_span_valid = .false.
        bf%trip_kindex = -1_C_INT
    end subroutine destroy_bodyforce

    subroutine enter_bodyforce_data(bf)
        type(bodyforce_type), intent(inout) :: bf

        if (.not. allocated(bf%f)) return

        !$omp target enter data map(to: bf)
        !$omp target enter data map(to: bf%f)
        if (allocated(bf%trip_blocks)) then
            !$omp target enter data map(to: bf%trip_blocks)
        end if
        if (allocated(bf%trip_gs)) then
            !$omp target enter data map(to: bf%trip_gs, bf%trip_g, bf%trip_gp1)
        end if
    end subroutine enter_bodyforce_data

    subroutine exit_bodyforce_data(bf)
        type(bodyforce_type), intent(inout) :: bf

        if (.not. allocated(bf%f)) return

        if (allocated(bf%trip_gs)) then
            !$omp target exit data map(delete: bf%trip_gs, bf%trip_g, bf%trip_gp1)
        end if
        if (allocated(bf%trip_blocks)) then
            !$omp target exit data map(delete: bf%trip_blocks)
        end if
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

    ! CaNS/SIMSON trip force on the wall-normal (v) component at time t:
    !   f_v = exp(-((x-x0)/lx)^2 - (y/ly)^2) * g(z,t)  (see the type comment).
    ! The random walk (redrawing the spanwise phases each ts) is advanced on
    ! the HOST -- cheap and rare (once per ts ~ hundreds of steps) -- and the
    ! FIELD is filled by a device kernel (fill_trip_kernel): nothing but a few
    ! scalars crosses to the device, not the whole force field every substage.
    ! No host-side full-field loop, no H2D copy of f. The spanwise signals are
    ! tabulated once per redraw (fill_trip_span), since they depend on z alone.
    subroutine fill_trip(bf, blk, t)
        type(bodyforce_type), intent(inout) :: bf
        type(block_set_type), intent(in) :: blk
        real(C_DOUBLE), intent(in) :: t

        integer(C_INT) :: kidx, kbefore
        real(C_DOUBLE) :: p, bstep

        ! Advance the walk until g_k / g_{k+1} bracket [k*ts, (k+1)*ts] ∋ t.
        kidx = int(floor(t/bf%trip_ts), C_INT)
        kbefore = bf%trip_kindex
        do while (bf%trip_kindex < kidx)
            bf%trip_phase_old = bf%trip_phase_new
            call trip_gen_phases(bf%trip_phase_new)
            bf%trip_kindex = bf%trip_kindex + 1_C_INT
        end do

        ! The spanwise signals depend only on the phases just (re)drawn, so
        ! they are retabulated here and not once per substage.
        if (bf%trip_kindex /= kbefore .or. .not. bf%trip_span_valid) then
            call fill_trip_span(bf, blk)
            bf%trip_span_valid = .true.
        end if

        p = t/bf%trip_ts - real(kidx, C_DOUBLE)
        bstep = p*p*(3.0d0 - 2.0d0*p)        ! 3p^2 - 2p^3, C^1 smooth step

        call fill_trip_kernel(bf, blk, bstep)
    end subroutine fill_trip

    ! Tabulate the CaNS spanwise signals g_s(z), g_k(z) and g_{k+1}(z) over the
    ! listed blocks' z lines. The sums are formed in the same order, from the
    ! same phases, as the per-cell version they replace, so the field they
    ! produce is bit-identical.
    subroutine fill_trip_span(bf, blk)
        type(bodyforce_type), intent(inout) :: bf
        type(block_set_type), intent(in) :: blk

        integer :: k, b, bb, m, nz, nTrip, nharm
        real(C_DOUBLE) :: z, arg, w, gs, gold, gnew, root2

        nTrip = int(bf%n_trip_blocks)
        if (nTrip == 0) return
        nz = int(blk%nb(3))
        nharm = int(bf%trip_nharm)
        w = 8.0d0*atan(1.0d0)/bf%trip_lz         ! 2*pi/Lz
        root2 = sqrt(2.0d0)

#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do collapse(2) &
        !$omp& map(to: w, root2, nharm, nz, nTrip, blk%z, bf%trip_blocks, &
        !$omp& bf%trip_phase_s(1:nharm+1), bf%trip_phase_old(1:nharm+1), &
        !$omp& bf%trip_phase_new(1:nharm+1)) &
        !$omp& map(to: bf%trip_gs, bf%trip_g, bf%trip_gp1) &
        !$omp& private(k,b,bb,m,z,arg,gs,gold,gnew)
#endif
        do bb = 1, nTrip
            do k = 1, nz
                b = int(bf%trip_blocks(bb))
                z = blk%z(k, VAR_V, b)
                ! CaNS mode 0 (spanwise-uniform) + nharm harmonics.
                gs   = root2*cos(bf%trip_phase_s(1))
                gold = root2*cos(bf%trip_phase_old(1))
                gnew = root2*cos(bf%trip_phase_new(1))
                do m = 1, nharm
                    arg = w*real(m, C_DOUBLE)*z
                    gs   = gs   + 2.0d0*cos(arg + bf%trip_phase_s(m+1))
                    gold = gold + 2.0d0*cos(arg + bf%trip_phase_old(m+1))
                    gnew = gnew + 2.0d0*cos(arg + bf%trip_phase_new(m+1))
                end do
                bf%trip_gs(k,bb)  = gs
                bf%trip_g(k,bb)   = gold
                bf%trip_gp1(k,bb) = gnew
            end do
        end do
#ifdef USE_OPENMP_OFFLOAD
        !$omp end target teams distribute parallel do
#endif
    end subroutine fill_trip_span

    ! Device kernel: fill bf%f's v-component from the tabulated spanwise
    ! signals, evaluating the Gaussian envelope per cell. bf%f and blk%{x,y} are
    ! already device-resident; only the scalars cross. On the CPU build this is
    ! a plain loop.
    !   f_v = env * [ (amp_s/nmodes) g_s + (amp_t/nmodes)((1-b) g_old + b g_new) ]
    !
    ! It runs over the LISTED blocks only (see trip_blocks), and writes the v
    ! component only: f_u and f_w are zero for a trip and were already zero from
    ! bf%f = 0.0d0 at allocation, so rewriting them every substage stored two
    ! thirds of a domain-sized array to no effect. Both restrictions write
    ! exactly the values that were there before -- this is bit-exact.
    subroutine fill_trip_kernel(bf, blk, bstep)
        type(bodyforce_type), intent(inout) :: bf
        type(block_set_type), intent(in) :: blk
        real(C_DOUBLE), intent(in) :: bstep

        integer :: i, j, k, b, bb, nx, ny, nz, nTrip
        real(C_DOUBLE) :: x, y, ex, env
        real(C_DOUBLE) :: amps, ampt, x0, lx, ly

        nTrip = int(bf%n_trip_blocks)
        if (nTrip == 0) return

        nx = int(blk%nb(1)); ny = int(blk%nb(2)); nz = int(blk%nb(3))
        ! CaNS scaling: divide the (non-normalized) signal by nmodes = 2*nharm.
        ampt = bf%trip_amp  /real(bf%trip_nmodes, C_DOUBLE)
        amps = bf%trip_amp_s/real(bf%trip_nmodes, C_DOUBLE)
        x0 = bf%trip_x0; lx = bf%trip_lx; ly = bf%trip_ly

#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do collapse(4) &
        !$omp& map(to: bstep, ampt, amps, x0, lx, ly, nx, ny, nz, nTrip, &
        !$omp& blk%x, blk%y, bf%trip_blocks, bf%trip_gs, bf%trip_g, bf%trip_gp1) &
        !$omp& map(to: bf%f) &
        !$omp& private(i,j,k,b,bb,x,y,ex,env)
#endif
        do bb = 1, nTrip
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    b = int(bf%trip_blocks(bb))
                    ! v (VAR_V) lives at blk%{x,y}(:,VAR_V,b); the z dependence
                    ! is entirely in the tabulated spanwise signals.
                    x = blk%x(i, VAR_V, b)
                    y = blk%y(j, VAR_V, b)
                    ex = -((x - x0)/lx)**2 - (y/ly)**2
                    if (ex < -50.0d0) then
                        bf%f(i,j,k,VAR_V,b) = 0.0d0
                    else
                        env = exp(ex)
                        bf%f(i,j,k,VAR_V,b) = env* &
                            (amps*bf%trip_gs(k,bb) + ampt*((1.0d0 - bstep)*bf%trip_g(k,bb) &
                            + bstep*bf%trip_gp1(k,bb)))
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
