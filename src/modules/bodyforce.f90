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

module bodyforce
    use, intrinsic :: iso_c_binding
    use :: init, only: dns_type, grid_type, VAR_U, VAR_V, VAR_W, NVAR, NVEL
    use :: blocks, only: block_set_type
    use :: io, only: read_force_file
    use :: volume_force, only: fill_volume_force
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
    ! Steady spatially-varying force defined by the student hook in
    ! volume_force.f90; filled once at init like a profile.
    integer(C_INT), parameter, public :: SRC_STEADY  = 4_C_INT

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

        if (bf%source == SRC_NONE) error stop "[force] type must be profile, steady, file or custom"
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
        case (SRC_STEADY)
            ! Steady student-defined force: fill once from volume_force.f90.
            call fill_volume_force(bf%f, blk)
        case (SRC_FILE)
            call read_force_file(bf%f, blk, dns, dns%force_file, c_has_terminal)
        case (SRC_CUSTOM)
            ! Left zeroed; the user fills it via update_bodyforce in the loop.
        end select
    end subroutine init_bodyforce

    subroutine destroy_bodyforce(bf)
        type(bodyforce_type), intent(inout) :: bf

        if (allocated(bf%f)) deallocate(bf%f)
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
    !------------------------------------------------------------------!
    subroutine update_bodyforce(bf, blk, dns, g, t)
        type(bodyforce_type), intent(inout) :: bf
        type(block_set_type), intent(in) :: blk
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        real(C_DOUBLE), intent(in) :: t

        ! Only the custom source is driven from the time loop; profile/file
        ! are already filled at init.
        if (.not. bodyforce_is_enabled(bf)) return
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
        case ("steady")
            id = SRC_STEADY
        case ("file")
            id = SRC_FILE
        case ("custom")
            id = SRC_CUSTOM
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
