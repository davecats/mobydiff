module boundarylayer_flow
    ! Spatially developing flat-plate boundary layer ([case] name =
    ! boundaryLayer). Moves the boundary-layer defaults into the case (as the
    ! channel/airfoil cases do): a Blasius similarity inlet (u + entrainment v),
    ! a Dirichlet-pressure outlet downstream and at the top, a no-slip wall, a
    ! spanwise-periodic direction, and a one-sided wall-clustered natural grid.
    ! The Schlatter & Orlu (2012) random wall-normal TRIP forcing is enabled by
    ! default from here (the [force] section still overrides). Statistics are
    ! profiles of (x, y) averaged in span + time (bl_stats), and the runtime
    ! line reports the divergence residual, Linf velocity, CFL, Peclet and
    ! wall-clock per step. All defaults are set-if-unset -- explicit [boundary]/
    ! [grid]/[force] keys in the ini win (parsed after apply_defaults).
    use, intrinsic :: iso_c_binding
    use :: flow_case_base, only: case_type
    use :: generic_flow, only: set_generic_defaults
    use :: init, only: dns_type, grid_type, VAR_U, VAR_V, VAR_W, GRID_BLAYER
    use :: blocks, only: block_set_type
    use :: boundary, only: boundary_type, boundary_face_id, &
        PATCH_INLET, PATCH_OUTLET, PATCH_WALL, PROFILE_BLASIUS
    use :: pressure_solver, only: pressure_solver_type
    use :: ibmm, only: ibm_type
    use :: comm, only: comm_type
    use :: case_config_helpers, only: next_config_entry, to_lower, clean_config_string
    use :: boundarylayer_profile, only: initialise_boundarylayer_field
    use :: boundarylayer_stats, only: bl_stats_type
    implicit none

    private

    character(len=*), parameter :: BL_CASE_NAME = "boundarylayer"
    integer, parameter :: STREAM_DIR = 1
    integer, parameter :: WALL_DIR = 2
    integer, parameter :: SPAN_DIR = 3
    real(C_DOUBLE), parameter :: BLASIUS_H = 2.591007d0     ! delta*/theta

    type, extends(case_type), public :: boundarylayer_case_type
        real(C_DOUBLE) :: u_inf = 1.0d0
        ! Inflow momentum thickness. Default = 1/H so that with lengths
        ! nondimensionalized by the inflow displacement thickness delta*_0 = 1
        ! (the Skote convention), [flow] re is Re_delta*,0.
        real(C_DOUBLE) :: theta_in = 1.0d0/BLASIUS_H
        ! One-sided wall-clustered grid. GRID_BLAYER resolves the turbulent
        ! layer [0, resolved_height] with natural clustering and coarsens the
        ! freestream above it (the domain is far taller than the boundary
        ! layer, so a plain natural line would waste most points aloft).
        ! resolved_height should be ~2-3x the outlet delta99.
        real(C_DOUBLE) :: natural_blend_index = 45.0d0
        real(C_DOUBLE) :: dyw_plus = 0.15d0
        real(C_DOUBLE) :: resolved_height = 30.0d0
        ! Trip forcing defaults (mirror the [force] trip_* keys; written into
        ! dns in apply_defaults, still overridable by an explicit [force]).
        logical :: trip_enabled = .true.
        real(C_DOUBLE) :: trip_x0 = 15.0d0, trip_lx = 4.0d0, trip_ly = 1.0d0
        real(C_DOUBLE) :: trip_amp = 0.15d0, trip_ts = 4.0d0
        integer(C_INT) :: trip_nmodes = 24_C_INT, trip_seed = 1_C_INT
        type(bl_stats_type) :: stats
    contains
        procedure :: read_config => bl_read_config
        procedure :: apply_defaults => bl_apply_defaults
        procedure :: setup_after_grid => bl_setup_after_grid
        procedure :: initialise_fields => bl_initialise_fields
        procedure :: after_step => bl_after_step
        procedure :: finalize => bl_finalize
    end type boundarylayer_case_type

    public :: create_boundarylayer_case, BL_CASE_NAME

contains

    subroutine create_boundarylayer_case(flow)
        class(case_type), allocatable, intent(out) :: flow

        allocate(boundarylayer_case_type :: flow)
        flow%name = BL_CASE_NAME
    end subroutine create_boundarylayer_case

    subroutine bl_apply_defaults(this, dns, g, bc, c, ps)
        class(boundarylayer_case_type), intent(inout) :: this
        type(dns_type), intent(inout) :: dns
        type(grid_type), intent(inout) :: g
        type(boundary_type), intent(inout) :: bc
        type(comm_type), intent(inout) :: c
        type(pressure_solver_type), intent(inout) :: ps

        integer :: inlet, outlet, wall, top

        call set_generic_defaults(dns, g, bc, c, ps)
        this%name = BL_CASE_NAME

        dns%ibm_enabled = .false.
        dns%forcing = 0.0d0                    ! no constant forcing (trip does the tripping)
        dns%initial_velocity = [this%u_inf, 0.0d0, 0.0d0]

        ! Faces: streamwise inflow/outflow + wall/top, spanwise periodic.
        bc%isPeriodic = .false.
        bc%isPeriodic(SPAN_DIR) = .true.

        inlet  = boundary_face_id(STREAM_DIR, 0)
        outlet = boundary_face_id(STREAM_DIR, 1)
        wall   = boundary_face_id(WALL_DIR, 0)
        top    = boundary_face_id(WALL_DIR, 1)

        ! Blasius inlet: u and the entrainment v carry the similarity profile
        ! (blasius value profile on the x_min face). blasius_theta = theta_in.
        bc%facePatchType(inlet) = PATCH_INLET
        bc%faceBcDefaultValue(VAR_U, inlet) = this%u_inf
        bc%faceBcDefaultValue(VAR_V, inlet) = this%u_inf
        bc%faceBcProfile(VAR_U, inlet) = PROFILE_BLASIUS
        bc%faceBcProfile(VAR_V, inlet) = PROFILE_BLASIUS
        bc%blasiusTheta = this%theta_in

        bc%facePatchType(outlet) = PATCH_OUTLET   ! zero-gradient velocity + Dirichlet p
        bc%facePatchType(wall)   = PATCH_WALL     ! no-slip plate
        bc%facePatchType(top)    = PATCH_OUTLET   ! displacement entrainment leaves

        ! Boundary-layer grid: natural wall clustering resolving
        ! [0, resolved_height], geometric coarsening in the freestream above.
        g%distribution(WALL_DIR) = GRID_BLAYER
        g%natural_one_sided(WALL_DIR) = .true.
        g%stretch(WALL_DIR) = this%natural_blend_index
        g%natural_dyw_plus(WALL_DIR) = this%dyw_plus
        g%natural_outer_height(WALL_DIR) = this%resolved_height

        ! Trip forcing (Schlatter & Orlu): write the case defaults into the
        ! [force] dns fields; an explicit [force] section still overrides.
        if (this%trip_enabled) then
            dns%force_enabled = .true.
            dns%force_type = "trip"
            dns%trip_x0 = this%trip_x0
            dns%trip_lx = this%trip_lx
            dns%trip_ly = this%trip_ly
            dns%trip_amp = this%trip_amp
            dns%trip_ts = this%trip_ts
            dns%trip_nmodes = this%trip_nmodes
            dns%trip_seed = this%trip_seed
        end if
    end subroutine bl_apply_defaults

    subroutine bl_setup_after_grid(this, blk, dns, g, bc, c)
        class(boundarylayer_case_type), intent(inout) :: this
        type(block_set_type), intent(inout) :: blk
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        type(boundary_type), intent(in) :: bc
        type(comm_type), intent(in) :: c

        call this%stats%setup(blk, dns, g, c)
    end subroutine bl_setup_after_grid

    subroutine bl_initialise_fields(this, blk, dns, g, bc, c)
        class(boundarylayer_case_type), intent(inout) :: this
        type(block_set_type), intent(inout) :: blk
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        type(boundary_type), intent(in) :: bc
        type(comm_type), intent(in) :: c

        call initialise_boundarylayer_field(blk, dns, this%u_inf, this%theta_in)
    end subroutine bl_initialise_fields

    subroutine bl_after_step(this, blk, dns, g, c, ibm)
        class(boundarylayer_case_type), intent(inout) :: this
        type(block_set_type), intent(inout) :: blk
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        type(comm_type), intent(in) :: c
        type(ibm_type), intent(in) :: ibm

        call this%stats%after_step(blk, dns, g, c)
    end subroutine bl_after_step

    subroutine bl_finalize(this, dns, g, c)
        class(boundarylayer_case_type), intent(inout) :: this
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        type(comm_type), intent(in) :: c

        call this%stats%finalize(dns, c)
    end subroutine bl_finalize

    subroutine bl_read_config(this, input_file, has_terminal)
        class(boundarylayer_case_type), intent(inout) :: this
        character(len=*), intent(in) :: input_file
        logical, intent(in), optional :: has_terminal

        integer :: unit, stat, line_no
        character(len=512) :: key, value
        character(len=64) :: section
        logical :: exists, ok

        section = ""
        line_no = 0
        inquire(file=trim(input_file), exist=exists)
        if (.not. exists) return
        open(newunit=unit, file=trim(input_file), status="old", action="read", iostat=stat)
        if (stat /= 0) return

        do
            call next_config_entry(unit, section, key, value, line_no, ok)
            if (.not. ok) exit
            if (trim(to_lower(section)) /= "case.boundarylayer") cycle
            call apply_bl_case_value(this, to_lower(key), value, line_no, has_terminal)
        end do
        close(unit)
    end subroutine bl_read_config

    subroutine apply_bl_case_value(this, key, value, line_no, has_terminal)
        type(boundarylayer_case_type), intent(inout) :: this
        character(len=*), intent(in) :: key, value
        integer, intent(in) :: line_no
        logical, intent(in), optional :: has_terminal

        integer :: int_value, stat
        real(C_DOUBLE) :: real_value
        logical :: bool_value, terminal

        terminal = .true.
        if (present(has_terminal)) terminal = has_terminal

        select case (trim(key))
        case ("u_inf")
            read(value, *, iostat=stat) real_value
            if (stat == 0) this%u_inf = real_value
        case ("theta_in", "blasius_theta")
            read(value, *, iostat=stat) real_value
            if (stat == 0) this%theta_in = real_value
        case ("natural_blend_index", "jb")
            read(value, *, iostat=stat) real_value
            if (stat == 0) this%natural_blend_index = real_value
        case ("dyw_plus")
            read(value, *, iostat=stat) real_value
            if (stat == 0) this%dyw_plus = real_value
        case ("resolved_height", "outer_height")
            read(value, *, iostat=stat) real_value
            if (stat == 0) this%resolved_height = real_value
        case ("trip_enabled")
            if (read_bool(value, bool_value)) this%trip_enabled = bool_value
        case ("trip_x0")
            read(value, *, iostat=stat) real_value
            if (stat == 0) this%trip_x0 = real_value
        case ("trip_lx")
            read(value, *, iostat=stat) real_value
            if (stat == 0) this%trip_lx = real_value
        case ("trip_ly")
            read(value, *, iostat=stat) real_value
            if (stat == 0) this%trip_ly = real_value
        case ("trip_amp")
            read(value, *, iostat=stat) real_value
            if (stat == 0) this%trip_amp = real_value
        case ("trip_ts")
            read(value, *, iostat=stat) real_value
            if (stat == 0) this%trip_ts = real_value
        case ("trip_nmodes")
            read(value, *, iostat=stat) int_value
            if (stat == 0) this%trip_nmodes = int(int_value, C_INT)
        case ("trip_seed")
            read(value, *, iostat=stat) int_value
            if (stat == 0) this%trip_seed = int(int_value, C_INT)
        case ("stats_sample_interval")
            read(value, *, iostat=stat) int_value
            if (stat == 0) this%stats%sample_interval = int_value
        case ("stats_write_interval")
            read(value, *, iostat=stat) int_value
            if (stat == 0) this%stats%write_interval = int_value
        case ("runtime_interval")
            read(value, *, iostat=stat) int_value
            if (stat == 0) this%stats%runtime_interval = int_value
        case ("stats_file")
            this%stats%file = clean_config_string(value)
        case ("runtime_file")
            this%stats%runtime_file = clean_config_string(value)
        case default
            if (terminal) print *, "warning: unknown boundaryLayer case key on input line", line_no, ": ", trim(key)
        end select
    end subroutine apply_bl_case_value

    logical function read_bool(value, out) result(ok)
        character(len=*), intent(in) :: value
        logical, intent(out) :: out
        character(len=:), allocatable :: v
        v = to_lower(clean_config_string(value))
        ok = .true.
        select case (trim(v))
        case ("true", ".true.", "1", "yes", "on")
            out = .true.
        case ("false", ".false.", "0", "no", "off")
            out = .false.
        case default
            ok = .false.
        end select
    end function read_bool

end module boundarylayer_flow
