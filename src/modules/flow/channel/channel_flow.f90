module channel_flow
    use, intrinsic :: iso_c_binding
    use :: flow_case_base, only: case_type
    use :: generic_flow, only: set_generic_defaults
    use :: init, only: dns_type, grid_type, field_type, VAR_U, VAR_V, VAR_W, VAR_P, GRID_NATURAL
    use :: boundary, only: boundary_type, boundary_face_id
    use :: pressure_solver, only: pressure_solver_type
    use :: comm, only: comm_type
    use :: case_config_helpers, only: strip_config_comment, parse_config_section, &
        split_config_key_value, to_lower, clean_config_string
    use :: channel_profile, only: initialise_channel_fields, channel_span_dir
    use :: channel_stats, only: channel_stats_type
    implicit none

    private

    character(len=*), parameter :: CHANNEL_CASE_NAME = "channel"

    type, extends(case_type), public :: channel_case_type
        integer(C_INT) :: stream_dir = 1_C_INT
        integer(C_INT) :: wall_dir = 0_C_INT
        integer(C_INT) :: span_dir = 3_C_INT
        integer :: n_walls = 2
        real(C_DOUBLE) :: natural_blend_index = 40.0d0
        real(C_DOUBLE) :: large_disturbance_amplitude = 1.0d-2
        real(C_DOUBLE) :: small_noise_amplitude = 1.0d-3
        type(channel_stats_type) :: stats
    contains
        procedure :: read_config => channel_read_config
        procedure :: apply_defaults => channel_apply_defaults
        procedure :: setup_after_grid => channel_setup_after_grid
        procedure :: initialise_fields => channel_initialise_fields
        procedure :: after_step => channel_after_step
        procedure :: finalize => channel_finalize
    end type channel_case_type

    public :: create_channel_case, CHANNEL_CASE_NAME

contains

    subroutine create_channel_case(flow)
        class(case_type), allocatable, intent(out) :: flow

        allocate(channel_case_type :: flow)
        flow%name = CHANNEL_CASE_NAME
    end subroutine create_channel_case

    subroutine channel_apply_defaults(this, dns, g, bc, c, ps)
        class(channel_case_type), intent(inout) :: this
        type(dns_type), intent(inout) :: dns
        type(grid_type), intent(inout) :: g
        type(boundary_type), intent(inout) :: bc
        type(comm_type), intent(inout) :: c
        type(pressure_solver_type), intent(inout) :: ps

        call resolve_channel_directions(this)
        call set_generic_defaults(dns, g, bc, c, ps)
        this%name = CHANNEL_CASE_NAME

        dns%leng(this%wall_dir) = merge(2.0d0, 1.0d0, this%n_walls == 2)
        dns%forcing = 0.0d0
        dns%forcing(this%stream_dir) = 1.0d0
        dns%ibm_enabled = .false.

        bc%isPeriodic(this%wall_dir) = .false.
        g%distribution(this%wall_dir) = GRID_NATURAL
        g%stretch(this%wall_dir) = this%natural_blend_index

        call set_channel_wall_bcs(this, bc)
    end subroutine channel_apply_defaults

    subroutine channel_setup_after_grid(this, f, dns, g, bc, c)
        class(channel_case_type), intent(inout) :: this
        type(field_type), intent(inout) :: f
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        type(boundary_type), intent(in) :: bc
        type(comm_type), intent(in) :: c

        call this%stats%setup(dns, g, this%wall_dir, c)
    end subroutine channel_setup_after_grid

    subroutine channel_initialise_fields(this, f, dns, g, bc, c)
        class(channel_case_type), intent(inout) :: this
        type(field_type), intent(inout) :: f
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        type(boundary_type), intent(in) :: bc
        type(comm_type), intent(in) :: c

        call initialise_channel_fields(f, dns, g, this%stream_dir, this%wall_dir, this%n_walls, &
            this%large_disturbance_amplitude, this%small_noise_amplitude)
    end subroutine channel_initialise_fields

    subroutine channel_after_step(this, f, dns, g, c)
        class(channel_case_type), intent(inout) :: this
        type(field_type), intent(inout) :: f
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        type(comm_type), intent(in) :: c

        if (this%stats%interval <= 0) return

        call this%stats%accumulate(f, dns, g, this%wall_dir)
        if (modulo(int(dns%step_current), this%stats%interval) == 0) then
            call this%stats%write(f, dns, g, c, this%stream_dir, this%wall_dir, this%span_dir)
        end if
    end subroutine channel_after_step

    subroutine channel_finalize(this, dns, g, c)
        class(channel_case_type), intent(inout) :: this
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        type(comm_type), intent(in) :: c

        call this%stats%finalize()
    end subroutine channel_finalize

    subroutine set_channel_wall_bcs(this, bc)
        class(channel_case_type), intent(in) :: this
        type(boundary_type), intent(inout) :: bc

        integer :: min_face, max_face

        min_face = boundary_face_id(int(this%wall_dir), 0)
        max_face = boundary_face_id(int(this%wall_dir), 1)

        bc%faceBcType(VAR_U:VAR_W,min_face) = 0_C_INT
        bc%faceBcDefaultValue(VAR_U:VAR_W,min_face) = 0.0d0
        bc%faceBcType(VAR_P,min_face) = 1_C_INT

        if (this%n_walls == 2) then
            bc%faceBcType(VAR_U:VAR_W,max_face) = 0_C_INT
        else
            bc%faceBcType(VAR_U:VAR_W,max_face) = 1_C_INT
            bc%faceBcType(this%wall_dir,max_face) = 0_C_INT
        end if
        bc%faceBcDefaultValue(VAR_U:VAR_W,max_face) = 0.0d0
        bc%faceBcType(VAR_P,max_face) = 1_C_INT
    end subroutine set_channel_wall_bcs

    subroutine resolve_channel_directions(this)
        class(channel_case_type), intent(inout) :: this

        if (this%wall_dir == 0_C_INT) then
            this%wall_dir = merge(2_C_INT, 1_C_INT, this%n_walls == 2)
        end if
        if (this%stream_dir == this%wall_dir) then
            this%stream_dir = merge(1_C_INT, 2_C_INT, this%wall_dir /= 1_C_INT)
        end if
        this%span_dir = channel_span_dir(this%stream_dir, this%wall_dir)
    end subroutine resolve_channel_directions

    subroutine channel_read_config(this, input_file, has_terminal)
        class(channel_case_type), intent(inout) :: this
        character(len=*), intent(in) :: input_file
        logical, intent(in), optional :: has_terminal

        integer :: unit, stat, line_no
        character(len=512) :: line, key, value
        character(len=64) :: section
        logical :: exists

        section = ""
        inquire(file=trim(input_file), exist=exists)
        if (.not. exists) return

        open(newunit=unit, file=trim(input_file), status="old", action="read", iostat=stat)
        if (stat /= 0) return

        line_no = 0
        do
            read(unit, '(A)', iostat=stat) line
            if (stat /= 0) exit
            line_no = line_no + 1
            call strip_config_comment(line)
            line = adjustl(line)
            if (len_trim(line) == 0) cycle
            if (line(1:1) == "[") then
                call parse_config_section(line, section)
                cycle
            end if
            if (trim(section) /= "case.channel") cycle
            call split_config_key_value(line, key, value)
            call apply_channel_case_value(this, to_lower(key), value, line_no, has_terminal)
        end do

        close(unit)
        call resolve_channel_directions(this)
    end subroutine channel_read_config

    subroutine apply_channel_case_value(this, key, value, line_no, has_terminal)
        type(channel_case_type), intent(inout) :: this
        character(len=*), intent(in) :: key, value
        integer, intent(in) :: line_no
        logical, intent(in), optional :: has_terminal

        integer :: int_value, stat
        real(C_DOUBLE) :: real_value
        logical :: terminal

        terminal = .true.
        if (present(has_terminal)) terminal = has_terminal

        select case (trim(key))
        case ("n_walls")
            read(value, *, iostat=stat) int_value
            if (stat == 0 .and. (int_value == 1 .or. int_value == 2)) this%n_walls = int_value
        case ("stream_dir", "stream_direction")
            read(value, *, iostat=stat) int_value
            if (stat == 0 .and. int_value >= 1 .and. int_value <= 3) this%stream_dir = int(int_value, C_INT)
        case ("wall_dir", "wall_direction")
            read(value, *, iostat=stat) int_value
            if (stat == 0 .and. int_value >= 1 .and. int_value <= 3) this%wall_dir = int(int_value, C_INT)
        case ("stats_interval")
            read(value, *, iostat=stat) int_value
            if (stat == 0) this%stats%interval = int_value
        case ("stats_file")
            this%stats%file = clean_config_string(value)
        case ("natural_blend_index", "jb")
            read(value, *, iostat=stat) real_value
            if (stat == 0) this%natural_blend_index = real_value
        case ("large_disturbance_amplitude")
            read(value, *, iostat=stat) real_value
            if (stat == 0) this%large_disturbance_amplitude = real_value
        case ("small_noise_amplitude")
            read(value, *, iostat=stat) real_value
            if (stat == 0) this%small_noise_amplitude = real_value
        case default
            if (terminal) print *, "warning: unknown channel case key on input line", line_no, ": ", trim(key)
        end select
    end subroutine apply_channel_case_value

end module channel_flow
