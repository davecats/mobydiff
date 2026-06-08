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
    use :: channel_profile, only: initialise_channel_fields
    use :: channel_stats, only: channel_stats_type
    implicit none

    private

    character(len=*), parameter :: CHANNEL_CASE_NAME = "channel"
    integer(C_INT), parameter :: CHANNEL_STREAM_DIR = 1_C_INT
    integer(C_INT), parameter :: CHANNEL_WALL_DIR = 2_C_INT

    type, extends(case_type), public :: channel_case_type
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

        call set_generic_defaults(dns, g, bc, c, ps)
        this%name = CHANNEL_CASE_NAME

        dns%leng(CHANNEL_WALL_DIR) = merge(2.0d0, 1.0d0, this%n_walls == 2)
        dns%forcing = 0.0d0
        dns%forcing(CHANNEL_STREAM_DIR) = 1.0d0
        dns%ibm_enabled = .false.

        bc%isPeriodic(CHANNEL_WALL_DIR) = .false.
        g%distribution(CHANNEL_WALL_DIR) = GRID_NATURAL
        g%stretch(CHANNEL_WALL_DIR) = this%natural_blend_index

        call set_channel_wall_bcs(this, bc)
    end subroutine channel_apply_defaults

    subroutine channel_setup_after_grid(this, f, dns, g, bc, c)
        class(channel_case_type), intent(inout) :: this
        type(field_type), intent(inout) :: f
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        type(boundary_type), intent(in) :: bc
        type(comm_type), intent(in) :: c

        call this%stats%setup(dns, g, c)
    end subroutine channel_setup_after_grid

    subroutine channel_initialise_fields(this, f, dns, g, bc, c)
        class(channel_case_type), intent(inout) :: this
        type(field_type), intent(inout) :: f
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        type(boundary_type), intent(in) :: bc
        type(comm_type), intent(in) :: c

        call initialise_channel_fields(f, dns, g, this%n_walls, &
            this%large_disturbance_amplitude, this%small_noise_amplitude)
    end subroutine channel_initialise_fields

    subroutine channel_after_step(this, f, dns, g, c)
        class(channel_case_type), intent(inout) :: this
        type(field_type), intent(inout) :: f
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        type(comm_type), intent(in) :: c

        integer :: runtime_interval, sample_interval, write_interval
        logical :: sample_stats, write_hdf5, write_runtime

        runtime_interval = this%stats%runtime_interval
        sample_interval = this%stats%sample_interval
        write_interval = this%stats%write_interval

        sample_stats = sample_interval > 0 .and. modulo(int(dns%step_current), sample_interval) == 0
        write_hdf5 = write_interval > 0 .and. modulo(int(dns%step_current), write_interval) == 0
        write_runtime = runtime_interval > 0 .and. modulo(int(dns%step_current), runtime_interval) == 0

        if (.not. (sample_stats .or. write_hdf5 .or. write_runtime)) return

        if (sample_stats) call this%stats%accumulate(f, dns, g)
        if (write_hdf5 .or. write_runtime) then
            call this%stats%write(f, dns, g, c, write_hdf5, write_runtime)
        end if
    end subroutine channel_after_step

    subroutine channel_finalize(this, dns, g, c)
        class(channel_case_type), intent(inout) :: this
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        type(comm_type), intent(in) :: c

        call this%stats%finalize(dns, c)
    end subroutine channel_finalize

    subroutine set_channel_wall_bcs(this, bc)
        class(channel_case_type), intent(in) :: this
        type(boundary_type), intent(inout) :: bc

        integer :: min_face, max_face

        min_face = boundary_face_id(int(CHANNEL_WALL_DIR), 0)
        max_face = boundary_face_id(int(CHANNEL_WALL_DIR), 1)

        bc%faceBcType(VAR_U:VAR_W,min_face) = 0_C_INT
        bc%faceBcDefaultValue(VAR_U:VAR_W,min_face) = 0.0d0
        bc%faceBcType(VAR_P,min_face) = 1_C_INT

        if (this%n_walls == 2) then
            bc%faceBcType(VAR_U:VAR_W,max_face) = 0_C_INT
        else
            bc%faceBcType(VAR_U:VAR_W,max_face) = 1_C_INT
            bc%faceBcType(VAR_V,max_face) = 0_C_INT
        end if
        bc%faceBcDefaultValue(VAR_U:VAR_W,max_face) = 0.0d0
        bc%faceBcType(VAR_P,max_face) = 1_C_INT
    end subroutine set_channel_wall_bcs

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
        case ("stats_sample_interval")
            read(value, *, iostat=stat) int_value
            if (stat == 0) this%stats%sample_interval = int_value
        case ("stats_write_interval")
            read(value, *, iostat=stat) int_value
            if (stat == 0) this%stats%write_interval = int_value
        case ("stats_file")
            this%stats%file = clean_config_string(value)
        case ("runtime_interval")
            read(value, *, iostat=stat) int_value
            if (stat == 0) this%stats%runtime_interval = int_value
        case ("runtime_file")
            this%stats%runtime_file = clean_config_string(value)
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
