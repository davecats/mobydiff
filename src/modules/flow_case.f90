module flow_case
    use :: flow_case_base, only: case_type
    use :: generic_flow, only: create_generic_case, GENERIC_CASE_NAME
    use :: channel_flow, only: create_channel_case, CHANNEL_CASE_NAME
    use :: case_config_helpers, only: strip_config_comment, parse_config_section, &
        split_config_key_value, to_lower, clean_config_string
    implicit none

    private

    public :: case_type, create_flow_case

contains

    subroutine create_flow_case(flow, input_file, has_terminal)
        class(case_type), allocatable, intent(out) :: flow
        character(len=*), intent(in) :: input_file
        logical, intent(in), optional :: has_terminal

        character(len=64) :: case_name

        case_name = "generic"
        call read_case_name(input_file, case_name)

        select case (trim(case_name))
        case (CHANNEL_CASE_NAME)
            call create_channel_case(flow)
        case (GENERIC_CASE_NAME, "")
            call create_generic_case(flow)
        case default
            if (present(has_terminal)) then
                if (has_terminal) print *, "warning: unknown case name, using generic case: ", trim(case_name)
            end if
            call create_generic_case(flow)
        end select

        call flow%read_config(input_file, has_terminal)
    end subroutine create_flow_case

    subroutine read_case_name(input_file, case_name)
        character(len=*), intent(in) :: input_file
        character(len=*), intent(inout) :: case_name

        integer :: unit, stat
        character(len=512) :: line, key, value
        character(len=64) :: section
        logical :: exists

        section = ""
        inquire(file=trim(input_file), exist=exists)
        if (.not. exists) return

        open(newunit=unit, file=trim(input_file), status="old", action="read", iostat=stat)
        if (stat /= 0) return

        do
            read(unit, '(A)', iostat=stat) line
            if (stat /= 0) exit
            call strip_config_comment(line)
            line = adjustl(line)
            if (len_trim(line) == 0) cycle
            if (line(1:1) == "[") then
                call parse_config_section(line, section)
                cycle
            end if
            if (trim(section) /= "case") cycle
            call split_config_key_value(line, key, value)
            if (trim(to_lower(key)) == "name") case_name = to_lower(clean_config_string(value))
        end do

        close(unit)
    end subroutine read_case_name

end module flow_case
