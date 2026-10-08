module flow_case
    use :: flow_case_base, only: case_type
    use :: generic_flow, only: create_generic_case, GENERIC_CASE_NAME
    use :: channel_flow, only: create_channel_case, CHANNEL_CASE_NAME
    use :: airfoil_flow, only: create_airfoil_case, AIRFOIL_CASE_NAME
    use :: boundarylayer_flow, only: create_boundarylayer_case, BL_CASE_NAME
    use :: case_config_helpers, only: next_config_entry, to_lower, clean_config_string
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
        case (AIRFOIL_CASE_NAME)
            call create_airfoil_case(flow)
        case (BL_CASE_NAME)
            call create_boundarylayer_case(flow)
        case (GENERIC_CASE_NAME, "")
            call create_generic_case(flow)
        case default
            if (present(has_terminal)) then
                if (has_terminal) print *, "error: unknown [case] name: ", trim(case_name), &
                    " (channel, airfoil, boundaryLayer or generic)"
            end if
            error stop
        end select

        call flow%read_config(input_file, has_terminal)
    end subroutine create_flow_case

    subroutine read_case_name(input_file, case_name)
        character(len=*), intent(in) :: input_file
        character(len=*), intent(inout) :: case_name

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
            if (trim(section) /= "case") cycle
            if (trim(to_lower(key)) == "name") case_name = to_lower(clean_config_string(value))
        end do
        close(unit)

        ! A [case.<other>] section is a configuration error, not a no-op: an
        ! ini carrying [case.boundaryLayer] without `[case] name = boundaryLayer`
        ! ran the GENERIC case with every trip/BC key silently dropped
        ! (2026-10-06, the gate that found it was itself the victim).
        open(newunit=unit, file=trim(input_file), status="old", action="read", iostat=stat)
        if (stat /= 0) return
        section = ""
        line_no = 0
        do
            call next_config_entry(unit, section, key, value, line_no, ok)
            if (.not. ok) exit
            if (len_trim(section) > 5 .and. to_lower(section(1:5)) == "case.") then
                if (to_lower(trim(section(6:))) /= trim(case_name)) then
                    write(*, '(A,A,A,A,A,I0,A)') "error: section [", trim(section), &
                        "] belongs to another case; [case] name = ", trim(case_name), &
                        " (input line ", line_no, ")"
                    error stop
                end if
            end if
        end do
        close(unit)
    end subroutine read_case_name

end module flow_case
