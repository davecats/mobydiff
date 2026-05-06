module config
    use, intrinsic :: iso_c_binding
    use :: init, only: grid_type, finalize_grid
    implicit none

contains

subroutine read_runtime_config(g, output_interval, output_prefix, input_file)
    type(grid_type), intent(inout) :: g
    integer, intent(out) :: output_interval
    character(len=*), intent(out) :: output_prefix
    character(len=*), intent(in) :: input_file

    integer :: unit, stat, line_no
    character(len=512) :: line, key, value
    character(len=64) :: section
    logical :: exists, nsteps_seen, t_final_seen

    output_interval = 1
    output_prefix = "data"
    section = ""
    nsteps_seen = .false.
    t_final_seen = .false.

    inquire(file=trim(input_file), exist=exists)
    if (.not. exists) then
        print *, "input file not found; using defaults: ", trim(input_file)
        call finalize_grid(g)
        return
    end if

    open(newunit=unit, file=trim(input_file), status="old", action="read", iostat=stat)
    if (stat /= 0) then
        print *, "could not open input file; using defaults: ", trim(input_file)
        call finalize_grid(g)
        return
    end if

    line_no = 0
    do
        read(unit, '(A)', iostat=stat) line
        if (stat /= 0) exit
        line_no = line_no + 1

        call strip_comment(line)
        line = adjustl(line)
        if (len_trim(line) == 0) cycle

        if (line(1:1) == "[") then
            call parse_section(line, section)
            cycle
        end if

        call split_key_value(line, key, value)
        if (len_trim(key) == 0) cycle
        call apply_config_value(section, key, value, g, output_interval, output_prefix, &
            nsteps_seen, t_final_seen, line_no)
    end do

    close(unit)

    if (t_final_seen .and. .not. nsteps_seen) then
        g%nsteps = max(1, ceiling(g%t_final/g%dt))
    end if
    call finalize_grid(g)
end subroutine read_runtime_config

subroutine apply_config_value(section, key, value, g, output_interval, output_prefix, &
        nsteps_seen, t_final_seen, line_no)
    character(len=*), intent(in) :: section, key, value
    type(grid_type), intent(inout) :: g
    integer, intent(inout) :: output_interval
    character(len=*), intent(inout) :: output_prefix
    logical, intent(inout) :: nsteps_seen, t_final_seen
    integer, intent(in) :: line_no

    character(len=:), allocatable :: section_l, key_l

    section_l = lower(trim(section))
    key_l = lower(trim(key))

    select case (section_l)
    case ("grid")
        select case (key_l)
        case ("nx")
            call read_integer(value, g%nx, line_no)
        case ("ny")
            call read_integer(value, g%ny, line_no)
        case ("nz")
            call read_integer(value, g%nz, line_no)
        case ("lx")
            call read_real(value, g%lx, line_no)
        case ("ly")
            call read_real(value, g%ly, line_no)
        case ("lz")
            call read_real(value, g%lz, line_no)
        end select
    case ("flow")
        if (key_l == "re") call read_real(value, g%re, line_no)
    case ("time")
        select case (key_l)
        case ("dt")
            call read_real(value, g%dt, line_no)
        case ("nsteps", "steps")
            call read_integer(value, g%nsteps, line_no)
            nsteps_seen = .true.
        case ("t_final")
            call read_real(value, g%t_final, line_no)
            t_final_seen = .true.
        case ("cflmax")
            call read_real(value, g%cflmax, line_no)
        case ("dtmax")
            call read_real(value, g%dtmax, line_no)
        end select
    case ("output")
        select case (key_l)
        case ("interval", "output_interval")
            call read_integer(value, output_interval, line_no)
        case ("prefix")
            output_prefix = clean_string(value)
        end select
    end select
end subroutine apply_config_value

subroutine parse_section(line, section)
    character(len=*), intent(in) :: line
    character(len=*), intent(out) :: section
    integer :: last

    last = index(line, "]")
    if (last > 2) then
        section = lower(trim(line(2:last-1)))
    end if
end subroutine parse_section

subroutine split_key_value(line, key, value)
    character(len=*), intent(in) :: line
    character(len=*), intent(out) :: key, value
    integer :: eq, sep

    key = ""
    value = ""
    eq = index(line, "=")

    if (eq > 0) then
        key = adjustl(line(:eq-1))
        value = adjustl(line(eq+1:))
    else
        sep = scan(line, " "//char(9))
        if (sep > 0) then
            key = adjustl(line(:sep-1))
            value = adjustl(line(sep+1:))
        else
            key = adjustl(line)
        end if
    end if

    key = trim(key)
    value = trim(value)
end subroutine split_key_value

subroutine strip_comment(line)
    character(len=*), intent(inout) :: line
    integer :: semicolon, hash, cut

    semicolon = index(line, ";")
    hash = index(line, "#")
    cut = 0

    if (semicolon > 0) cut = semicolon
    if (hash > 0 .and. (cut == 0 .or. hash < cut)) cut = hash
    if (cut > 0) line(cut:) = ""
end subroutine strip_comment

subroutine read_integer(value, target, line_no)
    character(len=*), intent(in) :: value
    integer, intent(inout) :: target
    integer, intent(in) :: line_no
    integer :: stat, parsed

    read(value, *, iostat=stat) parsed
    if (stat == 0) then
        target = parsed
    else
        print *, "warning: could not parse integer on input line", line_no
    end if
end subroutine read_integer

subroutine read_real(value, target, line_no)
    character(len=*), intent(in) :: value
    real(C_DOUBLE), intent(inout) :: target
    integer, intent(in) :: line_no
    integer :: stat
    real(C_DOUBLE) :: parsed

    read(value, *, iostat=stat) parsed
    if (stat == 0) then
        target = parsed
    else
        print *, "warning: could not parse real value on input line", line_no
    end if
end subroutine read_real

function lower(text) result(out)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: out
    integer :: i, c

    do i = 1, len(text)
        c = iachar(text(i:i))
        if (c >= iachar("A") .and. c <= iachar("Z")) then
            out(i:i) = achar(c + iachar("a") - iachar("A"))
        else
            out(i:i) = text(i:i)
        end if
    end do
end function lower

function clean_string(text) result(out)
    character(len=*), intent(in) :: text
    character(len=len_trim(text)) :: out
    character(len=len_trim(text)) :: tmp
    integer :: n

    tmp = trim(adjustl(text))
    n = len_trim(tmp)
    if (n >= 2) then
        if ((tmp(1:1) == '"' .and. tmp(n:n) == '"') .or. &
            (tmp(1:1) == "'" .and. tmp(n:n) == "'")) then
            out = tmp(2:n-1)
            return
        end if
    end if
    out = tmp
end function clean_string

end module config