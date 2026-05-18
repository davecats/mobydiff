module config
    use, intrinsic :: iso_c_binding
    use :: init, only: grid_type, finalize_grid
    use :: pressure_workspace, only: pressure_solver_type
    use :: boundary, only: boundary_type
    implicit none

contains

subroutine read_runtime_config(g, ps, bc, output_interval, output_prefix, field_interval, field_prefix, input_file)
    type(grid_type), intent(inout) :: g
    type(pressure_solver_type), intent(inout) :: ps
    type(boundary_type), intent(inout) :: bc
    integer, intent(out) :: output_interval
    integer, intent(out) :: field_interval
    character(len=*), intent(out) :: output_prefix
    character(len=*), intent(out) :: field_prefix
    character(len=*), intent(in) :: input_file

    integer :: unit, stat, line_no
    character(len=512) :: line, key, value
    character(len=64) :: section
    logical :: exists, nsteps_seen, t_final_seen

    output_interval = 1
    output_prefix = "data"
    field_interval = 0
    field_prefix = "field"
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
        call apply_config_value(section, key, value, g, ps, bc, output_interval, output_prefix, &
            field_interval, field_prefix, &
            nsteps_seen, t_final_seen, line_no)
    end do

    close(unit)

    if (t_final_seen .and. .not. nsteps_seen) then
        g%nsteps = max(1, ceiling(g%t_final/g%dt))
    end if
    call finalize_grid(g)
end subroutine read_runtime_config

subroutine apply_config_value(section, key, value, g, ps, bc, output_interval, output_prefix, &
        field_interval, field_prefix, &
        nsteps_seen, t_final_seen, line_no)
    character(len=*), intent(in) :: section, key, value
    type(grid_type), intent(inout) :: g
    type(pressure_solver_type), intent(inout) :: ps
    type(boundary_type), intent(inout) :: bc
    integer, intent(inout) :: output_interval
    integer, intent(inout) :: field_interval
    character(len=*), intent(inout) :: output_prefix
    character(len=*), intent(inout) :: field_prefix
    logical, intent(inout) :: nsteps_seen, t_final_seen
    integer, intent(in) :: line_no

    character(len=:), allocatable :: section_l, key_l
#ifdef USE_REDBLACK
    integer :: niter_value
#endif

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
        select case (key_l)
        case ("re")
            call read_real(value, g%re, line_no)
        case ("forcing_x", "force_x", "body_force_x", "meanpx")
            call read_real(value, g%forcing_x, line_no)
        case ("forcing_y", "force_y", "body_force_y", "meanpy")
            call read_real(value, g%forcing_y, line_no)
        case ("forcing_z", "force_z", "body_force_z", "meanpz")
            call read_real(value, g%forcing_z, line_no)
        end select
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
        case ("field_interval", "hdf5_interval", "h5_interval")
            call read_integer(value, field_interval, line_no)
        case ("field_prefix", "hdf5_prefix", "h5_prefix")
            field_prefix = clean_string(value)
        end select
    case ("pressure")
#ifdef USE_REDBLACK
        select case (key_l)
        case ("niter", "n_iter", "iterations")
            niter_value = int(ps%nIter)
            call read_integer(value, niter_value, line_no)
            if (niter_value > 0) then
                ps%nIter = int(niter_value, C_INT)
            else
                print *, "warning: pressure nIter must be positive on input line", line_no
            end if
        case ("sor", "omega")
            call read_real(value, ps%sor, line_no)
        end select
#endif
    case ("boundary")
        call apply_boundary_value(key_l, value, bc, line_no)
    end select
end subroutine apply_config_value

subroutine apply_boundary_value(key, value, bc, line_no)
    character(len=*), intent(in) :: key, value
    type(boundary_type), intent(inout) :: bc
    integer, intent(in) :: line_no

    integer :: dir, side, var
    character(len=16) :: field

    select case (trim(key))
    case ("periodic_x", "x_periodic", "isperiodic_x", "x_is_periodic")
        call read_bool(value, bc%isPeriodic(1), line_no)
    case ("periodic_y", "y_periodic", "isperiodic_y", "y_is_periodic")
        call read_bool(value, bc%isPeriodic(2), line_no)
    case ("periodic_z", "z_periodic", "isperiodic_z", "z_is_periodic")
        call read_bool(value, bc%isPeriodic(3), line_no)
    case default
        call parse_boundary_key(key, dir, side, var, field)
        if (dir == 0 .or. side < 0 .or. var < 0) then
            print *, "warning: unknown boundary key on input line", line_no, ": ", trim(key)
            return
        end if

        select case (trim(field))
        case ("type")
            call read_bc_type(value, bc%bcType(dir,side,var), line_no)
        case ("value")
            call read_real(value, bc%bcValue(dir,side,var), line_no)
        case default
            print *, "warning: boundary key must end in _type or _value on input line", line_no
        end select
    end select
end subroutine apply_boundary_value

subroutine parse_boundary_key(key, dir, side, var, field)
    character(len=*), intent(in) :: key
    integer, intent(out) :: dir, side, var
    character(len=*), intent(out) :: field

    integer :: p1, p2, p3

    dir = 0
    side = -1
    var = -1
    field = ""

    p1 = index(key, "_")
    if (p1 <= 1) return

    p2 = index(key(p1+1:), "_")
    if (p2 <= 1) return
    p2 = p1 + p2

    p3 = index(key(p2+1:), "_")
    if (p3 <= 1) return
    p3 = p2 + p3

    dir = boundary_direction_index(key(:p1-1))
    side = boundary_side_index(key(p1+1:p2-1))
    var = boundary_variable_index(key(p2+1:p3-1))
    field = trim(key(p3+1:))
end subroutine parse_boundary_key

integer function boundary_direction_index(token) result(idx)
    character(len=*), intent(in) :: token

    select case (trim(token))
    case ("x", "i", "streamwise")
        idx = 1
    case ("y", "j", "wallnormal", "wall_normal")
        idx = 2
    case ("z", "k", "spanwise")
        idx = 3
    case default
        idx = 0
    end select
end function boundary_direction_index

integer function boundary_side_index(token) result(idx)
    character(len=*), intent(in) :: token

    select case (trim(token))
    case ("min", "lo", "low", "lower", "left", "bottom", "front")
        idx = 0
    case ("max", "hi", "high", "upper", "right", "top", "back")
        idx = 1
    case default
        idx = -1
    end select
end function boundary_side_index

integer function boundary_variable_index(token) result(idx)
    character(len=*), intent(in) :: token

    select case (trim(token))
    case ("p", "pressure")
        idx = 0
    case ("u", "un")
        idx = 1
    case ("v", "vn")
        idx = 2
    case ("w", "wn")
        idx = 3
    case default
        idx = -1
    end select
end function boundary_variable_index

subroutine read_bool(value, target, line_no)
    character(len=*), intent(in) :: value
    logical(C_BOOL), intent(inout) :: target
    integer, intent(in) :: line_no

    character(len=:), allocatable :: value_l

    value_l = lower(clean_string(value))
    select case (trim(value_l))
    case ("true", "t", ".true.", "yes", "y", "1", "periodic")
        target = .true.
    case ("false", "f", ".false.", "no", "n", "0", "nonperiodic", "non-periodic")
        target = .false.
    case default
        print *, "warning: could not parse logical value on input line", line_no
    end select
end subroutine read_bool

subroutine read_bc_type(value, target, line_no)
    character(len=*), intent(in) :: value
    integer(C_INT), intent(inout) :: target
    integer, intent(in) :: line_no

    integer :: stat, parsed
    character(len=:), allocatable :: value_l

    value_l = lower(clean_string(value))
    select case (trim(value_l))
    case ("dirichlet", "d", "fixed", "value", "0")
        target = 0_C_INT
    case ("neumann", "n", "gradient", "normal_gradient", "1")
        target = 1_C_INT
    case default
        read(value_l, *, iostat=stat) parsed
        if (stat == 0 .and. (parsed == 0 .or. parsed == 1)) then
            target = int(parsed, C_INT)
        else
            print *, "warning: boundary type must be dirichlet/0 or neumann/1 on input line", line_no
        end if
    end select
end subroutine read_bc_type

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
