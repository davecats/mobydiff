module config
    use, intrinsic :: iso_c_binding
    use :: init, only: dns_type, grid_type, VAR_U, VAR_V, VAR_W, VAR_P, &
        GRID_UNIFORM, GRID_COSINE, GRID_TANH, GRID_NATURAL, config_seen_type
    use :: turbulence, only: turb_type, TURB_NONE, TURB_LES, TURB_RANS
    use :: les_model, only: les_type, LES_NONE, LES_SMAGORINSKY, LES_WALE
    use :: pressure_solver, only: pressure_solver_type
    use :: boundary, only: boundary_type, boundary_face_id, &
        PATCH_GENERIC, PATCH_WALL
    use :: comm, only: comm_type
    implicit none

    logical, save :: terminal_output = .true.

contains

subroutine read_runtime_config(dns, g, turb, les, ps, bc, c, input_file, has_terminal, seen_config)
    type(dns_type), intent(inout) :: dns
    type(grid_type), intent(inout) :: g
    type(turb_type), intent(inout) :: turb
    type(les_type), intent(inout) :: les
    type(pressure_solver_type), intent(inout) :: ps
    type(boundary_type), intent(inout) :: bc
    type(comm_type), intent(inout) :: c
    character(len=*), intent(in) :: input_file
    logical, intent(in), optional :: has_terminal
    type(config_seen_type), intent(out), optional :: seen_config

    integer :: unit, stat, line_no
    character(len=512) :: line, key, value
    character(len=64) :: section
    logical :: exists
    type(config_seen_type) :: seen

    if (present(has_terminal)) terminal_output = has_terminal

    section = ""

    inquire(file=trim(input_file), exist=exists)
    if (.not. exists) then
        if (terminal_output) print *, "error: input file not found: ", trim(input_file)
        error stop "missing input file"
    end if

    open(newunit=unit, file=trim(input_file), status="old", action="read", iostat=stat)
    if (stat /= 0) then
        if (terminal_output) print *, "error: could not open input file: ", trim(input_file)
        error stop "could not open input file"
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
            ! T1 hook (docs/next_session_iddes.md): the [rans] section's mere
            ! presence builds the SST geometry state at init, even with no
            ! keys, so the wall distance can be validated before the RANS
            ! transport phases land.
            if (trim(section) == "rans") dns%rans_configured = .true.
            cycle
        end if

        call split_key_value(line, key, value)
        if (len_trim(key) == 0) cycle
        call apply_config_value(section, key, value, dns, g, turb, les, ps, bc, c, seen, line_no)
    end do

    close(unit)
    if (present(seen_config)) seen_config = seen

    ! [turbulence] model selects the FAMILY (none|les|rans|iddes); when the
    ! key is absent, a configured [les] model implies the LES family, so inis
    ! predating the [turbulence] section run unchanged. An explicit
    ! model = none wins over a configured [les] section (an off switch).
    if (.not. seen%turbulence_model) then
        turb%model = merge(TURB_LES, TURB_NONE, les%model /= LES_NONE)
    end if

    if (has_restart_file(dns)) then
        if (dns%field_interval < 0) error stop "field interval must be non-negative"
    else
        call validate_runtime_config(dns, g, c, seen)
    end if
    call validate_turbulence_values(turb, les, dns)
end subroutine read_runtime_config

subroutine apply_config_value(section, key, value, dns, g, turb, les, ps, bc, c, seen, line_no)
    character(len=*), intent(in) :: section, key, value
    type(dns_type), intent(inout) :: dns
    type(grid_type), intent(inout) :: g
    type(turb_type), intent(inout) :: turb
    type(les_type), intent(inout) :: les
    type(pressure_solver_type), intent(inout) :: ps
    type(boundary_type), intent(inout) :: bc
    type(comm_type), intent(inout) :: c
    type(config_seen_type), intent(inout) :: seen
    integer, intent(in) :: line_no

    character(len=:), allocatable :: section_l, key_l
    integer(C_INT) :: niter_value

    section_l = lower(trim(section))
    key_l = lower(trim(key))

    select case (section_l)
    case ("grid")
        select case (key_l)
        case ("nx")
            call read_int(value, dns%globalSize(1), line_no)
            seen%size(1) = .true.
        case ("ny")
            call read_int(value, dns%globalSize(2), line_no)
            seen%size(2) = .true.
        case ("nz")
            call read_int(value, dns%globalSize(3), line_no)
            seen%size(3) = .true.
        case ("lx")
            call read_real(value, dns%leng(1), line_no)
            seen%length(1) = .true.
        case ("ly")
            call read_real(value, dns%leng(2), line_no)
            seen%length(2) = .true.
        case ("lz")
            call read_real(value, dns%leng(3), line_no)
            seen%length(3) = .true.
        end select
    case ("grid.x", "grid.y", "grid.z")
        call apply_grid_axis_value(section_l, key_l, value, dns, g, seen, line_no)
    case ("flow")
        select case (key_l)
        case ("re")
            call read_real(value, dns%re, line_no)
            seen%re = .true.
        case ("forcing_x")
            call read_real(value, dns%forcing(1), line_no)
            seen%forcing(1) = .true.
        case ("forcing_y")
            call read_real(value, dns%forcing(2), line_no)
            seen%forcing(2) = .true.
        case ("forcing_z")
            call read_real(value, dns%forcing(3), line_no)
            seen%forcing(3) = .true.
        case ("initial_u")
            call read_real(value, dns%initial_velocity(1), line_no)
        case ("initial_v")
            call read_real(value, dns%initial_velocity(2), line_no)
        case ("initial_w")
            call read_real(value, dns%initial_velocity(3), line_no)
        case ("initial_noise")
            call read_real(value, dns%initial_noise, line_no)
        case ("initial")
            dns%initial = clean_string(value)
        end select
    case ("time")
        select case (key_l)
        case ("dt")
            call read_real(value, dns%dt, line_no)
            seen%dt = .true.
        case ("nsteps")
            call read_int(value, dns%nsteps, line_no)
            seen%nsteps = .true.
        case ("t_final")
            call read_real(value, dns%t_final, line_no)
            seen%t_final = .true.
        case ("cflmax")
            call read_real(value, dns%cflmax, line_no)
            seen%cflmax = .true.
        case ("pecletmax")
            call read_real(value, dns%pecletmax, line_no)
            seen%pecletmax = .true.
        case ("dtmax")
            call read_real(value, dns%dtmax, line_no)
            seen%dtmax = .true.
        end select
    case ("output")
        select case (key_l)
        case ("field_interval")
            call read_int(value, dns%field_interval, line_no)
        case ("field_prefix")
            dns%field_prefix = clean_string(value)
        end select
    case ("restart")
        select case (key_l)
        case ("file")
            dns%restart_file = clean_string(value)
        end select
    case ("pressure")
        select case (key_l)
        case ("niter")
            niter_value = ps%nIter
            call read_int(value, niter_value, line_no)
            if (niter_value >= 0_C_INT) then
                ps%nIter = niter_value
                seen%pressure_niter = .true.
            else
                if (terminal_output) print *, "warning: pressure nIter must be non-negative on input line", line_no
            end if
        case ("sor")
            call read_real(value, ps%omega, line_no)
            seen%pressure_sor = .true.
        case ("accel")
            ! none | jacobi (default) or chebyshev (Chebyshev-Jacobi).
            select case (trim(lower(clean_string(value))))
            case ("chebyshev", "cheb")
                ps%cheb = .true.
            case ("none", "jacobi", "")
                ps%cheb = .false.
            case default
                if (terminal_output) print *, &
                    "warning: unknown pressure accel on input line", line_no, ": ", trim(value)
            end select
        case ("cheb_lmin")
            call read_real(value, ps%chebLmin, line_no)
        case ("cheb_lmax")
            call read_real(value, ps%chebLmax, line_no)
        end select
    case ("ibm")
        select case (key_l)
        case ("enabled")
            call read_bool(value, dns%ibm_enabled, line_no)
        case ("coeff_file")
            dns%ibm_coeff_file = clean_string(value)
        end select
    case ("blocks")
        select case (key_l)
        case ("nb")
            call read_int(value, dns%block_nb, line_no)
        case ("remove_solid")
            call read_bool(value, dns%block_remove_solid, line_no)
        case ("refine")
            ! Repeatable: each occurrence adds one refinement box.
            if (dns%block_refine_nboxes >= int(size(dns%block_refine_box, 2), C_INT)) then
                print *, "too many [blocks] refine boxes at line", line_no
                error stop "too many refinement boxes"
            end if
            dns%block_refine_nboxes = dns%block_refine_nboxes + 1_C_INT
            call read_real6(value, dns%block_refine_box(:, dns%block_refine_nboxes), line_no)
        case ("refine_levels")
            call read_int(value, dns%block_refine_levels, line_no)
        case ("refine_body")
            call read_bool(value, dns%block_refine_body, line_no)
        end select
    case ("force")
        select case (key_l)
        case ("enabled")
            call read_bool(value, dns%force_enabled, line_no)
        case ("type")
            dns%force_type = clean_string(value)
        case ("profile")
            dns%force_profile = clean_string(value)
        case ("amp_x")
            call read_real(value, dns%force_amp(1), line_no)
        case ("amp_y")
            call read_real(value, dns%force_amp(2), line_no)
        case ("amp_z")
            call read_real(value, dns%force_amp(3), line_no)
        case ("wavenumber_x", "k_x")
            call read_real(value, dns%force_wavenumber(1), line_no)
        case ("wavenumber_y", "k_y")
            call read_real(value, dns%force_wavenumber(2), line_no)
        case ("wavenumber_z", "k_z")
            call read_real(value, dns%force_wavenumber(3), line_no)
        case ("dir")
            call read_int(value, dns%force_dir, line_no)
        case ("file")
            dns%force_file = clean_string(value)
        end select
    case ("turbulence")
        call apply_turbulence_value(key_l, value, turb, seen, line_no)
    case ("les")
        call apply_les_value(key_l, value, les, line_no)
    case ("rans")
        call apply_rans_value(key_l, value, dns, line_no)
    case ("mpi")
        call apply_mpi_value(key_l, value, c, line_no)
    case ("boundary")
        call apply_boundary_value(key_l, value, bc, line_no)
    end select
end subroutine apply_config_value

subroutine validate_runtime_config(dns, g, c, seen)
    type(dns_type), intent(in) :: dns
    type(grid_type), intent(in) :: g
    type(comm_type), intent(in) :: c
    type(config_seen_type), intent(in) :: seen

    call validate_dns_values(dns, g)
    if (any(c%dims < 0)) error stop "MPI dimensions must be non-negative"
end subroutine validate_runtime_config

subroutine validate_turbulence_values(turb, les, dns)
    type(turb_type), intent(in) :: turb
    type(les_type), intent(in) :: les
    type(dns_type), intent(in) :: dns

    select case (turb%model)
    case (TURB_NONE, TURB_LES, TURB_RANS)
    case default
        error stop "invalid turbulence model"
    end select
    select case (les%model)
    case (LES_NONE, LES_SMAGORINSKY, LES_WALE)
    case default
        error stop "invalid LES model"
    end select
    if (turb%model == TURB_NONE) return

    if (turb%model == TURB_LES) then
        ! The LES family (and later iddes) needs an SGS kernel from [les].
        if (les%model == LES_NONE) then
            error stop "[turbulence] model = les requires an SGS model in [les]"
        end if
        if (les%model == LES_SMAGORINSKY .and. les%cs < 0.0d0) then
            error stop "LES Smagorinsky constant must be non-negative"
        end if
        if (les%model == LES_WALE .and. les%cw < 0.0d0) then
            error stop "LES WALE constant must be non-negative"
        end if
        if (les%delta_scale <= 0.0d0) error stop "LES delta_scale must be positive"
    end if

    if (turb%model == TURB_RANS) then
        ! T2/T3/T4 (docs/next_session_iddes.md): SST transport with resolved
        ! walls or wall functions, plus the gamma-Re_thetat transition
        ! variant on resolved walls.
        if (.not. dns%rans_configured .or. dns%rans_model == 0_C_INT) then
            error stop "[turbulence] model = rans requires [rans] model = sst"
        end if
        ! gamma-Re_theta needs y+ <~ 1; running it through log wall
        ! functions is meaningless.
        if (dns%rans_transition .and. dns%rans_wall_treatment /= 0_C_INT) then
            error stop "[rans] transition requires wall_treatment = resolved"
        end if
        if (dns%rans_tu <= 0.0d0) error stop "[rans] tu must be positive (percent)"
        if (dns%rans_nut_ratio <= 0.0d0) error stop "[rans] nut_ratio must be positive"
    end if
end subroutine validate_turbulence_values

logical function has_restart_file(dns)
    type(dns_type), intent(in) :: dns

    has_restart_file = len_trim(dns%restart_file) > 0
end function has_restart_file

subroutine validate_dns_values(dns, g)
    type(dns_type), intent(in) :: dns
    type(grid_type), intent(in) :: g

    if (any(dns%globalSize <= 0_C_INT)) error stop "grid sizes must be positive"
    if (any(dns%leng <= 0.0d0)) error stop "domain lengths must be positive"
    if (dns%re <= 0.0d0) error stop "Reynolds number must be positive"
    if (dns%dt <= 0.0d0) error stop "time step must be positive"
    if (dns%nsteps < 0_C_INT) error stop "number of time steps must be non-negative"
    if (dns%t_final < 0.0d0) error stop "final time must be non-negative"
    if (dns%nsteps <= 0_C_INT .and. dns%t_final <= 0.0d0) then
        error stop "either nsteps or t_final must be positive"
    end if
    if (dns%cflmax < 0.0d0) error stop "cflmax must be non-negative"
    if (dns%pecletmax < 0.0d0) error stop "pecletmax must be non-negative"
    if (dns%dtmax <= 0.0d0) error stop "dtmax must be positive"
    if (dns%field_interval < 0) error stop "field interval must be non-negative"
    if (dns%block_nb < 0_C_INT) error stop "block size nb must be non-negative"
    if (dns%block_nb > 0_C_INT) then
        if (dns%block_nb < 4_C_INT) error stop "block size nb must be at least 4"
        if (mod(dns%block_nb, 2_C_INT) /= 0_C_INT) error stop "block size nb must be even (red-black)"
    end if
    if (any(g%distribution < GRID_UNIFORM) .or. any(g%distribution > GRID_NATURAL)) then
        error stop "invalid grid distribution"
    end if
    if (any(g%stretch < 0.0d0)) error stop "grid stretch values must be non-negative"
    if (any(g%natural_dyw_plus <= 0.0d0)) error stop "natural grid dyw_plus values must be positive"
end subroutine validate_dns_values

subroutine apply_mpi_value(key, value, c, line_no)
    character(len=*), intent(in) :: key, value
    type(comm_type), intent(inout) :: c
    integer, intent(in) :: line_no

    select case (trim(key))
    case ("dims")
        call read_integer3(value, c%dims, line_no)
    end select
end subroutine apply_mpi_value

! [turbulence] holds the model FAMILY (and, later, the IDDES blend options);
! the family sub-models keep their own sections: [les] now, [rans] from T2.
subroutine apply_turbulence_value(key, value, turb, seen, line_no)
    character(len=*), intent(in) :: key, value
    type(turb_type), intent(inout) :: turb
    type(config_seen_type), intent(inout) :: seen
    integer, intent(in) :: line_no

    select case (trim(key))
    case ("model")
        call read_turbulence_model(value, turb%model, line_no)
        seen%turbulence_model = .true.
    end select
end subroutine apply_turbulence_value

! [rans] — the k-omega SST sub-model section (docs/next_session_iddes.md).
! T2: model = sst enables the transport equations under
! [turbulence] model = rans; T3: wall_treatment = wall_function;
! T4: transition = true adds the gamma-Re_thetat scalars (resolved walls
! only — transition with wall_function is a hard config error).
subroutine apply_rans_value(key, value, dns, line_no)
    character(len=*), intent(in) :: key, value
    type(dns_type), intent(inout) :: dns
    integer, intent(in) :: line_no

    character(len=:), allocatable :: value_l

    select case (trim(key))
    case ("model")
        value_l = lower(clean_string(value))
        select case (trim(value_l))
        case ("sst")
            dns%rans_model = 1_C_INT
        case ("none")
            dns%rans_model = 0_C_INT
        case default
            if (terminal_output) print *, "error: unknown [rans] model on input line", &
                line_no, ": ", trim(value_l)
            error stop "unknown [rans] model"
        end select
    case ("wall_treatment")
        value_l = lower(clean_string(value))
        select case (trim(value_l))
        case ("resolved")
            dns%rans_wall_treatment = 0_C_INT
        case ("wall_function")
            dns%rans_wall_treatment = 1_C_INT
        case default
            if (terminal_output) print *, "error: unknown [rans] wall_treatment on input line", &
                line_no, ": ", trim(value_l)
            error stop "unknown [rans] wall_treatment"
        end select
    case ("transition")
        call read_bool(value, dns%rans_transition, line_no)
    case ("tu")
        call read_real(value, dns%rans_tu, line_no)
    case ("nut_ratio")
        call read_real(value, dns%rans_nut_ratio, line_no)
    case ("dump_geometry")
        call read_bool(value, dns%rans_dump_geometry, line_no)
    case ("dwall_tol")
        call read_real(value, dns%rans_dwall_tol, line_no)
        if (dns%rans_dwall_tol <= 0.0d0) then
            if (terminal_output) print *, "error: [rans] dwall_tol must be positive"
            error stop
        end if
    case default
        if (terminal_output) print *, "warning: unknown [rans] key on input line", &
            line_no, ": ", trim(key)
    end select
end subroutine apply_rans_value

subroutine apply_les_value(key, value, les, line_no)
    character(len=*), intent(in) :: key, value
    type(les_type), intent(inout) :: les
    integer, intent(in) :: line_no

    select case (trim(key))
    case ("model")
        call read_les_model(value, les%model, line_no)
    case ("cs")
        call read_real(value, les%cs, line_no)
    case ("cw")
        call read_real(value, les%cw, line_no)
    case ("delta_scale")
        call read_real(value, les%delta_scale, line_no)
    case ("ibm_aware")
        call read_bool(value, les%ibm_aware, line_no)
    end select
end subroutine apply_les_value

subroutine apply_grid_axis_value(section, key, value, dns, g, seen, line_no)
    character(len=*), intent(in) :: section, key, value
    type(dns_type), intent(inout) :: dns
    type(grid_type), intent(inout) :: g
    type(config_seen_type), intent(inout) :: seen
    integer, intent(in) :: line_no

    integer :: dir

    dir = grid_axis_index(section)
    if (dir < 1) return

    select case (trim(key))
    case ("distribution")
        call read_grid_distribution(value, g%distribution(dir), line_no)
    case ("stretch")
        call read_real(value, g%stretch(dir), line_no)
    case ("natural_dyw_plus", "dyw_plus", "dy_wall_plus", "dyw+")
        call read_real(value, g%natural_dyw_plus(dir), line_no)
    case ("subdivided")
        call read_bool(value, g%subdivided(dir), line_no)
    case ("n")
        call read_int(value, dns%globalSize(dir), line_no)
        seen%size(dir) = .true.
    case ("length")
        call read_real(value, dns%leng(dir), line_no)
        seen%length(dir) = .true.
    end select
end subroutine apply_grid_axis_value

integer function grid_axis_index(section) result(dir)
    character(len=*), intent(in) :: section

    select case (trim(section))
    case ("grid.x")
        dir = 1
    case ("grid.y")
        dir = 2
    case ("grid.z")
        dir = 3
    case default
        dir = 0
    end select
end function grid_axis_index

subroutine apply_boundary_value(key, value, bc, line_no)
    character(len=*), intent(in) :: key, value
    type(boundary_type), intent(inout) :: bc
    integer, intent(in) :: line_no

    integer :: dir, side, var, face_id
    character(len=16) :: field

    select case (trim(key))
    case ("periodic_x")
        call read_bool(value, bc%isPeriodic(1), line_no)
    case ("periodic_y")
        call read_bool(value, bc%isPeriodic(2), line_no)
    case ("periodic_z")
        call read_bool(value, bc%isPeriodic(3), line_no)
    case ("x_min_patch", "x_max_patch", "y_min_patch", "y_max_patch", &
          "z_min_patch", "z_max_patch")
        ! Domain-face patch type (T5 STEP 0): wall | patch. Meaningful on
        ! non-periodic faces only (validated after parsing, when periodic_*
        ! is final); absent = the historical tangential-Dirichlet inference.
        dir = boundary_direction_index(key(1:1))
        side = boundary_side_index(key(3:5))
        call read_patch_type(value, bc%facePatchType(boundary_face_id(dir, side)), line_no)
    case default
        call parse_boundary_key(key, dir, side, var, field)
        if (dir == 0 .or. side < 0 .or. var < 0) then
            if (terminal_output) print *, "warning: unknown boundary key on input line", line_no, ": ", trim(key)
            return
        end if
        face_id = boundary_face_id(dir, side)

        select case (trim(field))
        case ("type")
            call read_bc_type(value, bc%faceBcType(var,face_id), line_no)
        case ("value")
            call read_real(value, bc%faceBcDefaultValue(var,face_id), line_no)
        case default
            if (terminal_output) print *, "warning: boundary key must end in _type or _value on input line", line_no
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
    case ("x")
        idx = 1
    case ("y")
        idx = 2
    case ("z")
        idx = 3
    case default
        idx = 0
    end select
end function boundary_direction_index

integer function boundary_side_index(token) result(idx)
    character(len=*), intent(in) :: token

    select case (trim(token))
    case ("min")
        idx = 0
    case ("max")
        idx = 1
    case default
        idx = -1
    end select
end function boundary_side_index

integer function boundary_variable_index(token) result(idx)
    character(len=*), intent(in) :: token

    select case (trim(token))
    case ("u")
        idx = int(VAR_U)
    case ("v")
        idx = int(VAR_V)
    case ("w")
        idx = int(VAR_W)
    case ("p")
        idx = int(VAR_P)
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
    case ("true", ".true.", "1", "yes")
        target = .true.
    case ("false", ".false.", "0", "no")
        target = .false.
    case default
        if (terminal_output) print *, "warning: could not parse logical value on input line", line_no
    end select
end subroutine read_bool

subroutine read_patch_type(value, target, line_no)
    character(len=*), intent(in) :: value
    integer(C_INT), intent(inout) :: target
    integer, intent(in) :: line_no

    character(len=:), allocatable :: value_l

    value_l = lower(clean_string(value))
    select case (trim(value_l))
    case ("wall")
        target = PATCH_WALL
    case ("patch", "generic")
        target = PATCH_GENERIC
    case default
        if (terminal_output) print *, "error: patch type must be wall or patch on input line", &
            line_no, ": ", trim(value_l)
        error stop "unknown [boundary] patch type"
    end select
end subroutine read_patch_type

subroutine read_bc_type(value, target, line_no)
    character(len=*), intent(in) :: value
    integer(C_INT), intent(inout) :: target
    integer, intent(in) :: line_no

    integer :: stat, parsed
    character(len=:), allocatable :: value_l

    value_l = lower(clean_string(value))
    select case (trim(value_l))
    case ("dirichlet", "0")
        target = 0_C_INT
    case ("neumann", "1")
        target = 1_C_INT
    case default
        read(value_l, *, iostat=stat) parsed
        if (stat == 0 .and. (parsed == 0 .or. parsed == 1)) then
            target = int(parsed, C_INT)
        else
            if (terminal_output) then
                print *, "warning: boundary type must be dirichlet/0 or neumann/1 on input line", line_no
            end if
        end if
    end select
end subroutine read_bc_type

subroutine read_grid_distribution(value, target, line_no)
    character(len=*), intent(in) :: value
    integer(C_INT), intent(inout) :: target
    integer, intent(in) :: line_no

    character(len=:), allocatable :: value_l

    value_l = lower(clean_string(value))
    select case (trim(value_l))
    case ("uniform")
        target = GRID_UNIFORM
    case ("cosine")
        target = GRID_COSINE
    case ("tanh")
        target = GRID_TANH
    case ("natural", "pirozzoli_orlandi", "pirozzoli-orlandi", "po")
        target = GRID_NATURAL
    case default
        if (terminal_output) print *, "warning: unknown grid distribution on input line", line_no, ": ", trim(value_l)
    end select
end subroutine read_grid_distribution

! Turbulence model FAMILY. rans / iddes are recognized but not implemented
! yet (docs/next_session_iddes.md); the SGS kernel choice lives in [les].
subroutine read_turbulence_model(value, target, line_no)
    character(len=*), intent(in) :: value
    integer(C_INT), intent(inout) :: target
    integer, intent(in) :: line_no

    character(len=:), allocatable :: value_l

    value_l = lower(clean_string(value))
    select case (trim(value_l))
    case ("none", "off", "disabled")
        target = TURB_NONE
    case ("les")
        target = TURB_LES
    case ("rans")
        target = TURB_RANS
    case ("sst")
        if (terminal_output) print *, "error: [turbulence] model selects the family;", &
            " use model = rans here and set model = sst in [rans] (input line", line_no, ")"
        error stop "RANS model belongs in the [rans] section"
    case ("iddes", "sst-iddes")
        if (terminal_output) print *, "error: turbulence model not implemented yet on input line", &
            line_no, ": ", trim(value_l)
        error stop "turbulence model not implemented yet"
    case ("smagorinsky", "smag", "wale")
        if (terminal_output) print *, "error: [turbulence] model selects the family;", &
            " use model = les here and set the SGS model in [les] (input line", line_no, ")"
        error stop "SGS model belongs in the [les] section"
    case default
        if (terminal_output) print *, "warning: unknown turbulence model on input line", &
            line_no, ": ", trim(value_l)
    end select
end subroutine read_turbulence_model

subroutine read_les_model(value, target, line_no)
    character(len=*), intent(in) :: value
    integer(C_INT), intent(inout) :: target
    integer, intent(in) :: line_no

    integer :: stat, parsed
    character(len=:), allocatable :: value_l

    value_l = lower(clean_string(value))
    select case (trim(value_l))
    case ("none", "off", "disabled", "0")
        target = LES_NONE
    case ("smagorinsky", "smag", "1")
        target = LES_SMAGORINSKY
    case ("wale", "2")
        target = LES_WALE
    case default
        read(value_l, *, iostat=stat) parsed
        if (stat == 0 .and. parsed >= LES_NONE .and. parsed <= LES_WALE) then
            target = int(parsed, C_INT)
        else
            if (terminal_output) print *, "warning: unknown LES model on input line", line_no, ": ", trim(value_l)
        end if
    end select
end subroutine read_les_model

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

! Parse a single integer config value. The C_INT target also accepts default
! integer actuals (same kind on this platform), so one routine serves both.
subroutine read_int(value, target, line_no)
    character(len=*), intent(in) :: value
    integer(C_INT), intent(inout) :: target
    integer, intent(in) :: line_no
    integer(C_INT) :: stat, parsed

    read(value, *, iostat=stat) parsed
    if (stat == 0) then
        target = parsed
    else
        if (terminal_output) print *, "warning: could not parse integer on input line", line_no
    end if
end subroutine read_int

subroutine read_real6(value, target, line_no)
    character(len=*), intent(in) :: value
    real(C_DOUBLE), intent(inout) :: target(1:6)
    integer, intent(in) :: line_no
    integer :: stat
    real(C_DOUBLE) :: parsed(6)

    read(value, *, iostat=stat) parsed
    if (stat == 0) then
        target = parsed
    else
        if (terminal_output) print *, "warning: could not parse six real values on input line", line_no
    end if
end subroutine read_real6

subroutine read_integer3(value, target, line_no)
    character(len=*), intent(in) :: value
    integer, intent(inout) :: target(1:3)
    integer, intent(in) :: line_no
    integer :: stat, parsed(3)

    read(value, *, iostat=stat) parsed
    if (stat == 0) then
        target = parsed
    else
        if (terminal_output) print *, "warning: could not parse three integer values on input line", line_no
    end if
end subroutine read_integer3

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
        if (terminal_output) print *, "warning: could not parse real value on input line", line_no
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
