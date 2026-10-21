! The [scalar] / [scalar.N] configuration (split out of scalar.f90 on
! 2026-10-21): the per-key parser config.f90 dispatches to, the one-shot
! array reservation it calls after its pre-scan, and the validation that
! derives dns%nScalar and resolves the conjugate defaults.
module scalar_config
    use, intrinsic :: iso_c_binding
    use :: init, only: dns_type, NVAR
    use :: boundary, only: NFACES, boundary_face_id, BC_DIRICHLET, BC_NEUMANN, PATCH_WALL, &
        PATCH_OUTLET
    use :: scalar_types, only: SC_CONV_ADV, SC_CONV_DIV, SC_CONV_SKEW, SC_IBM_ADIABATIC, &
        SC_IBM_CONJUGATE, SC_IBM_DIRICHLET, SC_INIT_LINEAR_Y, SC_INIT_UNIFORM, &
        SC_LAYOUT_PLANE, SC_LAYOUT_PROFILE, SC_PRT_CONSTANT, SC_PRT_KAYS, &
        SC_SRC_UNIFORM, SC_SRC_VELOCITY, scalar_type
    implicit none
    private
    public :: apply_scalar_config, deallocate_config, scalar_reserve, scalar_section_index, validate_conjugate_config
    public :: validate_scalar_config

contains

    ! Index N of a [scalar.N] section, 0 for anything else (the
    ! grid_axis_index pattern).
    integer function scalar_section_index(section) result(idx)
        character(len=*), intent(in) :: section

        integer :: stat

        idx = 0
        if (len_trim(section) <= 7) return
        if (section(1:7) /= "scalar.") return
        read(section(8:len_trim(section)), *, iostat=stat) idx
        if (stat /= 0) idx = 0
    end function scalar_section_index

    ! One [scalar] (index 0) or [scalar.N] key. The arrays grow on demand, so
    ! there is no maximum scalar count anywhere.
    subroutine apply_scalar_config(index, key, value, sc, line_no, terminal)
        integer, intent(in) :: index
        character(len=*), intent(in) :: key, value
        type(scalar_type), intent(inout) :: sc
        integer, intent(in) :: line_no
        logical, intent(in) :: terminal

        integer :: is

        if (index == 0) then
            select case (trim(key))
            case ("convection")
                sc%convSet = .true.
                select case (trim(value))
                case ("divergence", "div", "")
                    sc%convMode = SC_CONV_DIV
                case ("skew", "skew-symmetric", "skewsymmetric")
                    sc%convMode = SC_CONV_SKEW
                case ("advective", "advection", "nonconservative")
                    sc%convMode = SC_CONV_ADV
                case default
                    if (terminal) print *, "error: [scalar] convection must be", &
                        " divergence, skew or advective, on line", line_no
                    error stop "[scalar] convection"
                end select
            case ("count")
                is = int(read_int_value(value, key, line_no))
                if (is < 0) error stop "[scalar] count must be non-negative"
                sc%countKey = int(is, C_INT)
                call scalar_reserve(sc, is)
            ! Statistics (increment S4, scalar_stats.f90). Off unless asked.
            case ("stats_sample_interval")
                sc%statsSample = read_int_value(value, key, line_no)
            case ("stats_write_interval")
                sc%statsWrite = read_int_value(value, key, line_no)
            case ("stats_file")
                sc%statsFile = trim(value)
            case ("stats_layout")
                select case (trim(value))
                case ("profile", "channel")
                    sc%statsLayout = SC_LAYOUT_PROFILE
                case ("plane", "boundary_layer", "xy")
                    sc%statsLayout = SC_LAYOUT_PLANE
                case default
                    if (terminal) print *, "error: [scalar] stats_layout must be profile or", &
                        " plane, input line", line_no
                    error stop "unknown [scalar] stats_layout"
                end select
            case ("heat_interval")
                sc%heatInterval = read_int_value(value, key, line_no)
            ! Conjugate cut-face error indicator (increment C2). Needs a
            ! conjugate scalar; harmless (never called) otherwise.
            case ("indicator_interval")
                sc%indicatorInterval = read_int_value(value, key, line_no)
            case ("heat_file")
                sc%heatFile = trim(value)
            case default
                if (terminal) print *, "error: unknown [scalar] key on input line", &
                    line_no, ": ", trim(key)
                error stop
            end select
            return
        end if

        if (index < 1) error stop "[scalar.N] section index must be >= 1"
        call scalar_reserve(sc, index)
        sc%sectionSeen(index) = .true.
        is = index

        select case (trim(key))
        case ("name")
            sc%name(is) = trim(value)
        case ("pr", "prandtl", "schmidt", "sc")
            sc%pr(is) = read_real_value(value, key, line_no)
        case ("prt")
            sc%prt(is) = read_real_value(value, key, line_no)
        case ("prt_model")
            select case (trim(value))
            case ("constant")
                sc%prtModel(is) = SC_PRT_CONSTANT
            case ("kays", "kays_crawford", "kays-crawford")
                sc%prtModel(is) = SC_PRT_KAYS
            case default
                if (terminal) print *, "error: [scalar.N] prt_model must be constant or kays,", &
                    " input line", line_no
                error stop "unknown [scalar] prt_model"
            end select
        case ("initial", "init")
            sc%initValue(is) = read_real_value(value, key, line_no)
        case ("init_profile")
            select case (trim(value))
            case ("uniform")
                sc%initProfile(is) = SC_INIT_UNIFORM
            case ("linear_y", "linear-y")
                sc%initProfile(is) = SC_INIT_LINEAR_Y
            case default
                if (terminal) print *, "error: [scalar.N] init_profile must be uniform or", &
                    " linear_y, input line", line_no
                error stop "unknown [scalar] init_profile"
            end select
        case ("source")
            sc%source(is) = read_real_value(value, key, line_no)
        case ("source_type")
            select case (trim(value))
            case ("uniform")
                sc%srcType(is) = SC_SRC_UNIFORM
            case ("velocity", "kasagi")
                sc%srcType(is) = SC_SRC_VELOCITY
            case default
                if (terminal) print *, "error: [scalar.N] source_type must be", &
                    " uniform or velocity (= kasagi), input line", line_no
                error stop "unknown [scalar] source_type"
            end select
        case ("source_dir", "source_direction")
            select case (trim(value))
            case ("x")
                sc%srcDir(is) = 1_C_INT
            case ("y")
                sc%srcDir(is) = 2_C_INT
            case ("z")
                sc%srcDir(is) = 3_C_INT
            case default
                if (terminal) print *, "error: [scalar.N] source_dir must be", &
                    " x, y or z, input line", line_no
                error stop "unknown [scalar] source_dir"
            end select
            sc%srcDirSet(is) = .true.
        case ("ibm_wall")
            select case (trim(value))
            case ("dirichlet")
                sc%ibmMode(is) = SC_IBM_DIRICHLET
            case ("adiabatic", "neumann")
                sc%ibmMode(is) = SC_IBM_ADIABATIC
            case ("conjugate")
                sc%ibmMode(is) = SC_IBM_CONJUGATE
            case default
                if (terminal) print *, "error: [scalar.N] ibm_wall must be dirichlet,", &
                    " adiabatic or conjugate, input line", line_no
                error stop "unknown [scalar] ibm_wall"
            end select
        case ("ibm_value")
            sc%ibmValue(is) = read_real_value(value, key, line_no)
            sc%ibmValueSet(is) = .true.
        ! Conjugate solid properties (increment C1). All rejected outright on
        ! a non-conjugate scalar -- see validate_scalar_config.
        case ("solid_k", "solid_kappa")
            sc%solidK(is) = read_real_value(value, key, line_no)
            sc%solidKeySet(is) = .true.
        case ("solid_rhocp", "solid_capacity")
            sc%solidC(is) = read_real_value(value, key, line_no)
            sc%solidKeySet(is) = .true.
        case ("solid_init", "solid_initial")
            sc%solidInit(is) = read_real_value(value, key, line_no)
            sc%solidInitSet(is) = .true.
            sc%solidKeySet(is) = .true.
        case ("solid_source")
            sc%solidSource(is) = read_real_value(value, key, line_no)
            sc%solidKeySet(is) = .true.
        case ("contact_resistance")
            sc%contactR(is) = read_real_value(value, key, line_no)
            sc%solidKeySet(is) = .true.
        ! The solid's second material band (F2). Also rejected outright on a
        ! non-conjugate scalar, by the same solidKeySet guard.
        case ("solid_thickness")
            sc%solidDepth(is) = read_real_value(value, key, line_no)
            sc%solidKeySet(is) = .true.
        case ("solid_outer_k", "solid_outer_kappa")
            sc%outerK(is) = read_real_value(value, key, line_no)
            sc%solidKeySet(is) = .true.
        case ("solid_outer_rhocp", "solid_outer_capacity")
            sc%outerC(is) = read_real_value(value, key, line_no)
            sc%solidKeySet(is) = .true.
        ! The C2 escalation term. It rides solidKeySet, so it is rejected on
        ! a non-conjugate scalar by the same guard as the solid properties.
        case ("tangential_correction")
            sc%tangCorr(is) = merge(1_C_INT, 0_C_INT, &
                read_logical_value(value, key, line_no))
            sc%solidKeySet(is) = .true.
        case ("inlet")
            sc%inlet(is) = read_real_value(value, key, line_no)
        case default
            call apply_scalar_face_key(is, key, value, sc, line_no, terminal)
        end select
    end subroutine apply_scalar_config

    ! Per-face rows: <dir>_<side>_type / _value, the [boundary] key style.
    subroutine apply_scalar_face_key(is, key, value, sc, line_no, terminal)
        integer, intent(in) :: is
        character(len=*), intent(in) :: key, value
        type(scalar_type), intent(inout) :: sc
        integer, intent(in) :: line_no
        logical, intent(in) :: terminal

        integer :: p1, p2, dir, side, face_id

        p1 = index(key, "_")
        p2 = 0
        if (p1 > 1) p2 = index(key(p1+1:), "_")
        if (p1 <= 1 .or. p2 <= 1) then
            if (terminal) print *, "error: unknown [scalar.N] key on input line", &
                line_no, ": ", trim(key)
            error stop
        end if
        p2 = p1 + p2

        select case (key(:p1-1))
        case ("x"); dir = 1
        case ("y"); dir = 2
        case ("z"); dir = 3
        case default; dir = 0
        end select
        select case (key(p1+1:p2-1))
        case ("min"); side = 0
        case ("max"); side = 1
        case default; side = -1
        end select
        if (dir == 0 .or. side < 0) then
            if (terminal) print *, "error: unknown [scalar.N] key on input line", &
                line_no, ": ", trim(key)
            error stop
        end if
        face_id = boundary_face_id(dir, side)

        select case (trim(key(p2+1:)))
        case ("type")
            select case (trim(value))
            case ("dirichlet", "0")
                sc%bcType(is,face_id) = BC_DIRICHLET
            case ("neumann", "1")
                sc%bcType(is,face_id) = BC_NEUMANN
            case default
                if (terminal) print *, "error: [scalar.N] boundary type must be dirichlet", &
                    " or neumann, input line", line_no
                error stop "unknown [scalar] boundary type"
            end select
            sc%bcTypeSet(is,face_id) = .true.
        case ("value")
            sc%bcValue(is,face_id) = read_real_value(value, key, line_no)
            sc%bcValueSet(is,face_id) = .true.
        case default
            if (terminal) print *, "error: [scalar.N] boundary key must end in _type or", &
                " _value, input line", line_no
            error stop
        end select
    end subroutine apply_scalar_face_key

    ! Grow the per-scalar arrays to at least n entries, preserving values
    ! (host-only ini-parse cost; the count is unknown until the file ends).
    ! Size the per-scalar arrays ONCE for n scalars, with every default
    ! filled in. config.f90 calls it before the parse with the count its
    ! pre-scan found (scalar_count_hint), so the two calls the parser still
    ! makes (at [scalar] count and at each [scalar.N] header) are no-ops; a
    ! larger n after that is a pre-scan bug and stops.
    subroutine scalar_reserve(sc, n)
        type(scalar_type), intent(inout) :: sc
        integer, intent(in) :: n

        if (n <= int(sc%n)) return
        if (sc%n > 0_C_INT) error stop "[scalar]: more scalars than the config pre-scan reserved"

        allocate(sc%pr(n), sc%prt(n), sc%prtModel(n))
        allocate(sc%source(n), sc%initValue(n), sc%ibmValue(n), sc%inlet(n))
        allocate(sc%ibmMode(n), sc%initProfile(n), sc%srcType(n), sc%srcDir(n))
        allocate(sc%srcDirSet(n))
        allocate(sc%solidK(n), sc%solidC(n), sc%solidInit(n), sc%solidSource(n))
        allocate(sc%contactR(n), sc%tangCorr(n))
        allocate(sc%solidDepth(n), sc%outerK(n), sc%outerC(n))
        allocate(sc%solidKeySet(n), sc%solidInitSet(n), sc%ibmValueSet(n))
        allocate(sc%bcType(n,NFACES), sc%bcValue(n,NFACES))
        allocate(sc%bcTypeSet(n,NFACES), sc%bcValueSet(n,NFACES))
        allocate(sc%name(n), sc%sectionSeen(n))

        ! Defaults. Neumann 0 on every face = adiabatic / zero gradient, which
        ! is also what the PATCH_WALL / PATCH_OUTLET / generic faces resolve to.
        sc%pr = 1.0d0
        sc%prt = 0.85d0
        sc%prtModel = SC_PRT_CONSTANT
        sc%source = 0.0d0
        sc%srcType = SC_SRC_UNIFORM
        sc%srcDir = 1_C_INT
        sc%srcDirSet = .false.
        sc%initValue = 0.0d0
        sc%ibmValue = 0.0d0
        sc%inlet = 0.0d0
        sc%ibmMode = SC_IBM_DIRICHLET
        sc%initProfile = SC_INIT_UNIFORM
        ! kappa_s = C_s = 1 is the inert conjugate solid (same material as the
        ! fluid), which is also the k_s = k_f case where the scheme's dropped
        ! tangential term vanishes identically -- a useful control.
        sc%solidK = 1.0d0
        sc%solidC = 1.0d0
        sc%solidInit = 0.0d0
        sc%solidSource = 0.0d0
        sc%contactR = 0.0d0
        sc%tangCorr = 0_C_INT
        ! No band: one solid material, the C1/C3 arithmetic unchanged. When a
        ! band IS asked for, its default material is an insulator with a
        ! fluid-like capacity -- the adiabatic outer surface of F3, and a
        ! capacity that never multiplies a non-zero rhs there.
        sc%solidDepth = 0.0d0
        sc%outerK = 0.0d0
        sc%outerC = 1.0d0
        sc%solidKeySet = .false.
        sc%solidInitSet = .false.
        sc%ibmValueSet = .false.
        sc%bcType = BC_NEUMANN
        sc%bcValue = 0.0d0
        sc%bcTypeSet = .false.
        sc%bcValueSet = .false.
        sc%name = ""
        sc%sectionSeen = .false.

        sc%n = int(n, C_INT)
    end subroutine scalar_reserve

    ! Post-parse validation: contiguous sections, a count that agrees with
    ! the highest section index, positive Prandtl numbers and dataset names
    ! that do not collide with the reserved ones. Also derives dns%nVar.
    subroutine validate_scalar_config(sc, dns, terminal)
        type(scalar_type), intent(inout) :: sc
        type(dns_type), intent(inout) :: dns
        logical, intent(in) :: terminal

        integer :: is, js

        if (sc%countKey >= 0_C_INT .and. int(sc%countKey) /= int(sc%n)) then
            if (terminal) print *, "error: [scalar] count =", sc%countKey, &
                " disagrees with the highest [scalar.N] section", sc%n
            error stop "[scalar] count does not match the [scalar.N] sections"
        end if
        do is = 1, int(sc%n)
            if (.not. sc%sectionSeen(is)) then
                if (terminal) print *, "error: [scalar.", is, "] section missing;", &
                    " scalar sections must be contiguous from 1"
                error stop "missing [scalar.N] section"
            end if
            if (sc%pr(is) <= 0.0d0) error stop "[scalar.N] pr must be positive"
            if (sc%prt(is) <= 0.0d0) error stop "[scalar.N] prt must be positive"
            ! source_dir only means anything to the VELOCITY source. Rejected
            ! rather than ignored: naming a direction and getting the default
            ! one is exactly the silent wrong answer this campaign cannot
            ! afford (a pipe along z driven by u_x heats nothing).
            if (sc%srcDirSet(is) .and. sc%srcType(is) /= SC_SRC_VELOCITY) then
                if (terminal) print '(a,i0,a)', " error: [scalar.", is, &
                    "] source_dir needs source_type = velocity (= kasagi)"
                error stop "[scalar.N] source_dir without source_type = velocity"
            end if
            if (len_trim(sc%name(is)) == 0) write(sc%name(is), '(a,i0)') "s", is
            if (name_is_reserved(sc%name(is))) then
                if (terminal) print *, "error: [scalar.", is, "] name '", &
                    trim(sc%name(is)), "' is a reserved dataset name"
                error stop "[scalar.N] name collides with a reserved dataset"
            end if
            do js = 1, is - 1
                if (sc%name(js) == sc%name(is)) error stop "[scalar.N] duplicate name"
            end do
        end do

        call validate_conjugate_config(sc, dns, terminal)

        dns%nScalar = sc%n
        dns%nVar = NVAR + sc%n
    end subroutine validate_scalar_config

    ! Conjugate body mode (increment C1): the per-scalar solid properties and
    ! the case-level preconditions of docs/history/next_session_conjugate.md Section
    ! 12. Every one of these is a HARD ERROR rather than a warning: each
    ! silently produces a physically wrong answer that no gate in the passive
    ! set would flag.
    subroutine validate_conjugate_config(sc, dns, terminal)
        type(scalar_type), intent(inout) :: sc
        type(dns_type), intent(in) :: dns
        logical, intent(in) :: terminal

        integer :: is, nConj, nTang, nBand, nPlain

        nConj = 0
        nTang = 0
        nBand = 0
        do is = 1, int(sc%n)
            if (sc%ibmMode(is) /= SC_IBM_CONJUGATE) then
                if (sc%solidKeySet(is)) then
                    if (terminal) print '(a,i0,a)', " error: [scalar.", is, &
                        "] solid_k / solid_rhocp / solid_init / solid_source /" // &
                        " contact_resistance / solid_thickness / solid_outer_* /" // &
                        " tangential_correction need ibm_wall = conjugate"
                    error stop "[scalar.N] solid property without ibm_wall = conjugate"
                end if
                cycle
            end if
            nConj = nConj + 1
            if (sc%tangCorr(is) /= 0_C_INT) nTang = nTang + 1
            ! The body temperature is an OUTCOME of the conjugate problem,
            ! not an input -- an ibm_value would be silently ignored.
            if (sc%ibmValueSet(is)) then
                if (terminal) print '(a,i0,a)', " error: [scalar.", is, &
                    "] ibm_value is meaningless with ibm_wall = conjugate:" // &
                    " the interface temperature is solved for, not imposed"
                error stop "[scalar.N] ibm_value with ibm_wall = conjugate"
            end if
            if (sc%solidK(is) <= 0.0d0) error stop "[scalar.N] solid_k must be positive"
            if (sc%solidC(is) <= 0.0d0) error stop "[scalar.N] solid_rhocp must be positive"
            if (sc%contactR(is) < 0.0d0) &
                error stop "[scalar.N] contact_resistance must be non-negative"
            ! The second material band (F2). solid_thickness = 0 is "no band",
            ! so a negative one is a typo rather than a disable.
            if (sc%solidDepth(is) < 0.0d0) &
                error stop "[scalar.N] solid_thickness must be non-negative"
            if (sc%outerK(is) < 0.0d0) &
                error stop "[scalar.N] solid_outer_k must be non-negative"
            if (sc%outerC(is) <= 0.0d0) &
                error stop "[scalar.N] solid_outer_rhocp must be positive"
            if (sc%solidDepth(is) > 0.0d0) then
                nBand = nBand + 1
                ! The C2 correction's stencil, its Gershgorin rate and the
                ! interface indicator are all written against the single
                ! kappa_s of the C1 baseline. Combining them with a banded
                ! solid is untested, and C2 ships disabled by measurement, so
                ! the combination is rejected rather than half-implemented.
                if (sc%tangCorr(is) /= 0_C_INT) then
                    if (terminal) print '(a,i0,a)', " error: [scalar.", is, &
                        "] tangential_correction with solid_thickness is not" // &
                        " supported (the C2 correction assumes one solid material)"
                    error stop "[scalar.N] tangential_correction with solid_thickness"
                end if
            end if
            ! Default solid initial value = the scalar's own `initial`, so a
            ! conjugate run with no solid_init starts uniform.
            if (.not. sc%solidInitSet(is)) sc%solidInit(is) = sc%initValue(is)
        end do
        sc%nConjugate = int(nConj, C_INT)
        sc%nTangential = int(nTang, C_INT)
        sc%nBanded = int(nBand, C_INT)
        ! The kernel lists, in scalar order (the order is immaterial to the
        ! result: each scalar writes only its own slots).
        sc%nPlain = sc%n - sc%nConjugate
        if (allocated(sc%plainList)) deallocate(sc%plainList)
        if (allocated(sc%conjList)) deallocate(sc%conjList)
        allocate(sc%plainList(max(1, int(sc%nPlain))), sc%conjList(max(1, nConj)))
        sc%plainList = 0_C_INT
        sc%conjList = 0_C_INT
        nConj = 0
        nPlain = 0
        do is = 1, int(sc%n)
            if (sc%ibmMode(is) == SC_IBM_CONJUGATE) then
                nConj = nConj + 1
                sc%conjList(nConj) = int(is, C_INT)
            else
                nPlain = nPlain + 1
                sc%plainList(nPlain) = int(is, C_INT)
            end if
        end do
        if (nConj == 0) return

        ! The convective form (review finding F2, 2026-09-27). In conjugate
        ! mode convection is hard-masked on every cut face, including the
        ! fluid-centred ones that carry a real fluid velocity, so the five
        ! remaining fluxes of a fluid cut cell do not sum to zero: in
        ! divergence form a UNIFORM scalar drifts there at the rate
        ! s*u_face/h (skew halves it). Only the advective form, whose
        ! div(u) is built from the same masked velocities, preserves it
        ! exactly, and the same quantity then appears as the conservation
        ! defect a cut cell honestly has. So an unset key resolves to
        ! advective, and an explicit divergence is refused rather than run
        ! silently wrong. Non-conjugate scalars keep the divergence default
        ! (the early return above), which is their gated behaviour.
        if (.not. sc%convSet) then
            sc%convMode = SC_CONV_ADV
            if (terminal) print *, " [scalar] convection unset with a conjugate scalar:", &
                " resolved to advective (review finding F2)"
        else if (sc%convMode == SC_CONV_DIV) then
            if (terminal) print *, "error: [scalar] convection = divergence with", &
                " ibm_wall = conjugate: the masked cut-face flux makes a uniform", &
                " scalar drift in every fluid cut cell (numerics review F2)."
            if (terminal) print *, "   Use convection = advective (the conjugate default", &
                " when the key is absent)."
            error stop "[scalar] convection = divergence with a conjugate scalar"
        else if (sc%convMode == SC_CONV_SKEW .and. terminal) then
            print *, " warning: [scalar] convection = skew with a conjugate scalar", &
                " leaves half of the F2 cut-cell drift; advective removes it"
        end if

        if (.not. dns%ibm_enabled) then
            if (terminal) print *, "error: [scalar.N] ibm_wall = conjugate needs an", &
                " immersed body ([ibm] enabled = true)"
            error stop "[scalar] conjugate without an immersed body"
        end if
        ! The solid now carries a real temperature field, so a block buried
        ! inside the body is part of the SOLUTION domain. Removing it deletes
        ! the solid -- the same class of trap as the A3 penalization-force
        ! finding (validation/naca0012/README.md, in git history).
        if (all(dns%block_nb > 0_C_INT) .and. dns%block_remove_solid &
                .and. .not. dns%block_refine_body) then
            if (terminal) print *, "error: ibm_wall = conjugate needs [blocks]", &
                " remove_solid = false -- buried blocks carry the solid temperature field"
            error stop "[scalar] conjugate with [blocks] remove_solid"
        end if
        if (dns%block_refine_body .and. .not. dns%block_keep_buried) then
            if (terminal) print *, "error: ibm_wall = conjugate needs [blocks]", &
                " keep_buried = true under refine_body -- buried leaves carry", &
                " the solid temperature field (and the case file must be", &
                " prepared with it too)"
            error stop "[scalar] conjugate without [blocks] keep_buried"
        end if
        ! Wall functions model an unresolved fluid sublayer; a conjugate
        ! interface needs the near-wall layer RESOLVED on both sides
        ! (strategy doc Section 6). Not validated -- so, an error.
        if (dns%rans_wall_treatment /= 0_C_INT) then
            if (terminal) print *, "error: ibm_wall = conjugate with [rans]", &
                " wall_treatment = wall_function is not supported"
            error stop "[scalar] conjugate with wall functions"
        end if
        ! The C3 interface-heat diagnostic reports the C1 baseline cut-face
        ! flux. With the C2 tangential correction ON the kernel applies a
        ! DIFFERENT flux at those faces, and a diagnostic that does not report
        ! the flux the kernel applied is worse than no diagnostic (the
        ! invariant the S4 accumulators exist for). The correction ships
        ! disabled by measurement (C2's verdict), so rather than carry a
        ! second, ungated copy of its six-face stencil in the statistics, the
        ! combination is rejected and says so.
        if (sc%nTangential > 0_C_INT .and. sc%heatInterval > 0_C_INT) then
            if (terminal) print *, "error: [scalar] heat_interval with", &
                " [scalar.N] tangential_correction = true: the interface-heat", &
                " diagnostic reports the baseline cut-face flux, which is not", &
                " the one the kernel applies with the correction on"
            error stop "[scalar] heat_interval with tangential_correction"
        end if
    end subroutine validate_conjugate_config

    logical function name_is_reserved(name)
        character(len=*), intent(in) :: name

        select case (trim(name))
        case ("un", "vn", "wn", "pn", "nut", "k", "omega", "gamma", "rethetat", &
              "fd", "vfrac", "x", "y", "z", "blocks")
            name_is_reserved = .true.
        case default
            name_is_reserved = .false.
        end select
    end function name_is_reserved

    real(C_DOUBLE) function read_real_value(value, key, line_no) result(x)
        character(len=*), intent(in) :: value, key
        integer, intent(in) :: line_no

        integer :: stat

        read(value, *, iostat=stat) x
        if (stat /= 0) then
            print *, "error: [scalar] ", trim(key), " needs a real value, input line", line_no
            error stop "could not parse a [scalar] value"
        end if
    end function read_real_value

    integer(C_INT) function read_int_value(value, key, line_no) result(n)
        character(len=*), intent(in) :: value, key
        integer, intent(in) :: line_no

        integer :: stat

        read(value, *, iostat=stat) n
        if (stat /= 0) then
            print *, "error: [scalar] ", trim(key), " needs an integer value, input line", line_no
            error stop "could not parse a [scalar] value"
        end if
    end function read_int_value

    ! The [scalar] spellings of a logical, config.f90's read_bool verbatim --
    ! but a hard error rather than a warning, like the two readers above.
    logical function read_logical_value(value, key, line_no) result(flag)
        character(len=*), intent(in) :: value, key
        integer, intent(in) :: line_no

        select case (trim(adjustl(value)))
        case ("true", ".true.", "1", "yes")
            flag = .true.
        case ("false", ".false.", "0", "no")
            flag = .false.
        case default
            print *, "error: [scalar] ", trim(key), " needs a logical value, input line", line_no
            error stop "could not parse a [scalar] value"
        end select
    end function read_logical_value

    subroutine deallocate_config(sc)
        type(scalar_type), intent(inout) :: sc

        if (allocated(sc%pr)) deallocate(sc%pr)
        if (allocated(sc%prt)) deallocate(sc%prt)
        if (allocated(sc%prtModel)) deallocate(sc%prtModel)
        if (allocated(sc%source)) deallocate(sc%source)
        if (allocated(sc%srcType)) deallocate(sc%srcType)
        if (allocated(sc%srcDir)) deallocate(sc%srcDir)
        if (allocated(sc%srcDirSet)) deallocate(sc%srcDirSet)
        if (allocated(sc%initValue)) deallocate(sc%initValue)
        if (allocated(sc%ibmValue)) deallocate(sc%ibmValue)
        if (allocated(sc%inlet)) deallocate(sc%inlet)
        if (allocated(sc%ibmMode)) deallocate(sc%ibmMode)
        if (allocated(sc%initProfile)) deallocate(sc%initProfile)
        if (allocated(sc%solidK)) deallocate(sc%solidK)
        if (allocated(sc%solidC)) deallocate(sc%solidC)
        if (allocated(sc%solidInit)) deallocate(sc%solidInit)
        if (allocated(sc%solidSource)) deallocate(sc%solidSource)
        if (allocated(sc%contactR)) deallocate(sc%contactR)
        if (allocated(sc%tangCorr)) deallocate(sc%tangCorr)
        if (allocated(sc%plainList)) deallocate(sc%plainList)
        if (allocated(sc%conjList)) deallocate(sc%conjList)
        if (allocated(sc%solidDepth)) deallocate(sc%solidDepth)
        if (allocated(sc%outerK)) deallocate(sc%outerK)
        if (allocated(sc%outerC)) deallocate(sc%outerC)
        if (allocated(sc%solidKeySet)) deallocate(sc%solidKeySet)
        if (allocated(sc%solidInitSet)) deallocate(sc%solidInitSet)
        if (allocated(sc%ibmValueSet)) deallocate(sc%ibmValueSet)
        if (allocated(sc%bcType)) deallocate(sc%bcType)
        if (allocated(sc%bcValue)) deallocate(sc%bcValue)
        if (allocated(sc%bcTypeSet)) deallocate(sc%bcTypeSet)
        if (allocated(sc%bcValueSet)) deallocate(sc%bcValueSet)
        if (allocated(sc%name)) deallocate(sc%name)
        if (allocated(sc%sectionSeen)) deallocate(sc%sectionSeen)
    end subroutine deallocate_config

end module scalar_config
