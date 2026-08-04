!--------------------------!
!                          !
!     Passive scalars      !
!                          !
!--------------------------!
!
! Advected/diffused cell-centred scalars s_1..s_N (docs/next_session_scalar.md).
!
! DESIGN (increment S0/S1): a scalar is an EXTRA VARIABLE of blk%q, stored at
! the PRESSURE point, so
!   * the metric arrays never grow -- every scalar stencil reads the VAR_P
!     column of blk%x/d1x/...;
!   * the 2:1-interface halo exchange serves scalars for free: comm.f90 treats
!     any variable index > 3 as cell-centred (8-cell restrict average, injected
!     prolong, the (2 coarse + fine)/3 ghost blend), i.e. the pressure
!     treatment, and scalars ride the SAME MPI message as u,v,w,p;
!   * output/restart is one file, one collective write.
! `[scalar] count = 0` is bit-exact with a scalar-free build BY CONSTRUCTION:
! dns%nVar == NVAR reproduces every allocation shape, no kernel here is called,
! apply_bc is untouched and no extra exchange is issued.
!
! LANDMINE: the scalar index DIFFERS between the arrays -- q slot VAR_S0+is,
! qs/oldrhs slot SCR_S0+is (pressure has no scratch plane).
!
! Per-scalar state lives in scalar_type as flat ALLOCATABLE arrays mapped to
! the device once (enter_scalar_data), the turb_type/sst_type/bodyforce_type
! pattern. There is deliberately NO fixed bound on the scalar count.
module scalar
    use, intrinsic :: iso_c_binding
    use :: init, only: dns_type, VAR_U, VAR_V, VAR_W, VAR_P, &
        VAR_S0, SCR_S0, NVAR
    use :: blocks, only: block_set_type, FACE_CLOSED
    use :: boundary, only: boundary_type, NFACES, boundary_face_id, &
        apply_scalar_bc_q, BC_DIRICHLET, BC_NEUMANN, &
        PATCH_GENERIC, PATCH_WALL, PATCH_INLET, PATCH_OUTLET
    use :: comm, only: comm_type, exchange_halos
    use :: turbulence, only: turb_type, turbulence_is_enabled
    implicit none

    private

    ! Initial-condition profiles ([scalar.N] init_profile).
    integer(C_INT), parameter, public :: SC_INIT_UNIFORM = 0_C_INT
    integer(C_INT), parameter, public :: SC_INIT_LINEAR_Y = 1_C_INT
    ! Turbulent-Prandtl models ([scalar.N] prt_model): the per-scalar
    ! constant prt, or the Kays-Crawford correlation (prt_kays below).
    integer(C_INT), parameter, public :: SC_PRT_CONSTANT = 0_C_INT
    integer(C_INT), parameter, public :: SC_PRT_KAYS = 1_C_INT
    ! Immersed-body wall modes ([scalar.N] ibm_wall); used from S3.
    integer(C_INT), parameter, public :: SC_IBM_DIRICHLET = 0_C_INT
    integer(C_INT), parameter, public :: SC_IBM_ADIABATIC = 1_C_INT
    ! Statistics layouts ([scalar] stats_layout, increment S4): wall-normal
    ! rows x-z averaged (the channel form) or rows of the global (x,y) plane
    ! z averaged (the boundary-layer form). See scalar_stats.f90.
    integer(C_INT), parameter, public :: SC_LAYOUT_PROFILE = 0_C_INT
    integer(C_INT), parameter, public :: SC_LAYOUT_PLANE = 1_C_INT

    integer, parameter :: SC_NAME_LEN = 32

    type, public :: scalar_type
        integer(C_INT) :: n = 0_C_INT
        ! Per-scalar configuration (all sized n; host + device).
        real(C_DOUBLE), allocatable :: pr(:), prt(:)
        integer(C_INT), allocatable :: prtModel(:)
        real(C_DOUBLE), allocatable :: source(:), initValue(:), ibmValue(:), inlet(:)
        integer(C_INT), allocatable :: ibmMode(:), initProfile(:)
        ! Per-face boundary rows (n, NFACES): BC_DIRICHLET / BC_NEUMANN with
        ! the face value / normal derivative. They live HERE, not in
        ! boundary_type, so bc keeps its fixed VAR_U:VAR_P shape and the
        ! restart bc_type/bc_value attributes keep their 24-entry layout.
        integer(C_INT), allocatable :: bcType(:,:)
        real(C_DOUBLE), allocatable :: bcValue(:,:)
        logical, allocatable :: bcTypeSet(:,:), bcValueSet(:,:)   ! host only
        ! q variable indices of the scalars, for the batched halo exchange.
        integer(C_INT), allocatable :: varList(:)
        ! Inverse centre-to-centre distances at the p position (the scalar
        ! diffusion's face metric). Mirrors turb%inv_d?(:,VAR_P,:), but the
        ! turbulence module returns early when the model is none, so a DNS
        ! run has no tables -- the scalar module owns its own.
        real(C_DOUBLE), allocatable :: invDx(:,:), invDy(:,:), invDz(:,:)  ! (0:nb+1,nBlocks)
        ! 1-cell dummy standing in for turb%nut when no turbulence model is
        ! active (that array is then not allocated at all): it gives the
        ! transport kernel something to map, and every access to it is
        ! guarded by useNut -- the turb_type "1-cell dummies, uniform device
        ! maps, all accesses model-guarded" idiom.
        real(C_DOUBLE), allocatable :: nutNone(:,:,:,:)
        ! Host-only: dataset names and the config bookkeeping.
        character(len=SC_NAME_LEN), allocatable :: name(:)
        logical, allocatable :: sectionSeen(:)
        integer(C_INT) :: countKey = -1_C_INT    ! [scalar] count, -1 = unset
        ! Host-only statistics configuration ([scalar], increment S4). It is
        ! parsed here so ONE section configures everything about the scalars;
        ! the accumulator itself lives in scalar_stats.f90 (which uses this
        ! module, hence no state of its own can live in scalar_type). All
        ! intervals off by default: nothing in scalar_stats runs unless asked.
        integer(C_INT) :: statsLayout = SC_LAYOUT_PROFILE
        integer(C_INT) :: statsSample = -1_C_INT
        integer(C_INT) :: statsWrite = -1_C_INT
        integer(C_INT) :: heatInterval = -1_C_INT
        character(len=256) :: statsFile = "scalar_stats.h5"
        character(len=256) :: heatFile = "scalar_heat.txt"
    end type scalar_type

    public :: scalars_enabled, scalar_count
    public :: apply_scalar_config, validate_scalar_config, destroy_scalar
    public :: init_scalar, init_scalar_fields, enter_scalar_data, exit_scalar_data
    public :: scalar_sync, scalar_finish, scalar_transport
    public :: scalar_names, scalar_min_pr, scalar_min_prt
    public :: scalar_section_index, prt_kays, eddy_diffusivity

contains

    logical function scalars_enabled(sc)
        type(scalar_type), intent(in) :: sc

        scalars_enabled = sc%n > 0_C_INT
    end function scalars_enabled

    integer(C_INT) function scalar_count(sc) result(n)
        type(scalar_type), intent(in) :: sc

        n = sc%n
    end function scalar_count

    ! Host-only dataset-name table (io labels); empty when no scalars.
    function scalar_names(sc) result(names)
        type(scalar_type), intent(in) :: sc
        character(len=SC_NAME_LEN), allocatable :: names(:)

        if (allocated(sc%name)) then
            names = sc%name(1:int(sc%n))
        else
            allocate(names(0))
        end if
    end function scalar_names

    ! Smallest molecular Prandtl/Schmidt number (the binding Peclet limit).
    real(C_DOUBLE) function scalar_min_pr(sc) result(prmin)
        type(scalar_type), intent(in) :: sc

        prmin = 1.0d0
        if (sc%n > 0_C_INT) prmin = minval(sc%pr(1:int(sc%n)))
    end function scalar_min_pr

    ! Smallest turbulent Prandtl number any scalar can reach (the binding
    ! eddy-diffusivity Peclet limit). Under prt_model = kays the correlation
    ! stays in [prt, 2 prt] -- prt is its Pe_t -> infinity asymptote -- so the
    ! configured prt is the floor for both models.
    real(C_DOUBLE) function scalar_min_prt(sc) result(prtmin)
        type(scalar_type), intent(in) :: sc

        prtmin = 1.0d0
        if (sc%n > 0_C_INT) prtmin = minval(sc%prt(1:int(sc%n)))
    end function scalar_min_prt

    ! Kays-Crawford turbulent Prandtl number (Kays 1994; Weigand, Ferguson &
    ! Crawford 1997), a pointwise function of the turbulent Peclet number
    ! Pe_t = (nu_t/nu) Pr:
    !
    !   1/Pr_t = 1/(2 Prt_inf) + C Pe_t/sqrt(Prt_inf)
    !            - (C Pe_t)^2 [1 - exp(-1/(C Pe_t sqrt(Prt_inf)))],   C = 0.3
    !
    ! Limits: Pe_t -> 0 gives Pr_t = 2 Prt_inf (the molecular limit, inert --
    ! it multiplies nu_t = 0); Pe_t -> infinity gives Pr_t -> Prt_inf, which
    ! is the per-scalar [scalar.N] prt (0.85 by default), so `constant` and
    ! `kays` share one configuration key.
    !
    ! NUMERICS: written with x = 1/(C Pe_t sqrt(Prt_inf)), the bracket is
    ! (C Pe_t)^2 [1 - exp(-x)] = a/b - 1/(2b^2) + ... and the leading a/b
    ! cancels the second term of the sum exactly -- catastrophically so once
    ! Pe_t is large. The small-x branch therefore evaluates the analytically
    ! equivalent series
    !
    !   1/Pr_t = (1/Prt_inf) [ 1 - x/3! + x^2/4! - x^3/5! + ... ]
    !
    ! which is also where the Pe_t -> infinity limit is exact by inspection.
    ! The two branches agree to round-off at the x = 1/2 crossover.
    real(C_DOUBLE) function prt_kays(pet, prtinf) result(prt)
        !$omp declare target
        real(C_DOUBLE), intent(in) :: pet      ! turbulent Peclet number
        real(C_DOUBLE), intent(in) :: prtinf   ! Pr_t at Pe_t -> infinity

        real(C_DOUBLE), parameter :: CKC = 0.3d0
        real(C_DOUBLE) :: a, b, x, u, sum, invprt
        integer :: m

        b = sqrt(prtinf)
        a = CKC*max(pet, 0.0d0)
        if (a*b < 1.0d-300) then
            prt = 2.0d0*prtinf              ! the molecular limit
            return
        end if
        x = 1.0d0/(a*b)

        if (x < 0.5d0) then
            sum = 1.0d0
            u = -x/6.0d0                    ! the m = 1 term, -x/3!
            do m = 1, 16                    ! machine precision for x < 1/2
                sum = sum + u
                u = -u*x/real(m + 3, C_DOUBLE)
            end do
            invprt = sum/prtinf
        else
            invprt = 0.5d0/prtinf + a/b - a*a*(1.0d0 - exp(-x))
        end if

        prt = 1.0d0/invprt
    end function prt_kays

    ! Face eddy diffusivity nu_t/Pr_t for one scalar. Called per face inside
    ! the transport kernel, so it is a declare-target pure function.
    real(C_DOUBLE) function eddy_diffusivity(nutFace, pr, prt, prtModel, re) result(d)
        !$omp declare target
        real(C_DOUBLE), intent(in) :: nutFace, pr, prt, re
        integer(C_INT), intent(in) :: prtModel

        real(C_DOUBLE) :: nt

        nt = max(nutFace, 0.0d0)
        if (prtModel == SC_PRT_KAYS) then
            ! Pe_t = (nu_t/nu) Pr = nu_t Re Pr.
            d = nt/prt_kays(nt*re*pr, prt)
        else
            d = nt/prt
        end if
    end function eddy_diffusivity

    !--------------------------------------------------------------------
    ! Configuration
    !--------------------------------------------------------------------

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
            case ("count")
                is = int(read_int_value(value, key, line_no))
                if (is < 0) error stop "[scalar] count must be non-negative"
                sc%countKey = int(is, C_INT)
                call scalar_grow(sc, is)
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
            case ("heat_file")
                sc%heatFile = trim(value)
            case default
                if (terminal) print *, "warning: unknown [scalar] key on input line", &
                    line_no, ": ", trim(key)
            end select
            return
        end if

        if (index < 1) error stop "[scalar.N] section index must be >= 1"
        call scalar_grow(sc, index)
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
        case ("ibm_wall")
            select case (trim(value))
            case ("dirichlet")
                sc%ibmMode(is) = SC_IBM_DIRICHLET
            case ("adiabatic", "neumann")
                sc%ibmMode(is) = SC_IBM_ADIABATIC
            case default
                if (terminal) print *, "error: [scalar.N] ibm_wall must be dirichlet or", &
                    " adiabatic, input line", line_no
                error stop "unknown [scalar] ibm_wall"
            end select
        case ("ibm_value")
            sc%ibmValue(is) = read_real_value(value, key, line_no)
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
            if (terminal) print *, "warning: unknown [scalar.N] key on input line", &
                line_no, ": ", trim(key)
            return
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
            if (terminal) print *, "warning: unknown [scalar.N] key on input line", &
                line_no, ": ", trim(key)
            return
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
            if (terminal) print *, "warning: [scalar.N] boundary key must end in _type or", &
                " _value, input line", line_no
        end select
    end subroutine apply_scalar_face_key

    ! Grow the per-scalar arrays to at least n entries, preserving values
    ! (host-only ini-parse cost; the count is unknown until the file ends).
    subroutine scalar_grow(sc, n)
        type(scalar_type), intent(inout) :: sc
        integer, intent(in) :: n

        type(scalar_type) :: old
        integer :: nOld

        if (n <= int(sc%n)) return
        nOld = int(sc%n)
        if (nOld > 0) then
            old%pr = sc%pr; old%prt = sc%prt; old%prtModel = sc%prtModel
            old%source = sc%source; old%initValue = sc%initValue
            old%ibmValue = sc%ibmValue; old%inlet = sc%inlet
            old%ibmMode = sc%ibmMode; old%initProfile = sc%initProfile
            old%bcType = sc%bcType; old%bcValue = sc%bcValue
            old%bcTypeSet = sc%bcTypeSet; old%bcValueSet = sc%bcValueSet
            old%name = sc%name; old%sectionSeen = sc%sectionSeen
            call deallocate_config(sc)
        end if

        allocate(sc%pr(n), sc%prt(n), sc%prtModel(n))
        allocate(sc%source(n), sc%initValue(n), sc%ibmValue(n), sc%inlet(n))
        allocate(sc%ibmMode(n), sc%initProfile(n))
        allocate(sc%bcType(n,NFACES), sc%bcValue(n,NFACES))
        allocate(sc%bcTypeSet(n,NFACES), sc%bcValueSet(n,NFACES))
        allocate(sc%name(n), sc%sectionSeen(n))

        ! Defaults. Neumann 0 on every face = adiabatic / zero gradient, which
        ! is also what the PATCH_WALL / PATCH_OUTLET / generic faces resolve to.
        sc%pr = 1.0d0
        sc%prt = 0.85d0
        sc%prtModel = SC_PRT_CONSTANT
        sc%source = 0.0d0
        sc%initValue = 0.0d0
        sc%ibmValue = 0.0d0
        sc%inlet = 0.0d0
        sc%ibmMode = SC_IBM_DIRICHLET
        sc%initProfile = SC_INIT_UNIFORM
        sc%bcType = BC_NEUMANN
        sc%bcValue = 0.0d0
        sc%bcTypeSet = .false.
        sc%bcValueSet = .false.
        sc%name = ""
        sc%sectionSeen = .false.

        if (nOld > 0) then
            sc%pr(1:nOld) = old%pr; sc%prt(1:nOld) = old%prt
            sc%prtModel(1:nOld) = old%prtModel
            sc%source(1:nOld) = old%source; sc%initValue(1:nOld) = old%initValue
            sc%ibmValue(1:nOld) = old%ibmValue; sc%inlet(1:nOld) = old%inlet
            sc%ibmMode(1:nOld) = old%ibmMode; sc%initProfile(1:nOld) = old%initProfile
            sc%bcType(1:nOld,:) = old%bcType; sc%bcValue(1:nOld,:) = old%bcValue
            sc%bcTypeSet(1:nOld,:) = old%bcTypeSet
            sc%bcValueSet(1:nOld,:) = old%bcValueSet
            sc%name(1:nOld) = old%name; sc%sectionSeen(1:nOld) = old%sectionSeen
            call deallocate_config(old)
        end if
        sc%n = int(n, C_INT)
    end subroutine scalar_grow

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

        dns%nScalar = sc%n
        dns%nVar = NVAR + sc%n
    end subroutine validate_scalar_config

    logical function name_is_reserved(name)
        character(len=*), intent(in) :: name

        select case (trim(name))
        case ("un", "vn", "wn", "pn", "nut", "k", "omega", "gamma", "rethetat", &
              "fd", "x", "y", "z", "blocks")
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

    !--------------------------------------------------------------------
    ! Initialisation
    !--------------------------------------------------------------------

    ! Resolve the patch-derived boundary rows and build the per-block metric
    ! tables. Runs after the config/restart metadata is final (the patch types
    ! and the block set must both exist).
    subroutine init_scalar(sc, blk, bc, has_terminal)
        type(scalar_type), intent(inout) :: sc
        type(block_set_type), intent(in) :: blk
        type(boundary_type), intent(in) :: bc
        logical, intent(in) :: has_terminal

        integer :: is, dir, side, face_id, nx, ny, nz, i, b

        if (.not. scalars_enabled(sc)) return

        ! Patch-derived defaults, set-if-unset (the A0 face concept): a wall
        ! is adiabatic, an inlet holds the scalar's inlet value, an outlet
        ! and a generic patch are zero-gradient. BC_OUTFLOW is a
        ! face-staggered concept and does not apply to a cell-centred scalar.
        do dir = 1, 3
            do side = 0, 1
                face_id = boundary_face_id(dir, side)
                select case (bc%facePatchType(face_id))
                case (PATCH_INLET)
                    do is = 1, int(sc%n)
                        call resolve_scalar_row(sc, is, face_id, BC_DIRICHLET, &
                            sc%inlet(is), has_terminal, strict=.true.)
                    end do
                case (PATCH_OUTLET, PATCH_GENERIC)
                    do is = 1, int(sc%n)
                        call resolve_scalar_row(sc, is, face_id, BC_NEUMANN, 0.0d0, &
                            has_terminal, strict=.true.)
                    end do
                case (PATCH_WALL)
                    ! A wall is the ONE face kind where both scalar types are
                    ! physical: adiabatic (Neumann 0, the default) and
                    ! isothermal / fixed-concentration (Dirichlet) are equally
                    ! valid, so an explicit _type key is HONOURED here rather
                    ! than treated as a contradiction (Section 5 of
                    ! docs/next_session_scalar.md: "unless the ini gives a
                    ! type/value"). Inlet and outlet stay strict: an inlet must
                    ! impose a value, an outlet must be zero-gradient.
                    do is = 1, int(sc%n)
                        call resolve_scalar_row(sc, is, face_id, BC_NEUMANN, 0.0d0, &
                            has_terminal, strict=.false.)
                    end do
                end select
            end do
        end do

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))
        if (allocated(sc%invDx)) deallocate(sc%invDx, sc%invDy, sc%invDz)
        allocate(sc%invDx(0:nx+1,blk%nBlocks), sc%invDy(0:ny+1,blk%nBlocks), &
                 sc%invDz(0:nz+1,blk%nBlocks))
        do b = 1, int(blk%nBlocks)
            do i = 0, nx+1
                sc%invDx(i,b) = inv_delta(blk%x(i,VAR_P,b) - blk%x(i-1,VAR_P,b))
            end do
            do i = 0, ny+1
                sc%invDy(i,b) = inv_delta(blk%y(i,VAR_P,b) - blk%y(i-1,VAR_P,b))
            end do
            do i = 0, nz+1
                sc%invDz(i,b) = inv_delta(blk%z(i,VAR_P,b) - blk%z(i-1,VAR_P,b))
            end do
        end do

        if (allocated(sc%nutNone)) deallocate(sc%nutNone)
        allocate(sc%nutNone(0:0,0:0,0:0,1))
        sc%nutNone = 0.0d0

        if (allocated(sc%varList)) deallocate(sc%varList)
        allocate(sc%varList(int(sc%n)))
        do is = 1, int(sc%n)
            sc%varList(is) = VAR_S0 + int(is, C_INT)
        end do

        if (has_terminal) then
            print *, "passive scalars:", sc%n
            do is = 1, int(sc%n)
                print '(a,i0,a,a,a,f8.4,a,es10.3)', "   s", is, " '", trim(sc%name(is)), &
                    "'  pr =", sc%pr(is), "  source =", sc%source(is)
            end do
        end if
    end subroutine init_scalar

    ! Set-if-unset with the A0 contradiction rule: on a strict face an
    ! explicit _type key that disagrees with the declared patch is a hard
    ! config error; on a permissive one (a wall) it simply wins.
    subroutine resolve_scalar_row(sc, is, face_id, want, value, terminal, strict)
        type(scalar_type), intent(inout) :: sc
        integer, intent(in) :: is, face_id
        integer(C_INT), intent(in) :: want
        real(C_DOUBLE), intent(in) :: value
        logical, intent(in) :: terminal
        logical, intent(in) :: strict

        if (sc%bcTypeSet(is,face_id)) then
            if (strict .and. sc%bcType(is,face_id) /= want) then
                if (terminal) print '(a,i0,a,i0)', &
                    " error: explicit [scalar.N] _type key contradicts the declared " // &
                    "patch type on face ", face_id, ", scalar ", is
                error stop "[scalar] BC type contradicts the declared patch type"
            end if
        else
            sc%bcType(is,face_id) = want
        end if
        if (.not. sc%bcValueSet(is,face_id)) sc%bcValue(is,face_id) = value
    end subroutine resolve_scalar_row

    real(C_DOUBLE) function inv_delta(d) result(inv)
        real(C_DOUBLE), intent(in) :: d

        inv = 1.0d0/sign(max(abs(d), 1.0d-300), d)
    end function inv_delta

    ! Initial condition in the q scalar slots (host side, before the device
    ! map). Always called, including on restart: read_field then overwrites
    ! the scalars whose dataset the file carries, so an absent dataset (a
    ! newly added scalar, an older restart) keeps these values -- the RANS
    ! named-scalar restart precedent.
    subroutine init_scalar_fields(sc, dns, blk)
        type(scalar_type), intent(in) :: sc
        type(dns_type), intent(in) :: dns
        type(block_set_type), intent(inout) :: blk

        integer :: is, i, j, k, b, nx, ny, nz
        real(C_DOUBLE) :: s0, s1, y

        if (.not. scalars_enabled(sc)) return

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))

        do is = 1, int(sc%n)
            blk%q(:,:,:,VAR_S0+is,:) = sc%initValue(is)
            if (sc%initProfile(is) /= SC_INIT_LINEAR_Y) cycle
            ! linear_y: the straight line between the y faces' Dirichlet
            ! values (the conduction/heated-channel convenience IC).
            s0 = sc%bcValue(is, boundary_face_id(2,0))
            s1 = sc%bcValue(is, boundary_face_id(2,1))
            do b = 1, int(blk%nBlocks)
                do k = 0, nz+1
                    do j = 0, ny+1
                        y = blk%y(j,VAR_P,b)
                        do i = 0, nx+1
                            blk%q(i,j,k,VAR_S0+is,b) = s0 + (s1 - s0)*y/dns%leng(2)
                        end do
                    end do
                end do
            end do
        end do
    end subroutine init_scalar_fields

    subroutine enter_scalar_data(sc)
        type(scalar_type), intent(inout) :: sc

        if (.not. scalars_enabled(sc)) return

#ifdef USE_OPENMP_OFFLOAD
        !$omp target enter data map(to: sc)
        !$omp target enter data map(to: sc%pr, sc%prt, sc%prtModel, sc%source, &
        !$omp& sc%initValue, sc%ibmValue, sc%inlet, sc%ibmMode, sc%initProfile, &
        !$omp& sc%bcType, sc%bcValue, sc%invDx, sc%invDy, sc%invDz, sc%nutNone)
#endif
    end subroutine enter_scalar_data

    subroutine exit_scalar_data(sc)
        type(scalar_type), intent(inout) :: sc

        if (.not. scalars_enabled(sc)) return

#ifdef USE_OPENMP_OFFLOAD
        !$omp target exit data map(delete: sc%pr, sc%prt, sc%prtModel, sc%source, &
        !$omp& sc%initValue, sc%ibmValue, sc%inlet, sc%ibmMode, sc%initProfile, &
        !$omp& sc%bcType, sc%bcValue, sc%invDx, sc%invDy, sc%invDz, sc%nutNone)
        !$omp target exit data map(delete: sc)
#endif
    end subroutine exit_scalar_data

    subroutine deallocate_config(sc)
        type(scalar_type), intent(inout) :: sc

        if (allocated(sc%pr)) deallocate(sc%pr)
        if (allocated(sc%prt)) deallocate(sc%prt)
        if (allocated(sc%prtModel)) deallocate(sc%prtModel)
        if (allocated(sc%source)) deallocate(sc%source)
        if (allocated(sc%initValue)) deallocate(sc%initValue)
        if (allocated(sc%ibmValue)) deallocate(sc%ibmValue)
        if (allocated(sc%inlet)) deallocate(sc%inlet)
        if (allocated(sc%ibmMode)) deallocate(sc%ibmMode)
        if (allocated(sc%initProfile)) deallocate(sc%initProfile)
        if (allocated(sc%bcType)) deallocate(sc%bcType)
        if (allocated(sc%bcValue)) deallocate(sc%bcValue)
        if (allocated(sc%bcTypeSet)) deallocate(sc%bcTypeSet)
        if (allocated(sc%bcValueSet)) deallocate(sc%bcValueSet)
        if (allocated(sc%name)) deallocate(sc%name)
        if (allocated(sc%sectionSeen)) deallocate(sc%sectionSeen)
    end subroutine deallocate_config

    subroutine destroy_scalar(sc)
        type(scalar_type), intent(inout) :: sc

        call deallocate_config(sc)
        if (allocated(sc%varList)) deallocate(sc%varList)
        if (allocated(sc%nutNone)) deallocate(sc%nutNone)
        if (allocated(sc%invDx)) deallocate(sc%invDx)
        if (allocated(sc%invDy)) deallocate(sc%invDy)
        if (allocated(sc%invDz)) deallocate(sc%invDz)
        sc%n = 0_C_INT
    end subroutine destroy_scalar

    !--------------------------------------------------------------------
    ! Transport
    !--------------------------------------------------------------------

    ! One RK3 substage of every scalar, in ONE fused kernel (collapse(4) over
    ! b,k,j,i with an inner scalar loop, so the GPU sees a single launch
    ! regardless of the scalar count). Reads the START-of-substage q and
    ! writes only the qs/oldrhs scalar slots.
    !
    ! Convection: 2nd-order central in DIVERGENCE form on the p cell's own six
    ! face velocities -- no velocity interpolation is needed on the staggered
    ! mesh, and the flux form telescopes, so sum(s dV) changes only by the
    ! boundary flux.
    !
    ! [flow] convection = skew subtracts s*(div u)|stencil built from the SAME
    ! face velocities (docs/next_session_scalar.md Section 2). NOTE what that
    ! buys, since it differs from the momentum kernel's skew term: the FULL
    ! subtraction is the ADVECTIVE form u.grad s, which preserves a UNIFORM
    ! scalar exactly for ANY advecting field (div u never enters), at the cost
    ! of the exact global conservation the divergence form has. Subtracting
    ! HALF would instead be the skew-symmetric form (the momentum kernel's
    ! choice, energy-neutral in sum s^2 for any advecting field) but leaves
    ! 1/2 s div u on a uniform field. The projection makes div u small either
    ! way; the trade is uniform-preservation vs. sum s^2 neutrality.
    !
    ! There is deliberately NO upwind option (module header of rans.f90 and
    ! docs/next_session_scalar.md Section 2): central/skew only until the
    ! shared TVD increment lands.
    !
    ! Diffusion: face-flux form on the module's own inverse centre-to-centre
    ! distances, with the face diffusivity
    !
    !   D_face = 1/(Re Pr) + 1/2 (nu_t,L + nu_t,R)/Pr_t(face)
    !
    ! reading turb%nut -- the ONE blended eddy viscosity that LES (WALE /
    ! Smagorinsky), RANS (SST) and IDDES all write, so a single code path
    ! serves all three models (increment S2). Pr_t is the per-scalar constant
    ! or the Kays-Crawford correlation, selected by [scalar.N] prt_model.
    !
    ! Immersed body (increment S3), per scalar via [scalar.N] ibm_wall:
    !   * dirichlet -- volume penalization toward [scalar.N] ibm_value with
    !     the implicit factor mu_s = 1/(1 + dt_gamma coef_p/Pr) formed inline
    !     (ibm.f90's implicit form, so no dt restriction; no mu array per
    !     scalar). The stored coefficient carries the 1/Re scaling, so
    !     coef_p/Pr is exactly the scalar's own 1/(Re Pr) -- ONE cell-centred
    !     coefficient array serves every scalar, each with its own Pr. Inside
    !     the body coef_p = SOLID/Re, so mu_s -> 0 and the cell holds
    !     ibm_value. Diffusive fluxes are NOT masked: the solid cell holds
    !     the wall value and delivers the flux. First-order/staircase, like
    !     the velocity penalization.
    !   * adiabatic -- no penalization; instead BOTH the convective and the
    !     diffusive flux are masked on each of the six p-cell faces whose
    !     staggered velocity coefficient is solid (the solw/sole/... test
    !     rans.f90 uses for k). The mask is symmetric across a face, so the
    !     flux form still telescopes and the fluid conserves sum(s dV)
    !     exactly.
    ! With no immersed body (useIbm off) neither branch is taken and the
    ! arithmetic is the S2 kernel's, byte for byte.
    !
    ! The nut ghost cells follow the momentum SGS convention: the halo
    ! exchange fills them across block and rank boundaries, and they stay at
    ! zero on a physical face (a resolved wall has nu_t -> 0 there anyway;
    ! the T3 wall functions write the wall-cell value into the no-slip
    ! ghosts, which this reads for free).
    subroutine scalar_transport(sc, blk, dns, turb, coef, dt_alpha, dt_beta, dt_gamma)
        type(scalar_type), intent(in) :: sc
        type(block_set_type), intent(inout) :: blk
        type(dns_type), intent(in) :: dns
        type(turb_type), intent(in) :: turb
        ! ibm%coef: with scalars on it always carries the cell-centred VAR_P
        ! column (init_ibm), so there is no dummy-array case here -- the
        ! body is switched off by dns%ibm_enabled alone.
        real(C_DOUBLE), intent(in) :: coef(0:,0:,0:,1:,1:)
        real(C_DOUBLE), intent(in) :: dt_alpha, dt_beta, dt_gamma

        if (.not. scalars_enabled(sc)) return

        ! With no turbulence model turb%nut does not exist; hand the kernel
        ! the 1-cell dummy and switch the eddy term off, which makes the
        ! arithmetic byte-identical to the molecular-only (S1) kernel.
        if (turbulence_is_enabled(turb) .and. allocated(turb%nut)) then
            call scalar_transport_kernel(sc, blk, dns, dt_alpha, dt_beta, dt_gamma, &
                turb%nut, .true., coef)
        else
            call scalar_transport_kernel(sc, blk, dns, dt_alpha, dt_beta, dt_gamma, &
                sc%nutNone, .false., coef)
        end if
    end subroutine scalar_transport

    subroutine scalar_transport_kernel(sc, blk, dns, dt_alpha, dt_beta, dt_gamma, &
            nut, useNut, coef)
        type(scalar_type), intent(in) :: sc
        type(block_set_type), intent(inout) :: blk
        type(dns_type), intent(in) :: dns
        real(C_DOUBLE), intent(in) :: dt_alpha, dt_beta, dt_gamma
        real(C_DOUBLE), intent(in) :: nut(0:,0:,0:,1:)
        logical, intent(in) :: useNut
        real(C_DOUBLE), intent(in) :: coef(0:,0:,0:,1:,1:)

        ! The solid-coefficient test the SGS and RANS kernels use (les.f90,
        ! rans.f90's SOLID_FACE_THRESHOLD): SOLID/Re is 1e30/Re.
        real(C_DOUBLE), parameter :: SOLID_FACE_THRESHOLD = 1.0d20

        integer :: i, j, k, b, is, nx, ny, nz, nBlocks, nScal, var, scr
        real(C_DOUBLE) :: ire, re, uw, ue, vs, vn, wb, wt, divu
        real(C_DOUBLE) :: s0, conv, diff, rhs, dm, fw, fe, ss, mus, ipr
        real(C_DOUBLE) :: ntw, nte, nts, ntn, ntb, ntt
        real(C_DOUBLE) :: dxw, dxe, dys, dyn, dzb, dzt
        logical :: skew, useIbm, adiab
        logical :: clw, cle, cls, cln, clb, clt
        logical :: solw, sole, sols, soln, solb, solt
        logical :: mw, me, ms, mn, mb, mt

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))
        nBlocks = int(blk%nBlocks)
        nScal = int(sc%n)
        re = dns%re
        ire = 1.0d0/re
        skew = logical(dns%conv_skew)
        useIbm = logical(dns%ibm_enabled)

        !$omp target teams distribute parallel do collapse(4) &
        !$omp& map(to: ire, re, dt_alpha, dt_beta, dt_gamma, skew, useNut, useIbm, &
        !$omp& nScal, nx, ny, nz, &
        !$omp& blk%q, blk%d1x, blk%d1y, blk%d1z, blk%physLow, blk%physHigh, nut, coef, &
        !$omp& sc%pr, sc%prt, sc%prtModel, sc%source, sc%invDx, sc%invDy, sc%invDz, &
        !$omp& sc%ibmMode, sc%ibmValue) &
        !$omp& map(tofrom: blk%qs, blk%oldrhs) &
        !$omp& private(i,j,k,b,is,var,scr,uw,ue,vs,vn,wb,wt,divu, &
        !$omp& s0,conv,diff,rhs,dm,fw,fe,ss,mus,ipr,adiab, &
        !$omp& clw,cle,cls,cln,clb,clt,solw,sole,sols,soln,solb,solt, &
        !$omp& mw,me,ms,mn,mb,mt, &
        !$omp& ntw,nte,nts,ntn,ntb,ntt,dxw,dxe,dys,dyn,dzb,dzt)
        do b = 1, nBlocks
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    ! The p cell's six face velocities: exactly the six
                    ! staggered components, no interpolation.
                    uw = blk%q(i,  j,  k,  VAR_U,b)
                    ue = blk%q(i+1,j,  k,  VAR_U,b)
                    vs = blk%q(i,  j,  k,  VAR_V,b)
                    vn = blk%q(i,  j+1,k,  VAR_V,b)
                    wb = blk%q(i,  j,  k,  VAR_W,b)
                    wt = blk%q(i,  j,  k+1,VAR_W,b)
                    divu = (ue - uw)*blk%d1x(i,VAR_P,b) &
                         + (vn - vs)*blk%d1y(j,VAR_P,b) &
                         + (wt - wb)*blk%d1z(k,VAR_P,b)

                    ! A FACE_CLOSED face borders a block removed inside the
                    ! immersed body: its halo is held at zero, so the
                    ! DIFFUSIVE flux must be masked (the convective one
                    ! already vanishes -- the face velocity is pinned to 0).
                    ! FACE_PHYS faces are NOT masked: their ghosts carry the
                    ! boundary condition.
                    clw = i == 1  .and. blk%physLow(1,b)  == FACE_CLOSED
                    cle = i == nx .and. blk%physHigh(1,b) == FACE_CLOSED
                    cls = j == 1  .and. blk%physLow(2,b)  == FACE_CLOSED
                    cln = j == ny .and. blk%physHigh(2,b) == FACE_CLOSED
                    clb = k == 1  .and. blk%physLow(3,b)  == FACE_CLOSED
                    clt = k == nz .and. blk%physHigh(3,b) == FACE_CLOSED

                    ! Solid staggered faces of this p cell (adiabatic body
                    ! mode only; scalar-independent, so hoisted out of the
                    ! scalar loop). rans.f90's test for k, verbatim.
                    solw = .false.; sole = .false.; sols = .false.
                    soln = .false.; solb = .false.; solt = .false.
                    if (useIbm) then
                        solw = abs(coef(i,  j,  k,  VAR_U,b)) > SOLID_FACE_THRESHOLD
                        sole = abs(coef(i+1,j,  k,  VAR_U,b)) > SOLID_FACE_THRESHOLD
                        sols = abs(coef(i,  j,  k,  VAR_V,b)) > SOLID_FACE_THRESHOLD
                        soln = abs(coef(i,  j+1,k,  VAR_V,b)) > SOLID_FACE_THRESHOLD
                        solb = abs(coef(i,  j,  k,  VAR_W,b)) > SOLID_FACE_THRESHOLD
                        solt = abs(coef(i,  j,  k+1,VAR_W,b)) > SOLID_FACE_THRESHOLD
                    end if

                    ! Face eddy viscosities: the same 1/2 (L + R) average the
                    ! momentum SGS correction uses for its face-normal terms.
                    ! Scalar-independent, so they are hoisted out of the
                    ! scalar loop; Pr_t is not (it may be per-scalar).
                    ntw = 0.0d0; nte = 0.0d0; nts = 0.0d0
                    ntn = 0.0d0; ntb = 0.0d0; ntt = 0.0d0
                    if (useNut) then
                        ntw = 0.5d0*(nut(i-1,j,k,b) + nut(i,j,k,b))
                        nte = 0.5d0*(nut(i,j,k,b) + nut(i+1,j,k,b))
                        nts = 0.5d0*(nut(i,j-1,k,b) + nut(i,j,k,b))
                        ntn = 0.5d0*(nut(i,j,k,b) + nut(i,j+1,k,b))
                        ntb = 0.5d0*(nut(i,j,k-1,b) + nut(i,j,k,b))
                        ntt = 0.5d0*(nut(i,j,k,b) + nut(i,j,k+1,b))
                    end if

                    do is = 1, nScal
                        var = VAR_S0 + is
                        scr = SCR_S0 + is
                        s0 = blk%q(i,j,k,var,b)

                        ! An adiabatic body seals this cell's solid faces:
                        ! convective AND diffusive flux masked. Every other
                        ! configuration leaves adiab .false., so the six
                        ! masks below collapse to the FACE_CLOSED flags and
                        ! the S2 expressions are reproduced exactly.
                        adiab = useIbm .and. sc%ibmMode(is) == SC_IBM_ADIABATIC
                        mw = clw .or. (adiab .and. solw)
                        me = cle .or. (adiab .and. sole)
                        ms = cls .or. (adiab .and. sols)
                        mn = cln .or. (adiab .and. soln)
                        mb = clb .or. (adiab .and. solb)
                        mt = clt .or. (adiab .and. solt)

                        if (adiab) then
                            conv = (merge(0.0d0, ue*0.5d0*(s0 + blk%q(i+1,j,k,var,b)), me) &
                                  - merge(0.0d0, uw*0.5d0*(blk%q(i-1,j,k,var,b) + s0), mw)) &
                                    *blk%d1x(i,VAR_P,b) &
                                 + (merge(0.0d0, vn*0.5d0*(s0 + blk%q(i,j+1,k,var,b)), mn) &
                                  - merge(0.0d0, vs*0.5d0*(blk%q(i,j-1,k,var,b) + s0), ms)) &
                                    *blk%d1y(j,VAR_P,b) &
                                 + (merge(0.0d0, wt*0.5d0*(s0 + blk%q(i,j,k+1,var,b)), mt) &
                                  - merge(0.0d0, wb*0.5d0*(blk%q(i,j,k-1,var,b) + s0), mb)) &
                                    *blk%d1z(k,VAR_P,b)
                        else
                            conv = (ue*0.5d0*(s0 + blk%q(i+1,j,k,var,b)) &
                                  - uw*0.5d0*(blk%q(i-1,j,k,var,b) + s0))*blk%d1x(i,VAR_P,b) &
                                 + (vn*0.5d0*(s0 + blk%q(i,j+1,k,var,b)) &
                                  - vs*0.5d0*(blk%q(i,j-1,k,var,b) + s0))*blk%d1y(j,VAR_P,b) &
                                 + (wt*0.5d0*(s0 + blk%q(i,j,k+1,var,b)) &
                                  - wb*0.5d0*(blk%q(i,j,k-1,var,b) + s0))*blk%d1z(k,VAR_P,b)
                        end if
                        if (skew) conv = conv - s0*divu

                        ! Molecular diffusivity + the eddy part on each face.
                        ! With useNut off every face keeps dm exactly, so the
                        ! S1 arithmetic is reproduced bit-for-bit.
                        dm = ire/sc%pr(is)
                        dxw = dm; dxe = dm; dys = dm
                        dyn = dm; dzb = dm; dzt = dm
                        if (useNut) then
                            dxw = dm + eddy_diffusivity(ntw, sc%pr(is), sc%prt(is), &
                                sc%prtModel(is), re)
                            dxe = dm + eddy_diffusivity(nte, sc%pr(is), sc%prt(is), &
                                sc%prtModel(is), re)
                            dys = dm + eddy_diffusivity(nts, sc%pr(is), sc%prt(is), &
                                sc%prtModel(is), re)
                            dyn = dm + eddy_diffusivity(ntn, sc%pr(is), sc%prt(is), &
                                sc%prtModel(is), re)
                            dzb = dm + eddy_diffusivity(ntb, sc%pr(is), sc%prt(is), &
                                sc%prtModel(is), re)
                            dzt = dm + eddy_diffusivity(ntt, sc%pr(is), sc%prt(is), &
                                sc%prtModel(is), re)
                        end if

                        fw = merge(0.0d0, dxw*(s0 - blk%q(i-1,j,k,var,b)) &
                            *sc%invDx(i,b), mw)
                        fe = merge(0.0d0, dxe*(blk%q(i+1,j,k,var,b) - s0) &
                            *sc%invDx(i+1,b), me)
                        diff = (fe - fw)*blk%d1x(i,VAR_P,b)
                        fw = merge(0.0d0, dys*(s0 - blk%q(i,j-1,k,var,b)) &
                            *sc%invDy(j,b), ms)
                        fe = merge(0.0d0, dyn*(blk%q(i,j+1,k,var,b) - s0) &
                            *sc%invDy(j+1,b), mn)
                        diff = diff + (fe - fw)*blk%d1y(j,VAR_P,b)
                        fw = merge(0.0d0, dzb*(s0 - blk%q(i,j,k-1,var,b)) &
                            *sc%invDz(k,b), mb)
                        fe = merge(0.0d0, dzt*(blk%q(i,j,k+1,var,b) - s0) &
                            *sc%invDz(k+1,b), mt)
                        diff = diff + (fe - fw)*blk%d1z(k,VAR_P,b)

                        rhs = -conv + diff + sc%source(is)
                        ss = s0 + dt_alpha*rhs + dt_beta*blk%oldrhs(i,j,k,scr,b)
                        ! Implicit volume penalization toward the body value
                        ! (dirichlet mode). oldrhs keeps the UNpenalized rhs,
                        ! exactly as the momentum predictor does with mu.
                        if (useIbm .and. .not. adiab) then
                            ipr = 1.0d0/sc%pr(is)
                            mus = 1.0d0/(1.0d0 + dt_gamma*coef(i,j,k,VAR_P,b)*ipr)
                            ss = ss*mus + (1.0d0 - mus)*sc%ibmValue(is)
                        end if
                        blk%qs(i,j,k,scr,b) = ss
                        blk%oldrhs(i,j,k,scr,b) = rhs
                    end do
                end do
            end do
        end do
        end do
        !$omp end target teams distribute parallel do
    end subroutine scalar_transport_kernel

    ! Substage tail: qs -> q for the scalar slots, then the ghost/halo sync.
    ! Runs AFTER the pressure projection returns -- the scalar never enters
    ! the projection loop.
    subroutine scalar_finish(sc, blk, bc, c)
        type(scalar_type), intent(in) :: sc
        type(block_set_type), intent(inout) :: blk
        type(boundary_type), intent(in) :: bc
        type(comm_type), intent(inout) :: c

        integer :: i, j, k, b, is, nx, ny, nz, nBlocks, nScal

        if (.not. scalars_enabled(sc)) return

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))
        nBlocks = int(blk%nBlocks)
        nScal = int(sc%n)

        !$omp target teams distribute parallel do collapse(4) &
        !$omp& map(to: nScal, blk%qs) map(tofrom: blk%q) &
        !$omp& private(i,j,k,b,is)
        do b = 1, nBlocks
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    do is = 1, nScal
                        blk%q(i,j,k,VAR_S0+is,b) = blk%qs(i,j,k,SCR_S0+is,b)
                    end do
                end do
            end do
        end do
        end do
        !$omp end target teams distribute parallel do

        call scalar_sync(sc, blk, bc, c)
    end subroutine scalar_finish

    ! Physical ghosts + one batched halo exchange over all scalar variables
    ! (they ride the same entries, and the same message, as u,v,w,p).
    subroutine scalar_sync(sc, blk, bc, c)
        type(scalar_type), intent(in) :: sc
        type(block_set_type), intent(inout) :: blk
        type(boundary_type), intent(in) :: bc
        type(comm_type), intent(inout) :: c

        if (.not. scalars_enabled(sc)) return

        call apply_scalar_bc_q(blk, bc, int(sc%n), sc%bcType, sc%bcValue)
        call exchange_halos(c, blk, sc%varList)
    end subroutine scalar_sync

end module scalar
