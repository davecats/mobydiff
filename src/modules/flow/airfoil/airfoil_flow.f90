module airfoil_flow
    ! Quasi-2D immersed-boundary airfoil/cylinder case ([case] name = airfoil,
    ! docs/next_session_airfoil.md, phases A1/A2). The angle of attack turns
    ! into a Dirichlet freestream on the inlet faces (x_min, y_min, y_max) and
    ! a Dirichlet-pressure outlet at x_max, composed purely of PATCH TYPES
    ! (resolve_face_bcs derives the per-variable rows); geometry comes from
    ! the standard file-based IBM path. Runtime statistics are C_L(t)/C_D(t)
    ! from the penalization integral F = int coef*u dV -- exact discrete
    ! bookkeeping of the momentum the body removes (pressure + friction +
    ! modeled turbulent stress all included, since every momentum write is
    ! *mu-masked).
    use, intrinsic :: iso_c_binding
    use :: flow_case_base, only: case_type
    use :: generic_flow, only: set_generic_defaults
    use :: init, only: dns_type, grid_type, VAR_U, VAR_V, VAR_W
    use :: blocks, only: block_set_type
    use :: boundary, only: boundary_type, boundary_face_id, PATCH_INLET, PATCH_OUTLET
    use :: pressure_solver, only: pressure_solver_type
    use :: ibmm, only: ibm_type
    use :: comm, only: comm_type, comm_allreduce_sum
    use :: case_config_helpers, only: next_config_entry, to_lower, clean_config_string
    implicit none

    private

    character(len=*), parameter :: AIRFOIL_CASE_NAME = "airfoil"

    type, extends(case_type), public :: airfoil_case_type
        real(C_DOUBLE) :: aoa = 0.0d0            ! angle of attack [deg]
        real(C_DOUBLE) :: u_inf = 1.0d0
        real(C_DOUBLE) :: chord = 1.0d0
        integer :: force_sample_interval = 10
        character(len=256) :: runtime_file = "forces.txt"
        logical :: header_written = .false.
    contains
        procedure :: read_config => airfoil_read_config
        procedure :: apply_defaults => airfoil_apply_defaults
        procedure :: setup_after_grid => airfoil_setup_after_grid
        procedure :: initialise_fields => airfoil_initialise_fields
        procedure :: after_step => airfoil_after_step
        procedure :: finalize => airfoil_finalize
    end type airfoil_case_type

    public :: create_airfoil_case, AIRFOIL_CASE_NAME

contains

    subroutine create_airfoil_case(flow)
        class(case_type), allocatable, intent(out) :: flow

        allocate(airfoil_case_type :: flow)
        flow%name = AIRFOIL_CASE_NAME
    end subroutine create_airfoil_case

    ! Freestream direction cosines from the angle of attack.
    subroutine freestream(this, ux, uy)
        class(airfoil_case_type), intent(in) :: this
        real(C_DOUBLE), intent(out) :: ux, uy
        real(C_DOUBLE) :: a

        a = this%aoa*(4.0d0*atan(1.0d0))/180.0d0
        ux = this%u_inf*cos(a)
        uy = this%u_inf*sin(a)
    end subroutine freestream

    subroutine airfoil_apply_defaults(this, dns, g, bc, c, ps)
        class(airfoil_case_type), intent(inout) :: this
        type(dns_type), intent(inout) :: dns
        type(grid_type), intent(inout) :: g
        type(boundary_type), intent(inout) :: bc
        type(comm_type), intent(inout) :: c
        type(pressure_solver_type), intent(inout) :: ps

        real(C_DOUBLE) :: ux, uy
        integer :: face_id, dir, side

        call set_generic_defaults(dns, g, bc, c, ps)
        this%name = AIRFOIL_CASE_NAME
        call freestream(this, ux, uy)

        dns%ibm_enabled = .true.
        dns%forcing = 0.0d0             ! no volume forcing in airfoil runs
        dns%initial_velocity = [ux, uy, 0.0d0]

        ! Freestream composition in ONE vocabulary -- patch types (the ini is
        ! parsed AFTER apply_defaults, so explicit [boundary] keys still win):
        ! x_min, y_min, y_max Dirichlet inlets carrying (U cos a, U sin a, 0),
        ! x_max the Dirichlet-pressure outlet; z periodic (quasi-2D, nz = nb).
        ! Large |aoa| needs a taller domain -- the y faces are far-field
        ! Dirichlet, the standard penalization freestream.
        bc%isPeriodic(1) = .false.
        bc%isPeriodic(2) = .false.
        bc%isPeriodic(3) = .true.
        do dir = 1, 2
            do side = 0, 1
                face_id = boundary_face_id(dir, side)
                if (dir == 1 .and. side == 1) then
                    bc%facePatchType(face_id) = PATCH_OUTLET
                else
                    bc%facePatchType(face_id) = PATCH_INLET
                    bc%faceBcDefaultValue(VAR_U,face_id) = ux
                    bc%faceBcDefaultValue(VAR_V,face_id) = uy
                    bc%faceBcDefaultValue(VAR_W,face_id) = 0.0d0
                end if
            end do
        end do
    end subroutine airfoil_apply_defaults

    subroutine airfoil_setup_after_grid(this, blk, dns, g, bc, c)
        class(airfoil_case_type), intent(inout) :: this
        type(block_set_type), intent(inout) :: blk
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        type(boundary_type), intent(in) :: bc
        type(comm_type), intent(in) :: c
    end subroutine airfoil_setup_after_grid

    ! Uniform freestream everywhere (impulsive start; the IBM damps the body
    ! interior within the first steps).
    subroutine airfoil_initialise_fields(this, blk, dns, g, bc, c)
        class(airfoil_case_type), intent(inout) :: this
        type(block_set_type), intent(inout) :: blk
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        type(boundary_type), intent(in) :: bc
        type(comm_type), intent(in) :: c

        real(C_DOUBLE) :: ux, uy

        call freestream(this, ux, uy)
        blk%q(:,:,:,VAR_U,:) = ux
        blk%q(:,:,:,VAR_V,:) = uy
        blk%q(:,:,:,VAR_W,:) = 0.0d0
    end subroutine airfoil_initialise_fields

    ! Penalization-integral force on the body, F_d = sum coef*q*V_d over the
    ! interior staggered faces (1..nb counts every face exactly once
    ! globally; the block-redundant layer lives in halos). Evaluated on the
    ! END-OF-STEP field -- after the last projection exchange -- where
    ! q_final = mu*q_unmasked makes coef*q_final the exact momentum sink of
    ! the model as solved (viscous + pressure + modeled turbulent stress).
    ! Reduction determinism: per-BLOCK partial sums are scattered into the
    ! global block table (one contributor per entry, so the allreduce is
    ! exact) and the final sum runs in global-id order on every rank -- the
    ! sampled force is independent of the rank count by construction.
    subroutine airfoil_after_step(this, blk, dns, g, c, ibm)
        class(airfoil_case_type), intent(inout) :: this
        type(block_set_type), intent(inout) :: blk
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        type(comm_type), intent(in) :: c
        type(ibm_type), intent(in) :: ibm

        real(C_DOUBLE) :: f(3), sx, sy, sz, a, cl, cd, qref
        real(C_DOUBLE), allocatable :: fb(:,:), fbg(:)
        integer(C_INT) :: i, j, k, b, nx, ny, nz, nBlocks

        if (this%force_sample_interval <= 0) return
        if (modulo(int(dns%step_current), this%force_sample_interval) /= 0) return

        nx = blk%nb(1); ny = blk%nb(2); nz = blk%nb(3)
        nBlocks = blk%nBlocks
        allocate(fb(3, nBlocks))

#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute &
        !$omp& map(to: nx, ny, nz, blk%d1x, blk%d1y, blk%d1z, blk%q, ibm%coef) &
        !$omp& map(from: fb(1:3,1:nBlocks)) private(i,j,k,b,sx,sy,sz)
#endif
        do b = 1_C_INT, nBlocks
            sx = 0.0d0; sy = 0.0d0; sz = 0.0d0
#ifdef USE_OPENMP_OFFLOAD
            !$omp parallel do collapse(3) reduction(+:sx,sy,sz) private(i,j,k)
#endif
            do k = 1_C_INT, nz
                do j = 1_C_INT, ny
                    do i = 1_C_INT, nx
                        sx = sx + ibm%coef(i,j,k,VAR_U,b)*blk%q(i,j,k,VAR_U,b) &
                            /(blk%d1x(i,VAR_U,b)*blk%d1y(j,VAR_U,b)*blk%d1z(k,VAR_U,b))
                        sy = sy + ibm%coef(i,j,k,VAR_V,b)*blk%q(i,j,k,VAR_V,b) &
                            /(blk%d1x(i,VAR_V,b)*blk%d1y(j,VAR_V,b)*blk%d1z(k,VAR_V,b))
                        sz = sz + ibm%coef(i,j,k,VAR_W,b)*blk%q(i,j,k,VAR_W,b) &
                            /(blk%d1x(i,VAR_W,b)*blk%d1y(j,VAR_W,b)*blk%d1z(k,VAR_W,b))
                    end do
                end do
            end do
#ifdef USE_OPENMP_OFFLOAD
            !$omp end parallel do
#endif
            fb(1,b) = sx; fb(2,b) = sy; fb(3,b) = sz
        end do
#ifdef USE_OPENMP_OFFLOAD
        !$omp end target teams distribute
#endif

        ! Scatter to global ids (0-based), exact allreduce, ordered final sum.
        allocate(fbg(3*blk%nBlocksGlobal))
        fbg = 0.0d0
        do b = 1, int(nBlocks)
            fbg(3*blk%globalId(b)+1:3*blk%globalId(b)+3) = fb(:,b)
        end do
        call comm_allreduce_sum(c, fbg)
        f = 0.0d0
        do b = 1, int(blk%nBlocksGlobal)
            f = f + fbg(3*(b-1)+1:3*(b-1)+3)
        end do

        ! Drag along the freestream, lift normal to it (rho = 1):
        ! C = 2 F.e / (U_inf^2 * chord * Lz).
        a = this%aoa*(4.0d0*atan(1.0d0))/180.0d0
        qref = this%u_inf**2*this%chord*dns%leng(3)
        cd = 2.0d0*( f(1)*cos(a) + f(2)*sin(a))/qref
        cl = 2.0d0*(-f(1)*sin(a) + f(2)*cos(a))/qref

        if (c%has_terminal) call append_forces(this, dns, cl, cd)
    end subroutine airfoil_after_step

    subroutine append_forces(this, dns, cl, cd)
        class(airfoil_case_type), intent(inout) :: this
        type(dns_type), intent(in) :: dns
        real(C_DOUBLE), intent(in) :: cl, cd

        integer :: unit, stat
        character(len=*), parameter :: header = "iteration time cl cd"

        if (.not. this%header_written) then
            write(*,'(A)') header
            open(newunit=unit, file=trim(this%runtime_file), status="replace", &
                action="write", iostat=stat)
            if (stat == 0) then
                write(unit,'(A)') header
                close(unit)
            else
                print *, "warning: could not open runtime data file: ", trim(this%runtime_file)
            end if
            this%header_written = .true.
        end if

        write(*,'(I10,3(1X,ES16.8))') int(dns%step_current), dns%t_current, cl, cd
        open(newunit=unit, file=trim(this%runtime_file), status="old", &
            position="append", action="write", iostat=stat)
        if (stat == 0) then
            write(unit,'(I10,3(1X,ES16.8))') int(dns%step_current), dns%t_current, cl, cd
            close(unit)
        end if
    end subroutine append_forces

    subroutine airfoil_finalize(this, dns, g, c)
        class(airfoil_case_type), intent(inout) :: this
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        type(comm_type), intent(in) :: c
    end subroutine airfoil_finalize

    subroutine airfoil_read_config(this, input_file, has_terminal)
        class(airfoil_case_type), intent(inout) :: this
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
            if (trim(section) /= "case.airfoil") cycle
            call apply_airfoil_case_value(this, to_lower(key), value, line_no, has_terminal)
        end do

        close(unit)
    end subroutine airfoil_read_config

    subroutine apply_airfoil_case_value(this, key, value, line_no, has_terminal)
        type(airfoil_case_type), intent(inout) :: this
        character(len=*), intent(in) :: key, value
        integer, intent(in) :: line_no
        logical, intent(in), optional :: has_terminal

        integer :: int_value, stat
        real(C_DOUBLE) :: real_value
        logical :: terminal

        terminal = .true.
        if (present(has_terminal)) terminal = has_terminal

        select case (trim(key))
        case ("aoa")
            read(value, *, iostat=stat) real_value
            if (stat == 0) this%aoa = real_value
        case ("u_inf")
            read(value, *, iostat=stat) real_value
            if (stat == 0) this%u_inf = real_value
        case ("chord")
            read(value, *, iostat=stat) real_value
            if (stat == 0) this%chord = real_value
        case ("force_sample_interval")
            read(value, *, iostat=stat) int_value
            if (stat == 0) this%force_sample_interval = int_value
        case ("runtime_file")
            this%runtime_file = clean_config_string(value)
        case default
            if (terminal) print *, "warning: unknown airfoil case key on input line", line_no, ": ", trim(key)
        end select
    end subroutine apply_airfoil_case_value

end module airfoil_flow
