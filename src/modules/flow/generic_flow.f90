module generic_flow
    use, intrinsic :: iso_c_binding
    use :: flow_case_base, only: case_type
    use :: init, only: dns_type, grid_type, field_type, GRID_UNIFORM
    use :: boundary, only: boundary_type, init_bc
    use :: pressure_solver, only: pressure_solver_type
    use :: comm, only: comm_type
    implicit none

    private

    character(len=*), parameter :: GENERIC_CASE_NAME = "generic"

    type, extends(case_type), public :: generic_case_type
    contains
        procedure :: read_config => generic_read_config
        procedure :: apply_defaults => generic_apply_defaults
        procedure :: setup_after_grid => generic_setup_after_grid
        procedure :: initialise_fields => generic_initialise_fields
        procedure :: after_step => generic_after_step
        procedure :: finalize => generic_finalize
    end type generic_case_type

    public :: set_generic_defaults
    public :: create_generic_case, GENERIC_CASE_NAME

contains

    subroutine create_generic_case(flow)
        class(case_type), allocatable, intent(out) :: flow

        allocate(generic_case_type :: flow)
        flow%name = GENERIC_CASE_NAME
    end subroutine create_generic_case

    subroutine generic_read_config(this, input_file, has_terminal)
        class(generic_case_type), intent(inout) :: this
        character(len=*), intent(in) :: input_file
        logical, intent(in), optional :: has_terminal
    end subroutine generic_read_config

    subroutine generic_apply_defaults(this, dns, g, bc, c, ps)
        class(generic_case_type), intent(inout) :: this
        type(dns_type), intent(inout) :: dns
        type(grid_type), intent(inout) :: g
        type(boundary_type), intent(inout) :: bc
        type(comm_type), intent(inout) :: c
        type(pressure_solver_type), intent(inout) :: ps

        call set_generic_defaults(dns, g, bc, c, ps)
        this%name = GENERIC_CASE_NAME
    end subroutine generic_apply_defaults

    subroutine generic_setup_after_grid(this, f, dns, g, bc, c)
        class(generic_case_type), intent(inout) :: this
        type(field_type), intent(inout) :: f
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        type(boundary_type), intent(in) :: bc
        type(comm_type), intent(in) :: c
    end subroutine generic_setup_after_grid

    subroutine generic_initialise_fields(this, f, dns, g, bc, c)
        class(generic_case_type), intent(inout) :: this
        type(field_type), intent(inout) :: f
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        type(boundary_type), intent(in) :: bc
        type(comm_type), intent(in) :: c
    end subroutine generic_initialise_fields

    subroutine generic_after_step(this, f, dns, g, c)
        class(generic_case_type), intent(inout) :: this
        type(field_type), intent(inout) :: f
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        type(comm_type), intent(in) :: c
    end subroutine generic_after_step

    subroutine generic_finalize(this, dns, g, c)
        class(generic_case_type), intent(inout) :: this
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        type(comm_type), intent(in) :: c
    end subroutine generic_finalize

    subroutine set_generic_defaults(dns, g, bc, c, ps)
        type(dns_type), intent(inout) :: dns
        type(grid_type), intent(inout) :: g
        type(boundary_type), intent(inout) :: bc
        type(comm_type), intent(inout) :: c
        type(pressure_solver_type), intent(inout) :: ps

        call init_bc(bc)

        dns%globalSize = 0_C_INT
        dns%localSize = 0_C_INT
        dns%step_current = 0_C_INT
        dns%nsteps = 0_C_INT
        dns%leng = 0.0d0
        dns%re = 0.0d0
        dns%dt = 0.0d0
        dns%t_final = 0.0d0
        dns%t_current = 0.0d0
        dns%cfl = 0.0d0
        dns%cflmax = 0.0d0
        dns%pecletmax = 0.0d0
        dns%dtmax = huge(1.0d0)
        dns%forcing = 0.0d0
        dns%ibm_enabled = .true.
        dns%field_prefix = "field"
        dns%field_interval = 0
        dns%restart_file = ""

        g%distribution = GRID_UNIFORM
        g%stretch = 0.0d0

        c%dims = 0
        ps%nIter = 3_C_INT
        ps%sor = 1.5d0
    end subroutine set_generic_defaults

end module generic_flow
