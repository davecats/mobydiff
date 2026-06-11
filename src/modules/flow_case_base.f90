module flow_case_base
    use :: init, only: dns_type, grid_type
    use :: blocks, only: block_set_type
    use :: boundary, only: boundary_type
    use :: pressure_solver, only: pressure_solver_type
    use :: comm, only: comm_type
    implicit none

    private

    type, abstract, public :: case_type
        character(len=64) :: name = "generic"
    contains
        procedure(case_read_config), deferred :: read_config
        procedure(case_apply_defaults), deferred :: apply_defaults
        procedure(case_setup_after_grid), deferred :: setup_after_grid
        procedure(case_initialise_fields), deferred :: initialise_fields
        procedure(case_after_step), deferred :: after_step
        procedure(case_finalize), deferred :: finalize
    end type case_type

    abstract interface
        subroutine case_read_config(this, input_file, has_terminal)
            import :: case_type
            class(case_type), intent(inout) :: this
            character(len=*), intent(in) :: input_file
            logical, intent(in), optional :: has_terminal
        end subroutine case_read_config

        subroutine case_apply_defaults(this, dns, g, bc, c, ps)
            import :: case_type, dns_type, grid_type, boundary_type, comm_type, pressure_solver_type
            class(case_type), intent(inout) :: this
            type(dns_type), intent(inout) :: dns
            type(grid_type), intent(inout) :: g
            type(boundary_type), intent(inout) :: bc
            type(comm_type), intent(inout) :: c
            type(pressure_solver_type), intent(inout) :: ps
        end subroutine case_apply_defaults

        subroutine case_setup_after_grid(this, blk, dns, g, bc, c)
            import :: case_type, block_set_type, dns_type, grid_type, boundary_type, comm_type
            class(case_type), intent(inout) :: this
            type(block_set_type), intent(inout) :: blk
            type(dns_type), intent(in) :: dns
            type(grid_type), intent(in) :: g
            type(boundary_type), intent(in) :: bc
            type(comm_type), intent(in) :: c
        end subroutine case_setup_after_grid

        subroutine case_initialise_fields(this, blk, dns, g, bc, c)
            import :: case_type, block_set_type, dns_type, grid_type, boundary_type, comm_type
            class(case_type), intent(inout) :: this
            type(block_set_type), intent(inout) :: blk
            type(dns_type), intent(in) :: dns
            type(grid_type), intent(in) :: g
            type(boundary_type), intent(in) :: bc
            type(comm_type), intent(in) :: c
        end subroutine case_initialise_fields

        subroutine case_after_step(this, blk, dns, g, c)
            import :: case_type, block_set_type, dns_type, grid_type, comm_type
            class(case_type), intent(inout) :: this
            type(block_set_type), intent(inout) :: blk
            type(dns_type), intent(in) :: dns
            type(grid_type), intent(in) :: g
            type(comm_type), intent(in) :: c
        end subroutine case_after_step

        subroutine case_finalize(this, dns, g, c)
            import :: case_type, dns_type, grid_type, comm_type
            class(case_type), intent(inout) :: this
            type(dns_type), intent(in) :: dns
            type(grid_type), intent(in) :: g
            type(comm_type), intent(in) :: c
        end subroutine case_finalize
    end interface

end module flow_case_base
