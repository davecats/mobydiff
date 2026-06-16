module generic_flow
    use, intrinsic :: iso_c_binding
    use :: flow_case_base, only: case_type
    use :: init, only: dns_type, grid_type, GRID_UNIFORM, VAR_U
    use :: blocks, only: block_set_type
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

    subroutine generic_setup_after_grid(this, blk, dns, g, bc, c)
        class(generic_case_type), intent(inout) :: this
        type(block_set_type), intent(inout) :: blk
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        type(boundary_type), intent(in) :: bc
        type(comm_type), intent(in) :: c
    end subroutine generic_setup_after_grid

    subroutine generic_initialise_fields(this, blk, dns, g, bc, c)
        class(generic_case_type), intent(inout) :: this
        type(block_set_type), intent(inout) :: blk
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        type(boundary_type), intent(in) :: bc
        type(comm_type), intent(in) :: c

        integer :: i, j, k, v, b, nSeed
        integer, allocatable :: seed(:)
        real(C_DOUBLE) :: r
        character(len=32) :: shearEnv
        real(C_DOUBLE) :: shearAmp, twopi, ly, ycoord
        integer :: gx, gz, ios

        blk%q(:,:,:,1,:) = dns%initial_velocity(1)
        blk%q(:,:,:,2,:) = dns%initial_velocity(2)
        blk%q(:,:,:,3,:) = dns%initial_velocity(3)

        ! Diagnostic shear-mode gate (interface_review §vi): a smooth periodic
        ! mean shear u(y)=sin(2*pi*y/Ly) across the interface, plus a structured
        ! checkerboard seed v = a*(-1)^(gx+gz) on the interface-normal velocity
        ! in the FINE region only -- the zero-coarse-average [+a,-a,+a,-a]
        ! component the coarse pressure cannot see. MOBY_SHEAR_SEED=<amplitude>.
        call get_environment_variable("MOBY_SHEAR_SEED", shearEnv)
        if (len_trim(shearEnv) > 0) then
            read(shearEnv, *, iostat=ios) shearAmp
            if (ios /= 0) shearAmp = 1.0d-3
            twopi = 8.0d0*atan(1.0d0); ly = dns%leng(2)
            do b = 1, int(blk%nBlocks)
                do k = 0, int(blk%nb(3))+1
                    do j = 0, int(blk%nb(2))+1
                        do i = 0, int(blk%nb(1))+1
                            ycoord = blk%y(j, VAR_U, b)
                            ! Periodic y: Kolmogorov-like sin shear (projection-
                            ! only gate). Walls in y: steady laminar Poiseuille
                            ! u = 4 y(Ly-y)/Ly^2 (Umax=1), sustained by a uniform
                            ! forcing_x = 8/Re set in the input - a base that is
                            ! linearly stable below Re~5772, so any full-step
                            ! blow-up is attributable to the 2:1 interface.
                            if (bc%isPeriodic(2)) then
                                blk%q(i,j,k,1,b) = sin(twopi*ycoord/ly)
                            else
                                blk%q(i,j,k,1,b) = 4.0d0*ycoord*(ly-ycoord)/(ly*ly)
                            end if
                            blk%q(i,j,k,2,b) = 0.0d0
                            blk%q(i,j,k,3,b) = 0.0d0
                            if (blk%level(b) > 0_C_INT) then
                                gx = int(blk%origin(1,b)) + i - 1
                                gz = int(blk%origin(3,b)) + k - 1
                                blk%q(i,j,k,2,b) = shearAmp*real(1 - 2*modulo(gx+gz,2), C_DOUBLE)
                            end if
                        end do
                    end do
                end do
            end do
            return
        end if

        ! Optional white-noise perturbation ([flow] initial_noise), used by
        ! the interface-decay gate. Deterministic per rank.
        if (dns%initial_noise > 0.0d0) then
            call random_seed(size=nSeed)
            allocate(seed(nSeed))
            seed = 20260612 + 7919*c%cart_rank
            call random_seed(put=seed)
            do b = 1, int(blk%nBlocks)
                do v = 1, 3
                    do k = 1, int(blk%nb(3))
                        do j = 1, int(blk%nb(2))
                            do i = 1, int(blk%nb(1))
                                call random_number(r)
                                blk%q(i,j,k,v,b) = blk%q(i,j,k,v,b) &
                                    + dns%initial_noise*(2.0d0*r - 1.0d0)
                            end do
                        end do
                    end do
                end do
            end do
        end if
    end subroutine generic_initialise_fields

    subroutine generic_after_step(this, blk, dns, g, c)
        class(generic_case_type), intent(inout) :: this
        type(block_set_type), intent(inout) :: blk
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
        dns%initial_velocity = 0.0d0
        dns%ibm_enabled = .true.
        dns%ibm_coeff_file = ""
        dns%field_prefix = "field"
        dns%field_interval = 0
        dns%restart_file = ""

        g%distribution = GRID_UNIFORM
        g%stretch = 0.0d0
        g%natural_dyw_plus = 0.05d0
        g%natural_one_sided = .false.

        c%dims = 0
        ps%nIter = 3_C_INT
        ps%sor = 1.5d0
    end subroutine set_generic_defaults

end module generic_flow
