module generic_flow
    use, intrinsic :: iso_c_binding
    use :: flow_case_base, only: case_type
    use :: init, only: dns_type, grid_type, GRID_UNIFORM
    use :: blocks, only: block_set_type
    use :: boundary, only: boundary_type, init_bc
    use :: pressure_solver, only: pressure_solver_type
    use :: ibmm, only: ibm_type
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
        real(C_DOUBLE) :: r, kk

        if (trim(dns%initial) == "beltrami" .or. trim(dns%initial) == "tgv" &
            .or. trim(dns%initial) == "tgv3d") then
            ! Analytic fields on a 2*pi-periodic cube, k=2*pi/Lx, each velocity
            ! component set at its own staggered coordinate. Beltrami / ABC
            ! (curl u = k u, fully 3D, decays exp(-nu k^2 t)):
            !   u = sin(kz)+cos(ky)  v = sin(kx)+cos(kz)  w = sin(ky)+cos(kx)
            ! Taylor-Green (2D in x-y, decays exp(-2 nu k^2 t)):
            !   u =  sin(kx)cos(ky)  v = -cos(kx)sin(ky)  w = 0
            ! tgv3d (manufactured momentum-operator test, NOT an NS solution):
            !   u = sin(kx)cos(ky)cos(kz)
            !   v = cos(kx)sin(ky)cos(kz)
            !   w = cos(kx)cos(ky)sin(kz)
            ! Every component varies in every direction, so the wall-normal
            ! velocity varies in the normal direction at all three interface
            ! orientations (Beltrami's dv/dy=0 makes that term vanish). The
            ! Laplacian of each component is -3k^2 times the component, so the
            ! analytic diffusion is trivial; the advection is closed-form.
            ! Set interior + the single halo layer; use blk%q's real bounds
            ! (0:nb+1 here) since blk%x/y/z extend further (-1:nb+2) than q.
            kk = 8.0d0*atan(1.0d0)/dns%leng(1)
            do b = 1, int(blk%nBlocks)
                do k = lbound(blk%q,3), ubound(blk%q,3)
                    do j = lbound(blk%q,2), ubound(blk%q,2)
                        do i = lbound(blk%q,1), ubound(blk%q,1)
                            if (trim(dns%initial) == "beltrami") then
                                blk%q(i,j,k,1,b) = sin(kk*blk%z(k,1,b)) + cos(kk*blk%y(j,1,b))
                                blk%q(i,j,k,2,b) = sin(kk*blk%x(i,2,b)) + cos(kk*blk%z(k,2,b))
                                blk%q(i,j,k,3,b) = sin(kk*blk%y(j,3,b)) + cos(kk*blk%x(i,3,b))
                            else if (trim(dns%initial) == "tgv3d") then
                                blk%q(i,j,k,1,b) = sin(kk*blk%x(i,1,b))*cos(kk*blk%y(j,1,b))*cos(kk*blk%z(k,1,b))
                                blk%q(i,j,k,2,b) = cos(kk*blk%x(i,2,b))*sin(kk*blk%y(j,2,b))*cos(kk*blk%z(k,2,b))
                                blk%q(i,j,k,3,b) = cos(kk*blk%x(i,3,b))*cos(kk*blk%y(j,3,b))*sin(kk*blk%z(k,3,b))
                            else
                                blk%q(i,j,k,1,b) =  sin(kk*blk%x(i,1,b))*cos(kk*blk%y(j,1,b))
                                blk%q(i,j,k,2,b) = -cos(kk*blk%x(i,2,b))*sin(kk*blk%y(j,2,b))
                                blk%q(i,j,k,3,b) = 0.0d0
                            end if
                        end do
                    end do
                end do
            end do
            return
        end if

        blk%q(:,:,:,1,:) = dns%initial_velocity(1)
        blk%q(:,:,:,2,:) = dns%initial_velocity(2)
        blk%q(:,:,:,3,:) = dns%initial_velocity(3)

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

    subroutine generic_after_step(this, blk, dns, g, c, ibm)
        class(generic_case_type), intent(inout) :: this
        type(block_set_type), intent(inout) :: blk
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        type(comm_type), intent(in) :: c
        type(ibm_type), intent(in) :: ibm
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
        ! Damped-Jacobi relaxation factor: the projection DIVERGES for sor > 0.8
        ! (the retired red-black SOR used 1.5). See pressure_solver_type%omega.
        ps%omega = 0.8d0
    end subroutine set_generic_defaults

end module generic_flow
