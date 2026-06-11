! TODO:
!
!      (*) second-order IBM in space, then in time


program main
    use :: init
    use :: blocks, only: block_set_type, init_block_set, destroy_block_set, &
        block_set_matches_grid
    use :: chron, only: chron_type, start_chron, stop_chron, write_chron
    use :: flow_case, only: case_type, create_flow_case
    use :: config
    use :: boundary
    use :: io
    use :: step
    use :: pressure_solver
    use :: gpu_runtime
    use :: ibmm
    use :: les_model
    use :: comm, only: comm_type, comm_init_world, comm_init, comm_finalize, exchange_halos, exchange_scalar_halos
    implicit none

    integer :: arg_status, rkStage
    integer(C_INT) :: loop_steps
    real(C_DOUBLE) :: dt_alpha, dt_beta, dt_gamma
    real(C_DOUBLE) :: les_profile_start
    character(len=256) :: input_file
    type(chron_type) :: loop_timer
    class(case_type), allocatable :: flow
    type(dns_type) :: dns
    type(grid_type) :: g
    type(block_set_type) :: blk
    type(field_type) :: f
    type(boundary_type) :: bc
    type(pressure_solver_type) :: ps
    type(ibm_type) :: ibm
    type(les_type) :: les
    type(les_profile_type) :: les_prof
    type(config_seen_type) :: config_seen
    type(comm_type) :: c

    call comm_init_world(c)
    call splash(c%has_terminal)

    ! The first command-line argument can override the default input file.
    call get_command_argument(1, input_file, status=arg_status)
    if (arg_status /= 0 .or. len_trim(input_file) == 0) input_file = "input.ini"

    if (c%has_terminal) print *, "reading input data..."
    call create_flow_case(flow, input_file, c%has_terminal)
    call flow%apply_defaults(dns, g, bc, c, ps)
    call read_runtime_config(dns, g, les, ps, bc, c, input_file, c%has_terminal, config_seen)
    if (has_restart_file(dns)) then
        if (c%has_terminal) print *, "reading restart metadata: ", trim(dns%restart_file)
        call read_restart_metadata(dns, g, bc, ps%nIter, ps%sor, dns%restart_file, c, &
            preserve_cflmax=config_seen%cflmax, preserve_pecletmax=config_seen%pecletmax, &
            preserve_dtmax=config_seen%dtmax, preserve_t_final=config_seen%t_final)
    end if
    call comm_init(c, dns, bc)

    if (c%has_terminal) print *, "initialising grid..."
    call init_grid(g, dns, bc%isPeriodic)
    call validate_dns_values(dns, g)

    ! Phase 0 of the block refactor (docs/block_refinement_strategy.md):
    ! build the one-block-per-rank block set and check that it reproduces the
    ! grid metrics bit-for-bit. Temporary scaffold; the block set replaces
    ! grid_type/field_type as the solver's data layout in the next phase.
    call init_block_set(blk, dns, g, bc%isPeriodic, global_id=int(c%cart_rank, C_INT))
    if (.not. block_set_matches_grid(blk, g)) then
        error stop "phase-0 block set does not reproduce the grid metrics"
    end if
    call destroy_block_set(blk)
    call precompute_peclet_rate(dns, g, c)
    call init_boundary_faces(bc, dns, g)
    call init_openmp_offload(c%has_terminal)
    call enter_grid_data(g, dns)
    call enter_boundary_data(bc)

    if (c%has_terminal) print *, "initialising fields..."
    call init_field(f, dns)
    if (.not. has_restart_file(dns)) then
        call flow%initialise_fields(f, dns, g, bc, c)
    end if
    call enter_field_data(f, dns)
    if (has_restart_file(dns)) then
        if (c%has_terminal) print *, "reading restart fields: ", trim(dns%restart_file)
        call read_field(f, dns, dns%restart_file, c)
    end if

    if (c%has_terminal) print *, "initialising pressure solver..."
    call init_pressure_solver(ps, dns, bc, c%has_terminal)

    if (c%has_terminal) print *, "initialising IBM..."
    call init_ibm(ibm, dns)
    if (dns%ibm_enabled .and. len_trim(dns%ibm_coeff_file) > 0) then
        call read_ibm_coeff_file(ibm, dns, c%has_terminal)
        call enter_ibm_data(ibm, dns)
    else
        call enter_ibm_data(ibm, dns)
        call set_ibm_coeff(dns, g, ibm, VAR_U)
        call set_ibm_coeff(dns, g, ibm, VAR_V)
        call set_ibm_coeff(dns, g, ibm, VAR_W)
    end if

    if (les_is_enabled(les)) then
        if (c%has_terminal) print *, "initialising LES model..."
        call init_les(les, dns, g)
        call enter_les_data(les, dns)
    end if

    call apply_bc(f, dns, g, bc)
    call exchange_halos(c, f, [VAR_U, VAR_V, VAR_W, VAR_P])
    call flow%setup_after_grid(f, dns, g, bc, c)
    if (les_is_enabled(les)) then
        call update_les_viscosity(les, f, dns, g, ibm)
        call exchange_scalar_halos(c, les%nut)
        call update_timestep_limits(f, dns, g, c, les)
    else
        call update_timestep_limits(f, dns, g, c)
    end if

    if (c%has_terminal) print *, "main loop starting..."
    loop_steps = 0_C_INT
    call reset_les_profile(les_prof)
    call start_chron(loop_timer)
    do while (run_should_continue(dns, loop_steps))
        call trim_dt_for_final_time(dns)
        if (dns%dt <= 0.0d0) exit

        loop_steps = loop_steps + 1_C_INT
        dns%step_current = dns%step_current + 1_C_INT
        dns%t_current = dns%t_current + dns%dt

        do rkStage = 1,3
            dt_alpha = dns%dt*rk_alpha(rkStage)
            dt_beta  = dns%dt*rk_beta(rkStage)
            dt_gamma = dns%dt*rk_gamma(rkStage)

            ! Predictor: advance tentative staggered velocities, then enforce solid/body constraints.
            call update_ibm_mu(ibm, dt_gamma)
            if (les_is_enabled(les)) then
                les_profile_start = les_wall_seconds()
                call update_les_viscosity(les, f, dns, g, ibm)
                call add_les_profile(les_prof, LES_PROF_NUT, les_wall_seconds() - les_profile_start)
                les_profile_start = les_wall_seconds()
                call exchange_scalar_halos(c, les%nut)
                call add_les_profile(les_prof, LES_PROF_EXCHANGE, les_wall_seconds() - les_profile_start)
                call momentum(f, dns, g, dt_alpha, dt_beta, dt_gamma, ibm, bc, les, les_prof)
            else
                call momentum(f, dns, g, dt_alpha, dt_beta, dt_gamma, ibm, bc)
            end if
            call apply_bc(f, dns, g, bc)
            call exchange_halos(c, f, [VAR_U, VAR_V, VAR_W])

            ! Projection: solve for pressure correction and project tentative velocities.
            call pressure_projection(ps, f, dns, g, dt_gamma, ibm, bc, c)

        end do

        if (les_is_enabled(les)) then
            call update_timestep_limits(f, dns, g, c, les)
        else
            call update_timestep_limits(f, dns, g, c)
        end if

        if (dns%field_interval > 0) then
            call maybe_write_field(f, dns, g, int(dns%step_current), c, bc, ps%nIter, ps%sor)
        end if
        call flow%after_step(f, dns, g, c)

    end do
    call stop_chron(loop_timer, loop_steps)

    if (c%has_terminal) then
        print *, "main loop ended..."
        call write_chron(loop_timer)
        if (les_is_enabled(les)) call write_les_profile(les_prof, loop_steps)
    end if

    call write_field(f, dns, g, int(dns%step_current), c, bc, ps%nIter, ps%sor)

    ! Release device-side data before the host allocatables go out of scope.
    call flow%finalize(dns, g, c)
    if (les_is_enabled(les)) call exit_les_data(les, dns)
    call destroy_les(les)
    call exit_ibm_data(ibm, dns)
    call exit_field_data(f, dns)
    call exit_boundary_data(bc)
    call exit_grid_data(g, dns)
    call destroy_grid(g)
    call destroy_boundary_faces(bc)
    call comm_finalize(c)
end program main
