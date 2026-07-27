! TODO:
!
!      (*) second-order IBM in space, then in time


program moby_solve
    use :: init
    use :: blocks, only: block_set_type, init_block_set, destroy_block_set, &
        enter_block_data, exit_block_data, zero_closed_halos, &
        FACE_FINE, FACE_PHYS, FACE_CLOSED
    use :: chron, only: chron_type, start_chron, stop_chron, write_chron, &
        profiler_type, wall_seconds, profiler_add, write_profiler
    use :: flow_case, only: case_type, create_flow_case
    use :: config
    use :: boundary
    use :: io
    use :: step
    use :: pressure_solver
    use :: gpu_runtime
    use :: ibmm
    use :: turbulence
    use :: les_model
    use :: rans
    use :: boostconv, only: boostconv_type
    use :: bodyforce
    use :: comm, only: comm_type, comm_init_world, comm_init, comm_finalize, &
        init_block_exchange, exchange_halos, exchange_scalar_halos, &
        comm_allreduce_sum, comm_allreduce_max
    implicit none

    integer :: arg_status, rkStage
    integer(C_INT) :: loop_steps
    real(C_DOUBLE) :: dt_alpha, dt_beta, dt_gamma
    real(C_DOUBLE) :: turb_profile_start
    character(len=256) :: input_file
    type(chron_type) :: loop_timer
    class(case_type), allocatable :: flow
    type(dns_type) :: dns
    type(grid_type) :: g
    type(block_set_type) :: blk
    type(boundary_type) :: bc
    type(pressure_solver_type) :: ps
    type(ibm_type) :: ibm
    type(turb_type) :: turb
    type(les_type) :: les
    type(sst_type) :: sst
    type(profiler_type) :: turb_prof
    type(bodyforce_type) :: bf
    type(boostconv_type) :: bcv
    type(config_seen_type) :: config_seen
    type(comm_type) :: c
    integer(C_INT), allocatable :: blockActive(:)
    integer(C_INT), allocatable :: blockTouch(:,:), blockBuried(:,:)
    integer(C_INT), allocatable :: blockMaskLo(:,:), blockMaskDims(:,:)

    call comm_init_world(c)
    call splash(c%has_terminal)

    ! The first command-line argument can override the default input file.
    call get_command_argument(1, input_file, status=arg_status)
    if (arg_status /= 0 .or. len_trim(input_file) == 0) input_file = "input.ini"

    if (c%has_terminal) print *, "reading input data..."
    call create_flow_case(flow, input_file, c%has_terminal)
    call flow%apply_defaults(dns, g, bc, c, ps)
    call read_runtime_config(dns, g, turb, les, ps, bc, c, input_file, c%has_terminal, config_seen)
    if (has_restart_file(dns)) then
        if (c%has_terminal) print *, "reading restart metadata: ", trim(dns%restart_file)
        call read_restart_metadata(dns, g, bc, ps%nIter, ps%omega, dns%restart_file, c, config_seen)
    end if
    call comm_init(c, dns, bc)

    ! STL geometry is prepared offline (docs/prepare_solve_strategy.md):
    ! the solver consumes only the case file it produced.
    if (dns%ibm_enabled .and. dns%ibm_stl_count > 0_C_INT &
            .and. len_trim(dns%ibm_coeff_file) == 0) then
        if (c%has_terminal) print *, "error: [ibm] stl_file is a moby_prepare input;", &
            " run moby_prepare and point [ibm] coeff_file at its output"
        error stop
    end if

    if (c%has_terminal) print *, "initialising grid..."
    call init_grid(g, dns, bc%isPeriodic)
    call validate_dns_values(dns, g)

    ! Block refactor (docs/block_refinement_strategy.md): the solver state
    ! lives in a block set tiling the grid ([blocks] nb per block). With an
    ! immersed boundary, blocks buried inside the body are removed from the
    ! global table before the set is built.
    if (dns%block_nb > 0_C_INT .and. dns%block_refine_body) then
        ! Geometry-driven refinement (analytic or file IBM): refine to the
        ! finest level at the surface with a one-block buffer, removing
        ! buried blocks at every level. ibmm produces the geometry masks.
        ! The analytic classification is rank-split and exactly merged
        ! (P2) -- identical masks on any rank count.
        call classify_refinement_masks(blockTouch, blockBuried, blockMaskLo, &
            blockMaskDims, dns, g, ibm, bc%isPeriodic, c%has_terminal, isInBody, c=c)
        call init_block_set(blk, dns, g, bc%isPeriodic, int(c%cart_size, C_INT), &
            int(c%cart_rank, C_INT), touch=blockTouch, buried=blockBuried, &
            maskLo=blockMaskLo, maskDims=blockMaskDims)
        deallocate(blockTouch, blockBuried, blockMaskLo, blockMaskDims)
    else if (dns%block_nb > 0_C_INT .and. dns%ibm_enabled .and. dns%block_remove_solid) then
        ! Solid-block removal: drop blocks buried inside the immersed body.
        call classify_active_mask(blockActive, dns, g, ibm, bc%isPeriodic, &
            c%has_terminal, isInBody, c=c)
        call init_block_set(blk, dns, g, bc%isPeriodic, int(c%cart_size, C_INT), &
            int(c%cart_rank, C_INT), blockActive)
        deallocate(blockActive)
    else
        call init_block_set(blk, dns, g, bc%isPeriodic, int(c%cart_size, C_INT), &
            int(c%cart_rank, C_INT))
    end if
    call init_block_exchange(c, blk, dns)
    call precompute_peclet_rate(dns, blk, c)
    call init_boundary_faces(bc, blk, dns)
    call init_openmp_offload(c%has_terminal)
    call enter_grid_data(g)
    call enter_boundary_data(bc)

    if (c%has_terminal) print *, "initialising fields..."
    if (.not. has_restart_file(dns)) then
        call flow%initialise_fields(blk, dns, g, bc, c)
    end if
    call enter_block_data(blk)
    if (has_restart_file(dns)) then
        if (c%has_terminal) print *, "reading restart fields: ", trim(dns%restart_file)
        call read_field(blk, dns, dns%restart_file, c)
    end if
    call zero_closed_halos(blk)

    if (c%has_terminal) print *, "initialising pressure solver..."
    call init_pressure_solver(ps, dns, bc, c%has_terminal)

    if (c%has_terminal) print *, "initialising IBM..."
    call init_ibm(ibm, blk)
    if (dns%ibm_enabled .and. len_trim(dns%ibm_coeff_file) > 0) then
        call read_ibm_coeff_file(ibm, dns, blk, c%has_terminal)
        call enter_ibm_data(ibm, dns)
    else
        call enter_ibm_data(ibm, dns)
        call set_ibm_coeff(dns, blk, ibm, VAR_U)
        call set_ibm_coeff(dns, blk, ibm, VAR_V)
        call set_ibm_coeff(dns, blk, ibm, VAR_W)
    end if
    ! [ibm] band_filter: near-body band list from the device coefficients
    ! (off: nothing is built, allocated, or mapped).
    if (dns%ibm_enabled .and. dns%ibm_band_filter) &
        call init_ibm_band(ibm, dns, blk, c, c%has_terminal)

    ! A configured [rans] section builds the SST geometry state (T1: wall
    ! distance + IBM wall cells); [turbulence] model = rans additionally
    ! builds and advances the k-omega transport state (T2).
    if (dns%rans_configured) then
        if (c%has_terminal) print *, "initialising RANS geometry..."
#ifdef USE_OPENMP_OFFLOAD
        ! The analytic IBM coefficients are computed on the device; the
        ! wall-cell classification reads the host copy.
        !$omp target update from(ibm%coef)
#endif
        call init_rans_geometry(sst, dns, g, blk, bc, ibm, c)
        if (turb%model == TURB_RANS .or. turb%model == TURB_IDDES) &
            call init_rans_transport(sst, dns, blk, bc, ibm, c%has_terminal)
        call enter_rans_data(sst)
        ! [rans] boostconv: steady-state accelerator (RANS transport only)
        if (dns%rans_boostconv) then
            if (turb%model /= TURB_RANS) &
                error stop "[rans] boostconv requires [turbulence] model = rans"
            call init_boostconv_rans(bcv, blk, dns, c%has_terminal)
        end if
        if (dns%rans_dump_geometry) call write_rans_geometry(sst, blk, dns, c)
    end if

    if (turbulence_is_enabled(turb)) then
        if (c%has_terminal) print *, "initialising turbulence model..."
        call init_turbulence(turb, blk)
        ! Full IDDES: the static geometric fields (f_B, f_e1, the mesh
        ! length) need the RANS wall distance; fill them on the host before
        ! the device maps.
        if (turb%model == TURB_IDDES) &
            call init_iddes_geometry(turb, blk, sst%dwall, sst%yeff)
        call enter_turbulence_data(turb)
    end if

    if (dns%force_enabled) then
        if (c%has_terminal) print *, "initialising body force..."
        call init_bodyforce(bf, dns, blk, g, c%has_terminal)
        call enter_bodyforce_data(bf)
    end if

    call apply_bc(blk, bc, outflow_copy=.true.)
    call exchange_halos(c, blk, [VAR_U, VAR_V, VAR_W, VAR_P])

    call flow%setup_after_grid(blk, dns, g, bc, c)
    ! Every producer fills the same turb%nut; the consumer chain (exchange,
    ! dt limits) is model-agnostic.
    if (turbulence_is_enabled(turb)) then
        select case (turb%model)
        case (TURB_RANS)
            call rans_prepare(sst, turb, blk, dns, bc, c)
        case (TURB_IDDES)
            ! IDDES: SGS viscosity into the scratch, the DDES shielding fd
            ! (reads the lagged blended nut), the RANS pass (l_hyb rides
            ! turb%fd inside the kernel), then the blend.
            call update_sgs_viscosity(les, turb, blk, dns, ibm, turb%nut_sgs)
            call compute_iddes_fd(turb, blk, dns, sst%yeff)
            call rans_prepare(sst, turb, blk, dns, bc, c)
            call blend_iddes_nut(turb, blk)
        case default
            call update_sgs_viscosity(les, turb, blk, dns, ibm, turb%nut)
        end select
        call exchange_scalar_halos(c, turb%nut, blk)
        call update_timestep_limits(blk, dns, c, turb)
    else
        call update_timestep_limits(blk, dns, c)
    end if

    if (c%has_terminal) print *, "main loop starting..."
    loop_steps = 0_C_INT
    call init_turbulence_profiler(turb_prof)
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
            ! Body force: refresh the (custom) force for this substage; profile
            ! and file sources are filled once at init and this is a no-op.
            if (dns%force_enabled) call update_bodyforce(bf, blk, dns, g, dns%t_current)
            if (turbulence_is_enabled(turb)) then
                turb_profile_start = wall_seconds()
                select case (turb%model)
                case (TURB_RANS)
                    call rans_substage(sst, turb, blk, dns, ibm, bc, c, dt_alpha, dt_beta)
                case (TURB_IDDES)
                    call update_sgs_viscosity(les, turb, blk, dns, ibm, turb%nut_sgs)
                    call compute_iddes_fd(turb, blk, dns, sst%yeff)
                    call rans_substage(sst, turb, blk, dns, ibm, bc, c, dt_alpha, dt_beta)
                    call blend_iddes_nut(turb, blk)
                case default
                    call update_sgs_viscosity(les, turb, blk, dns, ibm, turb%nut)
                end select
                call profiler_add(turb_prof, TURB_PROF_NUT, wall_seconds() - turb_profile_start)
                turb_profile_start = wall_seconds()
                call exchange_scalar_halos(c, turb%nut, blk)
                call profiler_add(turb_prof, TURB_PROF_EXCHANGE, wall_seconds() - turb_profile_start)
                call momentum(blk, dns, dt_alpha, dt_beta, dt_gamma, ibm, turb, turb_prof, bf=bf)
            else
                call momentum(blk, dns, dt_alpha, dt_beta, dt_gamma, ibm, bf=bf)
            end if
            call apply_bc(blk, bc, outflow_copy=.true.)
            ! Post-predictor exchange with the conservation SYNC: the cross-level
            ! PROLONG/RESTRICT write the shared 2:1 face so the two stored copies
            ! start the projection mean-consistent (avg(fine)=coarse). The
            ! projection then owns the face and the composite stencil keeps it
            ! conservative.
            call exchange_halos(c, blk, [VAR_U, VAR_V, VAR_W], syncface=.true.)

            ! Projection: solve for pressure correction and project tentative velocities.
            call pressure_projection(ps, blk, dns, dt_gamma, ibm, bc, c)

        end do

        ! BoostConv activation (end of step, every interval-th step):
        ! replaces the state by the recombined one; the following steps
        ! absorb it like a restart (halos/BCs refreshed inside).
        if (dns%rans_boostconv) then
            if (mod(int(dns%step_current), int(bcv%interval)) == 0) &
                call boostconv_rans(bcv, sst, blk, dns, bc, c, c%has_terminal)
        end if

        if (turbulence_is_enabled(turb)) then
            call update_timestep_limits(blk, dns, c, turb)
        else
            call update_timestep_limits(blk, dns, c)
        end if

        if (dns%field_interval > 0) then
            call maybe_write_field(blk, dns, g, int(dns%step_current), c, bc, ps%nIter, ps%omega, &
                turb%nut, sst%k, sst%omg, sst%gam, sst%ret, turb%fd)
        end if
        call flow%after_step(blk, dns, g, c, ibm)

    end do
    call stop_chron(loop_timer, loop_steps)

    if (c%has_terminal) then
        print *, "main loop ended..."
        call write_chron(loop_timer)
        if (turbulence_is_enabled(turb)) call write_profiler(turb_prof, loop_steps)
    end if

    call write_field(blk, dns, g, int(dns%step_current), c, bc, ps%nIter, ps%omega, &
        turb%nut, sst%k, sst%omg, sst%gam, sst%ret, turb%fd)

    ! Release device-side data before the host allocatables go out of scope.
    call flow%finalize(dns, g, c)
    if (dns%force_enabled) call exit_bodyforce_data(bf)
    call destroy_bodyforce(bf)
    if (turbulence_is_enabled(turb)) call exit_turbulence_data(turb)
    call destroy_turbulence(turb)
    if (dns%rans_configured) call exit_rans_data(sst)
    call destroy_rans_geometry(sst)
    call exit_ibm_data(ibm, dns)
    call exit_block_data(blk)
    call exit_boundary_data(bc)
    call exit_grid_data(g)
    call destroy_block_set(blk)
    call destroy_grid(g)
    call destroy_boundary_faces(bc)
    call comm_finalize(c)
end program moby_solve
