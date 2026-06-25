! TODO:
!
!      (*) second-order IBM in space, then in time


program main
    use :: init
    use :: blocks, only: block_set_type, init_block_set, destroy_block_set, &
        enter_block_data, exit_block_data, zero_closed_halos, &
        FACE_FINE, FACE_PHYS, FACE_CLOSED
    use :: chron, only: chron_type, start_chron, stop_chron, write_chron
    use :: flow_case, only: case_type, create_flow_case
    use :: config
    use :: boundary
    use :: io
    use :: step
    use :: pressure_solver
    use :: pressure_solver, only: proj_t_exch, proj_t_ker
    use :: gpu_runtime
    use :: ibmm
    use :: les_model
    use :: comm, only: comm_type, comm_init_world, comm_init, comm_finalize, &
        init_block_exchange, exchange_halos, exchange_scalar_halos, &
        comm_allreduce_sum, comm_allreduce_max
    implicit none

    integer :: arg_status, rkStage, refluxDir, refluxComp
    integer(C_INT) :: loop_steps
    real(C_DOUBLE) :: dt_alpha, dt_beta, dt_gamma
    real(C_DOUBLE) :: les_profile_start
    character(len=256) :: input_file
    type(chron_type) :: loop_timer
    class(case_type), allocatable :: flow
    type(dns_type) :: dns
    type(grid_type) :: g
    type(block_set_type) :: blk
    type(boundary_type) :: bc
    type(pressure_solver_type) :: ps
    type(ibm_type) :: ibm
    type(les_type) :: les
    type(les_profile_type) :: les_prof
    type(config_seen_type) :: config_seen
    type(comm_type) :: c
    integer(C_INT), allocatable :: blockActive(:)
    integer(C_INT), allocatable :: blockTouch(:,:), blockBuried(:,:)
    integer :: refineLevel, maskCount
    logical :: blockActiveFound
    ! Interface-diagnosis toggles (off by default). MOBY_PROJONLY: skip the
    ! momentum predictor so the projection acts on the exact div-free initial
    ! field -- any change it makes IS the correction/interface defect. MOBY_PREDONLY:
    ! skip the projection so the field is the pure predictor. Either dumps the
    ! initial field (900000) so the change can be diffed.
    logical :: projOnly, predOnly, divdump, rhsdump, noRecon, stepDiv, phaseTime, keBal, skewIface
    ! MOBY_PHASETIME: accumulate wall time per major loop phase (target regions are
    ! synchronous, so host timers capture GPU time) and print the breakdown at the
    ! end -- to see where the refinement cost goes. Off by default.
    real(C_DOUBLE) :: pt0, pt1, tRecon, tReflux, tMom, tExch, tProj, tTstep
    real(C_DOUBLE) :: tRfC, tRfE, tRfA
    ! MOBY_IFFILT=<alpha>: coefficient of the coarse-interface-band tangential
    ! filter (fix i). 0 (default) = off / bit-exact.
    real(C_DOUBLE) :: ifFiltAlpha
    character(len=16) :: diagEnv
    ! MOBY_STEPDIV: per-timestep divergence monitor -- after the full RK step, print
    ! the final field's divergence L2 (rms) and max, and the global mass residual
    ! (Sum vol*div). Tracks whether an under-converged projection lets divergence
    ! accumulate / blow up. Off by default.
    real(C_DOUBLE), allocatable :: divStepBuf(:,:,:,:)
    ! MOBY_DIVDUMP: write the discrete (raw staggered) divergence the solver sees,
    ! BEFORE the last projection (div_pre = D u*, e.g. of the exact field under
    ! MOBY_PROJONLY) and AFTER (div_post = D u_after), as companion field files.
    ! div_post tests operator CONSISTENCY directly: a consistent projection leaves
    ! ~0 residual divergence; an O(h) interface residual means D and G are not
    ! consistent. Off by default; the normal path is bit-exact.
    real(C_DOUBLE), allocatable :: divPre(:,:,:,:), divPost(:,:,:,:)
    character(len=256) :: divBasePrefix
    ! MOBY_RHSDUMP: dump the discrete momentum RHS L_h(u) = -div(uu) + (1/Re)lap(u)
    ! (the predictor's blk%oldrhs, captured at the first RK substage so it acts
    ! on the pristine manufactured field) as un/vn/wn of a "<prefix>_rhs" field.
    ! tools/rhsband.py compares it to the analytic L(u) of the tgv3d field at the
    ! staggered points and reports the interface-band convergence order -- the
    ! momentum-operator-truncation gate. Run with MOBY_PREDONLY (no projection)
    ! and initial=tgv3d; the Re attr lets the script split advection vs diffusion.
    real(C_DOUBLE), allocatable :: rhsBuf(:,:,:,:,:)
    character(len=256) :: rhsBasePrefix

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
            preserve_dtmax=config_seen%dtmax, preserve_t_final=config_seen%t_final, &
            preserve_pressure_niter=config_seen%pressure_niter, &
            preserve_pressure_sor=config_seen%pressure_sor)
    end if
    call comm_init(c, dns, bc)

    if (c%has_terminal) print *, "initialising grid..."
    call init_grid(g, dns, bc%isPeriodic)
    call validate_dns_values(dns, g)

    ! Block refactor (docs/block_refinement_strategy.md): the solver state
    ! lives in a block set tiling the grid ([blocks] nb per block). With an
    ! immersed boundary, blocks buried inside the body are removed from the
    ! global table before the set is built.
    if (dns%block_nb > 0_C_INT .and. dns%block_refine_body) then
        ! Geometry-driven refinement (analytic IBM): refine to the finest
        ! level at the surface with a one-block buffer, removing buried
        ! blocks at every level.
        if (.not. dns%ibm_enabled) then
            error stop "[blocks] refine_body needs the IBM enabled"
        end if
        if (any(mod(dns%globalSize, dns%block_nb) /= 0_C_INT)) then
            error stop "[blocks] nb must divide the global grid in every direction"
        end if
        allocate(blockTouch(product(dns%globalSize/dns%block_nb)*8**dns%block_refine_levels, &
            dns%block_refine_levels + 1))
        allocate(blockBuried(size(blockTouch,1), size(blockTouch,2)))
        if (len_trim(dns%ibm_coeff_file) > 0) then
            ! File-based geometry: masks computed by mobygeom block-table.
            do refineLevel = 0, int(dns%block_refine_levels)
                maskCount = int(product(dns%globalSize/dns%block_nb))*(8**refineLevel)
                call read_block_masks(blockTouch(1:maskCount, refineLevel+1), &
                    blockBuried(1:maskCount, refineLevel+1), refineLevel, maskCount, &
                    blockActiveFound, dns, c%has_terminal)
                if (.not. blockActiveFound) then
                    error stop "coefficient file has no refinement masks; run mobygeom block-table"
                end if
            end do
        else
            call classify_block_geometry(blockTouch, blockBuried, dns, g, ibm, bc%isPeriodic, &
                int(dns%block_refine_levels) + 1)
        end if
        call init_block_set(blk, dns, g, bc%isPeriodic, int(c%cart_size, C_INT), &
            int(c%cart_rank, C_INT), touch=blockTouch, buried=blockBuried)
        deallocate(blockTouch, blockBuried)
    else if (dns%block_nb > 0_C_INT .and. dns%ibm_enabled .and. dns%block_remove_solid) then
        if (any(mod(dns%globalSize, dns%block_nb) /= 0_C_INT)) then
            error stop "[blocks] nb must divide the global grid in every direction"
        end if
        if (dns%block_refine_nboxes > 0_C_INT) then
            error stop "solid-block removal with box refinement is unsupported; use refine_body"
        end if
        allocate(blockActive(product(dns%globalSize/dns%block_nb)))
        if (len_trim(dns%ibm_coeff_file) > 0) then
            call read_block_active(blockActive, blockActiveFound, dns, c%has_terminal)
            if (.not. blockActiveFound) then
                if (c%has_terminal) print *, &
                    "coefficient file has no block_active table; keeping all blocks"
                blockActive = 1_C_INT
            end if
        else
            call classify_active_blocks(blockActive, dns, g, ibm, bc%isPeriodic)
        end if
        call init_block_set(blk, dns, g, bc%isPeriodic, int(c%cart_size, C_INT), &
            int(c%cart_rank, C_INT), blockActive)
        deallocate(blockActive)
    else
        call init_block_set(blk, dns, g, bc%isPeriodic, int(c%cart_size, C_INT), &
            int(c%cart_rank, C_INT))
    end if
    call init_block_exchange(c, blk, dns)
    call precompute_peclet_rate(dns, blk, c)
    call init_boundary_faces(bc, blk)
    call init_openmp_offload(c%has_terminal)
    call enter_grid_data(g, dns)
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

    if (les_is_enabled(les)) then
        if (c%has_terminal) print *, "initialising LES model..."
        call init_les(les, dns, blk)
        call enter_les_data(les, dns)
    end if

    call apply_bc(blk, bc)
    call exchange_halos(c, blk, [VAR_U, VAR_V, VAR_W, VAR_P])

    ! Temporary debugging hook: manufactured-field halo coherence audit.
    block
        character(len=16) :: auditEnv
        call get_environment_variable("MOBY_HALO_AUDIT", auditEnv)
        if (len_trim(auditEnv) > 0) then
            call halo_audit(blk, dns, bc, c)
            stop
        end if
    end block

    ! MOBY_TERMDUMP=<var>: dump each momentum term of component <var> (1=u,2=v,
    ! 3=w) SEPARATELY, computed from the post-IC-exchanged field exactly as the
    ! predictor would -- the three advection flux-divergence terms (x/y/z) as
    ! un/vn/wn of "<prefix>_adv", the three Laplacian terms (d2/dx2, d2/dy2,
    ! d2/dz2, times 1/Re) as un/vn/wn of "<prefix>_dif". tools/rhsterms.py
    ! compares each term to its analytic value at the staggered point and splits
    ! by interface band/row, so a broken interface stencil is attributed to the
    ! exact term (e.g. the wall-normal d(vv)/dy across a 2:1 face).
    block
        character(len=16) :: termEnv
        integer :: termVar
        real(C_DOUBLE), allocatable :: advT(:,:,:,:,:), difT(:,:,:,:,:)
        character(len=256) :: termBase
        call get_environment_variable("MOBY_TERMDUMP", termEnv)
        if (len_trim(termEnv) > 0) then
            read(termEnv, *) termVar
            ! Match the predictor: reconstruct the interface deep halos so the
            ! dumped terms reflect what momentum actually evaluates.
            call reconstruct_interface_halos(blk)
            allocate(advT(1:blk%nb(1),1:blk%nb(2),1:blk%nb(3),NVEL,blk%nBlocks))
            allocate(difT(1:blk%nb(1),1:blk%nb(2),1:blk%nb(3),NVEL,blk%nBlocks))
            call compute_momentum_terms(blk, dns, termVar, advT, difT)
            termBase = dns%field_prefix
            call overwrite_velocity(blk, advT)
            dns%field_prefix = trim(termBase)//"_adv"
            call write_field(blk, dns, g, int(dns%step_current), c, bc, ps%nIter, ps%sor)
            call overwrite_velocity(blk, difT)
            dns%field_prefix = trim(termBase)//"_dif"
            call write_field(blk, dns, g, int(dns%step_current), c, bc, ps%nIter, ps%sor)
            dns%field_prefix = termBase
            stop
        end if
    end block

    call flow%setup_after_grid(blk, dns, g, bc, c)
    if (les_is_enabled(les)) then
        call update_les_viscosity(les, blk, dns, ibm)
        call exchange_scalar_halos(c, les%nut, blk)
        call update_timestep_limits(blk, dns, c, les)
    else
        call update_timestep_limits(blk, dns, c)
    end if

    ! Interface-diagnosis toggles. Dump the (exact) initial field as 900000 so a
    ! projection-only / predictor-only run can be diffed against it.
    call get_environment_variable("MOBY_PROJONLY", diagEnv)
    projOnly = len_trim(diagEnv) > 0
    call get_environment_variable("MOBY_PREDONLY", diagEnv)
    predOnly = len_trim(diagEnv) > 0
    call get_environment_variable("MOBY_DIVDUMP", diagEnv)
    divdump = len_trim(diagEnv) > 0
    call get_environment_variable("MOBY_RHSDUMP", diagEnv)
    rhsdump = len_trim(diagEnv) > 0
    call get_environment_variable("MOBY_NORECON", diagEnv)
    noRecon = len_trim(diagEnv) > 0   ! debug: skip reconstruct_interface_halos
    call get_environment_variable("MOBY_STEPDIV", diagEnv)
    stepDiv = len_trim(diagEnv) > 0   ! per-step divergence L2/max + global mass monitor
    call get_environment_variable("MOBY_PHASETIME", diagEnv)
    phaseTime = len_trim(diagEnv) > 0
    call get_environment_variable("MOBY_KEBAL", diagEnv)
    keBal = len_trim(diagEnv) > 0   ! per-step convective KE-balance: band vs interior
    call get_environment_variable("MOBY_KESKEW", diagEnv)
    skewIface = len_trim(diagEnv) > 0   ! skew-symmetric (energy-conserving) interface convection
    call get_environment_variable("MOBY_IFFILT", diagEnv)
    ifFiltAlpha = 0.0d0
    if (len_trim(diagEnv) > 0) read(diagEnv, *) ifFiltAlpha
    tRecon = 0.0d0; tReflux = 0.0d0; tMom = 0.0d0; tExch = 0.0d0; tProj = 0.0d0; tTstep = 0.0d0
    tRfC = 0.0d0; tRfE = 0.0d0; tRfA = 0.0d0
    ! MOBY_MANUF=<amp>: add a pure-gradient perturbation amp*grad(phi),
    ! phi = cos(kx)cos(ky)cos(kz), to the (exact, div-free) initial field.
    ! The result is globally mass-conserving (a periodic gradient integrates
    ! to zero) but NOT locally divergence-free. A consistent, conservative
    ! projection must remove grad(phi) exactly -- recovering u_exact -- while
    ! keeping the global mass error at round-off, INCLUDING at the 2:1
    ! interface. Combine with MOBY_PROJONLY (project the manufactured field)
    ! and MOBY_DIVDUMP (Sum vol*div before/after); diff the final field vs
    ! the 900000 dump's exact part with tools/check_beltrami.py.
    call get_environment_variable("MOBY_MANUF", diagEnv)
    if (len_trim(diagEnv) > 0) then
        block
            real(C_DOUBLE) :: amp
            read(diagEnv, *) amp
            call manufacture_gradient_perturbation(blk, dns, amp)
        end block
    end if
    if (projOnly .or. predOnly) &
        call write_field(blk, dns, g, 900000, c, bc, ps%nIter, ps%sor)

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
            ! Skipped under MOBY_PROJONLY (the projection then acts on the exact field).
            if (.not. projOnly) then
                ! Reconstruct the velocity deep halos across each 2:1 interface so
                ! the predictor's advection/diffusion reaching into them (normal
                ! and tangential) are 2nd order (inert without an interface).
                if (phaseTime) pt0 = les_wall_seconds()
                if (.not. noRecon) call reconstruct_interface_halos(blk)
                if (phaseTime) tRecon = tRecon + les_wall_seconds() - pt0
                ! Momentum reflux: from the start-of-substage velocity, capture
                ! each interface direction's normal advective flux, restrict the
                ! fine flux into the coarse across-interface halo, and accumulate
                ! the coarse RHS correction; applied to the predicted field below.
                if (phaseTime) pt0 = les_wall_seconds()
                if (dns%block_momentum_reflux) then
                    call reflux_zero(blk)
                    do refluxDir = 1, 3
                        do refluxComp = 1, 3
                            if (phaseTime) pt1 = les_wall_seconds()
                            call reflux_compute_flux(blk, refluxDir, refluxComp)
                            if (phaseTime) tRfC = tRfC + les_wall_seconds() - pt1
                            if (phaseTime) pt1 = les_wall_seconds()
                            call exchange_scalar_halos(c, blk%refluxF, blk, ifaceRow=.true.)
                            if (phaseTime) tRfE = tRfE + les_wall_seconds() - pt1
                            if (phaseTime) pt1 = les_wall_seconds()
                            call reflux_accumulate(blk, refluxDir, refluxComp)
                            if (phaseTime) tRfA = tRfA + les_wall_seconds() - pt1
                        end do
                    end do
                    ! Energy-conserving (V&V) skew-symmetric interface convection:
                    ! add 1/2 u (div u) at the interface-band cells INTO refluxCorr
                    ! so reflux_apply lands it with the reflux correction. Gated by
                    ! MOBY_KESKEW; inert without a 2:1 interface.
                    if (skewIface) call skew_interface_correction(blk)
                end if
                if (phaseTime) tReflux = tReflux + les_wall_seconds() - pt0
                call update_ibm_mu(ibm, dt_gamma)
                if (les_is_enabled(les)) then
                    les_profile_start = les_wall_seconds()
                    call update_les_viscosity(les, blk, dns, ibm)
                    call add_les_profile(les_prof, LES_PROF_NUT, les_wall_seconds() - les_profile_start)
                    les_profile_start = les_wall_seconds()
                    call exchange_scalar_halos(c, les%nut, blk)
                    call add_les_profile(les_prof, LES_PROF_EXCHANGE, les_wall_seconds() - les_profile_start)
                    call momentum(blk, dns, dt_alpha, dt_beta, dt_gamma, ibm, les, les_prof)
                else
                    if (phaseTime) pt0 = les_wall_seconds()
                    call momentum(blk, dns, dt_alpha, dt_beta, dt_gamma, ibm)
                    if (phaseTime) tMom = tMom + les_wall_seconds() - pt0
                end if
                if (phaseTime) pt0 = les_wall_seconds()
                if (dns%block_momentum_reflux) call reflux_apply(blk, dt_alpha)
                if (phaseTime) tReflux = tReflux + les_wall_seconds() - pt0
                call apply_bc(blk, bc)
                ! Post-predictor exchange with the conservation SYNC: the
                ! cross-level PROLONG/RESTRICT write the shared 2:1 face so the
                ! two stored copies start the projection mean-consistent
                ! (avg(fine)=coarse). The projection then owns the face and the
                ! composite stencil keeps it conservative.
                if (phaseTime) pt0 = les_wall_seconds()
                call exchange_halos(c, blk, [VAR_U, VAR_V, VAR_W], syncface=.true.)
                if (phaseTime) tExch = tExch + les_wall_seconds() - pt0
                ! FIX (i): damp the coarse-interface-band tangential velocity, then
                ! refresh halos so the projection sees the filtered field. Inert at
                ! alpha=0.
                if (ifFiltAlpha /= 0.0d0) then
                    call filter_interface_band(blk, ifFiltAlpha)
                    call exchange_halos(c, blk, [VAR_U, VAR_V, VAR_W])
                end if
            end if

            ! Projection: solve for pressure correction and project tentative velocities.
            ! Skipped under MOBY_PREDONLY (the field is then the pure predictor).
            ! On the last substage capture the divergence the solver sees before
            ! (D u*) and after (D u_after) the projection for MOBY_DIVDUMP.
            ! Capture at the FIRST substage: divPre = D of the field entering the
            ! first projection (the exact field under MOBY_PROJONLY), divPost = D
            ! after that one projection. (Capturing later would measure D of an
            ! already multiply-projected field.)
            ! MOBY_RHSDUMP: capture L_h(u) of the pristine field (first substage,
            ! before dt_beta couples a previous RHS or a later substage mutates q).
            if (rhsdump .and. rkStage == 1) call capture_rhs(blk, rhsBuf)
            if (divdump .and. rkStage == 1) call capture_divergence(blk, divPre)
            if (phaseTime) pt0 = les_wall_seconds()
            if (.not. predOnly) &
                call pressure_projection(ps, blk, dns, dt_gamma, ibm, bc, c)
            if (phaseTime) tProj = tProj + les_wall_seconds() - pt0
            if (divdump .and. rkStage == 1) call capture_divergence(blk, divPost)

        end do

        if (stepDiv) call print_step_divergence(blk, c, dns, divStepBuf)
        if (keBal) call print_step_ke_balance(blk, c, dns, noRecon)

        if (phaseTime) pt0 = les_wall_seconds()
        if (les_is_enabled(les)) then
            call update_timestep_limits(blk, dns, c, les)
        else
            call update_timestep_limits(blk, dns, c)
        end if
        if (phaseTime) tTstep = tTstep + les_wall_seconds() - pt0

        if (dns%field_interval > 0) then
            call maybe_write_field(blk, dns, g, int(dns%step_current), c, bc, ps%nIter, ps%sor)
        end if
        call flow%after_step(blk, dns, g, c)

    end do
    call stop_chron(loop_timer, loop_steps)

    if (c%has_terminal) then
        print *, "main loop ended..."
        call write_chron(loop_timer)
        if (les_is_enabled(les)) call write_les_profile(les_prof, loop_steps)
        if (phaseTime .and. loop_steps > 0) then
            pt0 = tRecon + tReflux + tMom + tExch + tProj + tTstep
            print '(a,i0,a)', "PHASETIME (", int(loop_steps), " steps, ms/step | % of timed):"
            print '(a,f9.3,a,f5.1,a)', "  reconstruct  ", 1.0d3*tRecon /loop_steps, " | ", 100*tRecon /max(pt0,1d-30), "%"
            print '(a,f9.3,a,f5.1,a)', "  reflux       ", 1.0d3*tReflux/loop_steps, " | ", 100*tReflux/max(pt0,1d-30), "%"
            print '(a,f9.3,a,f9.3,a,f9.3,a)', "    rf.compute ", 1.0d3*tRfC/loop_steps, "  rf.exch ", 1.0d3*tRfE/loop_steps, "  rf.accum ", 1.0d3*tRfA/loop_steps, " (ms/step)"
            print '(a,f9.3,a,f5.1,a)', "  momentum     ", 1.0d3*tMom   /loop_steps, " | ", 100*tMom   /max(pt0,1d-30), "%"
            print '(a,f9.3,a,f5.1,a)', "  exchange(vel)", 1.0d3*tExch  /loop_steps, " | ", 100*tExch  /max(pt0,1d-30), "%"
            print '(a,f9.3,a,f5.1,a)', "  projection   ", 1.0d3*tProj  /loop_steps, " | ", 100*tProj  /max(pt0,1d-30), "%"
            print '(a,f9.3,a,f9.3,a)', "    proj.exch  ", 1.0d3*proj_t_exch/loop_steps, "    proj.kernels ", 1.0d3*proj_t_ker/loop_steps, " (ms/step)"
            print '(a,f9.3,a,f5.1,a)', "  timestep_lim ", 1.0d3*tTstep /loop_steps, " | ", 100*tTstep /max(pt0,1d-30), "%"
            print '(a,f9.3)',          "  timed total  ", 1.0d3*pt0/loop_steps
        end if
    end if

    call write_field(blk, dns, g, int(dns%step_current), c, bc, ps%nIter, ps%sor)

    ! MOBY_DIVDUMP companions: park the divergence in the pressure slot and write
    ! it as separate field files (blk%q is not reused after the final dump).
    if (divdump .and. allocated(divPost)) then
        divBasePrefix = dns%field_prefix
        call overwrite_var_p(blk, divPre)
        dns%field_prefix = trim(divBasePrefix)//"_divpre"
        call write_field(blk, dns, g, int(dns%step_current), c, bc, ps%nIter, ps%sor)
        call overwrite_var_p(blk, divPost)
        dns%field_prefix = trim(divBasePrefix)//"_divpost"
        call write_field(blk, dns, g, int(dns%step_current), c, bc, ps%nIter, ps%sor)
        dns%field_prefix = divBasePrefix
    end if

    ! MOBY_RHSDUMP companion: park the captured momentum RHS in the velocity
    ! slots and write it as "<prefix>_rhs" (un/vn/wn = rhsu/rhsv/rhsw).
    if (rhsdump .and. allocated(rhsBuf)) then
        rhsBasePrefix = dns%field_prefix
        call overwrite_velocity(blk, rhsBuf)
        dns%field_prefix = trim(rhsBasePrefix)//"_rhs"
        call write_field(blk, dns, g, int(dns%step_current), c, bc, ps%nIter, ps%sor)
        dns%field_prefix = rhsBasePrefix
    end if

    ! Release device-side data before the host allocatables go out of scope.
    call flow%finalize(dns, g, c)
    if (les_is_enabled(les)) call exit_les_data(les, dns)
    call destroy_les(les)
    call exit_ibm_data(ibm, dns)
    call exit_block_data(blk)
    call exit_boundary_data(bc)
    call exit_grid_data(g, dns)
    call destroy_block_set(blk)
    call destroy_grid(g)
    call destroy_boundary_faces(bc)
    call comm_finalize(c)

contains

    ! Manufactured-field halo audit (temporary debug tool). Fills every
    ! block interior with a per-variable linear field, sentinels the halos,
    ! runs one exchange and checks every exchange-written halo cell against
    ! the value the transfer DESIGN should produce: copies and restrictions
    ! of a linear field reproduce the field at the halo location; prolonged
    ! halos reproduce the field at the covering coarse cell/face location
    ! (pressure faces: blended with the first interior cell, 2:1 uniform
    ! weight 2/3). Wrapped periodic regions are skipped (linear field is
    ! discontinuous across the wrap), as are regions with no occupant.
    subroutine halo_audit(blk, dns, bc, c)
        use :: blocks, only: leaf_at, level_cells
        type(block_set_type), intent(inout) :: blk
        type(dns_type), intent(in) :: dns
        type(boundary_type), intent(in) :: bc
        type(comm_type), intent(inout) :: c

        real(C_DOUBLE), parameter :: SENTINEL = 1.0d33
        real(C_DOUBLE), parameter :: CXV(4) = [1.0d0, 0.3d0, 0.7d0, 1.3d0]
        real(C_DOUBLE), parameter :: CYV(4) = [2.0d0, 1.1d0, 0.2d0, 0.5d0]
        real(C_DOUBLE), parameter :: CZV(4) = [0.9d0, 1.7d0, 2.3d0, 0.4d0]
        integer(C_INT) :: b, var, i, j, k, d, l, nBad, nChecked, nb(3)
        integer(C_INT) :: off(3), idx(3), to(3), cl(3), cc(3), gnl, gidx, cov
        integer(C_INT) :: idSame, idParent, sx, sy, sz, op
        logical :: skip, anyChild
        real(C_DOUBLE) :: got, want, posD(3), srcD(3), wBlend
        ! breakdown: bad count per variable and per transfer op (1=copy 2=prolong
        ! 3=restrict), max abs error per op, total checked per op.
        integer(C_INT) :: badVar(4), badOp(3), chkOp(3)
        integer(C_INT) :: nSentinel, sentByOff(3), nSentPrint
        real(C_DOUBLE) :: maxErrOp(3)
        badVar = 0; badOp = 0; chkOp = 0; maxErrOp = 0.0d0; nSentinel = 0
        sentByOff = 0; nSentPrint = 0

        nb = blk%nb

        do b = 1, int(blk%nBlocks)
            do var = 1, 4
                do k = 0, nb(3) + 1
                    do j = 0, nb(2) + 1
                        do i = 0, nb(1) + 1
                            if (i >= 1 .and. i <= nb(1) .and. j >= 1 .and. j <= nb(2) &
                                .and. k >= 1 .and. k <= nb(3)) then
                                blk%q(i,j,k,var,b) = CXV(var)*blk%x(i,var,b) &
                                    + CYV(var)*blk%y(j,var,b) + CZV(var)*blk%z(k,var,b)
                            else
                                blk%q(i,j,k,var,b) = SENTINEL
                            end if
                        end do
                    end do
                end do
            end do
        end do
#ifdef USE_OPENMP_OFFLOAD
        !$omp target update to(blk%q)
#endif
        call exchange_halos(c, blk, [VAR_U, VAR_V, VAR_W, VAR_P])
#ifdef USE_OPENMP_OFFLOAD
        !$omp target update from(blk%q)
#endif

        nBad = 0
        nChecked = 0
        do b = 1, int(blk%nBlocks)
            l = int(blk%level(b))
            do var = 1, 4
                do k = 0, nb(3) + 1
                    do j = 0, nb(2) + 1
                        do i = 0, nb(1) + 1
                            idx = [i, j, k]
                            do d = 1, 3
                                off(d) = 0
                                if (idx(d) == 0) off(d) = -1
                                if (idx(d) == nb(d) + 1) off(d) = 1
                            end do
                            if (all(off == 0)) cycle

                            ! Neighbour region; skip wraps and walls.
                            skip = .false.
                            do d = 1, 3
                                gnl = level_cells(dns, d, int(l, C_INT))
                                to(d) = int(blk%origin(d,b)) + off(d)*nb(d)
                                if (to(d) < 0 .or. to(d) >= gnl) skip = .true.
                            end do
                            if (skip) cycle

                            cl = to/nb
                            idSame = int(leaf_at(blk, int(l, C_INT), int(cl, C_INT)))
                            idParent = -1
                            if (l > 0) idParent = int(leaf_at(blk, int(l - 1, C_INT), &
                                int(cl/2, C_INT)))

                            ! The fine block OWNS its 2:1 interface normal-velocity
                            ! face: the momentum predictor writes it, not the
                            ! exchange (comm.f90 lNrm / interface_normal_dim skips
                            ! var == the +1 face-normal dim on an interface face).
                            ! So the audit must not expect it written either.
                            if (idSame < 0 .and. sum(abs(off)) == 1 .and. var <= 3) then
                                if (off(var) == 1) cycle
                            end if

                            posD = [blk%x(i,var,b), blk%y(j,var,b), blk%z(k,var,b)]
                            if (idSame >= 0) then
                                want = CXV(var)*posD(1) + CYV(var)*posD(2) + CZV(var)*posD(3)
                            else if (idParent >= 0 .and. sum(abs(off)) == 1) then
                                ! FACE prolong: tangential interpolation + normal
                                ! (2C+F)/3 blend reproduces a linear field at the
                                ! fine halo location EXACTLY (the new correct design).
                                want = CXV(var)*posD(1) + CYV(var)*posD(2) + CZV(var)*posD(3)
                            else if (idParent >= 0) then
                                ! EDGE/CORNER prolong (no blend): still injection.
                                do d = 1, 3
                                    gnl = level_cells(dns, d, int(l, C_INT))
                                    gidx = int(blk%origin(d,b)) + idx(d) - 1
                                    if (var == d) then
                                        cov = modulo(gidx, gnl)/2
                                        srcD(d) = line_at(blk, d, l - 1, cov)
                                    else
                                        cov = modulo(gidx, gnl)/2
                                        srcD(d) = 0.5d0*(line_at(blk, d, l - 1, cov) &
                                                       + line_at(blk, d, l - 1, cov + 1))
                                    end if
                                end do
                                want = CXV(var)*srcD(1) + CYV(var)*srcD(2) + CZV(var)*srcD(3)
                            else
                                ! Finer occupants: restriction of a linear field
                                ! reproduces it at the halo location; skip if the
                                ! region is empty.
                                anyChild = .false.
                                do sz = 0, 1
                                do sy = 0, 1
                                do sx = 0, 1
                                    cc = 2*cl + [sx, sy, sz]
                                    if (leaf_at(blk, int(l + 1, C_INT), int(cc, C_INT)) >= 0) &
                                        anyChild = .true.
                                end do
                                end do
                                end do
                                if (.not. anyChild) cycle
                                want = CXV(var)*posD(1) + CYV(var)*posD(2) + CZV(var)*posD(3)
                            end if

                            got = blk%q(i,j,k,var,b)
                            nChecked = nChecked + 1
                            op = 1
                            if (idSame < 0) op = merge(2, 3, idParent >= 0)
                            chkOp(op) = chkOp(op) + 1
                            if (abs(got - want) > 1.0d-10*max(1.0d0, abs(want))) then
                                nBad = nBad + 1
                                badVar(var) = badVar(var) + 1
                                badOp(op) = badOp(op) + 1
                                if (abs(got) > 1.0d30) then
                                    nSentinel = nSentinel + 1   ! halo never written
                                    sentByOff(sum(abs(off))) = sentByOff(sum(abs(off))) + 1
                                    if (nSentPrint < 16) then
                                        nSentPrint = nSentPrint + 1
                                        print '(a,i6,a,i2,a,i2,a,3i3,a,3i4,a,i2)', &
                                            " UNWRITTEN b=", b, " lvl=", l, " var=", var, &
                                            " off=", off, " ijk=", i, j, k, " op=", op
                                    end if
                                else
                                    maxErrOp(op) = max(maxErrOp(op), abs(got - want))
                                end if
                                if (nBad <= 40) then
                                    print '(a,i6,a,i2,a,i2,a,3i4,a,3i5,a,2es22.14)', &
                                        " AUDIT BAD b=", b, " lvl=", l, " var=", var, &
                                        " ijk=", i, j, k, " orig=", blk%origin(:,b), &
                                        " got/want=", got, want
                                end if
                            end if
                        end do
                    end do
                end do
            end do
        end do
        print *, "HALO AUDIT: checked", nChecked, " bad", nBad
        print '(a,4i9)', " HALO AUDIT bad per var (u,v,w,p): ", badVar
        print '(a,3i9)', " HALO AUDIT checked per op (copy,prolong,restrict): ", chkOp
        print '(a,3i9)', " HALO AUDIT bad     per op (copy,prolong,restrict): ", badOp
        print '(a,3es12.4)', " HALO AUDIT maxErr  per op (copy,prolong,restrict): ", maxErrOp
        print '(a,i9)', " HALO AUDIT of which UN-WRITTEN (sentinel): ", nSentinel
        print '(a,3i9)', " HALO AUDIT un-written by off-type (face,edge,corner): ", sentByOff
    end subroutine halo_audit

    ! MOBY_DIVDUMP: raw staggered divergence of the current velocity (the exact D
    ! the Jacobi sweep minimises) into a separate buffer; interior cells, halos
    ! zeroed. Allocates + maps on first use; leaves blk%q untouched.
    subroutine capture_divergence(blk, arr)
        type(block_set_type), intent(in) :: blk
        real(C_DOUBLE), allocatable, intent(inout) :: arr(:,:,:,:)
        integer(C_INT) :: i, j, k, b, nx, ny, nz
        nx = blk%nb(1); ny = blk%nb(2); nz = blk%nb(3)
        if (.not. allocated(arr)) then
            allocate(arr(lbound(blk%q,1):ubound(blk%q,1), &
                         lbound(blk%q,2):ubound(blk%q,2), &
                         lbound(blk%q,3):ubound(blk%q,3), blk%nBlocks))
#ifdef USE_OPENMP_OFFLOAD
            !$omp target enter data map(alloc: arr)
#endif
        end if
#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do collapse(4) &
        !$omp& map(to: blk%q, blk%d1x, blk%d1y, blk%d1z) map(tofrom: arr) private(i,j,k,b)
#endif
        do b = 1_C_INT, blk%nBlocks
            do k = lbound(arr,3), ubound(arr,3)
                do j = lbound(arr,2), ubound(arr,2)
                    do i = lbound(arr,1), ubound(arr,1)
                        if (i >= 1 .and. i <= nx .and. j >= 1 .and. j <= ny &
                            .and. k >= 1 .and. k <= nz) then
                            arr(i,j,k,b) = &
                                (blk%q(i+1,j,k,VAR_U,b)-blk%q(i,j,k,VAR_U,b))*blk%d1x(i,VAR_P,b) &
                              + (blk%q(i,j+1,k,VAR_V,b)-blk%q(i,j,k,VAR_V,b))*blk%d1y(j,VAR_P,b) &
                              + (blk%q(i,j,k+1,VAR_W,b)-blk%q(i,j,k,VAR_W,b))*blk%d1z(k,VAR_P,b)
                        else
                            arr(i,j,k,b) = 0.0d0
                        end if
                    end do
                end do
            end do
        end do
#ifdef USE_OPENMP_OFFLOAD
        !$omp end target teams distribute parallel do
#endif
    end subroutine capture_divergence

    ! Per-timestep divergence monitor (MOBY_STEPDIV). After the full RK step, form
    ! the discrete divergence of the final velocity and reduce it to: L2 = rms over
    ! all interior cells, MAX = max|div|, and MASS = Sum(vol*div) (the global mass
    ! residual; round-off for a conservative scheme on a periodic box). Reductions
    ! are MPI-global. cell volume = 1/(d1x*d1y*d1z) on the pressure metric (per
    ! block, so correct across levels).
    subroutine print_step_divergence(blk, c, dns, divBuf)
        type(block_set_type), intent(inout) :: blk
        type(comm_type), intent(in) :: c
        type(dns_type), intent(in) :: dns
        real(C_DOUBLE), allocatable, intent(inout) :: divBuf(:,:,:,:)
        integer(C_INT) :: i, j, k, b, nx, ny, nz
        real(C_DOUBLE) :: sumsq, masssum, ncells, maxabs, vol, d
        real(C_DOUBLE) :: red(3), redmax(1)

        call capture_divergence(blk, divBuf)
        nx = blk%nb(1); ny = blk%nb(2); nz = blk%nb(3)
        sumsq = 0.0d0; masssum = 0.0d0; ncells = 0.0d0; maxabs = 0.0d0
#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do collapse(4) &
        !$omp& map(to: divBuf, blk%d1x, blk%d1y, blk%d1z) &
        !$omp& reduction(+:sumsq,masssum,ncells) reduction(max:maxabs) &
        !$omp& private(i,j,k,b,vol,d)
#endif
        do b = 1_C_INT, blk%nBlocks
            do k = 1, nz
                do j = 1, ny
                    do i = 1, nx
                        d = divBuf(i,j,k,b)
                        vol = 1.0d0/(blk%d1x(i,VAR_P,b)*blk%d1y(j,VAR_P,b)*blk%d1z(k,VAR_P,b))
                        sumsq = sumsq + d*d
                        masssum = masssum + vol*d
                        ncells = ncells + 1.0d0
                        maxabs = max(maxabs, abs(d))
                    end do
                end do
            end do
        end do
#ifdef USE_OPENMP_OFFLOAD
        !$omp end target teams distribute parallel do
#endif
        red = [sumsq, masssum, ncells]
        call comm_allreduce_sum(c, red)
        redmax = [maxabs]
        call comm_allreduce_max(c, redmax)
        if (c%has_terminal) then
            if (red(3) > 0.0d0) red(1) = sqrt(red(1)/red(3))
            print '(a,i7,a,es12.5,a,es12.5,a,es12.5)', &
                "STEPDIV step", int(dns%step_current), &
                "  div_l2=", red(1), "  div_max=", redmax(1), "  mass=", red(2)
        end if
    end subroutine print_step_divergence

    ! Per-timestep discrete kinetic-energy balance from CONVECTION (MOBY_KEBAL).
    ! The convective operator is energy-conserving (V&V 2003) iff it is
    ! skew-symmetric, in which case the convective KE production Sum_DOF vol*u*C(u)
    ! is identically zero for ANY field. This monitor forms the convective RHS (no
    ! pressure/diffusion/forcing) exactly as `momentum` does and reports its KE
    ! production split COARSE interface band (cells in a FACE_FINE row) vs interior,
    ! in TWO forms:
    !
    !  * DIV  = vol*u*C_div, the divergence-form production the scheme actually uses.
    !    For the div form the production telescopes to -Sum KE*(div u), so it is
    !    contaminated by the projection's residual divergence (large when niter is
    !    small or the interface projection is under-converged) -- NOT a clean
    !    interface signal.
    !  * SKEW = vol*( u*C_div + 1/2 u^2 * Div_cv ), where Div_cv is the discrete
    !    velocity divergence on the component's own control volume (same face
    !    neighbours as the convective fluxes). This is the energy production of the
    !    skew-symmetric form C_skew = C_div + 1/2 u (div u); it is identically zero
    !    in the interior for ANY field (constant-1/2 telescoping, divergence-state
    !    INDEPENDENT), so the band SKEW value is the pure interface energy defect.
    !    This is the pass/fail: the energy-conserving interface fix must drive the
    !    band SKEW production to round-off.
    !
    ! Halos are refreshed (exchange + interface reconstruction, matching the
    ! predictor's view) before the stencil; this mutates only halo cells, which the
    ! next substage re-fills. cell volume = 1/(d1*d1*d1) on each component's metric.
    subroutine print_step_ke_balance(blk, c, dns, noRecon)
        type(block_set_type), intent(inout) :: blk
        type(comm_type), intent(inout) :: c
        type(dns_type), intent(in) :: dns
        logical, intent(in) :: noRecon
        integer(C_INT) :: i, j, k, b, ip, im, jp, jm, kp, km, nx, ny, nz
        integer(C_INT) :: uStartX, vStartY, wStartZ
        real(C_DOUBLE) :: uu_p,uu_m,uv_p,uv_m,uw_p,uw_m
        real(C_DOUBLE) :: vu_p,vu_m,vv_p,vv_m,vw_p,vw_m
        real(C_DOUBLE) :: wu_p,wu_m,ww_p,ww_m,wv_p,wv_m
        real(C_DOUBLE) :: convU, convV, convW, prodD, prodS, divcv, vol, qd
        real(C_DOUBLE) :: keBandD, keIntD, keBandS, keIntS, absBandS, absIntS, red(6)
        logical :: bandU, bandV, bandW

        ! Match the predictor's halo state so the interface stencil reads the same
        ! restricted/reconstructed neighbours the momentum predictor reads.
        call exchange_halos(c, blk, [VAR_U, VAR_V, VAR_W])
        if (.not. noRecon) call reconstruct_interface_halos(blk)

        nx = blk%nb(1); ny = blk%nb(2); nz = blk%nb(3)
        keBandD = 0.0d0; keIntD = 0.0d0; keBandS = 0.0d0
        keIntS = 0.0d0; absBandS = 0.0d0; absIntS = 0.0d0
#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do collapse(4) &
        !$omp& map(to: blk%q, blk%d1x, blk%d1y, blk%d1z, blk%physLow, blk%physHigh) &
        !$omp& reduction(+:keBandD,keIntD,keBandS,keIntS,absBandS,absIntS) &
        !$omp& private(i,j,k,b,ip,im,jp,jm,kp,km,uStartX,vStartY,wStartZ, &
        !$omp& uu_p,uu_m,uv_p,uv_m,uw_p,uw_m,vu_p,vu_m,vv_p,vv_m,vw_p,vw_m, &
        !$omp& wu_p,wu_m,ww_p,ww_m,wv_p,wv_m,convU,convV,convW,prodD,prodS,divcv,vol,qd, &
        !$omp& bandU,bandV,bandW)
#endif
        do b = 1, blk%nBlocks
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    ip = i+1; im = i-1; jp = j+1; jm = j-1; kp = k+1; km = k-1
                    uStartX = merge(2, 1, blk%physLow(1,b) == FACE_PHYS .or. blk%physLow(1,b) == FACE_CLOSED)
                    vStartY = merge(2, 1, blk%physLow(2,b) == FACE_PHYS .or. blk%physLow(2,b) == FACE_CLOSED)
                    wStartZ = merge(2, 1, blk%physLow(3,b) == FACE_PHYS .or. blk%physLow(3,b) == FACE_CLOSED)
                    ! "interface band" = a cell in the interface row of any 2:1 face,
                    ! BOTH sides (coarse FACE_FINE and fine FACE_COARSE): the energy
                    ! defect lives on both, and only excluding both leaves a truly
                    ! interface-free interior (SKEW -> round-off). The handout's
                    ! coarse-cell band is the FACE_FINE subset.
                    bandU = ((blk%physLow(1,b) == FACE_FINE .or. blk%physLow(1,b) == FACE_COARSE) .and. i == 1) .or. &
                            ((blk%physHigh(1,b) == FACE_FINE .or. blk%physHigh(1,b) == FACE_COARSE) .and. i == nx) .or. &
                            ((blk%physLow(2,b) == FACE_FINE .or. blk%physLow(2,b) == FACE_COARSE) .and. j == 1) .or. &
                            ((blk%physHigh(2,b) == FACE_FINE .or. blk%physHigh(2,b) == FACE_COARSE) .and. j == ny) .or. &
                            ((blk%physLow(3,b) == FACE_FINE .or. blk%physLow(3,b) == FACE_COARSE) .and. k == 1) .or. &
                            ((blk%physHigh(3,b) == FACE_FINE .or. blk%physHigh(3,b) == FACE_COARSE) .and. k == nz)
                    bandV = bandU; bandW = bandU

                    if (i >= uStartX) then
                        uu_p = (blk%q(i,j,k,VAR_U,b) + blk%q(ip,j,k,VAR_U,b))**2
                        uu_m = (blk%q(im,j,k,VAR_U,b) + blk%q(i,j,k,VAR_U,b))**2
                        uv_p = (blk%q(i,j,k,VAR_U,b) + blk%q(i,jp,k,VAR_U,b)) &
                             * (blk%q(im,jp,k,VAR_V,b) + blk%q(i,jp,k,VAR_V,b))
                        uv_m = (blk%q(i,jm,k,VAR_U,b) + blk%q(i,j,k,VAR_U,b)) &
                             * (blk%q(im,j,k,VAR_V,b) + blk%q(i,j,k,VAR_V,b))
                        uw_p = (blk%q(i,j,k,VAR_U,b) + blk%q(i,j,kp,VAR_U,b)) &
                             * (blk%q(im,j,kp,VAR_W,b) + blk%q(i,j,kp,VAR_W,b))
                        uw_m = (blk%q(i,j,km,VAR_U,b) + blk%q(i,j,k,VAR_U,b)) &
                             * (blk%q(im,j,k,VAR_W,b) + blk%q(i,j,k,VAR_W,b))
                        convU = -0.25d0*( (uu_p-uu_m)*blk%d1x(i,VAR_U,b) &
                                         +(uv_p-uv_m)*blk%d1y(j,VAR_U,b) &
                                         +(uw_p-uw_m)*blk%d1z(k,VAR_U,b))
                        ! Velocity divergence on the u-control-volume (same face
                        ! averages as the convective fluxes above): the skew-form
                        ! correction term that telescopes the interior to round-off.
                        divcv = 0.5d0*( (blk%q(ip,j,k,VAR_U,b)-blk%q(im,j,k,VAR_U,b))*blk%d1x(i,VAR_U,b) &
                            + ((blk%q(im,jp,k,VAR_V,b)+blk%q(i,jp,k,VAR_V,b)) &
                              -(blk%q(im,j,k,VAR_V,b)+blk%q(i,j,k,VAR_V,b)))*blk%d1y(j,VAR_U,b) &
                            + ((blk%q(im,j,kp,VAR_W,b)+blk%q(i,j,kp,VAR_W,b)) &
                              -(blk%q(im,j,k,VAR_W,b)+blk%q(i,j,k,VAR_W,b)))*blk%d1z(k,VAR_U,b))
                        vol = 1.0d0/(blk%d1x(i,VAR_U,b)*blk%d1y(j,VAR_U,b)*blk%d1z(k,VAR_U,b))
                        qd = blk%q(i,j,k,VAR_U,b)
                        prodD = vol*qd*convU
                        prodS = prodD + 0.5d0*vol*qd*qd*divcv
                        if (bandU) then
                            keBandD = keBandD + prodD
                            keBandS = keBandS + prodS; absBandS = absBandS + abs(prodS)
                        else
                            keIntD = keIntD + prodD
                            keIntS = keIntS + prodS; absIntS = absIntS + abs(prodS)
                        end if
                    end if

                    if (j >= vStartY) then
                        vu_p = (blk%q(i,j,k,VAR_V,b) + blk%q(ip,j,k,VAR_V,b)) &
                             * (blk%q(ip,jm,k,VAR_U,b) + blk%q(ip,j,k,VAR_U,b))
                        vu_m = (blk%q(im,j,k,VAR_V,b) + blk%q(i,j,k,VAR_V,b)) &
                             * (blk%q(i,jm,k,VAR_U,b) + blk%q(i,j,k,VAR_U,b))
                        vv_p = (blk%q(i,j,k,VAR_V,b) + blk%q(i,jp,k,VAR_V,b))**2
                        vv_m = (blk%q(i,jm,k,VAR_V,b) + blk%q(i,j,k,VAR_V,b))**2
                        vw_p = (blk%q(i,j,k,VAR_V,b) + blk%q(i,j,kp,VAR_V,b)) &
                             * (blk%q(i,jm,kp,VAR_W,b) + blk%q(i,j,kp,VAR_W,b))
                        vw_m = (blk%q(i,j,km,VAR_V,b) + blk%q(i,j,k,VAR_V,b)) &
                             * (blk%q(i,jm,k,VAR_W,b) + blk%q(i,j,k,VAR_W,b))
                        convV = -0.25d0*( (vu_p-vu_m)*blk%d1x(i,VAR_V,b) &
                                         +(vv_p-vv_m)*blk%d1y(j,VAR_V,b) &
                                         +(vw_p-vw_m)*blk%d1z(k,VAR_V,b))
                        divcv = 0.5d0*( ((blk%q(ip,jm,k,VAR_U,b)+blk%q(ip,j,k,VAR_U,b)) &
                              -(blk%q(i,jm,k,VAR_U,b)+blk%q(i,j,k,VAR_U,b)))*blk%d1x(i,VAR_V,b) &
                            + (blk%q(i,jp,k,VAR_V,b)-blk%q(i,jm,k,VAR_V,b))*blk%d1y(j,VAR_V,b) &
                            + ((blk%q(i,jm,kp,VAR_W,b)+blk%q(i,j,kp,VAR_W,b)) &
                              -(blk%q(i,jm,k,VAR_W,b)+blk%q(i,j,k,VAR_W,b)))*blk%d1z(k,VAR_V,b))
                        vol = 1.0d0/(blk%d1x(i,VAR_V,b)*blk%d1y(j,VAR_V,b)*blk%d1z(k,VAR_V,b))
                        qd = blk%q(i,j,k,VAR_V,b)
                        prodD = vol*qd*convV
                        prodS = prodD + 0.5d0*vol*qd*qd*divcv
                        if (bandV) then
                            keBandD = keBandD + prodD
                            keBandS = keBandS + prodS; absBandS = absBandS + abs(prodS)
                        else
                            keIntD = keIntD + prodD
                            keIntS = keIntS + prodS; absIntS = absIntS + abs(prodS)
                        end if
                    end if

                    if (k >= wStartZ) then
                        wu_p = (blk%q(i,j,k,VAR_W,b) + blk%q(ip,j,k,VAR_W,b)) &
                             * (blk%q(ip,j,km,VAR_U,b) + blk%q(ip,j,k,VAR_U,b))
                        wu_m = (blk%q(im,j,k,VAR_W,b) + blk%q(i,j,k,VAR_W,b)) &
                             * (blk%q(i,j,km,VAR_U,b) + blk%q(i,j,k,VAR_U,b))
                        ww_p = (blk%q(i,j,k,VAR_W,b) + blk%q(i,j,kp,VAR_W,b))**2
                        ww_m = (blk%q(i,j,km,VAR_W,b) + blk%q(i,j,k,VAR_W,b))**2
                        wv_p = (blk%q(i,j,k,VAR_W,b) + blk%q(i,jp,k,VAR_W,b)) &
                             * (blk%q(i,jp,km,VAR_V,b) + blk%q(i,jp,k,VAR_V,b))
                        wv_m = (blk%q(i,jm,k,VAR_W,b) + blk%q(i,j,k,VAR_W,b)) &
                             * (blk%q(i,j,km,VAR_V,b) + blk%q(i,j,k,VAR_V,b))
                        convW = -0.25d0*( (wu_p-wu_m)*blk%d1x(i,VAR_W,b) &
                                         +(wv_p-wv_m)*blk%d1y(j,VAR_W,b) &
                                         +(ww_p-ww_m)*blk%d1z(k,VAR_W,b))
                        divcv = 0.5d0*( ((blk%q(ip,j,km,VAR_U,b)+blk%q(ip,j,k,VAR_U,b)) &
                              -(blk%q(i,j,km,VAR_U,b)+blk%q(i,j,k,VAR_U,b)))*blk%d1x(i,VAR_W,b) &
                            + ((blk%q(i,jp,km,VAR_V,b)+blk%q(i,jp,k,VAR_V,b)) &
                              -(blk%q(i,j,km,VAR_V,b)+blk%q(i,j,k,VAR_V,b)))*blk%d1y(j,VAR_W,b) &
                            + (blk%q(i,j,kp,VAR_W,b)-blk%q(i,j,km,VAR_W,b))*blk%d1z(k,VAR_W,b))
                        vol = 1.0d0/(blk%d1x(i,VAR_W,b)*blk%d1y(j,VAR_W,b)*blk%d1z(k,VAR_W,b))
                        qd = blk%q(i,j,k,VAR_W,b)
                        prodD = vol*qd*convW
                        prodS = prodD + 0.5d0*vol*qd*qd*divcv
                        if (bandW) then
                            keBandD = keBandD + prodD
                            keBandS = keBandS + prodS; absBandS = absBandS + abs(prodS)
                        else
                            keIntD = keIntD + prodD
                            keIntS = keIntS + prodS; absIntS = absIntS + abs(prodS)
                        end if
                    end if
                end do
            end do
        end do
        end do
#ifdef USE_OPENMP_OFFLOAD
        !$omp end target teams distribute parallel do
#endif
        red = [keBandD, keIntD, keBandS, keIntS, absBandS, absIntS]
        call comm_allreduce_sum(c, red)
        if (c%has_terminal) then
            ! SKEW band_net is the pass/fail (pure interface energy defect, should
            ! go to round-off); DIV is the actual divergence-form budget.
            print '(a,i7,a,es12.5,a,es10.3,a,a,es12.5,a,es10.3,a)', &
                "KEBAL step", int(dns%step_current), &
                "  SKEW band=", red(3), " (|.|=", red(5), ")", &
                "  int=", red(4), " (|.|=", red(6), ")"
            print '(a,es12.5,a,es12.5,a,es12.5)', &
                "           DIV  band=", red(1), "  int=", red(2), &
                "  total=", red(1)+red(2)
        end if
    end subroutine print_step_ke_balance

    ! MOBY_DIVDUMP: park a scalar buffer in the pressure slot so write_field emits
    ! it as "pn" (its target update from(blk%q) then pulls these values).
    subroutine overwrite_var_p(blk, arr)
        type(block_set_type), intent(inout) :: blk
        real(C_DOUBLE), allocatable, intent(in) :: arr(:,:,:,:)
        integer(C_INT) :: i, j, k, b
#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do collapse(4) &
        !$omp& map(to: arr) map(tofrom: blk%q) private(i,j,k,b)
#endif
        do b = 1_C_INT, blk%nBlocks
            do k = lbound(arr,3), ubound(arr,3)
                do j = lbound(arr,2), ubound(arr,2)
                    do i = lbound(arr,1), ubound(arr,1)
                        blk%q(i,j,k,VAR_P,b) = arr(i,j,k,b)
                    end do
                end do
            end do
        end do
#ifdef USE_OPENMP_OFFLOAD
        !$omp end target teams distribute parallel do
#endif
    end subroutine overwrite_var_p

    ! MOBY_RHSDUMP: copy the predictor's discrete momentum RHS (blk%oldrhs, the
    ! interior-only L_h(u) the time integrator stores) into a separate buffer.
    ! Allocates + maps on first use; leaves blk%oldrhs untouched.
    subroutine capture_rhs(blk, arr)
        type(block_set_type), intent(in) :: blk
        real(C_DOUBLE), allocatable, intent(inout) :: arr(:,:,:,:,:)
        integer(C_INT) :: i, j, k, v, b, nx, ny, nz
        nx = blk%nb(1); ny = blk%nb(2); nz = blk%nb(3)
        if (.not. allocated(arr)) then
            allocate(arr(1:nx, 1:ny, 1:nz, NVEL, blk%nBlocks))
#ifdef USE_OPENMP_OFFLOAD
            !$omp target enter data map(alloc: arr)
#endif
        end if
#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do collapse(5) &
        !$omp& map(to: blk%oldrhs) map(tofrom: arr) private(i,j,k,v,b)
#endif
        do b = 1_C_INT, blk%nBlocks
            do v = 1, NVEL
                do k = 1, nz
                    do j = 1, ny
                        do i = 1, nx
                            arr(i,j,k,v,b) = blk%oldrhs(i,j,k,v,b)
                        end do
                    end do
                end do
            end do
        end do
#ifdef USE_OPENMP_OFFLOAD
        !$omp end target teams distribute parallel do
#endif
    end subroutine capture_rhs

    ! MOBY_RHSDUMP: park a 3-component interior buffer into the velocity slots of
    ! blk%q so write_field emits it as un/vn/wn (its target update from(blk%q)
    ! then pulls these values). blk%q is not reused after the final dump.
    subroutine overwrite_velocity(blk, arr)
        type(block_set_type), intent(inout) :: blk
        real(C_DOUBLE), allocatable, intent(in) :: arr(:,:,:,:,:)
        integer(C_INT) :: i, j, k, b, nx, ny, nz
        nx = blk%nb(1); ny = blk%nb(2); nz = blk%nb(3)
#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do collapse(4) &
        !$omp& map(to: arr) map(tofrom: blk%q) private(i,j,k,b)
#endif
        do b = 1_C_INT, blk%nBlocks
            do k = 1, nz
                do j = 1, ny
                    do i = 1, nx
                        blk%q(i,j,k,VAR_U,b) = arr(i,j,k,VAR_U,b)
                        blk%q(i,j,k,VAR_V,b) = arr(i,j,k,VAR_V,b)
                        blk%q(i,j,k,VAR_W,b) = arr(i,j,k,VAR_W,b)
                    end do
                end do
            end do
        end do
#ifdef USE_OPENMP_OFFLOAD
        !$omp end target teams distribute parallel do
#endif
    end subroutine overwrite_velocity

    ! MOBY_TERMDUMP: recompute the individual momentum terms of component `var`
    ! from the (exchanged) field, exactly as the predictor kernel does, but keep
    ! the three advection flux-divergence terms (x/y/z) and the three Laplacian
    ! terms (x/y/z, times 1/Re) SEPARATE -> adv(:,:,:,1:3), dif(:,:,:,1:3). A
    ! diagnostic; run on build_cpu (host). The leading device pull keeps a GPU
    ! build's host copy current for the read.
    subroutine compute_momentum_terms(blk, dns, var, adv, dif)
        type(block_set_type), intent(inout) :: blk
        type(dns_type), intent(in) :: dns
        integer, intent(in) :: var
        real(C_DOUBLE), intent(out) :: adv(:,:,:,:,:), dif(:,:,:,:,:)
        integer(C_INT) :: i, j, k, b, ip, im, jp, jm, kp, km, nx, ny, nz
        real(C_DOUBLE) :: ire, fxp, fxm, fyp, fym, fzp, fzm
        associate(q => blk%q)
        nx = blk%nb(1); ny = blk%nb(2); nz = blk%nb(3)
        ire = 1.0d0/dns%re
#ifdef USE_OPENMP_OFFLOAD
        !$omp target update from(blk%q)
#endif
        do b = 1, int(blk%nBlocks)
            do k = 1, nz
                do j = 1, ny
                    do i = 1, nx
                        ip=i+1; im=i-1; jp=j+1; jm=j-1; kp=k+1; km=k-1
                        select case (var)
                        case (VAR_U)
                            fxp = (q(i,j,k,VAR_U,b)+q(ip,j,k,VAR_U,b))**2
                            fxm = (q(im,j,k,VAR_U,b)+q(i,j,k,VAR_U,b))**2
                            fyp = (q(i,j,k,VAR_U,b)+q(i,jp,k,VAR_U,b))*(q(im,jp,k,VAR_V,b)+q(i,jp,k,VAR_V,b))
                            fym = (q(i,jm,k,VAR_U,b)+q(i,j,k,VAR_U,b))*(q(im,j,k,VAR_V,b)+q(i,j,k,VAR_V,b))
                            fzp = (q(i,j,k,VAR_U,b)+q(i,j,kp,VAR_U,b))*(q(im,j,kp,VAR_W,b)+q(i,j,kp,VAR_W,b))
                            fzm = (q(i,j,km,VAR_U,b)+q(i,j,k,VAR_U,b))*(q(im,j,k,VAR_W,b)+q(i,j,k,VAR_W,b))
                        case (VAR_V)
                            fxp = (q(i,j,k,VAR_V,b)+q(ip,j,k,VAR_V,b))*(q(ip,jm,k,VAR_U,b)+q(ip,j,k,VAR_U,b))
                            fxm = (q(im,j,k,VAR_V,b)+q(i,j,k,VAR_V,b))*(q(i,jm,k,VAR_U,b)+q(i,j,k,VAR_U,b))
                            fyp = (q(i,j,k,VAR_V,b)+q(i,jp,k,VAR_V,b))**2
                            fym = (q(i,jm,k,VAR_V,b)+q(i,j,k,VAR_V,b))**2
                            fzp = (q(i,j,k,VAR_V,b)+q(i,j,kp,VAR_V,b))*(q(i,jm,kp,VAR_W,b)+q(i,j,kp,VAR_W,b))
                            fzm = (q(i,j,km,VAR_V,b)+q(i,j,k,VAR_V,b))*(q(i,jm,k,VAR_W,b)+q(i,j,k,VAR_W,b))
                        case default
                            fxp = (q(i,j,k,VAR_W,b)+q(ip,j,k,VAR_W,b))*(q(ip,j,km,VAR_U,b)+q(ip,j,k,VAR_U,b))
                            fxm = (q(im,j,k,VAR_W,b)+q(i,j,k,VAR_W,b))*(q(i,j,km,VAR_U,b)+q(i,j,k,VAR_U,b))
                            fyp = (q(i,j,k,VAR_W,b)+q(i,jp,k,VAR_W,b))*(q(i,jp,km,VAR_V,b)+q(i,jp,k,VAR_V,b))
                            fym = (q(i,jm,k,VAR_W,b)+q(i,j,k,VAR_W,b))*(q(i,j,km,VAR_V,b)+q(i,j,k,VAR_V,b))
                            fzp = (q(i,j,k,VAR_W,b)+q(i,j,kp,VAR_W,b))**2
                            fzm = (q(i,j,km,VAR_W,b)+q(i,j,k,VAR_W,b))**2
                        end select
                        adv(i,j,k,1,b) = -0.25d0*(fxp-fxm)*blk%d1x(i,var,b)
                        adv(i,j,k,2,b) = -0.25d0*(fyp-fym)*blk%d1y(j,var,b)
                        adv(i,j,k,3,b) = -0.25d0*(fzp-fzm)*blk%d1z(k,var,b)
                        dif(i,j,k,1,b) = ire*(blk%lapXm(i,var,b)*q(im,j,k,var,b) &
                            + blk%lapX0(i,var,b)*q(i,j,k,var,b) + blk%lapXp(i,var,b)*q(ip,j,k,var,b))
                        dif(i,j,k,2,b) = ire*(blk%lapYm(j,var,b)*q(i,jm,k,var,b) &
                            + blk%lapY0(j,var,b)*q(i,j,k,var,b) + blk%lapYp(j,var,b)*q(i,jp,k,var,b))
                        dif(i,j,k,3,b) = ire*(blk%lapZm(k,var,b)*q(i,j,km,var,b) &
                            + blk%lapZ0(k,var,b)*q(i,j,k,var,b) + blk%lapZp(k,var,b)*q(i,j,kp,var,b))
                    end do
                end do
            end do
        end do
        end associate
    end subroutine compute_momentum_terms

    ! MOBY_MANUF: add amp*grad(phi) to the velocity, phi = cos(kx)cos(ky)cos(kz),
    ! at each component's staggered coordinate (full q bounds, halos included so
    ! the fine-owned interface face carries the analytic value). A pure gradient
    ! is irrotational: a consistent projection removes it entirely and recovers
    ! the exact div-free field; a periodic gradient is globally mass-conserving.
    subroutine manufacture_gradient_perturbation(blk, dns, amp)
        type(block_set_type), intent(inout) :: blk
        type(dns_type), intent(in) :: dns
        real(C_DOUBLE), intent(in) :: amp
        integer(C_INT) :: i, j, k, b
        real(C_DOUBLE) :: kk
        kk = 8.0d0*atan(1.0d0)/dns%leng(1)
        do b = 1, int(blk%nBlocks)
            do k = lbound(blk%q,3), ubound(blk%q,3)
                do j = lbound(blk%q,2), ubound(blk%q,2)
                    do i = lbound(blk%q,1), ubound(blk%q,1)
                        blk%q(i,j,k,1,b) = blk%q(i,j,k,1,b) - amp*kk &
                            *sin(kk*blk%x(i,1,b))*cos(kk*blk%y(j,1,b))*cos(kk*blk%z(k,1,b))
                        blk%q(i,j,k,2,b) = blk%q(i,j,k,2,b) - amp*kk &
                            *cos(kk*blk%x(i,2,b))*sin(kk*blk%y(j,2,b))*cos(kk*blk%z(k,2,b))
                        blk%q(i,j,k,3,b) = blk%q(i,j,k,3,b) - amp*kk &
                            *cos(kk*blk%x(i,3,b))*cos(kk*blk%y(j,3,b))*sin(kk*blk%z(k,3,b))
                    end do
                end do
            end do
        end do
#ifdef USE_OPENMP_OFFLOAD
        !$omp target update to(blk%q)
#endif
    end subroutine manufacture_gradient_perturbation

    function line_at(blk, d, level, idx) result(v)
        type(block_set_type), intent(in) :: blk
        integer(C_INT), intent(in) :: d, idx
        integer, intent(in) :: level
        real(C_DOUBLE) :: v

        select case (d)
        case (1)
            v = blk%lineX(idx, level + 1)
        case (2)
            v = blk%lineY(idx, level + 1)
        case default
            v = blk%lineZ(idx, level + 1)
        end select
    end function line_at
end program main
