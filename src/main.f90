! TODO:
!
!      (*) second-order IBM in space, then in time


program main
    use :: init
    use :: blocks, only: block_set_type, init_block_set, destroy_block_set, &
        enter_block_data, exit_block_data, zero_closed_halos
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
    use :: comm, only: comm_type, comm_init_world, comm_init, comm_finalize, &
        init_block_exchange, exchange_halos, exchange_scalar_halos
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
        ! Analytic geometry only; without an IBM coef stays 0 (set in init_ibm).
        if (dns%ibm_enabled) then
            call set_ibm_coeff(dns, blk, ibm, VAR_U)
            call set_ibm_coeff(dns, blk, ibm, VAR_V)
            call set_ibm_coeff(dns, blk, ibm, VAR_W)
        end if
    end if

    if (les_is_enabled(les)) then
        if (c%has_terminal) print *, "initialising LES model..."
        call init_les(les, dns, blk)
        call enter_les_data(les, dns)
    end if

    call apply_bc(blk, bc)
    call exchange_halos(c, blk, [VAR_U, VAR_V, VAR_W, VAR_P])

    ! Halo coherence audit (MOBY_HALO_AUDIT): check every exchange-written halo
    ! cell against manufactured linear fields on the real block layout. Run this
    ! first when a 2:1 interface case misbehaves.
    block
        character(len=16) :: auditEnv
        call get_environment_variable("MOBY_HALO_AUDIT", auditEnv)
        if (len_trim(auditEnv) > 0) then
            call halo_audit(blk, dns, bc, c)
            stop
        end if
    end block

    ! Truncation-error probe (MOBY_TRUNC): with the exact TGV IC, set the
    ! interface velocity halos to the scheme's two-phase state (Pass A injection
    ! at line 172 above + Pass B linear here) and print the u-momentum operator
    ! term balance along a fine-owns interface (see step.f90:truncation_probe).
    block
        character(len=16) :: truncEnv
        call get_environment_variable("MOBY_TRUNC", truncEnv)
        if (len_trim(truncEnv) > 0) then
            call exchange_halos(c, blk, [VAR_U, VAR_V, VAR_W], linear_prolong=.true.)
            call truncation_probe(blk, dns)
            stop
        end if
    end block

    call flow%setup_after_grid(blk, dns, g, bc, c)
    if (les_is_enabled(les)) then
        call update_les_viscosity(les, blk, dns, ibm)
        call exchange_scalar_halos(c, les%nut)
        call update_timestep_limits(blk, dns, c, les)
    else
        call update_timestep_limits(blk, dns, c)
    end if

    ! Projection-consistency probe (MOBY_PROJPROBE): feed the projection the exact
    ! (divergence-free) TGV field as its predictor and print how much it spuriously
    ! changes it at the interface (step.f90:proj_consistency_probe).
    block
        character(len=16) :: projEnv
        call get_environment_variable("MOBY_PROJPROBE", projEnv)
        if (len_trim(projEnv) > 0) then
            call copy_q_to_qs(blk)
            call pressure_projection(ps, blk, dns%dt*rk_gamma(1), ibm, bc, c)
            call proj_consistency_probe(blk, dns)
            stop
        end if
    end block

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

            ! Predictor: advance tentative staggered velocities, then enforce
            ! solid/body constraints and exchange halos (the interface high faces
            ! are filled with the owner velocity here).
            ! Without an immersed boundary the penalization coefficient is zero
            ! everywhere, so mu = 1/(1+dt_gamma*coef) stays at its init value of 1
            ! and need not be recomputed each substage.
            if (dns%ibm_enabled) call update_ibm_mu(ibm, dt_gamma)
            if (les_is_enabled(les)) then
                les_profile_start = les_wall_seconds()
                call update_les_viscosity(les, blk, dns, ibm)
                call add_les_profile(les_prof, LES_PROF_NUT, les_wall_seconds() - les_profile_start)
                les_profile_start = les_wall_seconds()
                call exchange_scalar_halos(c, les%nut)
                call add_les_profile(les_prof, LES_PROF_EXCHANGE, les_wall_seconds() - les_profile_start)
                call momentum(blk, dns, dt_alpha, dt_beta, dt_gamma, ibm, les, les_prof)
            else
                call momentum(blk, dns, dt_alpha, dt_beta, dt_gamma, ibm)
            end if
            call apply_bc(blk, bc)
            call exchange_halos(c, blk, [VAR_U, VAR_V, VAR_W])

            ! Projection: solve for the pressure correction and project the
            ! tentative velocities to a divergence-free field.
            call pressure_projection(ps, blk, dt_gamma, ibm, bc, c)

        end do

        if (les_is_enabled(les)) then
            call update_timestep_limits(blk, dns, c, les)
        else
            call update_timestep_limits(blk, dns, c)
        end if

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
    end if

    call write_field(blk, dns, g, int(dns%step_current), c, bc, ps%nIter, ps%sor)

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
        use :: blocks, only: leaf_at, level_cells, FACE_COARSE
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
        integer(C_INT) :: idSame, idParent, sx, sy, sz, ghalo, gint
        logical :: skip, anyChild
        real(C_DOUBLE) :: got, want, posD(3), srcD(3), wBlend, aHalf, bHalf, cHalf

        nb = blk%nb

        do b = 1, int(blk%nBlocks)
            do var = 1, 4
                ! High side runs to nb+2 (the redundant top-face momentum
                ! halo); the low side keeps its single ghost at 0.
                do k = 0, nb(3) + 2
                    do j = 0, nb(2) + 2
                        do i = 0, nb(1) + 2
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
                ! High side runs to nb+2 (the redundant top-face momentum
                ! halo); the low side keeps its single ghost at 0.
                do k = 0, nb(3) + 2
                    do j = 0, nb(2) + 2
                        do i = 0, nb(1) + 2
                            idx = [i, j, k]
                            do d = 1, 3
                                off(d) = 0
                                if (idx(d) == 0) off(d) = -1
                                if (idx(d) >= nb(d) + 1) off(d) = 1
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

                            ! The 2nd high-side halo layer (idx == nb+2) is
                            ! audited for same-level copies; the interface
                            ! PROLONG normal stencil there is not.
                            if (any(idx > nb + 1) .and. idSame < 0) cycle

                            ! Momentum-owned interface faces: a fine block
                            ! computes its top normal face v(nb+1) (var == d)
                            ! where the +d neighbour is coarser, so it is NOT
                            ! exchange-written and keeps the sentinel (doc 6a).
                            skip = .false.
                            do d = 1, 3
                                if (off(d) == 1 .and. var == d .and. &
                                    blk%physHigh(d,b) == FACE_COARSE) skip = .true.
                            end do
                            if (skip) cycle

                            posD = [blk%x(i,var,b), blk%y(j,var,b), blk%z(k,var,b)]
                            if (idSame >= 0) then
                                want = CXV(var)*posD(1) + CYV(var)*posD(2) + CZV(var)*posD(3)
                            else if (idParent >= 0) then
                                ! Injection: field at the covering coarse location.
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
                                if (var == VAR_P .and. sum(abs(off)) == 1) then
                                    ! Pressure PROLONG ghost = w*p_coarse +
                                    ! (1-w)*p_fine with the node-line geometric
                                    ! weight w = (bHalf+cHalf)/(aHalf+cHalf) of
                                    ! the fine-halo / coarse-cover / fine-interior
                                    ! half-widths (mirrors comm entry_blend). A
                                    ! fixed 2/3 is only exact on a uniform grid.
                                    do d = 1, 3
                                        if (off(d) == 0) cycle
                                        ghalo = int(blk%origin(d,b)) + idx(d) - 1
                                        gint = ghalo - off(d)
                                        cov = ghalo/2
                                        bHalf = 0.5d0*(line_at(blk,d,l,ghalo+1) - line_at(blk,d,l,ghalo))
                                        cHalf = 0.5d0*(line_at(blk,d,l,gint+1) - line_at(blk,d,l,gint))
                                        aHalf = 0.5d0*(line_at(blk,d,l-1,cov+1) - line_at(blk,d,l-1,cov))
                                        wBlend = (bHalf + cHalf)/(aHalf + cHalf)
                                    end do
                                    want = wBlend*want + (1.0d0 - wBlend) &
                                        *(CXV(var)*blk%x(i-off(1),var,b) &
                                        + CYV(var)*blk%y(j-off(2),var,b) &
                                        + CZV(var)*blk%z(k-off(3),var,b))
                                end if
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
                            if (abs(got - want) > 1.0d-10*max(1.0d0, abs(want))) then
                                nBad = nBad + 1
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
    end subroutine halo_audit

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
