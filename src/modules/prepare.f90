! The case BUILDER (docs/prepare_solve_strategy.md P0-P3; numerics review
! step 7-3 makes it the one builder with two entry points). prepare_case
! runs the solver's own init pipeline -- node lines, block size (the nb
! rule), geometry classification, leaf table, IBM coefficients, RANS wall
! distance -- and writes ONE case file, the single source of truth for the
! grid and the leaf table, which `moby_solve` only reads (preparing and
! solving are separate executables; the solver refuses a missing file and
! names the moby_prepare command). A body-free case gets attrs + node lines
! + leaf table and no coefficient dataset. Host code throughout: the
! coefficients come from ONE kernel over the indicator, so a file is the
! same from the CPU and the GPU build (step 7-4 retired the device twin and
! its libm ulps).
!
! MPI-parallel: the global leaf table is built identically on every rank
! (exactly as in the solver) and the per-leaf work -- coefficient tiles,
! wall-distance tiles, file rows -- is split over the world ranks by the
! same Z-order closed form the solver uses, so any rank count produces the
! same file (for an explicit nb; a DERIVED nb is one block per rank of the
! layout in force, see derive_block_nb).
module prepare
    use, intrinsic :: iso_c_binding
    use :: init, only: dns_type, grid_type, init_grid, destroy_grid, &
        VAR_U, VAR_V, VAR_W, VAR_P
    use :: blocks, only: block_set_type, init_block_set, destroy_block_set, derive_block_nb
    use :: flow_case, only: case_type, create_flow_case
    use :: config, only: config_seen_type, read_runtime_config, validate_dns_values, &
        case_input_echo, case_file_name, has_restart_file
    use :: scalar, only: scalar_type, destroy_scalar, scalars_enabled, &
        scalar_conjugate_enabled
    use :: boundary, only: boundary_type
    use :: io, only: write_case_file, read_restart_metadata
    use :: ibmm, only: ibm_type, init_ibm, set_ibm_geometry, set_ibm_coeff, &
        classify_refinement_masks, classify_active_mask, isInBody, body_indicator_i
    use :: geometry_stl, only: stl_geometry_load, stl_geometry_destroy, &
        stl_is_in_body, stl_fill_dwall, stl_cull_box
    use :: rans, only: fill_body_distance_analytic
    use :: pressure_solver, only: pressure_solver_type
    use :: turbulence, only: turb_type
    use :: les_model, only: les_type
    use :: comm, only: comm_type, comm_cart_dims
    implicit none

    private
    public :: prepare_case

contains

    ! Build the case named by input_file and write it to case_file -- or,
    ! when that is empty, to the name the solver will look for
    ! (case_file_name: [case] file, else <field_prefix>.case.h5). c is an
    ! initialised world communicator; its [mpi] dims (read from the ini
    ! here too) describe the SOLVE's rank layout for the nb rule.
    subroutine prepare_case(input_file, case_file, c)
        character(len=*), intent(in) :: input_file, case_file
        type(comm_type), intent(inout) :: c

        class(case_type), allocatable :: flow
        type(dns_type) :: dns
        type(grid_type) :: g
        type(block_set_type) :: blk
        type(boundary_type) :: bc
        type(pressure_solver_type) :: ps
        type(turb_type) :: turb
        type(les_type) :: les
        type(ibm_type) :: ibm
        type(config_seen_type) :: config_seen
        type(scalar_type) :: sc
        ! Unallocated optionals stay absent in the write_case_file call: only
        ! the pieces this case needs are computed and written.
        integer(C_INT), allocatable :: blockActive(:)
        integer(C_INT), allocatable :: blockTouch(:,:), blockBuried(:,:)
        integer(C_INT), allocatable :: blockMaskLo(:,:), blockMaskDims(:,:)
        real(C_DOUBLE), allocatable :: dwall(:,:,:,:)
        ! The one geometry switch: every downstream stage takes the indicator.
        procedure(body_indicator_i), pointer :: inside => null()
        logical :: use_stl, restart_exists
        ! Solid-possible box for classification culling (allocated for STL
        ! only; unallocated stays absent in the classify calls).
        real(C_DOUBLE), allocatable :: cullLo(:), cullHi(:)
        ! [scalar] declared => the coefficient array gains its VAR_P column and
        ! the case file gains coef_p_blocks.
        logical :: cell_centred
        integer(C_INT) :: nCoefComp
        character(len=:), allocatable :: inputs, out_file

        call create_flow_case(flow, input_file, c%has_terminal)
        call flow%apply_defaults(dns, g, bc, c, ps)
        ! A [scalar] section here means one thing: the case file must also carry
        ! the CELL-CENTRED coefficient tiles the scalars penalise with
        ! (coef_p_blocks, increment S3). Nothing else about the scalars matters
        ! to prepare -- Pr, the wall mode and the boundary rows are solve-time
        ! configuration, and one coefficient array serves every scalar.
        call read_runtime_config(dns, g, turb, les, ps, bc, sc, c, input_file, &
            c%has_terminal, config_seen)
        cell_centred = scalars_enabled(sc)
        out_file = trim(case_file)
        if (len(out_file) == 0) out_file = case_file_name(dns)
        if (c%has_terminal) print *, "preparing case: ", trim(input_file), " -> ", out_file
        ! A restart file supplies whatever the ini leaves unset (grid, re,
        ! periodicity, ibm_enabled), exactly as in the solver -- the file
        ! must describe the case the SOLVER will see. When it does not exist
        ! yet (a case whose IC is minted after the prepare, les_ibm/setup.sh)
        ! the ini alone defines the case.
        if (has_restart_file(dns)) then
            inquire(file=trim(dns%restart_file), exist=restart_exists)
            if (restart_exists) call read_restart_metadata(dns, g, bc, ps%nIter, ps%omega, &
                dns%restart_file, c, config_seen)
        end if

        call init_grid(g, dns, bc%isPeriodic)
        call validate_dns_values(dns, g)
        ! The nb rule (step 7-1): an unset [blocks] nb becomes one block per
        ! rank of the Cartesian layout the solver forms -- [mpi] dims, zeros
        ! filled over the SOLVE's rank count (the caller's world size in the
        ! solver; standalone moby_prepare's own, whose explicit dims then
        ! describe the solve) -- and is stored in the file. Refinement and
        ! buried-block removal keep needing an explicit nb (rule 4): with a
        ! derived one they are not silently ignored (the old nb-less solver
        ! path did that) but refused.
        call comm_cart_dims(c, c%world_size, for_solve=.false.)
        call derive_block_nb(dns, c%dims, product(c%dims), c%has_terminal)
        if (dns%block_nb_auto .and. (dns%block_refine_body .or. dns%block_refine_nboxes > 0_C_INT)) &
            error stop "[blocks] refine / refine_body need an explicit [blocks] nb"
        ! What this file is a function of, echoed into it for the staleness
        ! check on read (case_input_echo).
        inputs = case_input_echo(dns, bc, sc)

        ! Geometry source: the analytic isInBody, or an STL body loaded behind
        ! the same indicator signature ([ibm] stl_file, P1).
        use_stl = dns%ibm_stl_count > 0_C_INT
        if (use_stl) then
            call stl_geometry_load(dns%ibm_stl_file(1:dns%ibm_stl_count), &
                dns%ibm_stl_scale, dns%ibm_stl_translate, dns%leng, &
                logical(bc%isPeriodic), c%has_terminal)
            inside => stl_is_in_body
            allocate(cullLo(3), cullHi(3))
            call stl_cull_box(cullLo, cullHi)
        else
            inside => isInBody
        end if

        ! Geometry classification + block set. classify_* keep their masks
        ! here so they can go into the case file after the block set consumed
        ! them. A body-free case (step 7-1) is the plain lattice; so is a body
        ! with a DERIVED nb, where buried-block removal is not applied (rule
        ! 4: the old nb-less solver path removed nothing, and a file must
        ! reproduce it).
        if (dns%ibm_enabled .and. c%has_terminal) print *, "classifying geometry..."
        if (dns%block_refine_body) then
            call classify_refinement_masks(blockTouch, blockBuried, blockMaskLo, &
                blockMaskDims, dns, g, ibm, bc%isPeriodic, c%has_terminal, inside, &
                cullLo, cullHi, c)
            call init_block_set(blk, dns, g, bc%isPeriodic, int(c%world_size, C_INT), &
                int(c%world_rank, C_INT), touch=blockTouch, buried=blockBuried, &
                maskLo=blockMaskLo, maskDims=blockMaskDims)
        else if (dns%ibm_enabled .and. dns%block_remove_solid .and. .not. dns%block_nb_auto) then
            call classify_active_mask(blockActive, dns, g, ibm, bc%isPeriodic, &
                c%has_terminal, inside, cullLo, cullHi, c)
            call init_block_set(blk, dns, g, bc%isPeriodic, int(c%world_size, C_INT), &
                int(c%world_rank, C_INT), blockActive)
        else
            call init_block_set(blk, dns, g, bc%isPeriodic, int(c%world_size, C_INT), &
                int(c%world_rank, C_INT))
        end if

        ! IBM coefficients on this rank's leaves: one host kernel over the
        ! indicator, analytic or STL.
        nCoefComp = 0_C_INT
        if (dns%ibm_enabled) then
            if (c%has_terminal) print *, "computing IBM coefficients..."
            call init_ibm(ibm, blk, cell_centred)
            ! Analytic wall geometry from the config. Both binaries apply it --
            ! see set_ibm_geometry for why that is not optional.
            call set_ibm_geometry(ibm, dns)
            nCoefComp = int(ubound(ibm%coef,4) - lbound(ibm%coef,4) + 1, C_INT)
            call set_ibm_coeff(dns, blk, ibm, VAR_U, inside)
            call set_ibm_coeff(dns, blk, ibm, VAR_V, inside)
            call set_ibm_coeff(dns, blk, ibm, VAR_W, inside)
            if (cell_centred) call set_ibm_coeff(dns, blk, ibm, VAR_P, inside)
        end if

        ! RANS wall distance: the raw body distance (the domain-wall min and
        ! half-cell floor stay solve-time, applied after the file read).
        ! Analytic = the indicator-driven walldist machinery, exactly what
        ! init_rans_geometry computed inline; STL = the exact BVH
        ! point-triangle distance (the same query mobygeom's dwall_blocks
        ! uses -- indicator bisections would cost millions of parity casts).
        ! Conjugate heat transfer (increment C1) reads the very same tiles: the
        ! signed distance whose sign the IBM marker supplies IS this field, so
        ! `ibm_wall = conjugate` triggers the build exactly as a [rans] section
        ! does. A TRIGGER only -- no new dataset, and a file prepared either way
        ! is identical (docs/next_session_conjugate.md Section 8, arrangement 1).
        if (dns%ibm_enabled .and. (dns%rans_configured .or. scalar_conjugate_enabled(sc))) then
            if (c%has_terminal) print *, "computing wall distance..."
            allocate(dwall(0:int(blk%nb(1))+1, 0:int(blk%nb(2))+1, &
                0:int(blk%nb(3))+1, blk%nBlocks))
            if (use_stl) then
                call stl_fill_dwall(dwall, blk)
            else
                call fill_body_distance_analytic(dwall, dns, blk, bc, ibm, &
                    c%has_terminal, inside)
            end if
        end if

        if (c%has_terminal) print *, "writing case file: ", out_file
        if (dns%ibm_enabled) then
            call write_case_file(out_file, blk, dns, g, bc, c, nCoefComp, inputs, c%has_terminal, &
                coef=ibm%coef, touch=blockTouch, buried=blockBuried, maskDims=blockMaskDims, &
                active=blockActive, dwall=dwall, maskLo=blockMaskLo)
        else
            call write_case_file(out_file, blk, dns, g, bc, c, nCoefComp, inputs, c%has_terminal)
        end if
        if (c%has_terminal) then
            print *, "case file written:", blk%nBlocksGlobal, "leaves,", &
                int(blk%nLevels) - 1, "refinement level(s)"
            print *, "solve with: mpirun -n ", product(c%dims), " moby_solve ", trim(input_file)
        end if

        if (use_stl) call stl_geometry_destroy()
        call destroy_block_set(blk)
        call destroy_grid(g)
        call destroy_scalar(sc)
    end subroutine prepare_case

end module prepare
