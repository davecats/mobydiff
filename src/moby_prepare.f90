! moby_prepare: the offline preparation step (docs/prepare_solve_strategy.md,
! phase P0). Runs the solver's own init pipeline for an ANALYTIC-geometry
! case -- node lines, geometry classification, leaf table, IBM coefficients,
! RANS wall distance -- and writes ONE case file in the block-table
! coefficient-file format the solver already reads. The solver then runs the
! identical case from the file ([ibm] coeff_file = <case.h5>) on any rank
! count, bit-exact vs its inline analytic path: every quantity in the file
! comes from the very kernels the solver would run at init.
!
! MPI-parallel: the global leaf table is built identically on every rank
! (exactly as in the solver) and the per-leaf work -- coefficient tiles,
! wall-distance tiles, file rows -- is split over the ranks by the same
! Z-order closed form the solver uses, so any rank count produces the same
! file.
program moby_prepare
    use, intrinsic :: iso_c_binding
    use :: init, only: dns_type, grid_type, init_grid, destroy_grid, &
        set_serial_local_size, VAR_U, VAR_V, VAR_W
    use :: blocks, only: block_set_type, init_block_set, destroy_block_set
    use :: flow_case, only: case_type, create_flow_case
    use :: config, only: config_seen_type, read_runtime_config, validate_dns_values
    use :: boundary, only: boundary_type
    use :: io, only: write_case_file
    use :: ibmm, only: ibm_type, init_ibm, enter_ibm_data, exit_ibm_data, &
        set_ibm_coeff, set_ibm_coeff_host, classify_refinement_masks, &
        classify_active_mask, isInBody, body_indicator_i
    use :: geometry_stl, only: stl_geometry_load, stl_geometry_destroy, &
        stl_is_in_body, stl_fill_dwall, stl_cull_box
    use :: rans, only: fill_body_distance_analytic
    use :: pressure_solver, only: pressure_solver_type
    use :: turbulence, only: turb_type
    use :: les_model, only: les_type
    use :: comm, only: comm_type, comm_init_world, comm_finalize
    implicit none

    character(len=256) :: input_file, output_file
    logical :: show_help
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
    type(comm_type) :: c
    ! Unallocated optionals stay absent in the write_case_file call: only
    ! the pieces this case needs are computed and written.
    integer(C_INT), allocatable :: blockActive(:)
    integer(C_INT), allocatable :: blockTouch(:,:), blockBuried(:,:)
    integer(C_INT), allocatable :: blockMaskLo(:,:), blockMaskDims(:,:)
    real(C_DOUBLE), allocatable :: dwall(:,:,:,:)
    ! The one geometry switch: every downstream stage takes the indicator.
    procedure(body_indicator_i), pointer :: inside => null()
    logical :: use_stl
    ! Solid-possible box for classification culling (allocated for STL
    ! only; unallocated stays absent in the classify calls).
    real(C_DOUBLE), allocatable :: cullLo(:), cullHi(:)

    call comm_init_world(c)
    call parse_prepare_args(input_file, output_file, show_help)

    if (show_help) then
        if (c%has_terminal) call print_usage()
        call comm_finalize(c)
        stop
    end if

    if (c%has_terminal) print *, "reading input data: ", trim(input_file)
    call create_flow_case(flow, input_file, c%has_terminal)
    call flow%apply_defaults(dns, g, bc, c, ps)
    call read_runtime_config(dns, g, turb, les, ps, bc, c, input_file, &
        c%has_terminal, config_seen)

    if (dns%block_nb <= 0_C_INT) &
        error stop "moby_prepare needs [blocks] nb: the case file is a block-table file"
    if (.not. dns%ibm_enabled) &
        error stop "moby_prepare needs [ibm] enabled = true (analytic or stl_file geometry)"
    if (len_trim(dns%ibm_coeff_file) > 0) &
        error stop "moby_prepare computes the coefficient file; drop [ibm] coeff_file from its input"

    call set_serial_local_size(dns)
    call init_grid(g, dns, bc%isPeriodic)
    call validate_dns_values(dns, g)

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

    ! Geometry classification + block set: the solver's init dispatch
    ! (main.f90), analytic branch. classify_* keep their masks here so they
    ! can go into the case file after the block set consumed them.
    if (c%has_terminal) print *, "classifying geometry..."
    if (dns%block_refine_body) then
        call classify_refinement_masks(blockTouch, blockBuried, blockMaskLo, &
            blockMaskDims, dns, g, ibm, bc%isPeriodic, c%has_terminal, inside, &
            cullLo, cullHi, c)
        call init_block_set(blk, dns, g, bc%isPeriodic, int(c%world_size, C_INT), &
            int(c%world_rank, C_INT), touch=blockTouch, buried=blockBuried, &
            maskLo=blockMaskLo, maskDims=blockMaskDims)
    else if (dns%block_remove_solid) then
        call classify_active_mask(blockActive, dns, g, ibm, bc%isPeriodic, &
            c%has_terminal, inside, cullLo, cullHi, c)
        call init_block_set(blk, dns, g, bc%isPeriodic, int(c%world_size, C_INT), &
            int(c%world_rank, C_INT), blockActive)
    else
        call init_block_set(blk, dns, g, bc%isPeriodic, int(c%world_size, C_INT), &
            int(c%world_rank, C_INT))
    end if

    ! IBM coefficients on this rank's leaves. Analytic keeps the solver's
    ! inline kernel verbatim (the P0 bit-exactness gate; on the device in
    ! offload builds -- prepare with the CPU build for the gates); STL runs
    ! the host twin over the indicator.
    if (c%has_terminal) print *, "computing IBM coefficients..."
    call init_ibm(ibm, blk)
    if (use_stl) then
        call set_ibm_coeff_host(dns, blk, ibm, VAR_U, inside)
        call set_ibm_coeff_host(dns, blk, ibm, VAR_V, inside)
        call set_ibm_coeff_host(dns, blk, ibm, VAR_W, inside)
    else
        call enter_ibm_data(ibm, dns)
        call set_ibm_coeff(dns, blk, ibm, VAR_U)
        call set_ibm_coeff(dns, blk, ibm, VAR_V)
        call set_ibm_coeff(dns, blk, ibm, VAR_W)
#ifdef USE_OPENMP_OFFLOAD
        !$omp target update from(ibm%coef)
#endif
    end if

    ! RANS wall distance: the raw body distance (the domain-wall min and
    ! half-cell floor stay solve-time, applied after the file read).
    ! Analytic = the indicator-driven walldist machinery, exactly what
    ! init_rans_geometry computes inline; STL = the exact BVH
    ! point-triangle distance (the same query mobygeom's dwall_blocks
    ! uses -- indicator bisections would cost millions of parity casts).
    if (dns%rans_configured) then
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

    if (c%has_terminal) print *, "writing case file: ", trim(output_file)
    call write_case_file(output_file, blk, dns, g, bc, c, ibm%coef, c%has_terminal, &
        touch=blockTouch, buried=blockBuried, maskDims=blockMaskDims, &
        active=blockActive, dwall=dwall, maskLo=blockMaskLo)
    if (c%has_terminal) then
        print *, "case file written:", blk%nBlocksGlobal, "leaves,", &
            int(blk%nLevels) - 1, "refinement level(s)"
        print *, "run the solver with [ibm] coeff_file = ", trim(output_file)
    end if

    if (use_stl) then
        call stl_geometry_destroy()
    else
        call exit_ibm_data(ibm, dns)
    end if
    call destroy_block_set(blk)
    call destroy_grid(g)
    call comm_finalize(c)

contains

    subroutine parse_prepare_args(input_file, output_file, show_help)
        character(len=*), intent(out) :: input_file
        character(len=*), intent(out) :: output_file
        logical, intent(out) :: show_help

        character(len=256) :: arg
        integer :: argc, i, positional

        input_file = "input.ini"
        output_file = "moby_case.h5"
        show_help = .false.
        positional = 0
        argc = command_argument_count()
        i = 1
        do while (i <= argc)
            call get_command_argument(i, arg)
            select case (trim(arg))
            case ("--help", "-h")
                show_help = .true.
                return
            case ("--input", "-i")
                i = i + 1
                if (i > argc) error stop "missing value after --input"
                call get_command_argument(i, input_file)
            case ("--output", "-o")
                i = i + 1
                if (i > argc) error stop "missing value after --output"
                call get_command_argument(i, output_file)
            case default
                positional = positional + 1
                select case (positional)
                case (1)
                    input_file = arg
                case (2)
                    output_file = arg
                case default
                    error stop "usage: moby_prepare [input.ini] [case.h5]"
                end select
            end select
            i = i + 1
        end do
    end subroutine parse_prepare_args

    subroutine print_usage()
        print '(A)', "usage: moby_prepare [input.ini] [case.h5]"
        print '(A)', "       moby_prepare --input input.ini --output case.h5"
    end subroutine print_usage

end program moby_prepare
