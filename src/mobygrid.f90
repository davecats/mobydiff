program mobygrid
    use, intrinsic :: iso_c_binding, only: C_INT
    use :: init, only: dns_type, grid_type, init_grid, destroy_grid
    use :: blocks, only: block_set_type, init_block_set, destroy_block_set
    use :: flow_case, only: case_type, create_flow_case
    use :: config, only: config_seen_type, read_runtime_config, has_restart_file, validate_dns_values
    use :: boundary, only: boundary_type
    use :: io, only: write_grid_export, read_restart_metadata
    use :: pressure_solver, only: pressure_solver_type
    use :: les_model, only: les_type
    use :: comm, only: comm_type, comm_init_world, comm_finalize
    implicit none

    character(len=256) :: input_file
    character(len=256) :: output_file
    logical :: show_help
    class(case_type), allocatable :: flow
    type(dns_type) :: dns
    type(grid_type) :: g
    type(block_set_type) :: blk
    type(boundary_type) :: bc
    type(pressure_solver_type) :: ps
    type(les_type) :: les
    type(config_seen_type) :: config_seen
    type(comm_type) :: c

    call comm_init_world(c)
    call parse_mobygrid_args(input_file, output_file, show_help)

    if (show_help) then
        if (c%has_terminal) call print_usage()
        call comm_finalize(c)
        stop
    end if

    if (c%world_size /= 1) then
        if (c%has_terminal) print *, "error: mobygrid writes one global grid file; run it with one MPI rank"
        call comm_finalize(c)
        error stop "mobygrid requires one MPI rank"
    end if

    if (c%has_terminal) print *, "reading input data: ", trim(input_file)
    call create_flow_case(flow, input_file, c%has_terminal)
    call flow%apply_defaults(dns, g, bc, c, ps)
    call read_runtime_config(dns, g, les, ps, bc, c, input_file, c%has_terminal, config_seen)
    if (has_restart_file(dns)) then
        if (c%has_terminal) print *, "reading restart metadata: ", trim(dns%restart_file)
        call read_restart_metadata(dns, g, bc, ps%nIter, ps%sor, dns%restart_file, c, &
            preserve_cflmax=config_seen%cflmax, preserve_pecletmax=config_seen%pecletmax, &
            preserve_dtmax=config_seen%dtmax, preserve_t_final=config_seen%t_final)
    end if

    call set_serial_local_size(dns)
    call init_grid(g, dns, bc%isPeriodic)
    call validate_dns_values(dns, g)

    ! Serial: one block spanning the whole grid provides the staggered
    ! coordinates for the export, regardless of any [blocks] nb setting.
    dns%block_nb = 0_C_INT
    call init_block_set(blk, dns, g, bc%isPeriodic)

    if (c%has_terminal) print *, "writing grid file: ", trim(output_file)
    call write_grid_export(dns, g, blk, bc, output_file, c%has_terminal)

    call destroy_block_set(blk)
    call destroy_grid(g)
    call comm_finalize(c)

contains

    subroutine parse_mobygrid_args(input_file, output_file, show_help)
        character(len=*), intent(out) :: input_file
        character(len=*), intent(out) :: output_file
        logical, intent(out) :: show_help

        character(len=256) :: arg
        integer :: argc, i, positional

        input_file = "input.ini"
        output_file = "mobygrid.h5"
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
                    error stop "usage: mobygrid [input.ini] [grid.h5]"
                end select
            end select
            i = i + 1
        end do
    end subroutine parse_mobygrid_args

    subroutine print_usage()
        print '(A)', "usage: mobygrid [input.ini] [grid.h5]"
        print '(A)', "       mobygrid --input input.ini --output grid.h5"
    end subroutine print_usage

    subroutine set_serial_local_size(dns)
        type(dns_type), intent(inout) :: dns
        integer :: dir

        do dir = 1, 3
            dns%localSize(dir,0) = 1_C_INT
            dns%localSize(dir,1) = dns%globalSize(dir)
            dns%localSize(dir,2) = dns%globalSize(dir)
        end do
    end subroutine set_serial_local_size

end program mobygrid
