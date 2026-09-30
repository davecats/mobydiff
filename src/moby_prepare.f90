! moby_prepare: the case builder's executable (src/modules/prepare.f90,
! docs/prepare_solve_strategy.md). `mpirun -n N moby_prepare input.ini`
! writes the case file moby_solve will look for ([case] file, else
! <field_prefix>.case.h5) with N ranks' worth of blocks when [blocks] nb is
! unset; the solver never builds one. A second argument names the output
! explicitly (the gate drivers' twins).
program moby_prepare
    use, intrinsic :: iso_c_binding
    use :: prepare, only: prepare_case
    use :: comm, only: comm_type, comm_init_world, comm_finalize
    implicit none

    character(len=256) :: input_file, output_file
    logical :: show_help
    type(comm_type) :: c

    call comm_init_world(c)
    call parse_prepare_args(input_file, output_file, show_help)

    if (show_help) then
        if (c%has_terminal) call print_usage()
        call comm_finalize(c)
        stop
    end if

    call prepare_case(input_file, trim(output_file), c)
    call comm_finalize(c)

contains

    subroutine parse_prepare_args(input_file, output_file, show_help)
        character(len=*), intent(out) :: input_file
        character(len=*), intent(out) :: output_file
        logical, intent(out) :: show_help

        character(len=256) :: arg
        integer :: argc, i, positional

        input_file = "input.ini"
        output_file = ""
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
        print '(A)', "       (no case.h5: the ini's [case] file, else <field_prefix>.case.h5)"
    end subroutine print_usage

end program moby_prepare
