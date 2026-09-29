! moby_prepare: the standalone entry point of the case builder
! (src/modules/prepare.f90, docs/prepare_solve_strategy.md). The solver
! runs the identical case from the file ([case] file = <case.h5>) on any
! rank count, bit-exact vs its inline path, and prepares in-process itself
! when the file is absent -- this executable exists for preparing once,
! elsewhere, or on a different rank count than the solve.
program moby_prepare
    use, intrinsic :: iso_c_binding
    use :: prepare, only: prepare_case
    use :: comm, only: comm_type, comm_init_world, comm_finalize
    implicit none

    character(len=256) :: input_file, output_file, legacy_file
    logical :: show_help
    type(comm_type) :: c

    call comm_init_world(c)
    call parse_prepare_args(input_file, output_file, legacy_file, show_help)

    if (show_help) then
        if (c%has_terminal) call print_usage()
        call comm_finalize(c)
        stop
    end if

    if (len_trim(legacy_file) > 0) then
        call prepare_case(input_file, output_file, c, legacy=trim(legacy_file))
    else
        call prepare_case(input_file, output_file, c)
    end if
    if (c%has_terminal) print *, "run the solver with [case] file = ", trim(output_file)
    call comm_finalize(c)

contains

    subroutine parse_prepare_args(input_file, output_file, legacy_file, show_help)
        character(len=*), intent(out) :: input_file
        character(len=*), intent(out) :: output_file
        character(len=*), intent(out) :: legacy_file
        logical, intent(out) :: show_help

        character(len=256) :: arg
        integer :: argc, i, positional

        input_file = "input.ini"
        output_file = "moby_case.h5"
        legacy_file = ""
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
            case ("--convert-legacy")
                ! One-off (step 7-3): copy a retired-mobygeom global-layout
                ! coefficient file into the block-table case-file layout,
                ! number for number.
                i = i + 1
                if (i > argc) error stop "missing value after --convert-legacy"
                call get_command_argument(i, legacy_file)
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
        print '(A)', "       moby_prepare --convert-legacy old_coeff.h5 input.ini case.h5"
    end subroutine print_usage

end program moby_prepare
