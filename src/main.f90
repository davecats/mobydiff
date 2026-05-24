! TODO:
!
!      (*) second-order IBM in space, then in time


program main
    use, intrinsic :: iso_fortran_env, only: int64
    use :: init
    use :: config
    use :: boundary
    use :: io
    use :: step
    use :: pressure_solver
    use :: gpu_runtime
    use :: ibmm
    use :: comm, only: comm_type, comm_init_world, comm_init, comm_finalize, comm_allreduce_max, &
                       exchange_halos
    implicit none

    integer :: i, arg_status, rkStage
    logical :: need_cfl
    real(C_DOUBLE) :: dt_alpha, dt_beta, dt_gamma
    real(C_DOUBLE) :: cfl_reduce(1)
    real(C_DOUBLE) :: loop_seconds, seconds_per_step
    integer(int64) :: loop_clock_start, loop_clock_end, clock_rate
    character(len=256) :: input_file
    type(dns_type) :: dns
    type(grid_type) :: g
    type(field_type) :: f
    type(boundary_type) :: bc
    type(pressure_solver_type) :: ps
    type(ibm_type) :: ibm
    type(comm_type) :: c

    call comm_init_world(c)
    call splash(c%has_terminal)

    ! The first command-line argument can override the default input file.
    call get_command_argument(1, input_file, status=arg_status)
    if (arg_status /= 0 .or. len_trim(input_file) == 0) input_file = "input.ini"

    if (c%has_terminal) print *, "reading input data..."
    call init_bc(bc)
    call read_runtime_config(dns, ps, bc, input_file, c%has_terminal)
    if (has_restart_file(dns)) then
        if (c%has_terminal) print *, "reading restart metadata: ", trim(dns%restart_file)
        call read_restart_metadata(dns, bc, ps%nIter, ps%sor, dns%restart_file, c)
    end if
    call comm_init(c, dns, bc, int(dns%mpiDims))

    if (c%has_terminal) print *, "initialising grid..."
    call init_grid(g, dns, bc%isPeriodic)
    call validate_dns_values(dns)
    call init_boundary_faces(bc, dns)
    call init_openmp_offload(c%has_terminal)
    call enter_grid_data(g, dns)

    if (c%has_terminal) print *, "initialising fields..."
    call init_field(f, dns)
    call enter_field_data(f, dns)
    if (has_restart_file(dns)) then
        if (c%has_terminal) print *, "reading restart fields: ", trim(dns%restart_file)
        call read_field(f, dns, dns%restart_file, c)
    end if

    if (c%has_terminal) print *, "initialising pressure solver..."
    call init_pressure_solver(ps, dns, bc, c%has_terminal)

    if (c%has_terminal) print *, "initialising IBM..."
    call init_ibm(ibm, dns, g)
    call enter_ibm_data(ibm, dns)
    call set_ibm_coeff(dns, g, ibm, VAR_U)
    call set_ibm_coeff(dns, g, ibm, VAR_V)
    call set_ibm_coeff(dns, g, ibm, VAR_W)

    call apply_bc(f, dns, g, bc)
    call exchange_halos(c, f, [VAR_U, VAR_V, VAR_W, VAR_P])

    if (c%has_terminal) print *, "main loop starting..."
    call system_clock(count_rate=clock_rate)
    call system_clock(count=loop_clock_start)
    do i = 1, int(dns%nsteps)
        dns%step_current = dns%step_current + 1_C_INT
        dns%t_current = dns%t_current + dns%dt

        do rkStage = 1,3
            dt_alpha = dns%dt*rk_alpha(rkStage)
            dt_beta  = dns%dt*rk_beta(rkStage)
            dt_gamma = dns%dt*rk_gamma(rkStage)

            ! Predictor: advance tentative staggered velocities, then enforce solid/body constraints.
            call momentum(f, dns, g, dt_alpha, dt_beta, dt_gamma, ibm, bc)
            call apply_bc(f, dns, g, bc)
            call exchange_halos(c, f, [VAR_U, VAR_V, VAR_W])

            ! Projection: solve for pressure correction and project tentative velocities.
            call pressure_projection(ps, f, dns, g, dt_gamma, ibm, bc, c)

        end do

        ! Compute CFL only when it drives adaptive time stepping.
        need_cfl = (dns%cflmax > 0.0d0)
        if (need_cfl) then
            dns%cfl = get_cfl(f, dns, g)
            cfl_reduce(1) = dns%cfl
            call comm_allreduce_max(c, cfl_reduce)
            dns%cfl = cfl_reduce(1)
        end if

        if (dns%cflmax > 0.0d0 .and. dns%cfl > 0.0d0) then
            dns%dt = min(dns%cflmax/dns%cfl, dns%dtmax)
        end if

        if (dns%field_interval > 0) then
            call maybe_write_field(f, dns, g, int(dns%step_current), c, bc, ps%nIter, ps%sor)
        end if

    end do
    call system_clock(count=loop_clock_end)
    loop_seconds = real(loop_clock_end - loop_clock_start, C_DOUBLE) / real(clock_rate, C_DOUBLE)
    seconds_per_step = loop_seconds / real(dns%nsteps, C_DOUBLE)

    if (c%has_terminal) then
        print *, "main loop ended..."
        write(*,'(A,1X,I0,1X,A,1X,ES16.8,1X,A,1X,ES16.8)') &
            "timing: nsteps", dns%nsteps, "loop_seconds", loop_seconds, "seconds_per_step", seconds_per_step
    end if

    ! Release device-side data before the host allocatables go out of scope.
    call exit_ibm_data(ibm, dns)
    call exit_field_data(f, dns)
    call exit_grid_data(g, dns)
    call destroy_grid(g)
    call destroy_boundary_faces(bc)
    call comm_finalize(c)
end program main
