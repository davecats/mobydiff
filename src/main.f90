! TODO:
!
!      (*) second-order IBM in space, then in time: implement bisection with direction
!       -  IBM in Poisson equation for redblack
!       -  multiGPU parallelisation
!       -  various optimisation



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
    implicit none

    integer :: i, arg_status, rkStage
    integer :: output_interval, field_interval
    logical :: need_cfl
    real(C_DOUBLE) :: dt_alpha, dt_beta, dt_gamma
    real(C_DOUBLE) :: loop_seconds, seconds_per_step
    integer(int64) :: loop_clock_start, loop_clock_end, clock_rate
    character(len=256) :: input_file, output_prefix, field_prefix
    type(grid_type) :: g
    type(field_type) :: f
    type(boundary_type) :: bc
    type(pressure_solver_type) :: ps
    type(output_workspace_type) :: out
    type(ibm_type) :: ibm

    ! The first command-line argument can override the default input file.
    call get_command_argument(1, input_file, status=arg_status)
    if (arg_status /= 0 .or. len_trim(input_file) == 0) input_file = "input.ini"

    print *, "initialising grid..."
    call init_grid(g)

    print *, "initialising BCs..."
    call init_bc(bc)
    call read_runtime_config(g, ps, bc, output_interval, output_prefix, field_interval, field_prefix, input_file)
    call init_openmp_offload()

    print *, "initialising fields..."
    call init_field(f, g)
    call enter_field_data(f, g)
    call init_output_workspace(out, g, output_interval)

    print *, "initialising pressure solver..."
    call init_pressure_solver(ps, g)

    print *, "initialising IBM..."
    call init_ibm(ibm, g)
    call set_ibm_coeff(g, ibm, ibm%coef_u, 1, 0, 0)
    call set_ibm_coeff(g, ibm, ibm%coef_v, 0, 1, 0)
    call set_ibm_coeff(g, ibm, ibm%coef_w, 0, 0, 1)
    call enter_ibm_data(ibm, g)

    print *, "main loop starting..."
    call system_clock(count_rate=clock_rate)
    call system_clock(count=loop_clock_start)
    do i = 1, g%nsteps
        g%t_current = g%t_current + g%dt

        do rkStage = 1,3
            dt_alpha = g%dt*rk_alpha(rkStage)
            dt_beta  = g%dt*rk_beta(rkStage)
            dt_gamma = g%dt*rk_gamma(rkStage)

            ! Predictor: advance tentative staggered velocities, then enforce solid/body constraints.
            call momentum(f, g, dt_alpha, dt_beta, dt_gamma, ibm, bc)

#ifndef USE_REDBLACK
            call apply_bc(f, g)
#endif

            ! Projection: solve for pressure correction and project tentative velocities.
            call pressure_projection(ps, f, g, dt_gamma, ibm, bc)
#ifdef USE_REDBLACK
            call apply_bc_redblack(f, g, bc)
#else
            call apply_bc(f, g)
#endif

        end do

        ! Compute CFL only when it drives dt or is needed for output reporting.
        need_cfl = (g%cflmax > 0.0d0) .or. output_is_due(i, output_interval)
        if (need_cfl) g%cfl = get_cfl(f, g)

        if (g%cflmax > 0.0d0 .and. g%cfl > 0.0d0) then
            g%dt = min(g%cflmax/g%cfl, g%dtmax)
        end if

        if (output_interval > 0) then
            call maybe_write_vtk(out, f, g, i, output_interval, output_prefix)
        end if
        if (field_interval > 0) then
            call maybe_write_field(f, g, i, field_interval, field_prefix)
        end if

    end do
    call system_clock(count=loop_clock_end)
    loop_seconds = real(loop_clock_end - loop_clock_start, C_DOUBLE) / real(clock_rate, C_DOUBLE)
    seconds_per_step = loop_seconds / real(g%nsteps, C_DOUBLE)

    print *, "main loop ended..."
    write(*,'(A,1X,I0,1X,A,1X,ES16.8,1X,A,1X,ES16.8)') &
        "timing: nsteps", g%nsteps, "loop_seconds", loop_seconds, "seconds_per_step", seconds_per_step

    ! Release device-side data before the host allocatables go out of scope.
    call destroy_pressure_solver(ps)
    call exit_ibm_data(ibm, g)
    call destroy_output_workspace(out, g)
    call exit_field_data(f, g)
end program main
