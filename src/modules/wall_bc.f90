module wall_bc
    !=========================================================================
    !  Time- and space-varying velocity boundary conditions on the two walls.
    !
    !  This module is meant to be edited by students. There is exactly ONE
    !  place to change: the subroutine `wall_velocity` below. It prescribes
    !  the velocity of the bottom wall (y = y_min) and the top wall
    !  (y = y_max) as a function of the position along the wall (x, z) and of
    !  the simulation time t. The rest of the module (`update_wall_bc`) is the
    !  driver that evaluates your function at every wall point each time step
    !  and hands the values to the solver -- you do not need to touch it.
    !=========================================================================
    use, intrinsic :: iso_c_binding
    use :: init, only: VAR_U, VAR_V, VAR_W, VAR_P
    use :: blocks, only: block_set_type
    use :: boundary, only: boundary_type, boundary_face_dir, boundary_face_side, &
        DIR_Y, SIDE_MIN
    implicit none
    private

    ! Wall identifier passed to `wall_velocity`.
    integer, parameter, public :: WALL_BOTTOM = 0   ! the wall at y = y_min
    integer, parameter, public :: WALL_TOP    = 1   ! the wall at y = y_max

    public :: update_wall_bc

contains

    !-------------------------------------------------------------------------
    !  STUDENT-EDITABLE WALL BOUNDARY CONDITION
    !
    !  Return the prescribed wall velocity (u, v, w) at a point on one of the
    !  two channel walls.
    !
    !    wall : WALL_BOTTOM (y = y_min)  or  WALL_TOP (y = y_max)
    !    x, z : coordinates of the point ALONG the wall
    !    t    : current simulation time
    !
    !  Components:
    !    u -> streamwise (x)   v -> wall-normal (y)   w -> spanwise (z)
    !
    !  The default is a stationary, no-slip wall: u = v = w = 0.
    !
    !  A few examples (replace the body below to try them):
    !
    !    ! Moving bottom wall -> plane Couette flow
    !    if (wall == WALL_BOTTOM) u = 1.0d0
    !
    !    ! Wall oscillating in the spanwise direction (Stokes layer)
    !    w = 0.5d0 * sin(2.0d0 * t)
    !
    !    ! Spanwise travelling wave of wall velocity (drag control)
    !    w = 2.0d0 * sin(1.0d0 * x - 0.5d0 * t)
    !
    !    ! Steady blowing / suction (careful: keep the net mass flux zero!)
    !    if (wall == WALL_BOTTOM) v =  0.01d0
    !    if (wall == WALL_TOP)    v =  0.01d0
    !-------------------------------------------------------------------------
    pure subroutine wall_velocity(wall, x, z, t, u, v, w)
        integer,        intent(in)  :: wall
        real(C_DOUBLE), intent(in)  :: x, z, t
        real(C_DOUBLE), intent(out) :: u, v, w

        u = 0.0d0
        v = 0.0d0
        w = 0.0d0
    end subroutine wall_velocity

    !-------------------------------------------------------------------------
    !  Driver (no need to edit): evaluate `wall_velocity` at every point of
    !  the two y-walls and store it as the per-point Dirichlet boundary value
    !  that `apply_bc` enforces. Called once per Runge-Kutta substage so the
    !  boundary condition can depend on time. Each velocity component is
    !  sampled at its own staggered location.
    !-------------------------------------------------------------------------
    subroutine update_wall_bc(bc, blk, t)
        type(boundary_type),   intent(inout) :: bc
        type(block_set_type),  intent(in)    :: blk
        real(C_DOUBLE),        intent(in)    :: t

        integer :: n, npts, b, i, k, face_id, wall
        real(C_DOUBLE) :: u, v, w

        npts = int(bc%nTotal)
        if (npts <= 0) return

        do n = 1, npts
            face_id = int(bc%pointFace(n))
            ! Only the two walls are wall-normal in y; skip any other physical face.
            if (boundary_face_dir(face_id) /= int(DIR_Y)) cycle
            wall = merge(WALL_BOTTOM, WALL_TOP, boundary_face_side(face_id) == int(SIDE_MIN))

            b = int(bc%slot(n))
            i = int(bc%i(n))
            k = int(bc%k(n))

            call wall_velocity(wall, blk%x(i,VAR_U,b), blk%z(k,VAR_U,b), t, u, v, w)
            bc%pointBcValue(VAR_U,n) = u
            call wall_velocity(wall, blk%x(i,VAR_V,b), blk%z(k,VAR_V,b), t, u, v, w)
            bc%pointBcValue(VAR_V,n) = v
            call wall_velocity(wall, blk%x(i,VAR_W,b), blk%z(k,VAR_W,b), t, u, v, w)
            bc%pointBcValue(VAR_W,n) = w
        end do

#ifdef USE_OPENMP_OFFLOAD
        ! The values are consumed inside a target region in apply_bc, so push
        ! the freshly computed host values to the device.
        !$omp target update to(bc%pointBcValue(VAR_U:VAR_P,1:npts))
#endif
    end subroutine update_wall_bc

end module wall_bc
