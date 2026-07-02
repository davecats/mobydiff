module wall_bc
    !=========================================================================
    !  Time- and space-varying velocity boundary conditions on the two walls.
    !
    !  This module is meant to be edited by students. There is exactly ONE
    !  place to change: the subroutine `wall_velocity` below. It prescribes
    !  the velocity of the bottom wall (y = y_min) and the top wall
    !  (y = y_max) as a function of the position along the wall (x, z), the
    !  simulation time t, and -- for feedback controls such as opposition
    !  control -- the wall-normal velocity sensed a short distance into the
    !  flow. The rest of the module (`update_wall_bc`) is the driver that
    !  evaluates your function at every wall point each substage and hands the
    !  values to the solver -- you do not need to touch it.
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

    ! Opposition control senses the wall-normal velocity on a "detection plane"
    ! this many cells away from the wall. Classic opposition control (Choi, Moin
    ! & Kim 1994) places it around y+ ~ 10-15; move it by changing this number.
    integer, parameter, public :: SENSOR_OFFSET = 12

    public :: update_wall_bc

contains

    !-------------------------------------------------------------------------
    !  STUDENT-EDITABLE WALL BOUNDARY CONDITION
    !
    !  Return the prescribed wall velocity (u, v, w) at a point on one of the
    !  two channel walls.
    !
    !    wall     : WALL_BOTTOM (y = y_min)  or  WALL_TOP (y = y_max)
    !    x, z     : coordinates of the point ALONG the wall
    !    t        : current simulation time
    !    v_sensed : wall-normal velocity measured on the detection plane
    !               SENSOR_OFFSET cells into the flow, above this point.
    !               Only needed for feedback (opposition) control; ignore it
    !               otherwise.
    !
    !  Components:
    !    u -> streamwise (x)   v -> wall-normal (y)   w -> spanwise (z)
    !
    !  The default is a stationary, no-slip wall: u = v = w = 0.
    !
    !  Examples (replace the body below to try them):
    !
    !    ! Spanwise oscillating wall (Stokes layer)
    !    w = 0.5d0 * sin(2.0d0 * t)
    !
    !    ! Streamwise travelling wave of spanwise wall velocity
    !    w = 2.0d0 * sin(1.0d0 * x - 0.5d0 * t)
    !
    !    ! Steady blowing / suction (careful: keep the net mass flux zero!)
    !    if (wall == WALL_BOTTOM) v =  0.01d0
    !    if (wall == WALL_TOP)    v =  0.01d0
    !
    !    ! Opposition control: oppose the sensed wall-normal velocity
    !    v = -v_sensed
    !-------------------------------------------------------------------------
    pure subroutine wall_velocity(wall, x, z, t, v_sensed, u, v, w)
        integer,        intent(in)  :: wall
        real(C_DOUBLE), intent(in)  :: x, z, t, v_sensed
        real(C_DOUBLE), intent(out) :: u, v, w

        u = 0.0d0
        v = 0.0d0
        w = 0.0d0
    end subroutine wall_velocity

    !-------------------------------------------------------------------------
    !  Driver (no need to edit). For every point of the two y-walls it:
    !    1. senses the wall-normal velocity on the detection plane
    !       (done where the field lives, so it is correct on the GPU too);
    !    2. evaluates the student `wall_velocity` law and stores the result as
    !       the per-point Dirichlet value that `apply_bc` enforces.
    !  Called once per Runge-Kutta substage so the condition can depend on
    !  time and on the current flow. Each component is sampled at its own
    !  staggered location.
    !-------------------------------------------------------------------------
    subroutine update_wall_bc(bc, blk, t)
        type(boundary_type),   intent(inout) :: bc
        type(block_set_type),  intent(in)    :: blk
        real(C_DOUBLE),        intent(in)    :: t

        integer :: n, npts, b, i, k, face_id, wall, dir, side, ny, j_det
        real(C_DOUBLE) :: u, v, w
        real(C_DOUBLE), allocatable :: v_sensed(:)

        npts = int(bc%nTotal)
        if (npts <= 0) return
        ny = int(blk%nb(2))

        allocate(v_sensed(npts))
        v_sensed = 0.0d0

        ! (1) Sense the wall-normal velocity on the detection plane for every
        !     wall point. This loop runs on the device (where the velocity field
        !     lives after the predictor); on the CPU it is an ordinary loop.
        !$omp target teams distribute parallel do &
        !$omp& map(to: npts, ny, blk%q, bc%pointFace(1:npts), bc%slot(1:npts), &
        !$omp& bc%i(1:npts), bc%k(1:npts)) map(tofrom: v_sensed(1:npts)) &
        !$omp& private(n, face_id, dir, side, b, i, k, j_det)
        do n = 1, npts
            face_id = int(bc%pointFace(n))
            dir = (face_id + 1)/2                       ! wall-normal direction
            if (dir /= 2) cycle                         ! DIR_Y walls only
            side = modulo(face_id - 1, 2)               ! 0 = bottom, 1 = top
            b = int(bc%slot(n))
            i = int(bc%i(n))
            k = int(bc%k(n))
            ! Detection plane SENSOR_OFFSET cells into the flow from this wall.
            j_det = merge(1 + SENSOR_OFFSET, ny + 1 - SENSOR_OFFSET, side == 0)
            v_sensed(n) = blk%q(i, j_det, k, VAR_V, b)
        end do

        ! (2) Evaluate the student control law on the host and store the wall
        !     boundary values.
        do n = 1, npts
            face_id = int(bc%pointFace(n))
            ! Only the two walls are wall-normal in y; skip any other physical face.
            if (boundary_face_dir(face_id) /= int(DIR_Y)) cycle
            wall = merge(WALL_BOTTOM, WALL_TOP, boundary_face_side(face_id) == int(SIDE_MIN))

            b = int(bc%slot(n))
            i = int(bc%i(n))
            k = int(bc%k(n))

            call wall_velocity(wall, blk%x(i,VAR_U,b), blk%z(k,VAR_U,b), t, v_sensed(n), u, v, w)
            bc%pointBcValue(VAR_U,n) = u
            call wall_velocity(wall, blk%x(i,VAR_V,b), blk%z(k,VAR_V,b), t, v_sensed(n), u, v, w)
            bc%pointBcValue(VAR_V,n) = v
            call wall_velocity(wall, blk%x(i,VAR_W,b), blk%z(k,VAR_W,b), t, v_sensed(n), u, v, w)
            bc%pointBcValue(VAR_W,n) = w
        end do

#ifdef USE_OPENMP_OFFLOAD
        ! The values are consumed inside a target region in apply_bc, so push
        ! the freshly computed host values to the device.
        !$omp target update to(bc%pointBcValue(VAR_U:VAR_P,1:npts))
#endif
        deallocate(v_sensed)
    end subroutine update_wall_bc

end module wall_bc
