module volume_force
    !=========================================================================
    !  Steady (time-independent) spatially-varying volumetric body force.
    !
    !  This module is meant to be edited by students, exactly like wall_bc.
    !  There is ONE place to change: the subroutine `body_force` below. It
    !  returns the force (fx, fy, fz) at a point (x, y, z) and is evaluated
    !  ONCE, at start-up, because the force does not change in time.
    !
    !  The force is OFF unless the input file turns it on:
    !        [force]
    !        enabled = true
    !        type    = steady
    !  With those two lines the solver fills the force from `body_force` and
    !  adds it to the momentum equation every step. Without them the body
    !  force is inert and the run is unchanged.
    !=========================================================================
    use, intrinsic :: iso_c_binding
    use :: init, only: VAR_U, VAR_V, VAR_W, NVEL
    use :: blocks, only: block_set_type
    implicit none
    private

    real(C_DOUBLE), parameter, public :: PI = 3.14159265358979323846d0

    public :: fill_volume_force

contains

    !-------------------------------------------------------------------------
    !  STUDENT-EDITABLE STEADY BODY FORCE
    !
    !  Return the volumetric force (fx, fy, fz) at the point (x, y, z).
    !    x -> streamwise   y -> wall-normal (0 at the bottom wall)   z -> spanwise
    !  The default is no force.
    !
    !  Example -- Schlatter & Canton steady streamwise vortices for skin-friction
    !  control (doi:10.1007/s10494-016-9723-8). The force is invariant in x and
    !  organises a row of streamwise rolls; with the wall-normal coordinate yc
    !  measured from the channel centre it reads
    !        fx = 0
    !        fy = A*beta*cos(beta*z)*(1 + cos(pi*yc/h))
    !        fz = A*(pi/h)*sin(beta*z)*sin(pi*yc/h)
    !  where A is the amplitude, h the channel half-height (centre at y = h) and
    !  beta the spanwise wavenumber, which MUST be an integer multiple of
    !  2*pi/Lz so the force is periodic in z. To use it, replace the body below
    !  with:
    !
    !        real(C_DOUBLE) :: A, beta, h, yc
    !        A    = 1.0d0                 ! amplitude
    !        h    = 1.0d0                 ! channel half-height (Ly/2)
    !        beta = 2.0d0                 ! spanwise wavenumber (multiple of 2*pi/Lz)
    !        yc   = y - h                 ! wall-normal distance from the centre
    !        fx = 0.0d0
    !        fy = A*beta*cos(beta*z)*(1.0d0 + cos(PI*yc/h))
    !        fz = A*(PI/h)*sin(beta*z)*sin(PI*yc/h)
    !-------------------------------------------------------------------------
    pure subroutine body_force(x, y, z, fx, fy, fz)
        real(C_DOUBLE), intent(in)  :: x, y, z
        real(C_DOUBLE), intent(out) :: fx, fy, fz

        fx = 0.0d0
        fy = 0.0d0
        fz = 0.0d0
    end subroutine body_force

    !-------------------------------------------------------------------------
    !  Driver (no need to edit): evaluate `body_force` once at each cell and
    !  store it into the force array f. Each component is sampled at its own
    !  staggered face location, matching how the momentum predictor reads it.
    !-------------------------------------------------------------------------
    subroutine fill_volume_force(f, blk)
        real(C_DOUBLE),       intent(inout) :: f(:,:,:,:,:)   ! (1:nb,1:nb,1:nb,NVEL,nBlocks)
        type(block_set_type), intent(in)    :: blk

        integer :: i, j, k, v, b, nx, ny, nz
        real(C_DOUBLE) :: fx, fy, fz

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))

        do b = 1, int(blk%nBlocks)
            do v = 1, int(NVEL)
                do k = 1, nz
                    do j = 1, ny
                        do i = 1, nx
                            ! Component v lives at blk%{x,y,z}(:,v,b).
                            call body_force(blk%x(i,v,b), blk%y(j,v,b), blk%z(k,v,b), fx, fy, fz)
                            select case (v)
                            case (VAR_U); f(i,j,k,v,b) = fx
                            case (VAR_V); f(i,j,k,v,b) = fy
                            case (VAR_W); f(i,j,k,v,b) = fz
                            end select
                        end do
                    end do
                end do
            end do
        end do
    end subroutine fill_volume_force

end module volume_force
