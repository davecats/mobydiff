module turbulence
    use, intrinsic :: iso_c_binding
    use :: chron, only: profiler_type, init_profiler
    use :: init, only: VAR_U, VAR_V, VAR_W, VAR_P
    use :: blocks, only: block_set_type
    implicit none

    private

    ! Turbulence model families. Only NONE and LES are implemented; TURB_RANS
    ! and TURB_IDDES are reserved for the k-omega SST and IDDES phases
    ! (docs/next_session_iddes.md).
    integer(C_INT), parameter, public :: TURB_NONE = 0_C_INT
    integer(C_INT), parameter, public :: TURB_LES = 1_C_INT
    integer(C_INT), parameter, public :: TURB_RANS = 2_C_INT
    integer(C_INT), parameter, public :: TURB_IDDES = 3_C_INT
    ! Category indices into the turbulence timing profiler (see init_turbulence_profiler).
    integer, parameter, public :: TURB_PROF_NUT = 1
    integer, parameter, public :: TURB_PROF_EXCHANGE = 2
    integer, parameter, public :: TURB_PROF_SGS = 3

    ! Model-agnostic turbulence state. Every producer (the LES today; RANS and
    ! the IDDES blend later) fills the same cell-centred eddy viscosity nut, so
    ! everything downstream of nut (halo exchange, momentum correction, dt
    ! limit, io) never needs to know which model ran. The 1D tables are grid
    ! metrics shared by the producers' stencils, not LES physics.
    type, public :: turb_type
        integer(C_INT) :: model = TURB_NONE
        real(C_DOUBLE), allocatable :: nut(:,:,:,:)   ! (0:nb+1,...,nBlocks)
        ! Per-block 1D metric tables (trailing block index).
        real(C_DOUBLE), allocatable :: filter_x(:,:), filter_y(:,:), filter_z(:,:)
        real(C_DOUBLE), allocatable :: d1xm(:,:,:), d1x0(:,:,:), d1xp(:,:,:)
        real(C_DOUBLE), allocatable :: d1ym(:,:,:), d1y0(:,:,:), d1yp(:,:,:)
        real(C_DOUBLE), allocatable :: d1zm(:,:,:), d1z0(:,:,:), d1zp(:,:,:)
        real(C_DOUBLE), allocatable :: p_from_u_x(:,:), p_from_v_y(:,:), p_from_w_z(:,:)
        real(C_DOUBLE), allocatable :: u_from_p_x(:,:), v_from_p_y(:,:), w_from_p_z(:,:)
        real(C_DOUBLE), allocatable :: inv_dx(:,:,:), inv_dy(:,:,:), inv_dz(:,:,:)
    end type turb_type

    public :: init_turbulence, destroy_turbulence, enter_turbulence_data, exit_turbulence_data
    public :: turbulence_is_enabled
    public :: init_turbulence_profiler

contains

    ! Build the turbulence timing profiler: output tag "turb_timing" and the
    ! three phase labels, ordered to match the TURB_PROF_* category indices.
    subroutine init_turbulence_profiler(profile)
        type(profiler_type), intent(out) :: profile

        call init_profiler(profile, "turb_timing", &
            [character(len=24) :: "nut_update", "nut_exchange", "sgs_correction"])
    end subroutine init_turbulence_profiler

    logical function turbulence_is_enabled(turb)
        type(turb_type), intent(in) :: turb

        turbulence_is_enabled = turb%model /= TURB_NONE
    end function turbulence_is_enabled

    subroutine init_turbulence(turb, blk)
        type(turb_type), intent(inout) :: turb
        type(block_set_type), intent(in) :: blk

        integer :: nx, ny, nz

        call destroy_turbulence(turb)
        if (.not. turbulence_is_enabled(turb)) return

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))

        allocate(turb%nut(0:nx+1,0:ny+1,0:nz+1,blk%nBlocks))
        allocate(turb%filter_x(0:nx+1,blk%nBlocks), turb%filter_y(0:ny+1,blk%nBlocks), &
            turb%filter_z(0:nz+1,blk%nBlocks))
        allocate(turb%d1xm(0:nx+1,VAR_U:VAR_P,blk%nBlocks), turb%d1x0(0:nx+1,VAR_U:VAR_P,blk%nBlocks), &
            turb%d1xp(0:nx+1,VAR_U:VAR_P,blk%nBlocks))
        allocate(turb%d1ym(0:ny+1,VAR_U:VAR_P,blk%nBlocks), turb%d1y0(0:ny+1,VAR_U:VAR_P,blk%nBlocks), &
            turb%d1yp(0:ny+1,VAR_U:VAR_P,blk%nBlocks))
        allocate(turb%d1zm(0:nz+1,VAR_U:VAR_P,blk%nBlocks), turb%d1z0(0:nz+1,VAR_U:VAR_P,blk%nBlocks), &
            turb%d1zp(0:nz+1,VAR_U:VAR_P,blk%nBlocks))
        allocate(turb%p_from_u_x(0:nx+1,blk%nBlocks), turb%p_from_v_y(0:ny+1,blk%nBlocks), &
            turb%p_from_w_z(0:nz+1,blk%nBlocks))
        allocate(turb%u_from_p_x(0:nx+1,blk%nBlocks), turb%v_from_p_y(0:ny+1,blk%nBlocks), &
            turb%w_from_p_z(0:nz+1,blk%nBlocks))
        allocate(turb%inv_dx(0:nx+1,VAR_U:VAR_P,blk%nBlocks), turb%inv_dy(0:ny+1,VAR_U:VAR_P,blk%nBlocks), &
            turb%inv_dz(0:nz+1,VAR_U:VAR_P,blk%nBlocks))

        turb%nut = 0.0d0
        call precompute_turbulence_metrics(turb, blk, nx, ny, nz)
    end subroutine init_turbulence

    subroutine destroy_turbulence(turb)
        type(turb_type), intent(inout) :: turb

        if (allocated(turb%nut)) deallocate(turb%nut)
        if (allocated(turb%filter_x)) deallocate(turb%filter_x)
        if (allocated(turb%filter_y)) deallocate(turb%filter_y)
        if (allocated(turb%filter_z)) deallocate(turb%filter_z)
        if (allocated(turb%d1xm)) deallocate(turb%d1xm)
        if (allocated(turb%d1x0)) deallocate(turb%d1x0)
        if (allocated(turb%d1xp)) deallocate(turb%d1xp)
        if (allocated(turb%d1ym)) deallocate(turb%d1ym)
        if (allocated(turb%d1y0)) deallocate(turb%d1y0)
        if (allocated(turb%d1yp)) deallocate(turb%d1yp)
        if (allocated(turb%d1zm)) deallocate(turb%d1zm)
        if (allocated(turb%d1z0)) deallocate(turb%d1z0)
        if (allocated(turb%d1zp)) deallocate(turb%d1zp)
        if (allocated(turb%p_from_u_x)) deallocate(turb%p_from_u_x)
        if (allocated(turb%p_from_v_y)) deallocate(turb%p_from_v_y)
        if (allocated(turb%p_from_w_z)) deallocate(turb%p_from_w_z)
        if (allocated(turb%u_from_p_x)) deallocate(turb%u_from_p_x)
        if (allocated(turb%v_from_p_y)) deallocate(turb%v_from_p_y)
        if (allocated(turb%w_from_p_z)) deallocate(turb%w_from_p_z)
        if (allocated(turb%inv_dx)) deallocate(turb%inv_dx)
        if (allocated(turb%inv_dy)) deallocate(turb%inv_dy)
        if (allocated(turb%inv_dz)) deallocate(turb%inv_dz)
    end subroutine destroy_turbulence

    subroutine enter_turbulence_data(turb)
        type(turb_type), intent(inout) :: turb

        if (.not. allocated(turb%nut)) return

        !$omp target enter data map(to: turb)
        !$omp target enter data map(to: turb%nut)
        !$omp target enter data map(to: turb%filter_x, turb%filter_y, turb%filter_z)
        !$omp target enter data map(to: turb%d1xm, turb%d1x0, turb%d1xp)
        !$omp target enter data map(to: turb%d1ym, turb%d1y0, turb%d1yp)
        !$omp target enter data map(to: turb%d1zm, turb%d1z0, turb%d1zp)
        !$omp target enter data map(to: turb%p_from_u_x, turb%p_from_v_y, turb%p_from_w_z)
        !$omp target enter data map(to: turb%u_from_p_x, turb%v_from_p_y, turb%w_from_p_z)
        !$omp target enter data map(to: turb%inv_dx, turb%inv_dy, turb%inv_dz)
    end subroutine enter_turbulence_data

    subroutine exit_turbulence_data(turb)
        type(turb_type), intent(inout) :: turb

        if (.not. allocated(turb%nut)) return

        !$omp target exit data map(delete: turb%inv_dx, turb%inv_dy, turb%inv_dz)
        !$omp target exit data map(delete: turb%u_from_p_x, turb%v_from_p_y, turb%w_from_p_z)
        !$omp target exit data map(delete: turb%p_from_u_x, turb%p_from_v_y, turb%p_from_w_z)
        !$omp target exit data map(delete: turb%d1zm, turb%d1z0, turb%d1zp)
        !$omp target exit data map(delete: turb%d1ym, turb%d1y0, turb%d1yp)
        !$omp target exit data map(delete: turb%d1xm, turb%d1x0, turb%d1xp)
        !$omp target exit data map(delete: turb%filter_x, turb%filter_y, turb%filter_z)
        !$omp target exit data map(delete: turb%nut)
        !$omp target exit data map(delete: turb)
    end subroutine exit_turbulence_data

    subroutine precompute_turbulence_metrics(turb, blk, nx, ny, nz)
        type(turb_type), intent(inout) :: turb
        type(block_set_type), intent(in) :: blk
        integer, intent(in) :: nx, ny, nz

        integer :: i, var, b

        do b = 1, int(blk%nBlocks)
        do i = 0, nx+1
            turb%filter_x(i,b) = max(1.0d0/blk%d1x(i,VAR_P,b), 1.0d-30)**(1.0d0/3.0d0)
            turb%p_from_u_x(i,b) = linear_weight(blk%x(i,VAR_U,b), blk%x(i+1,VAR_U,b), blk%x(i,VAR_P,b))
            turb%u_from_p_x(i,b) = linear_weight(blk%x(i-1,VAR_P,b), blk%x(i,VAR_P,b), blk%x(i,VAR_U,b))
            do var = VAR_U, VAR_P
                call first_derivative_coeffs(blk%x(i-1,var,b), blk%x(i,var,b), blk%x(i+1,var,b), &
                    turb%d1xm(i,var,b), turb%d1x0(i,var,b), turb%d1xp(i,var,b))
                turb%inv_dx(i,var,b) = safe_inv_delta(blk%x(i,var,b) - blk%x(i-1,var,b))
            end do
        end do

        do i = 0, ny+1
            turb%filter_y(i,b) = max(1.0d0/blk%d1y(i,VAR_P,b), 1.0d-30)**(1.0d0/3.0d0)
            turb%p_from_v_y(i,b) = linear_weight(blk%y(i,VAR_V,b), blk%y(i+1,VAR_V,b), blk%y(i,VAR_P,b))
            turb%v_from_p_y(i,b) = linear_weight(blk%y(i-1,VAR_P,b), blk%y(i,VAR_P,b), blk%y(i,VAR_V,b))
            do var = VAR_U, VAR_P
                call first_derivative_coeffs(blk%y(i-1,var,b), blk%y(i,var,b), blk%y(i+1,var,b), &
                    turb%d1ym(i,var,b), turb%d1y0(i,var,b), turb%d1yp(i,var,b))
                turb%inv_dy(i,var,b) = safe_inv_delta(blk%y(i,var,b) - blk%y(i-1,var,b))
            end do
        end do

        do i = 0, nz+1
            turb%filter_z(i,b) = max(1.0d0/blk%d1z(i,VAR_P,b), 1.0d-30)**(1.0d0/3.0d0)
            turb%p_from_w_z(i,b) = linear_weight(blk%z(i,VAR_W,b), blk%z(i+1,VAR_W,b), blk%z(i,VAR_P,b))
            turb%w_from_p_z(i,b) = linear_weight(blk%z(i-1,VAR_P,b), blk%z(i,VAR_P,b), blk%z(i,VAR_W,b))
            do var = VAR_U, VAR_P
                call first_derivative_coeffs(blk%z(i-1,var,b), blk%z(i,var,b), blk%z(i+1,var,b), &
                    turb%d1zm(i,var,b), turb%d1z0(i,var,b), turb%d1zp(i,var,b))
                turb%inv_dz(i,var,b) = safe_inv_delta(blk%z(i,var,b) - blk%z(i-1,var,b))
            end do
        end do
        end do
    end subroutine precompute_turbulence_metrics

    pure subroutine first_derivative_coeffs(xm, x0, xp, cm, c0, cp)
        real(C_DOUBLE), intent(in) :: xm, x0, xp
        real(C_DOUBLE), intent(out) :: cm, c0, cp

        real(C_DOUBLE) :: hm, hp

        hm = max(x0 - xm, 1.0d-30)
        hp = max(xp - x0, 1.0d-30)
        cm = -hp/(hm*(hm + hp))
        c0 = (hp - hm)/(hm*hp)
        cp = hm/(hp*(hm + hp))
    end subroutine first_derivative_coeffs

    pure real(C_DOUBLE) function linear_weight(x0, x1, x) result(w)
        real(C_DOUBLE), intent(in) :: x0, x1, x

        if (abs(x1 - x0) <= 1.0d-30) then
            w = 0.5d0
        else
            w = (x - x0)/(x1 - x0)
        end if
    end function linear_weight

    pure real(C_DOUBLE) function safe_inv_delta(delta) result(inv)
        real(C_DOUBLE), intent(in) :: delta

        inv = 1.0d0/max(delta, 1.0d-30)
    end function safe_inv_delta

end module turbulence
