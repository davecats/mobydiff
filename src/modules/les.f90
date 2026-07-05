module les_model
    use, intrinsic :: iso_c_binding
    use :: chron, only: profiler_type, init_profiler
    use :: init, only: dns_type, VAR_U, VAR_V, VAR_W, VAR_P
    use :: blocks, only: block_set_type
    use :: ibmm, only: ibm_type
    implicit none

    private

    integer(C_INT), parameter, public :: LES_NONE = 0_C_INT
    integer(C_INT), parameter, public :: LES_SMAGORINSKY = 1_C_INT
    integer(C_INT), parameter, public :: LES_WALE = 2_C_INT
    ! Category indices into the LES timing profiler (see init_les_profiler).
    integer, parameter, public :: LES_PROF_NUT = 1
    integer, parameter, public :: LES_PROF_EXCHANGE = 2
    integer, parameter, public :: LES_PROF_SGS = 3

    type, public :: les_type
        integer(C_INT) :: model = LES_NONE
        real(C_DOUBLE) :: cs = 0.10d0
        real(C_DOUBLE) :: cw = 0.325d0
        real(C_DOUBLE) :: delta_scale = 1.0d0
        logical(C_BOOL) :: ibm_aware = .true.
        real(C_DOUBLE), allocatable :: nut(:,:,:,:)   ! (0:nb+1,...,nBlocks)
        ! Per-block 1D metric tables (trailing block index).
        real(C_DOUBLE), allocatable :: filter_x(:,:), filter_y(:,:), filter_z(:,:)
        real(C_DOUBLE), allocatable :: d1xm(:,:,:), d1x0(:,:,:), d1xp(:,:,:)
        real(C_DOUBLE), allocatable :: d1ym(:,:,:), d1y0(:,:,:), d1yp(:,:,:)
        real(C_DOUBLE), allocatable :: d1zm(:,:,:), d1z0(:,:,:), d1zp(:,:,:)
        real(C_DOUBLE), allocatable :: p_from_u_x(:,:), p_from_v_y(:,:), p_from_w_z(:,:)
        real(C_DOUBLE), allocatable :: u_from_p_x(:,:), v_from_p_y(:,:), w_from_p_z(:,:)
        real(C_DOUBLE), allocatable :: inv_dx(:,:,:), inv_dy(:,:,:), inv_dz(:,:,:)
    end type les_type

    public :: init_les, destroy_les, enter_les_data, exit_les_data
    public :: les_is_enabled, update_les_viscosity
    public :: init_les_profiler

contains

    ! Build the LES timing profiler: output tag "les_timing" and the three
    ! phase labels, ordered to match the LES_PROF_* category indices.
    subroutine init_les_profiler(profile)
        type(profiler_type), intent(out) :: profile

        call init_profiler(profile, "les_timing", &
            [character(len=24) :: "nut_update", "nut_exchange", "sgs_correction"])
    end subroutine init_les_profiler

    logical function les_is_enabled(les)
        type(les_type), intent(in) :: les

        les_is_enabled = les%model /= LES_NONE
    end function les_is_enabled

    subroutine init_les(les, dns, blk)
        type(les_type), intent(inout) :: les
        type(dns_type), intent(in) :: dns
        type(block_set_type), intent(in) :: blk

        integer :: nx, ny, nz

        call destroy_les(les)
        if (.not. les_is_enabled(les)) return

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))

        allocate(les%nut(0:nx+1,0:ny+1,0:nz+1,blk%nBlocks))
        allocate(les%filter_x(0:nx+1,blk%nBlocks), les%filter_y(0:ny+1,blk%nBlocks), &
            les%filter_z(0:nz+1,blk%nBlocks))
        allocate(les%d1xm(0:nx+1,VAR_U:VAR_P,blk%nBlocks), les%d1x0(0:nx+1,VAR_U:VAR_P,blk%nBlocks), &
            les%d1xp(0:nx+1,VAR_U:VAR_P,blk%nBlocks))
        allocate(les%d1ym(0:ny+1,VAR_U:VAR_P,blk%nBlocks), les%d1y0(0:ny+1,VAR_U:VAR_P,blk%nBlocks), &
            les%d1yp(0:ny+1,VAR_U:VAR_P,blk%nBlocks))
        allocate(les%d1zm(0:nz+1,VAR_U:VAR_P,blk%nBlocks), les%d1z0(0:nz+1,VAR_U:VAR_P,blk%nBlocks), &
            les%d1zp(0:nz+1,VAR_U:VAR_P,blk%nBlocks))
        allocate(les%p_from_u_x(0:nx+1,blk%nBlocks), les%p_from_v_y(0:ny+1,blk%nBlocks), &
            les%p_from_w_z(0:nz+1,blk%nBlocks))
        allocate(les%u_from_p_x(0:nx+1,blk%nBlocks), les%v_from_p_y(0:ny+1,blk%nBlocks), &
            les%w_from_p_z(0:nz+1,blk%nBlocks))
        allocate(les%inv_dx(0:nx+1,VAR_U:VAR_P,blk%nBlocks), les%inv_dy(0:ny+1,VAR_U:VAR_P,blk%nBlocks), &
            les%inv_dz(0:nz+1,VAR_U:VAR_P,blk%nBlocks))

        les%nut = 0.0d0
        call precompute_les_metrics(les, blk, nx, ny, nz)
    end subroutine init_les

    subroutine destroy_les(les)
        type(les_type), intent(inout) :: les

        if (allocated(les%nut)) deallocate(les%nut)
        if (allocated(les%filter_x)) deallocate(les%filter_x)
        if (allocated(les%filter_y)) deallocate(les%filter_y)
        if (allocated(les%filter_z)) deallocate(les%filter_z)
        if (allocated(les%d1xm)) deallocate(les%d1xm)
        if (allocated(les%d1x0)) deallocate(les%d1x0)
        if (allocated(les%d1xp)) deallocate(les%d1xp)
        if (allocated(les%d1ym)) deallocate(les%d1ym)
        if (allocated(les%d1y0)) deallocate(les%d1y0)
        if (allocated(les%d1yp)) deallocate(les%d1yp)
        if (allocated(les%d1zm)) deallocate(les%d1zm)
        if (allocated(les%d1z0)) deallocate(les%d1z0)
        if (allocated(les%d1zp)) deallocate(les%d1zp)
        if (allocated(les%p_from_u_x)) deallocate(les%p_from_u_x)
        if (allocated(les%p_from_v_y)) deallocate(les%p_from_v_y)
        if (allocated(les%p_from_w_z)) deallocate(les%p_from_w_z)
        if (allocated(les%u_from_p_x)) deallocate(les%u_from_p_x)
        if (allocated(les%v_from_p_y)) deallocate(les%v_from_p_y)
        if (allocated(les%w_from_p_z)) deallocate(les%w_from_p_z)
        if (allocated(les%inv_dx)) deallocate(les%inv_dx)
        if (allocated(les%inv_dy)) deallocate(les%inv_dy)
        if (allocated(les%inv_dz)) deallocate(les%inv_dz)
    end subroutine destroy_les

    subroutine enter_les_data(les, dns)
        type(les_type), intent(inout) :: les
        type(dns_type), intent(in) :: dns

        if (.not. allocated(les%nut)) return

        !$omp target enter data map(to: les)
        !$omp target enter data map(to: les%nut)
        !$omp target enter data map(to: les%filter_x, les%filter_y, les%filter_z)
        !$omp target enter data map(to: les%d1xm, les%d1x0, les%d1xp)
        !$omp target enter data map(to: les%d1ym, les%d1y0, les%d1yp)
        !$omp target enter data map(to: les%d1zm, les%d1z0, les%d1zp)
        !$omp target enter data map(to: les%p_from_u_x, les%p_from_v_y, les%p_from_w_z)
        !$omp target enter data map(to: les%u_from_p_x, les%v_from_p_y, les%w_from_p_z)
        !$omp target enter data map(to: les%inv_dx, les%inv_dy, les%inv_dz)
    end subroutine enter_les_data

    subroutine exit_les_data(les, dns)
        type(les_type), intent(inout) :: les
        type(dns_type), intent(in) :: dns

        if (.not. allocated(les%nut)) return

        !$omp target exit data map(delete: les%inv_dx, les%inv_dy, les%inv_dz)
        !$omp target exit data map(delete: les%u_from_p_x, les%v_from_p_y, les%w_from_p_z)
        !$omp target exit data map(delete: les%p_from_u_x, les%p_from_v_y, les%p_from_w_z)
        !$omp target exit data map(delete: les%d1zm, les%d1z0, les%d1zp)
        !$omp target exit data map(delete: les%d1ym, les%d1y0, les%d1yp)
        !$omp target exit data map(delete: les%d1xm, les%d1x0, les%d1xp)
        !$omp target exit data map(delete: les%filter_x, les%filter_y, les%filter_z)
        !$omp target exit data map(delete: les%nut)
        !$omp target exit data map(delete: les)
    end subroutine exit_les_data

    subroutine precompute_les_metrics(les, blk, nx, ny, nz)
        type(les_type), intent(inout) :: les
        type(block_set_type), intent(in) :: blk
        integer, intent(in) :: nx, ny, nz

        integer :: i, var, b

        do b = 1, int(blk%nBlocks)
        do i = 0, nx+1
            les%filter_x(i,b) = max(1.0d0/blk%d1x(i,VAR_P,b), 1.0d-30)**(1.0d0/3.0d0)
            les%p_from_u_x(i,b) = linear_weight(blk%x(i,VAR_U,b), blk%x(i+1,VAR_U,b), blk%x(i,VAR_P,b))
            les%u_from_p_x(i,b) = linear_weight(blk%x(i-1,VAR_P,b), blk%x(i,VAR_P,b), blk%x(i,VAR_U,b))
            do var = VAR_U, VAR_P
                call first_derivative_coeffs(blk%x(i-1,var,b), blk%x(i,var,b), blk%x(i+1,var,b), &
                    les%d1xm(i,var,b), les%d1x0(i,var,b), les%d1xp(i,var,b))
                les%inv_dx(i,var,b) = safe_inv_delta(blk%x(i,var,b) - blk%x(i-1,var,b))
            end do
        end do

        do i = 0, ny+1
            les%filter_y(i,b) = max(1.0d0/blk%d1y(i,VAR_P,b), 1.0d-30)**(1.0d0/3.0d0)
            les%p_from_v_y(i,b) = linear_weight(blk%y(i,VAR_V,b), blk%y(i+1,VAR_V,b), blk%y(i,VAR_P,b))
            les%v_from_p_y(i,b) = linear_weight(blk%y(i-1,VAR_P,b), blk%y(i,VAR_P,b), blk%y(i,VAR_V,b))
            do var = VAR_U, VAR_P
                call first_derivative_coeffs(blk%y(i-1,var,b), blk%y(i,var,b), blk%y(i+1,var,b), &
                    les%d1ym(i,var,b), les%d1y0(i,var,b), les%d1yp(i,var,b))
                les%inv_dy(i,var,b) = safe_inv_delta(blk%y(i,var,b) - blk%y(i-1,var,b))
            end do
        end do

        do i = 0, nz+1
            les%filter_z(i,b) = max(1.0d0/blk%d1z(i,VAR_P,b), 1.0d-30)**(1.0d0/3.0d0)
            les%p_from_w_z(i,b) = linear_weight(blk%z(i,VAR_W,b), blk%z(i+1,VAR_W,b), blk%z(i,VAR_P,b))
            les%w_from_p_z(i,b) = linear_weight(blk%z(i-1,VAR_P,b), blk%z(i,VAR_P,b), blk%z(i,VAR_W,b))
            do var = VAR_U, VAR_P
                call first_derivative_coeffs(blk%z(i-1,var,b), blk%z(i,var,b), blk%z(i+1,var,b), &
                    les%d1zm(i,var,b), les%d1z0(i,var,b), les%d1zp(i,var,b))
                les%inv_dz(i,var,b) = safe_inv_delta(blk%z(i,var,b) - blk%z(i-1,var,b))
            end do
        end do
        end do
    end subroutine precompute_les_metrics

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

    subroutine update_les_viscosity(les, blk, dns, ibm)
        type(les_type), intent(inout) :: les
        type(block_set_type), intent(inout) :: blk
        type(dns_type), intent(in) :: dns
        type(ibm_type), intent(in) :: ibm

        if (.not. allocated(les%nut)) return

        select case (les%model)
        case (LES_NONE)
            return
        case (LES_WALE)
            call update_wale_viscosity(les, blk, dns, ibm)
        case default
            call update_generic_les_viscosity(les, blk, dns, ibm)
        end select
    end subroutine update_les_viscosity

    ! Velocity-gradient tensor g(a,b) = du_a/dx_b at the centre of cell (i,j,k)
    ! in block b, on the staggered mesh: the diagonal terms are the face
    ! difference across the cell; each off-diagonal term interpolates the
    ! neighbouring face-difference to the cell centre with the p_from_* staggered
    ! weights and the les d1?? stencils. Shared verbatim by the Smagorinsky and
    ! WALE kernels; declared target so both device kernels can call it.
    subroutine velocity_gradient_tensor(blk, les, i, j, k, b, &
            g11, g12, g13, g21, g22, g23, g31, g32, g33)
!$omp declare target
        type(block_set_type), intent(in) :: blk
        type(les_type), intent(in) :: les
        integer, intent(in) :: i, j, k, b
        real(C_DOUBLE), intent(out) :: g11, g12, g13, g21, g22, g23, g31, g32, g33

        real(C_DOUBLE) :: d0, d1

        g11 = (blk%q(i+1,j,k,VAR_U,b) - blk%q(i,j,k,VAR_U,b))*blk%d1x(i,VAR_P,b)
        g22 = (blk%q(i,j+1,k,VAR_V,b) - blk%q(i,j,k,VAR_V,b))*blk%d1y(j,VAR_P,b)
        g33 = (blk%q(i,j,k+1,VAR_W,b) - blk%q(i,j,k,VAR_W,b))*blk%d1z(k,VAR_P,b)

        d0 = les%d1ym(j,VAR_U,b)*blk%q(i,j-1,k,VAR_U,b) &
           + les%d1y0(j,VAR_U,b)*blk%q(i,j,k,VAR_U,b) &
           + les%d1yp(j,VAR_U,b)*blk%q(i,j+1,k,VAR_U,b)
        d1 = les%d1ym(j,VAR_U,b)*blk%q(i+1,j-1,k,VAR_U,b) &
           + les%d1y0(j,VAR_U,b)*blk%q(i+1,j,k,VAR_U,b) &
           + les%d1yp(j,VAR_U,b)*blk%q(i+1,j+1,k,VAR_U,b)
        g12 = (1.0d0 - les%p_from_u_x(i,b))*d0 + les%p_from_u_x(i,b)*d1

        d0 = les%d1zm(k,VAR_U,b)*blk%q(i,j,k-1,VAR_U,b) &
           + les%d1z0(k,VAR_U,b)*blk%q(i,j,k,VAR_U,b) &
           + les%d1zp(k,VAR_U,b)*blk%q(i,j,k+1,VAR_U,b)
        d1 = les%d1zm(k,VAR_U,b)*blk%q(i+1,j,k-1,VAR_U,b) &
           + les%d1z0(k,VAR_U,b)*blk%q(i+1,j,k,VAR_U,b) &
           + les%d1zp(k,VAR_U,b)*blk%q(i+1,j,k+1,VAR_U,b)
        g13 = (1.0d0 - les%p_from_u_x(i,b))*d0 + les%p_from_u_x(i,b)*d1

        d0 = les%d1xm(i,VAR_V,b)*blk%q(i-1,j,k,VAR_V,b) &
           + les%d1x0(i,VAR_V,b)*blk%q(i,j,k,VAR_V,b) &
           + les%d1xp(i,VAR_V,b)*blk%q(i+1,j,k,VAR_V,b)
        d1 = les%d1xm(i,VAR_V,b)*blk%q(i-1,j+1,k,VAR_V,b) &
           + les%d1x0(i,VAR_V,b)*blk%q(i,j+1,k,VAR_V,b) &
           + les%d1xp(i,VAR_V,b)*blk%q(i+1,j+1,k,VAR_V,b)
        g21 = (1.0d0 - les%p_from_v_y(j,b))*d0 + les%p_from_v_y(j,b)*d1

        d0 = les%d1zm(k,VAR_V,b)*blk%q(i,j,k-1,VAR_V,b) &
           + les%d1z0(k,VAR_V,b)*blk%q(i,j,k,VAR_V,b) &
           + les%d1zp(k,VAR_V,b)*blk%q(i,j,k+1,VAR_V,b)
        d1 = les%d1zm(k,VAR_V,b)*blk%q(i,j+1,k-1,VAR_V,b) &
           + les%d1z0(k,VAR_V,b)*blk%q(i,j+1,k,VAR_V,b) &
           + les%d1zp(k,VAR_V,b)*blk%q(i,j+1,k+1,VAR_V,b)
        g23 = (1.0d0 - les%p_from_v_y(j,b))*d0 + les%p_from_v_y(j,b)*d1

        d0 = les%d1xm(i,VAR_W,b)*blk%q(i-1,j,k,VAR_W,b) &
           + les%d1x0(i,VAR_W,b)*blk%q(i,j,k,VAR_W,b) &
           + les%d1xp(i,VAR_W,b)*blk%q(i+1,j,k,VAR_W,b)
        d1 = les%d1xm(i,VAR_W,b)*blk%q(i-1,j,k+1,VAR_W,b) &
           + les%d1x0(i,VAR_W,b)*blk%q(i,j,k+1,VAR_W,b) &
           + les%d1xp(i,VAR_W,b)*blk%q(i+1,j,k+1,VAR_W,b)
        g31 = (1.0d0 - les%p_from_w_z(k,b))*d0 + les%p_from_w_z(k,b)*d1

        d0 = les%d1ym(j,VAR_W,b)*blk%q(i,j-1,k,VAR_W,b) &
           + les%d1y0(j,VAR_W,b)*blk%q(i,j,k,VAR_W,b) &
           + les%d1yp(j,VAR_W,b)*blk%q(i,j+1,k,VAR_W,b)
        d1 = les%d1ym(j,VAR_W,b)*blk%q(i,j-1,k+1,VAR_W,b) &
           + les%d1y0(j,VAR_W,b)*blk%q(i,j,k+1,VAR_W,b) &
           + les%d1yp(j,VAR_W,b)*blk%q(i,j+1,k+1,VAR_W,b)
        g32 = (1.0d0 - les%p_from_w_z(k,b))*d0 + les%p_from_w_z(k,b)*d1
    end subroutine velocity_gradient_tensor

    subroutine update_generic_les_viscosity(les, blk, dns, ibm)
        type(les_type), intent(inout) :: les
        type(block_set_type), intent(inout) :: blk
        type(dns_type), intent(in) :: dns
        type(ibm_type), intent(in) :: ibm

        integer :: i, j, k, b, nx, ny, nz, nBlocks
        integer(C_INT) :: model
        logical(C_BOOL) :: ibm_aware, ibm_enabled
        real(C_DOUBLE) :: cs, cs2, delta_scale
        real(C_DOUBLE) :: delta
        real(C_DOUBLE) :: g11, g12, g13, g21, g22, g23, g31, g32, g33
        real(C_DOUBLE) :: s11, s22, s33, s12, s13, s23, s2, strain_mag
        real(C_DOUBLE) :: solid_threshold
        logical :: solid_cell

        if (.not. allocated(les%nut)) return

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))
        nBlocks = int(blk%nBlocks)
        model = les%model
        cs = les%cs
        cs2 = cs*cs
        delta_scale = les%delta_scale
        ibm_aware = les%ibm_aware
        ibm_enabled = dns%ibm_enabled
        solid_threshold = 1.0d20

        if (model == LES_NONE) return

        !$omp target teams distribute parallel do collapse(4) &
        !$omp& map(to: model, cs2, delta_scale, ibm_aware, ibm_enabled, solid_threshold, &
        !$omp& blk%q, blk%d1x, blk%d1y, blk%d1z, ibm%coef, &
        !$omp& les%filter_x, les%filter_y, les%filter_z, &
        !$omp& les%d1xm, les%d1x0, les%d1xp, les%d1ym, les%d1y0, les%d1yp, &
        !$omp& les%d1zm, les%d1z0, les%d1zp, &
        !$omp& les%p_from_u_x, les%p_from_v_y, les%p_from_w_z) &
        !$omp& map(tofrom: les%nut) &
        !$omp& private(i,j,k,b,delta,g11,g12,g13,g21,g22,g23,g31,g32,g33, &
        !$omp& s11,s22,s33,s12,s13,s23,s2,strain_mag,solid_cell)
        do b = 1, nBlocks
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    les%nut(i,j,k,b) = 0.0d0

                    solid_cell = .false.
                    if (ibm_aware .and. ibm_enabled) then
                        solid_cell = abs(ibm%coef(i,j,k,VAR_U,b)) > solid_threshold .or. &
                                     abs(ibm%coef(i+1,j,k,VAR_U,b)) > solid_threshold .or. &
                                     abs(ibm%coef(i,j,k,VAR_V,b)) > solid_threshold .or. &
                                     abs(ibm%coef(i,j+1,k,VAR_V,b)) > solid_threshold .or. &
                                     abs(ibm%coef(i,j,k,VAR_W,b)) > solid_threshold .or. &
                                     abs(ibm%coef(i,j,k+1,VAR_W,b)) > solid_threshold
                    end if
                    if (solid_cell) cycle

                    delta = delta_scale*les%filter_x(i,b)*les%filter_y(j,b)*les%filter_z(k,b)

                    call velocity_gradient_tensor(blk, les, i, j, k, b, &
                        g11, g12, g13, g21, g22, g23, g31, g32, g33)

                    s11 = g11
                    s22 = g22
                    s33 = g33
                    s12 = 0.5d0*(g12 + g21)
                    s13 = 0.5d0*(g13 + g31)
                    s23 = 0.5d0*(g23 + g32)
                    s2 = s11*s11 + s22*s22 + s33*s33 + 2.0d0*(s12*s12 + s13*s13 + s23*s23)

                    ! WALE is dispatched to update_wale_viscosity; this generic
                    ! path handles the Smagorinsky (and any future algebraic
                    ! eddy-viscosity) model.
                    select case (model)
                    case (LES_SMAGORINSKY)
                        strain_mag = sqrt(max(0.0d0, 2.0d0*s2))
                        les%nut(i,j,k,b) = cs2*delta*delta*strain_mag
                    end select
                end do
            end do
        end do
        end do
        !$omp end target teams distribute parallel do
    end subroutine update_generic_les_viscosity

    subroutine update_wale_viscosity(les, blk, dns, ibm)
        type(les_type), intent(inout) :: les
        type(block_set_type), intent(inout) :: blk
        type(dns_type), intent(in) :: dns
        type(ibm_type), intent(in) :: ibm

        integer :: i, j, k, b, nx, ny, nz, nBlocks
        logical(C_BOOL) :: ibm_aware, ibm_enabled
        real(C_DOUBLE) :: cw2, delta_scale, delta
        real(C_DOUBLE) :: g11, g12, g13, g21, g22, g23, g31, g32, g33
        real(C_DOUBLE) :: s2, g2_11, g2_22, g2_33, trace_g2
        real(C_DOUBLE) :: sd11, sd22, sd33, sd12, sd13, sd23, sd2, denom
        real(C_DOUBLE) :: sqrt_s2, sqrt_sd2, s2_52, sd2_32, sd2_54
        real(C_DOUBLE) :: solid_threshold
        logical :: solid_cell

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))
        nBlocks = int(blk%nBlocks)
        cw2 = les%cw*les%cw
        delta_scale = les%delta_scale
        ibm_aware = les%ibm_aware
        ibm_enabled = dns%ibm_enabled
        solid_threshold = 1.0d20

        !$omp target teams distribute parallel do collapse(4) &
        !$omp& map(to: cw2, delta_scale, ibm_aware, ibm_enabled, solid_threshold, &
        !$omp& blk%q, blk%d1x, blk%d1y, blk%d1z, ibm%coef, &
        !$omp& les%filter_x, les%filter_y, les%filter_z, &
        !$omp& les%d1xm, les%d1x0, les%d1xp, les%d1ym, les%d1y0, les%d1yp, &
        !$omp& les%d1zm, les%d1z0, les%d1zp, &
        !$omp& les%p_from_u_x, les%p_from_v_y, les%p_from_w_z) &
        !$omp& map(tofrom: les%nut) &
        !$omp& private(i,j,k,b,delta,g11,g12,g13,g21,g22,g23,g31,g32,g33, &
        !$omp& s2,g2_11,g2_22,g2_33,trace_g2,sd11,sd22,sd33,sd12,sd13,sd23, &
        !$omp& sd2,denom,sqrt_s2,sqrt_sd2,s2_52,sd2_32,sd2_54,solid_cell)
        do b = 1, nBlocks
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    les%nut(i,j,k,b) = 0.0d0

                    solid_cell = .false.
                    if (ibm_aware .and. ibm_enabled) then
                        solid_cell = abs(ibm%coef(i,j,k,VAR_U,b)) > solid_threshold .or. &
                                     abs(ibm%coef(i+1,j,k,VAR_U,b)) > solid_threshold .or. &
                                     abs(ibm%coef(i,j,k,VAR_V,b)) > solid_threshold .or. &
                                     abs(ibm%coef(i,j+1,k,VAR_V,b)) > solid_threshold .or. &
                                     abs(ibm%coef(i,j,k,VAR_W,b)) > solid_threshold .or. &
                                     abs(ibm%coef(i,j,k+1,VAR_W,b)) > solid_threshold
                    end if
                    if (solid_cell) cycle

                    delta = delta_scale*les%filter_x(i,b)*les%filter_y(j,b)*les%filter_z(k,b)

                    call velocity_gradient_tensor(blk, les, i, j, k, b, &
                        g11, g12, g13, g21, g22, g23, g31, g32, g33)

                    s2 = g11*g11 + g22*g22 + g33*g33 &
                       + 0.5d0*((g12 + g21)*(g12 + g21) &
                               + (g13 + g31)*(g13 + g31) &
                               + (g23 + g32)*(g23 + g32))

                    g2_11 = g11*g11 + g12*g21 + g13*g31
                    g2_22 = g21*g12 + g22*g22 + g23*g32
                    g2_33 = g31*g13 + g32*g23 + g33*g33

                    trace_g2 = (g2_11 + g2_22 + g2_33)/3.0d0
                    sd11 = g2_11 - trace_g2
                    sd22 = g2_22 - trace_g2
                    sd33 = g2_33 - trace_g2
                    sd12 = 0.5d0*(g11*g12 + g12*g22 + g13*g32 &
                                 + g21*g11 + g22*g21 + g23*g31)
                    sd13 = 0.5d0*(g11*g13 + g12*g23 + g13*g33 &
                                 + g31*g11 + g32*g21 + g33*g31)
                    sd23 = 0.5d0*(g21*g13 + g22*g23 + g23*g33 &
                                 + g31*g12 + g32*g22 + g33*g32)
                    sd2 = sd11*sd11 + sd22*sd22 + sd33*sd33 + &
                          2.0d0*(sd12*sd12 + sd13*sd13 + sd23*sd23)

                    sqrt_s2 = sqrt(s2)
                    sqrt_sd2 = sqrt(sd2)
                    s2_52 = s2*s2*sqrt_s2
                    sd2_32 = sd2*sqrt_sd2
                    sd2_54 = sd2*sqrt(sqrt_sd2)
                    denom = s2_52 + sd2_54 + 1.0d-30
                    les%nut(i,j,k,b) = cw2*delta*delta*sd2_32/denom
                end do
            end do
        end do
        end do
        !$omp end target teams distribute parallel do
    end subroutine update_wale_viscosity

end module les_model
