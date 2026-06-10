module les_model
    use, intrinsic :: iso_c_binding
    use, intrinsic :: iso_fortran_env, only: int64
    use :: init, only: dns_type, grid_type, field_type, VAR_U, VAR_V, VAR_W, VAR_P
    use :: ibmm, only: ibm_type
    implicit none

    private

    integer(C_INT), parameter, public :: LES_NONE = 0_C_INT
    integer(C_INT), parameter, public :: LES_SMAGORINSKY = 1_C_INT
    integer(C_INT), parameter, public :: LES_WALE = 2_C_INT
    integer, parameter, public :: LES_PROF_NUT = 1
    integer, parameter, public :: LES_PROF_EXCHANGE = 2
    integer, parameter, public :: LES_PROF_SGS = 3
    integer, parameter, public :: LES_PROF_NCATS = 3

    type, public :: les_type
        integer(C_INT) :: model = LES_NONE
        real(C_DOUBLE) :: cs = 0.10d0
        real(C_DOUBLE) :: cw = 0.325d0
        real(C_DOUBLE) :: delta_scale = 1.0d0
        logical(C_BOOL) :: ibm_aware = .true.
        real(C_DOUBLE), allocatable :: nut(:,:,:)
        real(C_DOUBLE), allocatable :: filter_x(:), filter_y(:), filter_z(:)
        real(C_DOUBLE), allocatable :: d1xm(:,:), d1x0(:,:), d1xp(:,:)
        real(C_DOUBLE), allocatable :: d1ym(:,:), d1y0(:,:), d1yp(:,:)
        real(C_DOUBLE), allocatable :: d1zm(:,:), d1z0(:,:), d1zp(:,:)
        real(C_DOUBLE), allocatable :: p_from_u_x(:), p_from_v_y(:), p_from_w_z(:)
        real(C_DOUBLE), allocatable :: u_from_p_x(:), v_from_p_y(:), w_from_p_z(:)
        real(C_DOUBLE), allocatable :: inv_dx(:,:), inv_dy(:,:), inv_dz(:,:)
    end type les_type

    type, public :: les_profile_type
        real(C_DOUBLE) :: seconds(LES_PROF_NCATS) = 0.0d0
        integer(C_INT) :: calls(LES_PROF_NCATS) = 0_C_INT
    end type les_profile_type

    public :: init_les, destroy_les, enter_les_data, exit_les_data
    public :: les_is_enabled, update_les_viscosity
    public :: reset_les_profile, add_les_profile, write_les_profile, les_wall_seconds

contains

    real(C_DOUBLE) function les_wall_seconds() result(seconds)
        integer(int64) :: count, rate

        call system_clock(count=count, count_rate=rate)
        seconds = real(count, C_DOUBLE)/real(rate, C_DOUBLE)
    end function les_wall_seconds

    subroutine reset_les_profile(profile)
        type(les_profile_type), intent(inout) :: profile

        profile%seconds = 0.0d0
        profile%calls = 0_C_INT
    end subroutine reset_les_profile

    subroutine add_les_profile(profile, category, elapsed_seconds)
        type(les_profile_type), intent(inout) :: profile
        integer, intent(in) :: category
        real(C_DOUBLE), intent(in) :: elapsed_seconds

        if (category < 1 .or. category > LES_PROF_NCATS) return
        profile%seconds(category) = profile%seconds(category) + elapsed_seconds
        profile%calls(category) = profile%calls(category) + 1_C_INT
    end subroutine add_les_profile

    subroutine write_les_profile(profile, nsteps)
        type(les_profile_type), intent(in) :: profile
        integer(C_INT), intent(in) :: nsteps

        real(C_DOUBLE) :: total_seconds

        total_seconds = sum(profile%seconds)
        call write_les_profile_line("nut_update", profile%seconds(LES_PROF_NUT), &
            profile%calls(LES_PROF_NUT), nsteps)
        call write_les_profile_line("nut_exchange", profile%seconds(LES_PROF_EXCHANGE), &
            profile%calls(LES_PROF_EXCHANGE), nsteps)
        call write_les_profile_line("sgs_correction", profile%seconds(LES_PROF_SGS), &
            profile%calls(LES_PROF_SGS), nsteps)
        call write_les_profile_line("total_measured", total_seconds, sum(profile%calls), nsteps)
    end subroutine write_les_profile

    subroutine write_les_profile_line(label, seconds, calls, nsteps)
        character(len=*), intent(in) :: label
        real(C_DOUBLE), intent(in) :: seconds
        integer(C_INT), intent(in) :: calls, nsteps

        real(C_DOUBLE) :: seconds_per_step, seconds_per_call

        seconds_per_step = 0.0d0
        seconds_per_call = 0.0d0
        if (nsteps > 0_C_INT) seconds_per_step = seconds/real(nsteps, C_DOUBLE)
        if (calls > 0_C_INT) seconds_per_call = seconds/real(calls, C_DOUBLE)

        write(*,'(A,1X,A,1X,A,1X,I0,1X,A,1X,I0,1X,A,1X,ES16.8,1X,A,1X,ES16.8,1X,A,1X,ES16.8)') &
            "les_timing:", trim(label), "calls", calls, "nsteps", nsteps, &
            "seconds", seconds, "seconds_per_step", seconds_per_step, &
            "seconds_per_call", seconds_per_call
    end subroutine write_les_profile_line

    logical function les_is_enabled(les)
        type(les_type), intent(in) :: les

        les_is_enabled = les%model /= LES_NONE
    end function les_is_enabled

    subroutine init_les(les, dns, g)
        type(les_type), intent(inout) :: les
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g

        integer :: nx, ny, nz

        call destroy_les(les)
        if (.not. les_is_enabled(les)) return

        nx = int(dns%localSize(1,2))
        ny = int(dns%localSize(2,2))
        nz = int(dns%localSize(3,2))

        allocate(les%nut(0:nx+1,0:ny+1,0:nz+1))
        allocate(les%filter_x(0:nx+1), les%filter_y(0:ny+1), les%filter_z(0:nz+1))
        allocate(les%d1xm(0:nx+1,VAR_U:VAR_P), les%d1x0(0:nx+1,VAR_U:VAR_P), &
            les%d1xp(0:nx+1,VAR_U:VAR_P))
        allocate(les%d1ym(0:ny+1,VAR_U:VAR_P), les%d1y0(0:ny+1,VAR_U:VAR_P), &
            les%d1yp(0:ny+1,VAR_U:VAR_P))
        allocate(les%d1zm(0:nz+1,VAR_U:VAR_P), les%d1z0(0:nz+1,VAR_U:VAR_P), &
            les%d1zp(0:nz+1,VAR_U:VAR_P))
        allocate(les%p_from_u_x(0:nx+1), les%p_from_v_y(0:ny+1), les%p_from_w_z(0:nz+1))
        allocate(les%u_from_p_x(0:nx+1), les%v_from_p_y(0:ny+1), les%w_from_p_z(0:nz+1))
        allocate(les%inv_dx(0:nx+1,VAR_U:VAR_P), les%inv_dy(0:ny+1,VAR_U:VAR_P), &
            les%inv_dz(0:nz+1,VAR_U:VAR_P))

        les%nut = 0.0d0
        call precompute_les_metrics(les, g, nx, ny, nz)
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

    subroutine precompute_les_metrics(les, g, nx, ny, nz)
        type(les_type), intent(inout) :: les
        type(grid_type), intent(in) :: g
        integer, intent(in) :: nx, ny, nz

        integer :: i, var

        do i = 0, nx+1
            les%filter_x(i) = max(1.0d0/g%d1x(i,VAR_P), 1.0d-30)**(1.0d0/3.0d0)
            les%p_from_u_x(i) = linear_weight(g%x(i,VAR_U), g%x(i+1,VAR_U), g%x(i,VAR_P))
            les%u_from_p_x(i) = linear_weight(g%x(i-1,VAR_P), g%x(i,VAR_P), g%x(i,VAR_U))
            do var = VAR_U, VAR_P
                call first_derivative_coeffs(g%x(i-1,var), g%x(i,var), g%x(i+1,var), &
                    les%d1xm(i,var), les%d1x0(i,var), les%d1xp(i,var))
                les%inv_dx(i,var) = safe_inv_delta(g%x(i,var) - g%x(i-1,var))
            end do
        end do

        do i = 0, ny+1
            les%filter_y(i) = max(1.0d0/g%d1y(i,VAR_P), 1.0d-30)**(1.0d0/3.0d0)
            les%p_from_v_y(i) = linear_weight(g%y(i,VAR_V), g%y(i+1,VAR_V), g%y(i,VAR_P))
            les%v_from_p_y(i) = linear_weight(g%y(i-1,VAR_P), g%y(i,VAR_P), g%y(i,VAR_V))
            do var = VAR_U, VAR_P
                call first_derivative_coeffs(g%y(i-1,var), g%y(i,var), g%y(i+1,var), &
                    les%d1ym(i,var), les%d1y0(i,var), les%d1yp(i,var))
                les%inv_dy(i,var) = safe_inv_delta(g%y(i,var) - g%y(i-1,var))
            end do
        end do

        do i = 0, nz+1
            les%filter_z(i) = max(1.0d0/g%d1z(i,VAR_P), 1.0d-30)**(1.0d0/3.0d0)
            les%p_from_w_z(i) = linear_weight(g%z(i,VAR_W), g%z(i+1,VAR_W), g%z(i,VAR_P))
            les%w_from_p_z(i) = linear_weight(g%z(i-1,VAR_P), g%z(i,VAR_P), g%z(i,VAR_W))
            do var = VAR_U, VAR_P
                call first_derivative_coeffs(g%z(i-1,var), g%z(i,var), g%z(i+1,var), &
                    les%d1zm(i,var), les%d1z0(i,var), les%d1zp(i,var))
                les%inv_dz(i,var) = safe_inv_delta(g%z(i,var) - g%z(i-1,var))
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

    subroutine update_les_viscosity(les, f, dns, g, ibm)
        type(les_type), intent(inout) :: les
        type(field_type), intent(in) :: f
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        type(ibm_type), intent(in) :: ibm

        if (.not. allocated(les%nut)) return

        select case (les%model)
        case (LES_NONE)
            return
        case (LES_WALE)
            call update_wale_viscosity(les, f, dns, g, ibm)
        case default
            call update_generic_les_viscosity(les, f, dns, g, ibm)
        end select
    end subroutine update_les_viscosity

    subroutine update_generic_les_viscosity(les, f, dns, g, ibm)
        type(les_type), intent(inout) :: les
        type(field_type), intent(in) :: f
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        type(ibm_type), intent(in) :: ibm

        integer :: i, j, k, nx, ny, nz
        integer(C_INT) :: model
        logical(C_BOOL) :: ibm_aware, ibm_enabled
        real(C_DOUBLE) :: cs, cw, cs2, cw2, delta_scale
        real(C_DOUBLE) :: delta
        real(C_DOUBLE) :: g11, g12, g13, g21, g22, g23, g31, g32, g33
        real(C_DOUBLE) :: s11, s22, s33, s12, s13, s23, s2, strain_mag
        real(C_DOUBLE) :: g2_11, g2_12, g2_13, g2_21, g2_22, g2_23, g2_31, g2_32, g2_33
        real(C_DOUBLE) :: trace_g2, sd11, sd22, sd33, sd12, sd13, sd23, sd2, denom
        real(C_DOUBLE) :: sqrt_s2, sqrt_sd2, s2_52, sd2_32, sd2_54
        real(C_DOUBLE) :: d0, d1, solid_threshold
        logical :: solid_cell

        if (.not. allocated(les%nut)) return

        nx = int(dns%localSize(1,2))
        ny = int(dns%localSize(2,2))
        nz = int(dns%localSize(3,2))
        model = les%model
        cs = les%cs
        cw = les%cw
        cs2 = cs*cs
        cw2 = cw*cw
        delta_scale = les%delta_scale
        ibm_aware = les%ibm_aware
        ibm_enabled = dns%ibm_enabled
        solid_threshold = 1.0d20

        if (model == LES_NONE) return

        !$omp target teams distribute parallel do collapse(3) &
        !$omp& map(to: model, cs2, cw2, delta_scale, ibm_aware, ibm_enabled, solid_threshold, &
        !$omp& f%q, g%d1x, g%d1y, g%d1z, ibm%coef, &
        !$omp& les%filter_x, les%filter_y, les%filter_z, &
        !$omp& les%d1xm, les%d1x0, les%d1xp, les%d1ym, les%d1y0, les%d1yp, &
        !$omp& les%d1zm, les%d1z0, les%d1zp, &
        !$omp& les%p_from_u_x, les%p_from_v_y, les%p_from_w_z) &
        !$omp& map(tofrom: les%nut) &
        !$omp& private(i,j,k,delta,g11,g12,g13,g21,g22,g23,g31,g32,g33, &
        !$omp& s11,s22,s33,s12,s13,s23,s2,strain_mag,g2_11,g2_12,g2_13, &
        !$omp& g2_21,g2_22,g2_23,g2_31,g2_32,g2_33,trace_g2,sd11,sd22,sd33, &
        !$omp& sd12,sd13,sd23,sd2,denom,sqrt_s2,sqrt_sd2, &
        !$omp& s2_52,sd2_32,sd2_54,d0,d1,solid_cell)
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    les%nut(i,j,k) = 0.0d0

                    solid_cell = .false.
                    if (ibm_aware .and. ibm_enabled) then
                        solid_cell = abs(ibm%coef(i,j,k,VAR_U)) > solid_threshold .or. &
                                     abs(ibm%coef(i+1,j,k,VAR_U)) > solid_threshold .or. &
                                     abs(ibm%coef(i,j,k,VAR_V)) > solid_threshold .or. &
                                     abs(ibm%coef(i,j+1,k,VAR_V)) > solid_threshold .or. &
                                     abs(ibm%coef(i,j,k,VAR_W)) > solid_threshold .or. &
                                     abs(ibm%coef(i,j,k+1,VAR_W)) > solid_threshold
                    end if
                    if (solid_cell) cycle

                    delta = delta_scale*les%filter_x(i)*les%filter_y(j)*les%filter_z(k)

                    g11 = (f%q(i+1,j,k,VAR_U) - f%q(i,j,k,VAR_U))*g%d1x(i,VAR_P)
                    g22 = (f%q(i,j+1,k,VAR_V) - f%q(i,j,k,VAR_V))*g%d1y(j,VAR_P)
                    g33 = (f%q(i,j,k+1,VAR_W) - f%q(i,j,k,VAR_W))*g%d1z(k,VAR_P)

                    d0 = les%d1ym(j,VAR_U)*f%q(i,j-1,k,VAR_U) &
                       + les%d1y0(j,VAR_U)*f%q(i,j,k,VAR_U) &
                       + les%d1yp(j,VAR_U)*f%q(i,j+1,k,VAR_U)
                    d1 = les%d1ym(j,VAR_U)*f%q(i+1,j-1,k,VAR_U) &
                       + les%d1y0(j,VAR_U)*f%q(i+1,j,k,VAR_U) &
                       + les%d1yp(j,VAR_U)*f%q(i+1,j+1,k,VAR_U)
                    g12 = (1.0d0 - les%p_from_u_x(i))*d0 + les%p_from_u_x(i)*d1

                    d0 = les%d1zm(k,VAR_U)*f%q(i,j,k-1,VAR_U) &
                       + les%d1z0(k,VAR_U)*f%q(i,j,k,VAR_U) &
                       + les%d1zp(k,VAR_U)*f%q(i,j,k+1,VAR_U)
                    d1 = les%d1zm(k,VAR_U)*f%q(i+1,j,k-1,VAR_U) &
                       + les%d1z0(k,VAR_U)*f%q(i+1,j,k,VAR_U) &
                       + les%d1zp(k,VAR_U)*f%q(i+1,j,k+1,VAR_U)
                    g13 = (1.0d0 - les%p_from_u_x(i))*d0 + les%p_from_u_x(i)*d1

                    d0 = les%d1xm(i,VAR_V)*f%q(i-1,j,k,VAR_V) &
                       + les%d1x0(i,VAR_V)*f%q(i,j,k,VAR_V) &
                       + les%d1xp(i,VAR_V)*f%q(i+1,j,k,VAR_V)
                    d1 = les%d1xm(i,VAR_V)*f%q(i-1,j+1,k,VAR_V) &
                       + les%d1x0(i,VAR_V)*f%q(i,j+1,k,VAR_V) &
                       + les%d1xp(i,VAR_V)*f%q(i+1,j+1,k,VAR_V)
                    g21 = (1.0d0 - les%p_from_v_y(j))*d0 + les%p_from_v_y(j)*d1

                    d0 = les%d1zm(k,VAR_V)*f%q(i,j,k-1,VAR_V) &
                       + les%d1z0(k,VAR_V)*f%q(i,j,k,VAR_V) &
                       + les%d1zp(k,VAR_V)*f%q(i,j,k+1,VAR_V)
                    d1 = les%d1zm(k,VAR_V)*f%q(i,j+1,k-1,VAR_V) &
                       + les%d1z0(k,VAR_V)*f%q(i,j+1,k,VAR_V) &
                       + les%d1zp(k,VAR_V)*f%q(i,j+1,k+1,VAR_V)
                    g23 = (1.0d0 - les%p_from_v_y(j))*d0 + les%p_from_v_y(j)*d1

                    d0 = les%d1xm(i,VAR_W)*f%q(i-1,j,k,VAR_W) &
                       + les%d1x0(i,VAR_W)*f%q(i,j,k,VAR_W) &
                       + les%d1xp(i,VAR_W)*f%q(i+1,j,k,VAR_W)
                    d1 = les%d1xm(i,VAR_W)*f%q(i-1,j,k+1,VAR_W) &
                       + les%d1x0(i,VAR_W)*f%q(i,j,k+1,VAR_W) &
                       + les%d1xp(i,VAR_W)*f%q(i+1,j,k+1,VAR_W)
                    g31 = (1.0d0 - les%p_from_w_z(k))*d0 + les%p_from_w_z(k)*d1

                    d0 = les%d1ym(j,VAR_W)*f%q(i,j-1,k,VAR_W) &
                       + les%d1y0(j,VAR_W)*f%q(i,j,k,VAR_W) &
                       + les%d1yp(j,VAR_W)*f%q(i,j+1,k,VAR_W)
                    d1 = les%d1ym(j,VAR_W)*f%q(i,j-1,k+1,VAR_W) &
                       + les%d1y0(j,VAR_W)*f%q(i,j,k+1,VAR_W) &
                       + les%d1yp(j,VAR_W)*f%q(i,j+1,k+1,VAR_W)
                    g32 = (1.0d0 - les%p_from_w_z(k))*d0 + les%p_from_w_z(k)*d1

                    s11 = g11
                    s22 = g22
                    s33 = g33
                    s12 = 0.5d0*(g12 + g21)
                    s13 = 0.5d0*(g13 + g31)
                    s23 = 0.5d0*(g23 + g32)
                    s2 = s11*s11 + s22*s22 + s33*s33 + 2.0d0*(s12*s12 + s13*s13 + s23*s23)

                    select case (model)
                    case (LES_SMAGORINSKY)
                        strain_mag = sqrt(max(0.0d0, 2.0d0*s2))
                        les%nut(i,j,k) = cs2*delta*delta*strain_mag
                    case (LES_WALE)
                        g2_11 = g11*g11 + g12*g21 + g13*g31
                        g2_12 = g11*g12 + g12*g22 + g13*g32
                        g2_13 = g11*g13 + g12*g23 + g13*g33
                        g2_21 = g21*g11 + g22*g21 + g23*g31
                        g2_22 = g21*g12 + g22*g22 + g23*g32
                        g2_23 = g21*g13 + g22*g23 + g23*g33
                        g2_31 = g31*g11 + g32*g21 + g33*g31
                        g2_32 = g31*g12 + g32*g22 + g33*g32
                        g2_33 = g31*g13 + g32*g23 + g33*g33

                        trace_g2 = g2_11 + g2_22 + g2_33
                        sd11 = g2_11 - trace_g2/3.0d0
                        sd22 = g2_22 - trace_g2/3.0d0
                        sd33 = g2_33 - trace_g2/3.0d0
                        sd12 = 0.5d0*(g2_12 + g2_21)
                        sd13 = 0.5d0*(g2_13 + g2_31)
                        sd23 = 0.5d0*(g2_23 + g2_32)
                        sd2 = sd11*sd11 + sd22*sd22 + sd33*sd33 + &
                              2.0d0*(sd12*sd12 + sd13*sd13 + sd23*sd23)

                        sqrt_s2 = sqrt(s2)
                        sqrt_sd2 = sqrt(sd2)
                        s2_52 = s2*s2*sqrt_s2
                        sd2_32 = sd2*sqrt_sd2
                        sd2_54 = sd2*sqrt(sqrt_sd2)
                        denom = s2_52 + sd2_54 + 1.0d-30
                        les%nut(i,j,k) = cw2*delta*delta*sd2_32/denom
                    end select
                end do
            end do
        end do
        !$omp end target teams distribute parallel do
    end subroutine update_generic_les_viscosity

    subroutine update_wale_viscosity(les, f, dns, g, ibm)
        type(les_type), intent(inout) :: les
        type(field_type), intent(in) :: f
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        type(ibm_type), intent(in) :: ibm

        integer :: i, j, k, nx, ny, nz
        logical(C_BOOL) :: ibm_aware, ibm_enabled
        real(C_DOUBLE) :: cw2, delta_scale, delta
        real(C_DOUBLE) :: g11, g12, g13, g21, g22, g23, g31, g32, g33
        real(C_DOUBLE) :: s2, g2_11, g2_22, g2_33, trace_g2
        real(C_DOUBLE) :: sd11, sd22, sd33, sd12, sd13, sd23, sd2, denom
        real(C_DOUBLE) :: sqrt_s2, sqrt_sd2, s2_52, sd2_32, sd2_54
        real(C_DOUBLE) :: d0, d1, solid_threshold
        logical :: solid_cell

        nx = int(dns%localSize(1,2))
        ny = int(dns%localSize(2,2))
        nz = int(dns%localSize(3,2))
        cw2 = les%cw*les%cw
        delta_scale = les%delta_scale
        ibm_aware = les%ibm_aware
        ibm_enabled = dns%ibm_enabled
        solid_threshold = 1.0d20

        !$omp target teams distribute parallel do collapse(3) &
        !$omp& map(to: cw2, delta_scale, ibm_aware, ibm_enabled, solid_threshold, &
        !$omp& f%q, g%d1x, g%d1y, g%d1z, ibm%coef, &
        !$omp& les%filter_x, les%filter_y, les%filter_z, &
        !$omp& les%d1xm, les%d1x0, les%d1xp, les%d1ym, les%d1y0, les%d1yp, &
        !$omp& les%d1zm, les%d1z0, les%d1zp, &
        !$omp& les%p_from_u_x, les%p_from_v_y, les%p_from_w_z) &
        !$omp& map(tofrom: les%nut) &
        !$omp& private(i,j,k,delta,g11,g12,g13,g21,g22,g23,g31,g32,g33, &
        !$omp& s2,g2_11,g2_22,g2_33,trace_g2,sd11,sd22,sd33,sd12,sd13,sd23, &
        !$omp& sd2,denom,sqrt_s2,sqrt_sd2,s2_52,sd2_32,sd2_54,d0,d1,solid_cell)
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    les%nut(i,j,k) = 0.0d0

                    solid_cell = .false.
                    if (ibm_aware .and. ibm_enabled) then
                        solid_cell = abs(ibm%coef(i,j,k,VAR_U)) > solid_threshold .or. &
                                     abs(ibm%coef(i+1,j,k,VAR_U)) > solid_threshold .or. &
                                     abs(ibm%coef(i,j,k,VAR_V)) > solid_threshold .or. &
                                     abs(ibm%coef(i,j+1,k,VAR_V)) > solid_threshold .or. &
                                     abs(ibm%coef(i,j,k,VAR_W)) > solid_threshold .or. &
                                     abs(ibm%coef(i,j,k+1,VAR_W)) > solid_threshold
                    end if
                    if (solid_cell) cycle

                    delta = delta_scale*les%filter_x(i)*les%filter_y(j)*les%filter_z(k)

                    g11 = (f%q(i+1,j,k,VAR_U) - f%q(i,j,k,VAR_U))*g%d1x(i,VAR_P)
                    g22 = (f%q(i,j+1,k,VAR_V) - f%q(i,j,k,VAR_V))*g%d1y(j,VAR_P)
                    g33 = (f%q(i,j,k+1,VAR_W) - f%q(i,j,k,VAR_W))*g%d1z(k,VAR_P)

                    d0 = les%d1ym(j,VAR_U)*f%q(i,j-1,k,VAR_U) &
                       + les%d1y0(j,VAR_U)*f%q(i,j,k,VAR_U) &
                       + les%d1yp(j,VAR_U)*f%q(i,j+1,k,VAR_U)
                    d1 = les%d1ym(j,VAR_U)*f%q(i+1,j-1,k,VAR_U) &
                       + les%d1y0(j,VAR_U)*f%q(i+1,j,k,VAR_U) &
                       + les%d1yp(j,VAR_U)*f%q(i+1,j+1,k,VAR_U)
                    g12 = (1.0d0 - les%p_from_u_x(i))*d0 + les%p_from_u_x(i)*d1

                    d0 = les%d1zm(k,VAR_U)*f%q(i,j,k-1,VAR_U) &
                       + les%d1z0(k,VAR_U)*f%q(i,j,k,VAR_U) &
                       + les%d1zp(k,VAR_U)*f%q(i,j,k+1,VAR_U)
                    d1 = les%d1zm(k,VAR_U)*f%q(i+1,j,k-1,VAR_U) &
                       + les%d1z0(k,VAR_U)*f%q(i+1,j,k,VAR_U) &
                       + les%d1zp(k,VAR_U)*f%q(i+1,j,k+1,VAR_U)
                    g13 = (1.0d0 - les%p_from_u_x(i))*d0 + les%p_from_u_x(i)*d1

                    d0 = les%d1xm(i,VAR_V)*f%q(i-1,j,k,VAR_V) &
                       + les%d1x0(i,VAR_V)*f%q(i,j,k,VAR_V) &
                       + les%d1xp(i,VAR_V)*f%q(i+1,j,k,VAR_V)
                    d1 = les%d1xm(i,VAR_V)*f%q(i-1,j+1,k,VAR_V) &
                       + les%d1x0(i,VAR_V)*f%q(i,j+1,k,VAR_V) &
                       + les%d1xp(i,VAR_V)*f%q(i+1,j+1,k,VAR_V)
                    g21 = (1.0d0 - les%p_from_v_y(j))*d0 + les%p_from_v_y(j)*d1

                    d0 = les%d1zm(k,VAR_V)*f%q(i,j,k-1,VAR_V) &
                       + les%d1z0(k,VAR_V)*f%q(i,j,k,VAR_V) &
                       + les%d1zp(k,VAR_V)*f%q(i,j,k+1,VAR_V)
                    d1 = les%d1zm(k,VAR_V)*f%q(i,j+1,k-1,VAR_V) &
                       + les%d1z0(k,VAR_V)*f%q(i,j+1,k,VAR_V) &
                       + les%d1zp(k,VAR_V)*f%q(i,j+1,k+1,VAR_V)
                    g23 = (1.0d0 - les%p_from_v_y(j))*d0 + les%p_from_v_y(j)*d1

                    d0 = les%d1xm(i,VAR_W)*f%q(i-1,j,k,VAR_W) &
                       + les%d1x0(i,VAR_W)*f%q(i,j,k,VAR_W) &
                       + les%d1xp(i,VAR_W)*f%q(i+1,j,k,VAR_W)
                    d1 = les%d1xm(i,VAR_W)*f%q(i-1,j,k+1,VAR_W) &
                       + les%d1x0(i,VAR_W)*f%q(i,j,k+1,VAR_W) &
                       + les%d1xp(i,VAR_W)*f%q(i+1,j,k+1,VAR_W)
                    g31 = (1.0d0 - les%p_from_w_z(k))*d0 + les%p_from_w_z(k)*d1

                    d0 = les%d1ym(j,VAR_W)*f%q(i,j-1,k,VAR_W) &
                       + les%d1y0(j,VAR_W)*f%q(i,j,k,VAR_W) &
                       + les%d1yp(j,VAR_W)*f%q(i,j+1,k,VAR_W)
                    d1 = les%d1ym(j,VAR_W)*f%q(i,j-1,k+1,VAR_W) &
                       + les%d1y0(j,VAR_W)*f%q(i,j,k+1,VAR_W) &
                       + les%d1yp(j,VAR_W)*f%q(i,j+1,k+1,VAR_W)
                    g32 = (1.0d0 - les%p_from_w_z(k))*d0 + les%p_from_w_z(k)*d1

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
                    les%nut(i,j,k) = cw2*delta*delta*sd2_32/denom
                end do
            end do
        end do
        !$omp end target teams distribute parallel do
    end subroutine update_wale_viscosity

end module les_model
