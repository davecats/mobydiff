module turbulence
    use, intrinsic :: iso_c_binding
    use :: chron, only: profiler_type, init_profiler
    use :: init, only: dns_type, VAR_U, VAR_V, VAR_W, VAR_P
    use :: blocks, only: block_set_type
    implicit none

    private

    ! Turbulence model families (docs/next_session_iddes.md). TURB_IDDES is
    ! the SST-based hybrid: RANS transport everywhere, the DDES shielding
    ! function f_d switching the k-destruction length and the nut assembly
    ! between the RANS and the SGS producer.
    integer(C_INT), parameter, public :: TURB_NONE = 0_C_INT
    integer(C_INT), parameter, public :: TURB_LES = 1_C_INT
    integer(C_INT), parameter, public :: TURB_RANS = 2_C_INT
    integer(C_INT), parameter, public :: TURB_IDDES = 3_C_INT
    ! Category indices into the turbulence timing profiler (see init_turbulence_profiler).
    integer, parameter, public :: TURB_PROF_NUT = 1
    integer, parameter, public :: TURB_PROF_EXCHANGE = 2
    integer, parameter, public :: TURB_PROF_SGS = 3

    ! IDDES (DDES-shielding form; SST-IDDES calibration, Gritskevich et al.
    ! 2012): l_LES = C_DES Delta with C_DES the F1-blend of set 1/2; the
    ! stored fd = tanh((8 r_d)^3) is the RANS-RETENTION weight (see
    ! compute_iddes_fd) on r_d = (nut + nu)/(kappa^2 y_eff^2 |grad u|).
    real(C_DOUBLE), parameter, public :: IDDES_CDES1 = 0.78d0
    real(C_DOUBLE), parameter, public :: IDDES_CDES2 = 0.61d0
    real(C_DOUBLE), parameter :: IDDES_KAPPA = 0.41d0

    ! Model-agnostic turbulence state. Every producer (LES, RANS, the IDDES
    ! blend) fills the same cell-centred eddy viscosity nut, so everything
    ! downstream of nut (halo exchange, momentum correction, dt limit, io)
    ! never needs to know which model ran. The 1D tables are grid metrics
    ! shared by the producers' stencils, not LES physics.
    type, public :: turb_type
        integer(C_INT) :: model = TURB_NONE
        real(C_DOUBLE), allocatable :: nut(:,:,:,:)   ! (0:nb+1,...,nBlocks)
        ! IDDES state: the SGS producer writes nut_sgs, compute_iddes_fd the
        ! shielding function fd; blend_iddes_nut combines them with the RANS
        ! nut. Full arrays only under model = iddes -- 1-cell dummies
        ! otherwise (uniform device maps; every access is model-guarded).
        real(C_DOUBLE), allocatable :: nut_sgs(:,:,:,:), fd(:,:,:,:)
        ! [turbulence] fd_force: validation hook forcing fd to a constant
        ! (0 = pure-SGS limit, 1 = pure-RANS limit); < 0 = off (production).
        real(C_DOUBLE) :: fd_force = -1.0d0
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
    public :: velocity_gradient_tensor
    public :: compute_iddes_fd, blend_iddes_nut, iddes_k_sink_coeff

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
        if (turb%model == TURB_IDDES) then
            allocate(turb%nut_sgs(0:nx+1,0:ny+1,0:nz+1,blk%nBlocks))
            allocate(turb%fd(0:nx+1,0:ny+1,0:nz+1,blk%nBlocks))
        else
            ! 1-cell dummies: uniform device maps, all accesses model-guarded.
            allocate(turb%nut_sgs(0:0,0:0,0:0,1), turb%fd(0:0,0:0,0:0,1))
        end if
        turb%nut_sgs = 0.0d0
        turb%fd = 1.0d0   ! the RANS limit, until the first compute_iddes_fd
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
        if (allocated(turb%nut_sgs)) deallocate(turb%nut_sgs)
        if (allocated(turb%fd)) deallocate(turb%fd)
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
        !$omp target enter data map(to: turb%nut_sgs, turb%fd)
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
        !$omp target exit data map(delete: turb%nut_sgs, turb%fd)
        !$omp target exit data map(delete: turb%nut)
        !$omp target exit data map(delete: turb)
    end subroutine exit_turbulence_data

    ! Velocity-gradient tensor g(a,b) = du_a/dx_b at the centre of cell (i,j,k)
    ! in block b, on the staggered mesh: the diagonal terms are the face
    ! difference across the cell; each off-diagonal term interpolates the
    ! neighbouring face-difference to the cell centre with the p_from_* staggered
    ! weights and the turb d1?? stencils. Shared verbatim by the SGS kernels
    ! (les.f90), the SST producer (rans.f90) and the IDDES shielding function
    ! below -- it lives HERE because turbulence.f90 sits below both producers
    ! in the module graph; declared target so every device kernel can call it.
    subroutine velocity_gradient_tensor(blk, turb, i, j, k, b, &
            g11, g12, g13, g21, g22, g23, g31, g32, g33)
!$omp declare target
        type(block_set_type), intent(in) :: blk
        type(turb_type), intent(in) :: turb
        integer, intent(in) :: i, j, k, b
        real(C_DOUBLE), intent(out) :: g11, g12, g13, g21, g22, g23, g31, g32, g33

        real(C_DOUBLE) :: d0, d1

        g11 = (blk%q(i+1,j,k,VAR_U,b) - blk%q(i,j,k,VAR_U,b))*blk%d1x(i,VAR_P,b)
        g22 = (blk%q(i,j+1,k,VAR_V,b) - blk%q(i,j,k,VAR_V,b))*blk%d1y(j,VAR_P,b)
        g33 = (blk%q(i,j,k+1,VAR_W,b) - blk%q(i,j,k,VAR_W,b))*blk%d1z(k,VAR_P,b)

        d0 = turb%d1ym(j,VAR_U,b)*blk%q(i,j-1,k,VAR_U,b) &
           + turb%d1y0(j,VAR_U,b)*blk%q(i,j,k,VAR_U,b) &
           + turb%d1yp(j,VAR_U,b)*blk%q(i,j+1,k,VAR_U,b)
        d1 = turb%d1ym(j,VAR_U,b)*blk%q(i+1,j-1,k,VAR_U,b) &
           + turb%d1y0(j,VAR_U,b)*blk%q(i+1,j,k,VAR_U,b) &
           + turb%d1yp(j,VAR_U,b)*blk%q(i+1,j+1,k,VAR_U,b)
        g12 = (1.0d0 - turb%p_from_u_x(i,b))*d0 + turb%p_from_u_x(i,b)*d1

        d0 = turb%d1zm(k,VAR_U,b)*blk%q(i,j,k-1,VAR_U,b) &
           + turb%d1z0(k,VAR_U,b)*blk%q(i,j,k,VAR_U,b) &
           + turb%d1zp(k,VAR_U,b)*blk%q(i,j,k+1,VAR_U,b)
        d1 = turb%d1zm(k,VAR_U,b)*blk%q(i+1,j,k-1,VAR_U,b) &
           + turb%d1z0(k,VAR_U,b)*blk%q(i+1,j,k,VAR_U,b) &
           + turb%d1zp(k,VAR_U,b)*blk%q(i+1,j,k+1,VAR_U,b)
        g13 = (1.0d0 - turb%p_from_u_x(i,b))*d0 + turb%p_from_u_x(i,b)*d1

        d0 = turb%d1xm(i,VAR_V,b)*blk%q(i-1,j,k,VAR_V,b) &
           + turb%d1x0(i,VAR_V,b)*blk%q(i,j,k,VAR_V,b) &
           + turb%d1xp(i,VAR_V,b)*blk%q(i+1,j,k,VAR_V,b)
        d1 = turb%d1xm(i,VAR_V,b)*blk%q(i-1,j+1,k,VAR_V,b) &
           + turb%d1x0(i,VAR_V,b)*blk%q(i,j+1,k,VAR_V,b) &
           + turb%d1xp(i,VAR_V,b)*blk%q(i+1,j+1,k,VAR_V,b)
        g21 = (1.0d0 - turb%p_from_v_y(j,b))*d0 + turb%p_from_v_y(j,b)*d1

        d0 = turb%d1zm(k,VAR_V,b)*blk%q(i,j,k-1,VAR_V,b) &
           + turb%d1z0(k,VAR_V,b)*blk%q(i,j,k,VAR_V,b) &
           + turb%d1zp(k,VAR_V,b)*blk%q(i,j,k+1,VAR_V,b)
        d1 = turb%d1zm(k,VAR_V,b)*blk%q(i,j+1,k-1,VAR_V,b) &
           + turb%d1z0(k,VAR_V,b)*blk%q(i,j+1,k,VAR_V,b) &
           + turb%d1zp(k,VAR_V,b)*blk%q(i,j+1,k+1,VAR_V,b)
        g23 = (1.0d0 - turb%p_from_v_y(j,b))*d0 + turb%p_from_v_y(j,b)*d1

        d0 = turb%d1xm(i,VAR_W,b)*blk%q(i-1,j,k,VAR_W,b) &
           + turb%d1x0(i,VAR_W,b)*blk%q(i,j,k,VAR_W,b) &
           + turb%d1xp(i,VAR_W,b)*blk%q(i+1,j,k,VAR_W,b)
        d1 = turb%d1xm(i,VAR_W,b)*blk%q(i-1,j,k+1,VAR_W,b) &
           + turb%d1x0(i,VAR_W,b)*blk%q(i,j,k+1,VAR_W,b) &
           + turb%d1xp(i,VAR_W,b)*blk%q(i+1,j,k+1,VAR_W,b)
        g31 = (1.0d0 - turb%p_from_w_z(k,b))*d0 + turb%p_from_w_z(k,b)*d1

        d0 = turb%d1ym(j,VAR_W,b)*blk%q(i,j-1,k,VAR_W,b) &
           + turb%d1y0(j,VAR_W,b)*blk%q(i,j,k,VAR_W,b) &
           + turb%d1yp(j,VAR_W,b)*blk%q(i,j+1,k,VAR_W,b)
        d1 = turb%d1ym(j,VAR_W,b)*blk%q(i,j-1,k+1,VAR_W,b) &
           + turb%d1y0(j,VAR_W,b)*blk%q(i,j,k+1,VAR_W,b) &
           + turb%d1yp(j,VAR_W,b)*blk%q(i,j+1,k+1,VAR_W,b)
        g32 = (1.0d0 - turb%p_from_w_z(k,b))*d0 + turb%p_from_w_z(k,b)*d1
    end subroutine velocity_gradient_tensor

    ! IDDES piece (1): the DDES shielding field. We STORE the RANS-retention
    ! weight fd = tanh((8 r_d)^3) = 1 - f_d^Spalart(2006), so fd -> 1 deep in
    ! the (attached) boundary layer and -> 0 in the LES region -- the weight
    ! the blend and the l_hyb formula multiply nut_rans/l_RANS by directly
    ! (storing Spalart's LES-ward f_d with the same blend hands the WALL
    ! layer to the SGS model: measured +16% log-layer error before the flip).
    ! r_d = (nut + nu)/(kappa^2 y_eff^2 sqrt(sum g_ij g_ij)) reads the lagged
    ! BLENDED nut (the previous assembly), the explicit-scheme stance the
    ! RANS diffusivities already take. y_eff is passed as a plain array so
    ! this module stays independent of the RANS module that owns the
    ! wall-distance state. [turbulence] fd_force >= 0 overrides the field
    ! (validation hook: 0 = pure-SGS limit, 1 = pure-RANS limit).
    subroutine compute_iddes_fd(turb, blk, dns, yeff)
        type(turb_type), intent(inout) :: turb
        type(block_set_type), intent(inout) :: blk
        type(dns_type), intent(in) :: dns
        real(C_DOUBLE), intent(in) :: yeff(0:,0:,0:,1:)

        integer :: i, j, k, b, nx, ny, nz, nBlocks
        real(C_DOUBLE) :: nu, force, gg, rd, y
        real(C_DOUBLE) :: g11, g12, g13, g21, g22, g23, g31, g32, g33

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))
        nBlocks = int(blk%nBlocks)
        nu = 1.0d0/dns%re
        force = turb%fd_force

        !$omp target teams distribute parallel do collapse(4) &
        !$omp& map(to: nu, force, blk%q, blk%d1x, blk%d1y, blk%d1z, yeff, turb%nut, &
        !$omp& turb%d1xm, turb%d1x0, turb%d1xp, turb%d1ym, turb%d1y0, turb%d1yp, &
        !$omp& turb%d1zm, turb%d1z0, turb%d1zp, &
        !$omp& turb%p_from_u_x, turb%p_from_v_y, turb%p_from_w_z) &
        !$omp& map(tofrom: turb%fd) &
        !$omp& private(i,j,k,b,gg,rd,y,g11,g12,g13,g21,g22,g23,g31,g32,g33)
        do b = 1, nBlocks
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    if (force >= 0.0d0) then
                        turb%fd(i,j,k,b) = force
                        cycle
                    end if
                    call velocity_gradient_tensor(blk, turb, i, j, k, b, &
                        g11, g12, g13, g21, g22, g23, g31, g32, g33)
                    gg = sqrt(g11*g11 + g12*g12 + g13*g13 &
                            + g21*g21 + g22*g22 + g23*g23 &
                            + g31*g31 + g32*g32 + g33*g33)
                    y = yeff(i,j,k,b)
                    ! Stagnant cells (gg -> 0, e.g. inside the IBM solid) get
                    ! a huge r_d -> fd = 1 (RANS retention); both nut factors
                    ! are zero there, so the blend is unaffected. The cap
                    ! keeps (8 r_d)^3 finite.
                    rd = min((turb%nut(i,j,k,b) + nu) &
                        /max(IDDES_KAPPA*IDDES_KAPPA*y*y*gg, 1.0d-30), 1.0d10)
                    turb%fd(i,j,k,b) = tanh((8.0d0*rd)**3)
                end do
            end do
        end do
        end do
        !$omp end target teams distribute parallel do
    end subroutine compute_iddes_fd

    ! IDDES piece (2): the final eddy-viscosity assembly
    ! nut = f_d nut_rans + (1 - f_d) nut_sgs. Runs AFTER the RANS assembly
    ! wrote nut_rans into turb%nut and the SGS producer wrote nut_sgs. In
    ! solid cells both producers are zero, so the blend stays zero for any
    ! f_d; at IBM wall cells f_d -> 1 and WALE -> 0 keep the (1-f_d) leak
    ! negligible.
    subroutine blend_iddes_nut(turb, blk)
        type(turb_type), intent(inout) :: turb
        type(block_set_type), intent(in) :: blk

        integer :: i, j, k, b, nx, ny, nz, nBlocks

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))
        nBlocks = int(blk%nBlocks)

        !$omp target teams distribute parallel do collapse(4) &
        !$omp& map(to: turb%fd, turb%nut_sgs) &
        !$omp& map(tofrom: turb%nut) &
        !$omp& private(i,j,k,b)
        do b = 1, nBlocks
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    turb%nut(i,j,k,b) = turb%fd(i,j,k,b)*turb%nut(i,j,k,b) &
                        + (1.0d0 - turb%fd(i,j,k,b))*turb%nut_sgs(i,j,k,b)
                end do
            end do
        end do
        end do
        !$omp end target teams distribute parallel do
    end subroutine blend_iddes_nut

    ! IDDES piece (1b): the POINT-IMPLICIT k-destruction coefficient
    ! D_k/k = sqrt(k)/l_hyb with l_hyb = f_d l_RANS + (1 - f_d) C_DES Delta
    ! and l_RANS = sqrt(k)/(beta* omega) -- the hybrid length replaces
    ! l_RANS in D_k ONLY. Keeping the sink point-implicit is load-bearing:
    ! an explicit k^{3/2}/l_hyb sink reintroduces exactly the near-wall
    ! stiffness the Patankar treatment removed. With f_d = 1 this equals
    ! beta* omega up to ROUND-OFF only, so the pure-RANS caller must keep
    ! its original arithmetic on a separate branch for bit-exactness.
    pure real(C_DOUBLE) function iddes_k_sink_coeff(kv, wv, fdv, cdes_delta, &
            betastar) result(coef)
!$omp declare target
        real(C_DOUBLE), intent(in) :: kv, wv, fdv, cdes_delta, betastar

        real(C_DOUBLE) :: sqrtk, lhyb

        sqrtk = sqrt(max(kv, 0.0d0))
        lhyb = fdv*sqrtk/(betastar*wv) + (1.0d0 - fdv)*cdes_delta
        coef = sqrtk/max(lhyb, 1.0d-30)
    end function iddes_k_sink_coeff

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
