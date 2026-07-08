module les_model
    use, intrinsic :: iso_c_binding
    use :: init, only: dns_type, VAR_U, VAR_V, VAR_W, VAR_P
    use :: blocks, only: block_set_type
    use :: ibmm, only: ibm_type
    use :: turbulence, only: turb_type
    implicit none

    private

    integer(C_INT), parameter, public :: LES_NONE = 0_C_INT
    integer(C_INT), parameter, public :: LES_SMAGORINSKY = 1_C_INT
    integer(C_INT), parameter, public :: LES_WALE = 2_C_INT

    ! Algebraic SGS eddy-viscosity models. The state (nut, the grid-metric
    ! tables) lives in turb_type (turbulence.f90); les_type carries only the
    ! kernel constants, and the kernels write into a caller-supplied nut
    ! target so future hybrid models can direct the SGS viscosity into a
    ! scratch array.
    type, public :: les_type
        integer(C_INT) :: model = LES_NONE
        real(C_DOUBLE) :: cs = 0.10d0
        real(C_DOUBLE) :: cw = 0.325d0
        real(C_DOUBLE) :: delta_scale = 1.0d0
        logical(C_BOOL) :: ibm_aware = .true.
    end type les_type

    public :: update_sgs_viscosity
    ! Shared with the RANS producer (SST production/blending need the same
    ! cell-centred gradient tensor).
    public :: velocity_gradient_tensor

contains

    subroutine update_sgs_viscosity(les, turb, blk, dns, ibm, nut)
        type(les_type), intent(in) :: les
        type(turb_type), intent(in) :: turb
        type(block_set_type), intent(inout) :: blk
        type(dns_type), intent(in) :: dns
        type(ibm_type), intent(in) :: ibm
        real(C_DOUBLE), intent(inout) :: nut(0:,0:,0:,1:)

        select case (les%model)
        case (LES_NONE)
            return
        case (LES_WALE)
            call update_wale_viscosity(les, turb, blk, dns, ibm, nut)
        case default
            call update_generic_les_viscosity(les, turb, blk, dns, ibm, nut)
        end select
    end subroutine update_sgs_viscosity

    ! Velocity-gradient tensor g(a,b) = du_a/dx_b at the centre of cell (i,j,k)
    ! in block b, on the staggered mesh: the diagonal terms are the face
    ! difference across the cell; each off-diagonal term interpolates the
    ! neighbouring face-difference to the cell centre with the p_from_* staggered
    ! weights and the turb d1?? stencils. Shared verbatim by the Smagorinsky and
    ! WALE kernels; declared target so both device kernels can call it.
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

    subroutine update_generic_les_viscosity(les, turb, blk, dns, ibm, nut)
        type(les_type), intent(in) :: les
        type(turb_type), intent(in) :: turb
        type(block_set_type), intent(inout) :: blk
        type(dns_type), intent(in) :: dns
        type(ibm_type), intent(in) :: ibm
        real(C_DOUBLE), intent(inout) :: nut(0:,0:,0:,1:)

        integer :: i, j, k, b, nx, ny, nz, nBlocks
        integer(C_INT) :: model
        logical(C_BOOL) :: ibm_aware, ibm_enabled
        real(C_DOUBLE) :: cs, cs2, delta_scale
        real(C_DOUBLE) :: delta
        real(C_DOUBLE) :: g11, g12, g13, g21, g22, g23, g31, g32, g33
        real(C_DOUBLE) :: s11, s22, s33, s12, s13, s23, s2, strain_mag
        real(C_DOUBLE) :: solid_threshold
        logical :: solid_cell

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
        !$omp& turb%filter_x, turb%filter_y, turb%filter_z, &
        !$omp& turb%d1xm, turb%d1x0, turb%d1xp, turb%d1ym, turb%d1y0, turb%d1yp, &
        !$omp& turb%d1zm, turb%d1z0, turb%d1zp, &
        !$omp& turb%p_from_u_x, turb%p_from_v_y, turb%p_from_w_z) &
        !$omp& map(tofrom: nut) &
        !$omp& private(i,j,k,b,delta,g11,g12,g13,g21,g22,g23,g31,g32,g33, &
        !$omp& s11,s22,s33,s12,s13,s23,s2,strain_mag,solid_cell)
        do b = 1, nBlocks
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    nut(i,j,k,b) = 0.0d0

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

                    delta = delta_scale*turb%filter_x(i,b)*turb%filter_y(j,b)*turb%filter_z(k,b)

                    call velocity_gradient_tensor(blk, turb, i, j, k, b, &
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
                        nut(i,j,k,b) = cs2*delta*delta*strain_mag
                    end select
                end do
            end do
        end do
        end do
        !$omp end target teams distribute parallel do
    end subroutine update_generic_les_viscosity

    subroutine update_wale_viscosity(les, turb, blk, dns, ibm, nut)
        type(les_type), intent(in) :: les
        type(turb_type), intent(in) :: turb
        type(block_set_type), intent(inout) :: blk
        type(dns_type), intent(in) :: dns
        type(ibm_type), intent(in) :: ibm
        real(C_DOUBLE), intent(inout) :: nut(0:,0:,0:,1:)

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
        !$omp& turb%filter_x, turb%filter_y, turb%filter_z, &
        !$omp& turb%d1xm, turb%d1x0, turb%d1xp, turb%d1ym, turb%d1y0, turb%d1yp, &
        !$omp& turb%d1zm, turb%d1z0, turb%d1zp, &
        !$omp& turb%p_from_u_x, turb%p_from_v_y, turb%p_from_w_z) &
        !$omp& map(tofrom: nut) &
        !$omp& private(i,j,k,b,delta,g11,g12,g13,g21,g22,g23,g31,g32,g33, &
        !$omp& s2,g2_11,g2_22,g2_33,trace_g2,sd11,sd22,sd33,sd12,sd13,sd23, &
        !$omp& sd2,denom,sqrt_s2,sqrt_sd2,s2_52,sd2_32,sd2_54,solid_cell)
        do b = 1, nBlocks
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    nut(i,j,k,b) = 0.0d0

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

                    delta = delta_scale*turb%filter_x(i,b)*turb%filter_y(j,b)*turb%filter_z(k,b)

                    call velocity_gradient_tensor(blk, turb, i, j, k, b, &
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
                    nut(i,j,k,b) = cw2*delta*delta*sd2_32/denom
                end do
            end do
        end do
        end do
        !$omp end target teams distribute parallel do
    end subroutine update_wale_viscosity

end module les_model
