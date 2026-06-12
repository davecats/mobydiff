module pressure_solver
    use, intrinsic :: iso_c_binding
    use :: init, only: dns_type, VAR_U, VAR_V, VAR_W, VAR_P
    use :: blocks, only: block_set_type, FACE_OPEN, FACE_PHYS, FACE_CLOSED
    use :: ibmm, only: ibm_type
    use :: boundary, only: boundary_type, apply_bc
    use :: comm, only: comm_type, exchange_halos

    implicit none

    private
    public :: pressure_solver_type, init_pressure_solver, pressure_projection

    type :: pressure_solver_type
        integer(C_INT) :: nIter=3
        real(C_DOUBLE) :: sor=1.5d0
    end type pressure_solver_type

contains

    subroutine init_pressure_solver(ps, dns, bc, has_terminal)
        type(pressure_solver_type), intent(inout) :: ps
        type(dns_type), intent(in) :: dns
        type(boundary_type), intent(in) :: bc
        logical, intent(in), optional :: has_terminal
        integer :: dir
        logical :: terminal

        terminal = .true.
        if (present(has_terminal)) terminal = has_terminal


        do dir = 1, 3
            if (bc%isPeriodic(dir) .and. mod(dns%globalSize(dir), 2_C_INT) /= 0) then
                if (terminal) then
                    print *, "invalid red-black grid: periodic direction", dir, &
                             "has odd global size", dns%globalSize(dir)
                end if
                error stop "red-black pressure solver requires even global sizes in periodic directions"
            end if
        end do

    end subroutine init_pressure_solver

    subroutine pressure_projection(ps, blk, dns, dt_gamma, ibm, bc, c)
        type(pressure_solver_type), intent(in) :: ps
        type(block_set_type), intent(inout) :: blk
        type(dns_type), intent(in) :: dns
        real(C_DOUBLE), intent(in) :: dt_gamma
        type(ibm_type), intent(in) :: ibm
        type(boundary_type), intent(in) :: bc
        type(comm_type), intent(inout) :: c

        integer(C_INT) :: iIter, color

        do iIter = 1_C_INT, ps%nIter
            do color = 1_C_INT, 0_C_INT, -1_C_INT
                call redblack_sweep(ps, blk, dt_gamma, ibm, color)
                call apply_bc(blk, bc)
                if (iIter == ps%nIter .and. color == 0_C_INT) then
                    call exchange_halos(c, blk, [VAR_U, VAR_V, VAR_W, VAR_P])
                else
                    call exchange_halos(c, blk, [VAR_U, VAR_V, VAR_W])
                end if
            end do
        end do

    end subroutine pressure_projection

    subroutine redblack_sweep(ps, blk, dt_gamma, ibm, color)
        type(pressure_solver_type), intent(in) :: ps
        type(block_set_type), intent(inout) :: blk
        real(C_DOUBLE), intent(in) :: dt_gamma
        type(ibm_type), intent(in) :: ibm
        integer(C_INT), intent(in) :: color

        real(C_DOUBLE) :: phi,denom,idt,sor
        real(C_DOUBLE) :: div
        real(C_DOUBLE) :: mu_u_i,mu_u_ip,mu_v_j,mu_v_jp,mu_w_k,mu_w_kp
        integer(C_INT) :: i,ip,j,jp,k,kp,b,nBlocks,nLowerHaloDirections,iColor,nColorX
        integer(C_INT) :: iLo,jLo,kLo,colorOffset
        integer(C_INT) :: hi(1:3)

        ! Each block sweeps from 0 (its halo layer, redundantly with the
        ! neighbour that owns those cells) except on physical boundaries,
        ! exactly the rank-level scheme one level down. Red-black parity is
        ! anchored to the global index space through the block origin.
        hi = blk%nb(1:3)
        sor = ps%sor
        nColorX = (hi(1) + 2_C_INT)/2_C_INT
        nBlocks = blk%nBlocks

        idt = 1.0_C_DOUBLE/dt_gamma

#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do collapse(4) &
        !$omp& map(to: color, nColorX, sor, idt, dt_gamma, hi(1:3), &
        !$omp& blk%origin, blk%physLow, blk%physHigh, &
        !$omp& blk%d1x, blk%d1y, blk%d1z, ibm%mu) &
        !$omp& map(tofrom: blk%q) &
        !$omp& private(i,ip,j,jp,k,kp,b,iColor,iLo,jLo,kLo,colorOffset, &
        !$omp& phi,denom,div,nLowerHaloDirections, &
        !$omp& mu_u_i,mu_u_ip,mu_v_j,mu_v_jp,mu_w_k,mu_w_kp)
#endif
        DO b=1_C_INT,nBlocks
        DO k=0_C_INT,hi(3)
            DO j=0_C_INT,hi(2)
                DO iColor=0_C_INT,nColorX-1_C_INT
                    iLo = merge(1_C_INT, 0_C_INT, blk%physLow(1,b) /= 0_C_INT)
                    jLo = merge(1_C_INT, 0_C_INT, blk%physLow(2,b) /= 0_C_INT)
                    kLo = merge(1_C_INT, 0_C_INT, blk%physLow(3,b) /= 0_C_INT)
                    if (j < jLo .or. k < kLo) cycle
                    colorOffset = modulo(blk%origin(1,b) + blk%origin(2,b) + blk%origin(3,b), 2_C_INT)
                    i = iLo + modulo(color - modulo(iLo+j+k+colorOffset, 2_C_INT), 2_C_INT) &
                        + 2_C_INT*iColor
                    if (i > hi(1)) cycle

                    nLowerHaloDirections = 0_C_INT
                    if (i == 0_C_INT) nLowerHaloDirections = nLowerHaloDirections + 1_C_INT
                    if (j == 0_C_INT) nLowerHaloDirections = nLowerHaloDirections + 1_C_INT
                    if (k == 0_C_INT) nLowerHaloDirections = nLowerHaloDirections + 1_C_INT
                    if (nLowerHaloDirections > 1_C_INT) cycle

                    jp = j + 1
                    kp = k + 1
                    ip = i + 1

                    mu_u_i  = ibm%mu(i,j,k,VAR_U,b)
                    mu_u_ip = ibm%mu(ip,j,k,VAR_U,b)
                    mu_v_j  = ibm%mu(i,j,k,VAR_V,b)
                    mu_v_jp = ibm%mu(i,jp,k,VAR_V,b)
                    mu_w_k  = ibm%mu(i,j,k,VAR_W,b)
                    mu_w_kp = ibm%mu(i,j,kp,VAR_W,b)

                    denom = (merge(0.0d0, mu_u_i*blk%d1x(i,VAR_U,b), &
                                      noflux_low(blk%physLow(1,b)) .and. i == 1_C_INT) &
                           + merge(0.0d0, mu_u_ip*blk%d1x(ip,VAR_U,b), &
                                      noflux_high(blk%physHigh(1,b)) .and. i == hi(1)))*blk%d1x(i,VAR_P,b) &
                          + (merge(0.0d0, mu_v_j*blk%d1y(j,VAR_V,b), &
                                      noflux_low(blk%physLow(2,b)) .and. j == 1_C_INT) &
                           + merge(0.0d0, mu_v_jp*blk%d1y(jp,VAR_V,b), &
                                      noflux_high(blk%physHigh(2,b)) .and. j == hi(2)))*blk%d1y(j,VAR_P,b) &
                          + (merge(0.0d0, mu_w_k*blk%d1z(k,VAR_W,b), &
                                      noflux_low(blk%physLow(3,b)) .and. k == 1_C_INT) &
                           + merge(0.0d0, mu_w_kp*blk%d1z(kp,VAR_W,b), &
                                      noflux_high(blk%physHigh(3,b)) .and. k == hi(3)))*blk%d1z(k,VAR_P,b)

                    div = (blk%q(ip,j,k,VAR_U,b)-blk%q(i,j,k,VAR_U,b))*blk%d1x(i,VAR_P,b) &
                        + (blk%q(i,jp,k,VAR_V,b)-blk%q(i,j,k,VAR_V,b))*blk%d1y(j,VAR_P,b) &
                        + (blk%q(i,j,kp,VAR_W,b)-blk%q(i,j,k,VAR_W,b))*blk%d1z(k,VAR_P,b)

                    phi = -sor*div/denom

                    blk%q(i,j,k,VAR_P,b) = blk%q(i,j,k,VAR_P,b) + phi*idt

                    ! Face corrections are masked on no-flux faces with the
                    ! same conditions as denom. On physical walls the value
                    ! was dead anyway (apply_bc overwrites it before any
                    ! read); on FACE_CLOSED faces this keeps the pinned
                    ! interface velocity exactly zero.
                    blk%q(i,j,k,VAR_U,b) = blk%q(i,j,k,VAR_U,b) &
                        - merge(0.0d0, phi*blk%d1x(i,VAR_U,b)*mu_u_i, &
                                noflux_low(blk%physLow(1,b)) .and. i == 1_C_INT)
                    blk%q(ip,j,k,VAR_U,b) = blk%q(ip,j,k,VAR_U,b) &
                        + merge(0.0d0, phi*blk%d1x(ip,VAR_U,b)*mu_u_ip, &
                                noflux_high(blk%physHigh(1,b)) .and. i == hi(1))
                    blk%q(i,j,k,VAR_V,b) = blk%q(i,j,k,VAR_V,b) &
                        - merge(0.0d0, phi*blk%d1y(j,VAR_V,b)*mu_v_j, &
                                noflux_low(blk%physLow(2,b)) .and. j == 1_C_INT)
                    blk%q(i,jp,k,VAR_V,b) = blk%q(i,jp,k,VAR_V,b) &
                        + merge(0.0d0, phi*blk%d1y(jp,VAR_V,b)*mu_v_jp, &
                                noflux_high(blk%physHigh(2,b)) .and. j == hi(2))
                    blk%q(i,j,k,VAR_W,b) = blk%q(i,j,k,VAR_W,b) &
                        - merge(0.0d0, phi*blk%d1z(k,VAR_W,b)*mu_w_k, &
                                noflux_low(blk%physLow(3,b)) .and. k == 1_C_INT)
                    blk%q(i,j,kp,VAR_W,b) = blk%q(i,j,kp,VAR_W,b) &
                        + merge(0.0d0, phi*blk%d1z(kp,VAR_W,b)*mu_w_kp, &
                                noflux_high(blk%physHigh(3,b)) .and. k == hi(3))
                END DO
            END DO
        END DO
        END DO
#ifdef USE_OPENMP_OFFLOAD
        !$omp end target teams distribute parallel do
#endif

    end subroutine redblack_sweep

    ! Faces excluded from the sweep denominator and corrections. The
    ! LOW-side block owns a 2:1 shared face (strategy doc 6): low faces
    ! are masked only at physical walls and closed faces, while a block's
    ! HIGH face is masked at any interface too - its halo copy is
    ! refreshed by the exchange (restriction or injection) and only feeds
    ! the divergence.
    pure logical function noflux_low(fk)
!$omp declare target
        integer(C_INT), intent(in) :: fk

        noflux_low = fk == FACE_PHYS .or. fk == FACE_CLOSED
    end function noflux_low

    pure logical function noflux_high(fk)
!$omp declare target
        integer(C_INT), intent(in) :: fk

        noflux_high = fk /= FACE_OPEN
    end function noflux_high

end module pressure_solver
