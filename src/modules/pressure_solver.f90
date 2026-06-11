module pressure_solver
    use, intrinsic :: iso_c_binding
    use :: init, only: dns_type, VAR_U, VAR_V, VAR_W, VAR_P
    use :: blocks, only: block_set_type
    use :: ibmm, only: ibm_type
    use :: boundary, only: boundary_type, apply_bc
    use :: comm, only: comm_type, exchange_halos

    implicit none

    private
    public :: pressure_solver_type, init_pressure_solver, pressure_projection

    type :: pressure_solver_type
        integer(C_INT) :: nIter=3
        real(C_DOUBLE) :: sor=1.5d0
        integer(C_INT) :: sweepLo(1:3) = 1_C_INT
        logical(C_BOOL) :: pressureNeumannLow(1:3) = .false.
        logical(C_BOOL) :: pressureNeumannHigh(1:3) = .false.
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

        ps%sweepLo = 1_C_INT
        ps%pressureNeumannLow = .false.
        ps%pressureNeumannHigh = .false.

        do dir = 1, 3
            if (bc%isPeriodic(dir) .or. dns%localSize(dir,0) > 1_C_INT) then
                ps%sweepLo(dir) = 0_C_INT
            end if
            ps%pressureNeumannLow(dir) = (.not. bc%isPeriodic(dir)) .and. dns%localSize(dir,0) == 1_C_INT
            ps%pressureNeumannHigh(dir) = (.not. bc%isPeriodic(dir)) .and. dns%localSize(dir,1) == dns%globalSize(dir)
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
        integer(C_INT) :: local_n(1:3), colorOffset

        local_n = blk%nb(1:3)
        ! Red-black parity from the block's global cell origin. Phase 0: one block
        ! per rank, so a single offset; Phase 1 moves this inside the sweep, one
        ! offset per block.
        colorOffset = modulo(sum(blk%origin(1:3,1)), 2_C_INT)

        do iIter = 1_C_INT, ps%nIter
            do color = 1_C_INT, 0_C_INT, -1_C_INT
                call redblack_sweep(ps, blk, dt_gamma, ibm, local_n, color, colorOffset)
                call apply_bc(blk, bc)
                if (iIter == ps%nIter .and. color == 0_C_INT) then
                    call exchange_halos(c, blk, [VAR_U, VAR_V, VAR_W, VAR_P])
                else
                    call exchange_halos(c, blk, [VAR_U, VAR_V, VAR_W])
                end if
            end do
        end do

    end subroutine pressure_projection

    subroutine redblack_sweep(ps, blk, dt_gamma, ibm, local_n, color, colorOffset)
        type(pressure_solver_type), intent(in) :: ps
        type(block_set_type), intent(inout) :: blk
        real(C_DOUBLE), intent(in) :: dt_gamma
        type(ibm_type), intent(in) :: ibm
        integer(C_INT), intent(in) :: local_n(1:3), color, colorOffset

        real(C_DOUBLE) :: phi,denom,idt,sor
        real(C_DOUBLE) :: div
        real(C_DOUBLE) :: mu_u_i,mu_u_ip,mu_v_j,mu_v_jp,mu_w_k,mu_w_kp
        integer(C_INT) :: i,ip,j,jp,k,kp,b,nBlocks,nLowerHaloDirections,iColor,nColorX
        integer(C_INT) :: lo(1:3), hi(1:3)
        logical(C_BOOL) :: pressureNeumannLow(1:3), pressureNeumannHigh(1:3)

        lo = ps%sweepLo
        hi = local_n
        pressureNeumannLow = ps%pressureNeumannLow
        pressureNeumannHigh = ps%pressureNeumannHigh
        sor = ps%sor
        nColorX = (hi(1) - lo(1) + 2_C_INT)/2_C_INT
        nBlocks = blk%nBlocks

        idt = 1.0_C_DOUBLE/dt_gamma

#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do collapse(4) &
        !$omp& map(to: color, colorOffset, nColorX, sor, idt, dt_gamma, &
        !$omp& lo(1:3), hi(1:3), pressureNeumannLow(1:3), pressureNeumannHigh(1:3), &
        !$omp& blk%d1x, blk%d1y, blk%d1z, ibm%mu) &
        !$omp& map(tofrom: blk%q) &
        !$omp& private(i,ip,j,jp,k,kp,b,iColor,phi,denom,div,nLowerHaloDirections, &
        !$omp& mu_u_i,mu_u_ip,mu_v_j,mu_v_jp,mu_w_k,mu_w_kp)
#endif
        DO b=1_C_INT,nBlocks
        DO k=lo(3),hi(3)
            DO j=lo(2),hi(2)
                DO iColor=0_C_INT,nColorX-1_C_INT
                    i = lo(1) + modulo(color - modulo(lo(1)+j+k+colorOffset, 2_C_INT), 2_C_INT) &
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
                                      pressureNeumannLow(1) .and. i == 1_C_INT) &
                           + merge(0.0d0, mu_u_ip*blk%d1x(ip,VAR_U,b), &
                                      pressureNeumannHigh(1) .and. i == hi(1)))*blk%d1x(i,VAR_P,b) &
                          + (merge(0.0d0, mu_v_j*blk%d1y(j,VAR_V,b), &
                                      pressureNeumannLow(2) .and. j == 1_C_INT) &
                           + merge(0.0d0, mu_v_jp*blk%d1y(jp,VAR_V,b), &
                                      pressureNeumannHigh(2) .and. j == hi(2)))*blk%d1y(j,VAR_P,b) &
                          + (merge(0.0d0, mu_w_k*blk%d1z(k,VAR_W,b), &
                                      pressureNeumannLow(3) .and. k == 1_C_INT) &
                           + merge(0.0d0, mu_w_kp*blk%d1z(kp,VAR_W,b), &
                                      pressureNeumannHigh(3) .and. k == hi(3)))*blk%d1z(k,VAR_P,b)

                    div = (blk%q(ip,j,k,VAR_U,b)-blk%q(i,j,k,VAR_U,b))*blk%d1x(i,VAR_P,b) &
                        + (blk%q(i,jp,k,VAR_V,b)-blk%q(i,j,k,VAR_V,b))*blk%d1y(j,VAR_P,b) &
                        + (blk%q(i,j,kp,VAR_W,b)-blk%q(i,j,k,VAR_W,b))*blk%d1z(k,VAR_P,b)

                    phi = -sor*div/denom

                    blk%q(i,j,k,VAR_P,b) = blk%q(i,j,k,VAR_P,b) + phi*idt

                    blk%q(i,j,k,VAR_U,b) = blk%q(i,j,k,VAR_U,b) - phi*blk%d1x(i,VAR_U,b)*mu_u_i
                    blk%q(ip,j,k,VAR_U,b) = blk%q(ip,j,k,VAR_U,b) + phi*blk%d1x(ip,VAR_U,b)*mu_u_ip
                    blk%q(i,j,k,VAR_V,b) = blk%q(i,j,k,VAR_V,b) - phi*blk%d1y(j,VAR_V,b)*mu_v_j
                    blk%q(i,jp,k,VAR_V,b) = blk%q(i,jp,k,VAR_V,b) + phi*blk%d1y(jp,VAR_V,b)*mu_v_jp
                    blk%q(i,j,k,VAR_W,b) = blk%q(i,j,k,VAR_W,b) - phi*blk%d1z(k,VAR_W,b)*mu_w_k
                    blk%q(i,j,kp,VAR_W,b) = blk%q(i,j,kp,VAR_W,b) + phi*blk%d1z(kp,VAR_W,b)*mu_w_kp
                END DO
            END DO
        END DO
        END DO
#ifdef USE_OPENMP_OFFLOAD
        !$omp end target teams distribute parallel do
#endif

    end subroutine redblack_sweep

end module pressure_solver
