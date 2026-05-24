module pressure_solver
    use, intrinsic :: iso_c_binding
    use :: init, only: dns_type, grid_type, field_type, VAR_U, VAR_V, VAR_W, VAR_P
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

    subroutine pressure_projection(ps, f, dns, g, dt_gamma, ibm, bc, c)
        type(pressure_solver_type), intent(in) :: ps
        type(field_type), intent(inout) :: f
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        real(C_DOUBLE), intent(in) :: dt_gamma
        type(ibm_type), intent(in) :: ibm
        type(boundary_type), intent(in) :: bc
        type(comm_type), intent(inout) :: c

        integer(C_INT) :: iIter, color
        integer(C_INT) :: local_n(1:3), colorOffset

        local_n = dns%localSize(1:3,2)
        colorOffset = modulo(sum(dns%localSize(1:3,0) - 1_C_INT), 2_C_INT)

        do iIter = 1_C_INT, ps%nIter
            do color = 1_C_INT, 0_C_INT, -1_C_INT
                call redblack_sweep(ps, f, g, dt_gamma, ibm, local_n, color, colorOffset)
                call apply_bc(f, dns, g, bc)
                if (iIter == ps%nIter .and. color == 0_C_INT) then
                    call exchange_halos(c, f, [VAR_U, VAR_V, VAR_W, VAR_P])
                else
                    call exchange_halos(c, f, [VAR_U, VAR_V, VAR_W])
                end if
            end do
        end do

    end subroutine pressure_projection

    subroutine redblack_sweep(ps, f, g, dt_gamma, ibm, local_n, color, colorOffset)
        type(pressure_solver_type), intent(in) :: ps
        type(field_type), intent(inout) :: f
        type(grid_type), intent(in) :: g
        real(C_DOUBLE), intent(in) :: dt_gamma
        type(ibm_type), intent(in) :: ibm
        integer(C_INT), intent(in) :: local_n(1:3), color, colorOffset

        real(C_DOUBLE) :: phi,denom,idt,sor
        real(C_DOUBLE) :: div
        real(C_DOUBLE) :: mu_u_i,mu_u_ip,mu_v_j,mu_v_jp,mu_w_k,mu_w_kp
        integer(C_INT) :: i,ip,j,jp,k,kp,nLowerHaloDirections,iColor,nColorX
        integer(C_INT) :: lo(1:3), hi(1:3)
        logical(C_BOOL) :: pressureNeumannLow(1:3), pressureNeumannHigh(1:3)

        lo = ps%sweepLo
        hi = local_n
        pressureNeumannLow = ps%pressureNeumannLow
        pressureNeumannHigh = ps%pressureNeumannHigh
        sor = ps%sor
        nColorX = (hi(1) - lo(1) + 2_C_INT)/2_C_INT

        idt = 1.0_C_DOUBLE/dt_gamma

#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do collapse(3) &
        !$omp& map(to: color, colorOffset, nColorX, sor, idt, dt_gamma, &
        !$omp& lo(1:3), hi(1:3), pressureNeumannLow(1:3), pressureNeumannHigh(1:3), &
        !$omp& g%d1x, g%d1y, g%d1z, ibm%coef) &
        !$omp& map(tofrom: f%q) &
        !$omp& private(i,ip,j,jp,k,kp,iColor,phi,denom,div,nLowerHaloDirections, &
        !$omp& mu_u_i,mu_u_ip,mu_v_j,mu_v_jp,mu_w_k,mu_w_kp)
#endif
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

                    mu_u_i  = 1.0d0/(1.0d0+dt_gamma*ibm%coef(i,j,k,VAR_U))
                    mu_u_ip = 1.0d0/(1.0d0+dt_gamma*ibm%coef(ip,j,k,VAR_U))
                    mu_v_j  = 1.0d0/(1.0d0+dt_gamma*ibm%coef(i,j,k,VAR_V))
                    mu_v_jp = 1.0d0/(1.0d0+dt_gamma*ibm%coef(i,jp,k,VAR_V))
                    mu_w_k  = 1.0d0/(1.0d0+dt_gamma*ibm%coef(i,j,k,VAR_W))
                    mu_w_kp = 1.0d0/(1.0d0+dt_gamma*ibm%coef(i,j,kp,VAR_W))

                    denom = (merge(0.0d0, mu_u_i*g%d1x(i,VAR_U), &
                                      pressureNeumannLow(1) .and. i == 1_C_INT) &
                           + merge(0.0d0, mu_u_ip*g%d1x(ip,VAR_U), &
                                      pressureNeumannHigh(1) .and. i == hi(1)))*g%d1x(i,VAR_P) &
                          + (merge(0.0d0, mu_v_j*g%d1y(j,VAR_V), &
                                      pressureNeumannLow(2) .and. j == 1_C_INT) &
                           + merge(0.0d0, mu_v_jp*g%d1y(jp,VAR_V), &
                                      pressureNeumannHigh(2) .and. j == hi(2)))*g%d1y(j,VAR_P) &
                          + (merge(0.0d0, mu_w_k*g%d1z(k,VAR_W), &
                                      pressureNeumannLow(3) .and. k == 1_C_INT) &
                           + merge(0.0d0, mu_w_kp*g%d1z(kp,VAR_W), &
                                      pressureNeumannHigh(3) .and. k == hi(3)))*g%d1z(k,VAR_P)

                    div = (f%q(ip,j,k,VAR_U)-f%q(i,j,k,VAR_U))*g%d1x(i,VAR_P) &
                        + (f%q(i,jp,k,VAR_V)-f%q(i,j,k,VAR_V))*g%d1y(j,VAR_P) &
                        + (f%q(i,j,kp,VAR_W)-f%q(i,j,k,VAR_W))*g%d1z(k,VAR_P)

                    phi = -sor*div/denom

                    f%q(i,j,k,VAR_P) = f%q(i,j,k,VAR_P) + phi*idt

                    f%q(i,j,k,VAR_U) = f%q(i,j,k,VAR_U) - phi*g%d1x(i,VAR_U)*mu_u_i
                    f%q(ip,j,k,VAR_U) = f%q(ip,j,k,VAR_U) + phi*g%d1x(ip,VAR_U)*mu_u_ip
                    f%q(i,j,k,VAR_V) = f%q(i,j,k,VAR_V) - phi*g%d1y(j,VAR_V)*mu_v_j
                    f%q(i,jp,k,VAR_V) = f%q(i,jp,k,VAR_V) + phi*g%d1y(jp,VAR_V)*mu_v_jp
                    f%q(i,j,k,VAR_W) = f%q(i,j,k,VAR_W) - phi*g%d1z(k,VAR_W)*mu_w_k
                    f%q(i,j,kp,VAR_W) = f%q(i,j,kp,VAR_W) + phi*g%d1z(kp,VAR_W)*mu_w_kp
                END DO
            END DO
        END DO
#ifdef USE_OPENMP_OFFLOAD
        !$omp end target teams distribute parallel do
#endif

    end subroutine redblack_sweep

end module pressure_solver
