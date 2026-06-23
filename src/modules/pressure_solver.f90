module pressure_solver
    use, intrinsic :: iso_c_binding
    use :: init, only: dns_type, VAR_U, VAR_V, VAR_W, VAR_P
    use :: blocks, only: block_set_type, FACE_PHYS, FACE_CLOSED
    use :: ibmm, only: ibm_type
    use :: boundary, only: boundary_type, apply_bc
    use :: comm, only: comm_type, exchange_halos, exchange_scalar_halos

    implicit none

    private
    public :: pressure_solver_type, init_pressure_solver, pressure_projection

    type :: pressure_solver_type
        integer(C_INT) :: nIter=3
        ! Damped-Jacobi relaxation factor. The pressure projection is now a
        ! simple (un-coloured) damped Jacobi iteration; for the Poisson-like
        ! projection operator the high-frequency (checkerboard) mode forces
        ! the factor strictly below 1 (omega=1 is marginally unstable). 0.8 is
        ! a safe default; the config key is still "sor".
        real(C_DOUBLE) :: sor=0.8d0
    end type pressure_solver_type

    ! Pressure-increment buffer for the Jacobi sweep. Unlike the in-place
    ! red-black SOR, Jacobi computes every cell's increment from the SAME frozen
    ! divergence, then applies pressure and the velocity-face corrections in
    ! separate race-free passes. phi shares blk%q's spatial bounds (0:nb+1) so
    ! exchange_scalar_halos can fill the halo layer the face corrections read.
    real(C_DOUBLE), allocatable :: phi(:,:,:,:)

contains

    subroutine init_pressure_solver(ps, dns, bc, has_terminal)
        type(pressure_solver_type), intent(inout) :: ps
        type(dns_type), intent(in) :: dns
        type(boundary_type), intent(in) :: bc
        logical, intent(in), optional :: has_terminal

        ! Damped Jacobi needs no red-black colouring, so the even-global-size
        ! restriction of the red-black scheme no longer applies.
    end subroutine init_pressure_solver

    subroutine pressure_projection(ps, blk, dns, dt_gamma, ibm, bc, c)
        type(pressure_solver_type), intent(in) :: ps
        type(block_set_type), intent(inout) :: blk
        type(dns_type), intent(in) :: dns
        real(C_DOUBLE), intent(in) :: dt_gamma
        type(ibm_type), intent(in) :: ibm
        type(boundary_type), intent(in) :: bc
        type(comm_type), intent(inout) :: c

        integer(C_INT) :: iIter

        call allocate_phi(blk)

        ! Damped-Jacobi projection. Each iteration: (1) compute the pressure
        ! increment phi = -omega*div/denom from the frozen velocity, (2) exchange
        ! phi's halos so the face corrections see the neighbour increment, (3)
        ! apply phi to the pressure and the velocity faces, (4) refresh the
        ! velocity (and, on the last iteration, pressure) halos for the next
        ! divergence. The 2:1 interface is still handled by the exchange transfer
        ! (RESTRICT/PROLONG) and reconciled by the final full exchange.
        do iIter = 1_C_INT, ps%nIter
            call jacobi_compute_phi(ps, blk, dt_gamma, ibm)
            call exchange_scalar_halos(c, phi)
            call jacobi_apply(ps, blk, dt_gamma, ibm)
            call apply_bc(blk, bc)
            if (iIter == ps%nIter) then
                call exchange_halos(c, blk, [VAR_U, VAR_V, VAR_W, VAR_P])
            else
                call exchange_halos(c, blk, [VAR_U, VAR_V, VAR_W], interp=.false.)
            end if
        end do
    end subroutine pressure_projection

    ! Allocate the device-resident phi buffer once (same 0:nb+1 spatial bounds
    ! as blk%q so the scalar halo exchange and the face corrections line up).
    subroutine allocate_phi(blk)
        type(block_set_type), intent(in) :: blk
        if (.not. allocated(phi)) then
            allocate(phi(lbound(blk%q,1):ubound(blk%q,1), &
                         lbound(blk%q,2):ubound(blk%q,2), &
                         lbound(blk%q,3):ubound(blk%q,3), blk%nBlocks))
            phi = 0.0d0
#ifdef USE_OPENMP_OFFLOAD
            !$omp target enter data map(to: phi)
#endif
        end if
    end subroutine allocate_phi

    ! Jacobi step 1: phi(i,j,k) = -omega * div / denom for every interior cell,
    ! from the frozen velocity. denom is the projection diagonal (sum of the
    ! cell's non-pinned face metrics); pinned faces (walls, closed faces) leave
    ! the diagonal exactly as in the red-black scheme.
    subroutine jacobi_compute_phi(ps, blk, dt_gamma, ibm)
        type(pressure_solver_type), intent(in) :: ps
        type(block_set_type), intent(inout) :: blk
        real(C_DOUBLE), intent(in) :: dt_gamma
        type(ibm_type), intent(in) :: ibm

        real(C_DOUBLE) :: denom, div, omega
        real(C_DOUBLE) :: mu_u_i, mu_u_ip, mu_v_j, mu_v_jp, mu_w_k, mu_w_kp
        integer(C_INT) :: i, ip, j, jp, k, kp, b, nBlocks, nx, ny, nz

        nx = blk%nb(1); ny = blk%nb(2); nz = blk%nb(3)
        omega = ps%sor
        nBlocks = blk%nBlocks

#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do collapse(4) &
        !$omp& map(to: omega, nx, ny, nz, blk%physLow, blk%physHigh, &
        !$omp& blk%d1x, blk%d1y, blk%d1z, ibm%mu) map(tofrom: phi, blk%q) &
        !$omp& private(i,ip,j,jp,k,kp,b,denom,div, &
        !$omp& mu_u_i,mu_u_ip,mu_v_j,mu_v_jp,mu_w_k,mu_w_kp)
#endif
        do b = 1_C_INT, nBlocks
        do k = 1_C_INT, nz
            do j = 1_C_INT, ny
                do i = 1_C_INT, nx
                    ip = i + 1; jp = j + 1; kp = k + 1

                    mu_u_i  = ibm%mu(i,j,k,VAR_U,b);  mu_u_ip = ibm%mu(ip,j,k,VAR_U,b)
                    mu_v_j  = ibm%mu(i,j,k,VAR_V,b);  mu_v_jp = ibm%mu(i,jp,k,VAR_V,b)
                    mu_w_k  = ibm%mu(i,j,k,VAR_W,b);  mu_w_kp = ibm%mu(i,j,kp,VAR_W,b)

                    denom = (merge(0.0d0, mu_u_i*blk%d1x(i,VAR_U,b), &
                                      face_pinned(blk%physLow(1,b)) .and. i == 1_C_INT) &
                           + merge(0.0d0, mu_u_ip*blk%d1x(ip,VAR_U,b), &
                                      face_pinned(blk%physHigh(1,b)) .and. i == nx))*blk%d1x(i,VAR_P,b) &
                          + (merge(0.0d0, mu_v_j*blk%d1y(j,VAR_V,b), &
                                      face_pinned(blk%physLow(2,b)) .and. j == 1_C_INT) &
                           + merge(0.0d0, mu_v_jp*blk%d1y(jp,VAR_V,b), &
                                      face_pinned(blk%physHigh(2,b)) .and. j == ny))*blk%d1y(j,VAR_P,b) &
                          + (merge(0.0d0, mu_w_k*blk%d1z(k,VAR_W,b), &
                                      face_pinned(blk%physLow(3,b)) .and. k == 1_C_INT) &
                           + merge(0.0d0, mu_w_kp*blk%d1z(kp,VAR_W,b), &
                                      face_pinned(blk%physHigh(3,b)) .and. k == nz))*blk%d1z(k,VAR_P,b)

                    div = (blk%q(ip,j,k,VAR_U,b)-blk%q(i,j,k,VAR_U,b))*blk%d1x(i,VAR_P,b) &
                        + (blk%q(i,jp,k,VAR_V,b)-blk%q(i,j,k,VAR_V,b))*blk%d1y(j,VAR_P,b) &
                        + (blk%q(i,j,kp,VAR_W,b)-blk%q(i,j,k,VAR_W,b))*blk%d1z(k,VAR_P,b)

                    phi(i,j,k,b) = -omega*div/denom
                end do
            end do
        end do
        end do
#ifdef USE_OPENMP_OFFLOAD
        !$omp end target teams distribute parallel do
#endif
    end subroutine jacobi_compute_phi

    ! Jacobi step 3: apply the frozen increment. Pressure gets p += phi/dt_gamma;
    ! each velocity face is corrected once from the two cells it separates,
    !   q_face += (phi_below - phi_above) * d1(face) * mu(face),
    ! the exact sum of the two adjacent red-black corrections, but read from the
    ! frozen phi (with the halo layer filled by the scalar exchange) so there is
    ! no in-place race and no colouring. Only the cell's own LOW faces (1..nb)
    ! are written here; each block's high halo face is the neighbour's low face,
    ! filled by the velocity exchange. Pinned faces are left untouched.
    subroutine jacobi_apply(ps, blk, dt_gamma, ibm)
        type(pressure_solver_type), intent(in) :: ps
        type(block_set_type), intent(inout) :: blk
        real(C_DOUBLE), intent(in) :: dt_gamma
        type(ibm_type), intent(in) :: ibm

        real(C_DOUBLE) :: idt
        integer(C_INT) :: i, j, k, b, nBlocks, nx, ny, nz

        nx = blk%nb(1); ny = blk%nb(2); nz = blk%nb(3)
        nBlocks = blk%nBlocks
        idt = 1.0_C_DOUBLE/dt_gamma

        ! Pressure update (interior cells).
#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do collapse(4) &
        !$omp& map(to: idt, nx, ny, nz) map(tofrom: blk%q, phi) private(i,j,k,b)
#endif
        do b = 1_C_INT, nBlocks
        do k = 1_C_INT, nz
            do j = 1_C_INT, ny
                do i = 1_C_INT, nx
                    blk%q(i,j,k,VAR_P,b) = blk%q(i,j,k,VAR_P,b) + phi(i,j,k,b)*idt
                end do
            end do
        end do
        end do
#ifdef USE_OPENMP_OFFLOAD
        !$omp end target teams distribute parallel do
#endif

        ! Velocity-face corrections (each block's low faces 1..nb; the high halo
        ! face is the neighbour's low face, supplied by the exchange).
#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do collapse(4) &
        !$omp& map(to: nx, ny, nz, blk%physLow, blk%d1x, blk%d1y, blk%d1z, ibm%mu) &
        !$omp& map(tofrom: blk%q, phi) private(i,j,k,b)
#endif
        do b = 1_C_INT, nBlocks
        do k = 1_C_INT, nz
            do j = 1_C_INT, ny
                do i = 1_C_INT, nx
                    if (.not. (i == 1_C_INT .and. face_pinned(blk%physLow(1,b)))) &
                        blk%q(i,j,k,VAR_U,b) = blk%q(i,j,k,VAR_U,b) &
                            + (phi(i-1,j,k,b) - phi(i,j,k,b))*blk%d1x(i,VAR_U,b)*ibm%mu(i,j,k,VAR_U,b)
                    if (.not. (j == 1_C_INT .and. face_pinned(blk%physLow(2,b)))) &
                        blk%q(i,j,k,VAR_V,b) = blk%q(i,j,k,VAR_V,b) &
                            + (phi(i,j-1,k,b) - phi(i,j,k,b))*blk%d1y(j,VAR_V,b)*ibm%mu(i,j,k,VAR_V,b)
                    if (.not. (k == 1_C_INT .and. face_pinned(blk%physLow(3,b)))) &
                        blk%q(i,j,k,VAR_W,b) = blk%q(i,j,k,VAR_W,b) &
                            + (phi(i,j,k-1,b) - phi(i,j,k,b))*blk%d1z(k,VAR_W,b)*ibm%mu(i,j,k,VAR_W,b)
                end do
            end do
        end do
        end do
#ifdef USE_OPENMP_OFFLOAD
        !$omp end target teams distribute parallel do
#endif
    end subroutine jacobi_apply

    ! Pinned faces carry zero flux forever (physical walls, closed faces against
    ! removed blocks): they leave both the diagonal and the corrections.
    pure logical function face_pinned(fk)
!$omp declare target
        integer(C_INT), intent(in) :: fk

        face_pinned = fk == FACE_PHYS .or. fk == FACE_CLOSED
    end function face_pinned

end module pressure_solver
