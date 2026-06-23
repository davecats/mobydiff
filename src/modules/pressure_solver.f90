module pressure_solver
    use, intrinsic :: iso_c_binding
    use :: init, only: dns_type, VAR_U, VAR_V, VAR_W, VAR_P
    use :: blocks, only: block_set_type, FACE_PHYS, FACE_CLOSED, FACE_COARSE, FACE_FINE
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

                    ! Diagonal: each face's pressure-gradient metric. face_grad
                    ! returns 0 for a pinned wall face, the coarse-fine gradient
                    ! 1/d for a 2:1 interface face (the symmetric composite stencil),
                    ! and the regular 1/h otherwise.
                    denom = (face_grad(blk%physLow(1,b), i == 1_C_INT, blk%d1x(i,VAR_U,b))*mu_u_i &
                           + face_grad(blk%physHigh(1,b), i == nx, blk%d1x(ip,VAR_U,b))*mu_u_ip)*blk%d1x(i,VAR_P,b) &
                          + (face_grad(blk%physLow(2,b), j == 1_C_INT, blk%d1y(j,VAR_V,b))*mu_v_j &
                           + face_grad(blk%physHigh(2,b), j == ny, blk%d1y(jp,VAR_V,b))*mu_v_jp)*blk%d1y(j,VAR_P,b) &
                          + (face_grad(blk%physLow(3,b), k == 1_C_INT, blk%d1z(k,VAR_W,b))*mu_w_k &
                           + face_grad(blk%physHigh(3,b), k == nz, blk%d1z(kp,VAR_W,b))*mu_w_kp)*blk%d1z(k,VAR_P,b)

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

        real(C_DOUBLE) :: idt, cf
        integer(C_INT) :: i, ip, j, jp, k, kp, b, nBlocks, nx, ny, nz

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

        ! Velocity-face corrections. Each cell corrects its LOW faces (1..nb):
        ! q_face += (phi_below - phi_self) * grad(face) * mu, with grad(face) from
        ! face_grad (0 if pinned, ifGrad if a 2:1 interface, 1/h otherwise). The
        ! HIGH face is normally the neighbour's low face (filled by the exchange),
        ! so it is corrected here ONLY when it is itself a 2:1 interface face --
        ! then this block owns it, reconstructing it from the halo phi (the coarse
        ! halo = avg of the fine increments / the fine halo = the coarse increment,
        ! both from the scalar exchange), which keeps the two sides conservative.
#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do collapse(4) &
        !$omp& map(to: nx, ny, nz, blk%physLow, blk%physHigh, blk%d1x, blk%d1y, blk%d1z, ibm%mu) &
        !$omp& map(tofrom: blk%q, phi) private(i,ip,j,jp,k,kp,b,cf)
#endif
        do b = 1_C_INT, nBlocks
        do k = 1_C_INT, nz
            do j = 1_C_INT, ny
                do i = 1_C_INT, nx
                    ip = i + 1; jp = j + 1; kp = k + 1

                    cf = face_grad(blk%physLow(1,b), i == 1_C_INT, blk%d1x(i,VAR_U,b))
                    if (cf /= 0.0d0) blk%q(i,j,k,VAR_U,b) = blk%q(i,j,k,VAR_U,b) &
                        + (phi(i-1,j,k,b) - phi(i,j,k,b))*cf*ibm%mu(i,j,k,VAR_U,b)
                    cf = face_grad(blk%physLow(2,b), j == 1_C_INT, blk%d1y(j,VAR_V,b))
                    if (cf /= 0.0d0) blk%q(i,j,k,VAR_V,b) = blk%q(i,j,k,VAR_V,b) &
                        + (phi(i,j-1,k,b) - phi(i,j,k,b))*cf*ibm%mu(i,j,k,VAR_V,b)
                    cf = face_grad(blk%physLow(3,b), k == 1_C_INT, blk%d1z(k,VAR_W,b))
                    if (cf /= 0.0d0) blk%q(i,j,k,VAR_W,b) = blk%q(i,j,k,VAR_W,b) &
                        + (phi(i,j,k-1,b) - phi(i,j,k,b))*cf*ibm%mu(i,j,k,VAR_W,b)

                    ! High faces: only the owned 2:1 interface ones.
                    if (i == nx .and. is_interface(blk%physHigh(1,b))) &
                        blk%q(ip,j,k,VAR_U,b) = blk%q(ip,j,k,VAR_U,b) &
                            + (phi(i,j,k,b) - phi(ip,j,k,b)) &
                              *face_grad(blk%physHigh(1,b), .true., blk%d1x(ip,VAR_U,b))*ibm%mu(ip,j,k,VAR_U,b)
                    if (j == ny .and. is_interface(blk%physHigh(2,b))) &
                        blk%q(i,jp,k,VAR_V,b) = blk%q(i,jp,k,VAR_V,b) &
                            + (phi(i,j,k,b) - phi(i,jp,k,b)) &
                              *face_grad(blk%physHigh(2,b), .true., blk%d1y(jp,VAR_V,b))*ibm%mu(i,jp,k,VAR_V,b)
                    if (k == nz .and. is_interface(blk%physHigh(3,b))) &
                        blk%q(i,j,kp,VAR_W,b) = blk%q(i,j,kp,VAR_W,b) &
                            + (phi(i,j,k,b) - phi(i,j,kp,b)) &
                              *face_grad(blk%physHigh(3,b), .true., blk%d1z(kp,VAR_W,b))*ibm%mu(i,j,kp,VAR_W,b)
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

    ! A 2:1 coarse-fine interface face (either orientation).
    pure logical function is_interface(fk)
!$omp declare target
        integer(C_INT), intent(in) :: fk

        is_interface = fk == FACE_COARSE .or. fk == FACE_FINE
    end function is_interface

    ! Pressure-gradient metric for one cell face. atBnd marks a face on the
    ! block boundary (where the face kind fk applies); interior faces always use
    ! the regular metric d1f. A pinned wall/closed face contributes 0. A 2:1
    ! interface face uses the coarse-fine gradient 1/d over the coarse+fine
    ! half-cell gap: on a uniform dyadic grid d = 3/4 h_coarse = 3/2 h_fine, so
    ! 1/d = (2/3) d1_fine = (4/3) d1_coarse -- each cell forms it from its OWN d1
    ! (FACE_COARSE => this block is fine, FACE_FINE => this block is coarse),
    ! giving the same physical 1/d on both sides (the symmetric composite stencil).
    pure real(C_DOUBLE) function face_grad(fk, atBnd, d1f)
!$omp declare target
        integer(C_INT), intent(in) :: fk
        logical, intent(in) :: atBnd
        real(C_DOUBLE), intent(in) :: d1f

        if (.not. atBnd) then
            face_grad = d1f
        else if (fk == FACE_PHYS .or. fk == FACE_CLOSED) then
            face_grad = 0.0d0
        else if (fk == FACE_COARSE) then
            face_grad = (2.0d0/3.0d0)*d1f
        else if (fk == FACE_FINE) then
            face_grad = (4.0d0/3.0d0)*d1f
        else
            face_grad = d1f
        end if
    end function face_grad

end module pressure_solver
