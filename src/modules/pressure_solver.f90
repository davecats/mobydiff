module pressure_solver
    use, intrinsic :: iso_c_binding
    use :: init, only: dns_type, VAR_U, VAR_V, VAR_W, VAR_P
    use :: blocks, only: block_set_type, FACE_PHYS, FACE_CLOSED, FACE_COARSE, FACE_FINE
    use :: ibmm, only: ibm_type
    use :: boundary, only: boundary_type, apply_bc, apply_scalar_bc, &
        boundary_face_id, NFACES, PATCH_OUTLET, SCALAR_BC_NONE, SCALAR_BC_MIRROR
    use :: comm, only: comm_type, exchange_halos, exchange_scalar_halos, sync_divergence_halos
    use :: profiling, only: prof_tic, prof_toc, proj_prof, &
        PROF_SWEEP, PROF_APPLY, PROF_PHI_EXCHANGE, PROF_PROJ_VEL_EXCHANGE, &
        PROF_PROJ_BC, PROF_PROJ_SETUP

    implicit none

    private
    public :: pressure_solver_type, init_pressure_solver, pressure_projection
    public :: PRESSURE_JACOBI, PRESSURE_REDBLACK

    ! Pressure smoother family ([pressure] solver): the current damped-Jacobi /
    ! Chebyshev-Jacobi projection (default), or the original single-level
    ! red-black SOR restored as a selectable alternative. They are mutually
    ! exclusive with the 2:1 interface -- red-black is single-level only.
    integer(C_INT), parameter :: PRESSURE_JACOBI = 0_C_INT
    integer(C_INT), parameter :: PRESSURE_REDBLACK = 1_C_INT

    type :: pressure_solver_type
        integer(C_INT) :: nIter=3
        integer(C_INT) :: method=PRESSURE_JACOBI
        ! Damped-Jacobi relaxation factor. The pressure projection is now a
        ! simple (un-coloured) damped Jacobi iteration; for the Poisson-like
        ! projection operator the high-frequency (checkerboard) mode forces
        ! the factor strictly below 1 (omega=1 is marginally unstable). 0.8 is
        ! a safe default; the config key is still "sor".
        real(C_DOUBLE) :: omega=0.8d0
        ! Chebyshev-Jacobi acceleration ([pressure] accel = chebyshev):
        ! a Chebyshev semi-iteration over the diagonal(Jacobi)-
        ! preconditioned projection operator, whose spectrum is bounded in
        ! [chebLmin, chebLmax]. lmax~2 by Gershgorin regardless of grid
        ! stretching / 2:1 interface; lmin is the condition-number knob,
        ! auto-set to the lowest-mode eigenvalue (2/3)sin^2(pi/N) of the
        ! Jacobi-Poisson on the base grid. A non-positive bound means "auto".
        ! Off by default (plain damped Jacobi).
        logical :: cheb=.false.
        real(C_DOUBLE) :: chebLmin=-1.0d0, chebLmax=-1.0d0
    end type pressure_solver_type

    ! Pressure-increment buffer for the Jacobi sweep. Unlike the in-place
    ! red-black SOR, Jacobi computes every cell's increment from the SAME frozen
    ! divergence, then applies pressure and the velocity-face corrections in
    ! separate race-free passes. phi shares blk%q's spatial bounds (0:nb+1) so
    ! exchange_scalar_halos can fill the halo layer the face corrections read.
    real(C_DOUBLE), allocatable :: phi(:,:,:,:)

    ! Chebyshev search increment delta_k = alpha_k z_k + gamma_k delta_{k-1}
    ! (the actual per-iteration update applied to p and the velocity faces, so
    ! the apply step is identical to Jacobi's -- it just adds delta instead of
    ! phi). Allocated only when Chebyshev is enabled.
    real(C_DOUBLE), allocatable :: delta(:,:,:,:)

contains

    subroutine init_pressure_solver(ps, dns, bc, has_terminal)
        type(pressure_solver_type), intent(inout) :: ps
        type(dns_type), intent(in) :: dns
        type(boundary_type), intent(in) :: bc
        logical, intent(in), optional :: has_terminal

        character(len=32) :: env
        real(C_DOUBLE) :: pi
        integer(C_INT) :: nMax
        integer :: dir

        ! Damped Jacobi needs no red-black colouring, so the even-global-size
        ! restriction of the red-black scheme no longer applies.
        ! [pressure] accel = chebyshev sets ps%cheb in config; the eigenvalue
        ! bounds chebLmin/chebLmax are auto-derived below (config default -1).

        ! Auto eigenvalue bounds (non-positive => auto): lmax = 2 (Gershgorin,
        ! holds for any Jacobi-preconditioned Poisson with zero interior row
        ! sum -- stretched grids and the 2:1 interface included). lmin = the
        ! Jacobi-Poisson lowest-mode eigenvalue (2/3) sin^2(pi/N_max) on the
        ! base global grid; the lowest mode is domain-scale so refinement does
        ! not change it.
        if (ps%cheb) then
            if (ps%chebLmax <= 0.0d0) ps%chebLmax = 2.0d0
            if (ps%chebLmin <= 0.0d0) then
                pi = 4.0d0*atan(1.0d0)
                nMax = maxval(dns%globalSize)
                ps%chebLmin = (2.0d0/3.0d0)*sin(pi/real(nMax, C_DOUBLE))**2
            end if
            if (present(has_terminal)) then
                if (has_terminal) print '(a,es12.4,a,es12.4)', &
                    " Chebyshev-Jacobi projection: lmin=", ps%chebLmin, " lmax=", ps%chebLmax
            end if
        end if

        ! Red-black SOR: kept mutually exclusive with the 2:1 interface (which it
        ! predates and never handled) and with Chebyshev acceleration; and, like
        ! the original scheme, it needs even global sizes in periodic directions
        ! for a consistent colouring.
        if (ps%method == PRESSURE_REDBLACK) then
            if (ps%cheb) error stop &
                "[pressure] solver = redblack is incompatible with accel = chebyshev"
            if (dns%block_refine_levels > 1_C_INT .or. dns%block_refine_body .or. &
                dns%block_refine_nboxes > 0_C_INT) error stop &
                "[pressure] solver = redblack does not support 2:1 block refinement"
            do dir = 1, 3
                if (bc%isPeriodic(dir) .and. mod(dns%globalSize(dir), 2_C_INT) /= 0_C_INT) &
                    error stop "red-black pressure solver requires even global sizes in periodic directions"
            end do
        end if
    end subroutine init_pressure_solver

    subroutine pressure_projection(ps, blk, dns, dt_gamma, ibm, bc, c)
        type(pressure_solver_type), intent(in) :: ps
        type(block_set_type), intent(inout) :: blk
        type(dns_type), intent(in) :: dns
        real(C_DOUBLE), intent(in) :: dt_gamma
        type(ibm_type), intent(in) :: ibm
        type(boundary_type), intent(in) :: bc
        type(comm_type), intent(inout) :: c

        integer(C_INT) :: iIter, dir
        real(C_DOUBLE) :: omega
        real(C_DOUBLE) :: dd, cc, alpha, alphaPrev, beta, gamma
        logical(C_BOOL) :: outLow(3), outHigh(3), refd(3)
        logical :: anyOutlet
        integer(C_INT) :: phiMode(NFACES)
        real(C_DOUBLE) :: t0

        if (ps%method == PRESSURE_REDBLACK) then
            call redblack_projection(ps, blk, dt_gamma, ibm, bc, c)
            return
        end if

        t0 = prof_tic()
        call allocate_phi(blk)
        if (ps%cheb) call allocate_delta(blk)

        ! Per-direction refinement flags ([blocks] refine_dims): a 2:1 face
        ! normal to an unrefined direction uses the uniform gradient metric
        ! (face_grad), all true in xyz mode.
        refd = blk%refMask == 1_C_INT

        ! Dirichlet-pressure outlet faces (declared [boundary] _patch = outlet).
        ! The flags are per DOMAIN face; inside the kernels they act only where
        ! the block's own face is FACE_PHYS (i.e. the block face IS that domain
        ! face). phi's ghost at an outlet face is mirrored after every halo
        ! exchange (phi_face = 0: the increment of a HELD outlet pressure), so
        ! the projection owns and corrects the outlet face flux -- this is what
        ! keeps the Poisson problem solvable under net prescribed inflow.
        do dir = 1, 3
            outLow(dir) = bc%facePatchType(boundary_face_id(int(dir), 0)) == PATCH_OUTLET
            outHigh(dir) = bc%facePatchType(boundary_face_id(int(dir), 1)) == PATCH_OUTLET
            phiMode(boundary_face_id(int(dir), 0)) = &
                merge(SCALAR_BC_MIRROR, SCALAR_BC_NONE, outLow(dir))
            phiMode(boundary_face_id(int(dir), 1)) = &
                merge(SCALAR_BC_MIRROR, SCALAR_BC_NONE, outHigh(dir))
        end do
        anyOutlet = any(outLow) .or. any(outHigh)

        ! Chebyshev semi-iteration coefficients over [lmin, lmax].
        dd = 0.5d0*(ps%chebLmax + ps%chebLmin)
        cc = 0.5d0*(ps%chebLmax - ps%chebLmin)
        omega = merge(1.0d0, ps%omega, ps%cheb)
        alpha = 0.0d0

        ! Projection. Each iteration: (1) preconditioned residual into phi
        ! (omega*(-div/denom); Jacobi omega=sor folds the damping in, Chebyshev
        ! omega=1 leaves z and the Chebyshev coefficients do the damping), (2)
        ! for Chebyshev combine into the increment delta_k = alpha z + gamma
        ! delta_{k-1} (phi := delta_k), (3) exchange phi's halos (the 2:1
        ! interface-row restrict keeps the interface corrections conservative),
        ! (4) apply phi to the pressure and velocity faces, (5) refresh the
        ! velocity (+ pressure on the last iteration) halos for the next
        ! divergence.
        call prof_toc(proj_prof, PROF_PROJ_SETUP, t0)

        do iIter = 1_C_INT, ps%nIter
            t0 = prof_tic()
            call jacobi_compute_phi(blk, ibm, omega, outLow, outHigh, refd)
            if (ps%cheb) then
                if (iIter == 1_C_INT) then
                    alpha = 1.0d0/dd
                    gamma = 0.0d0
                else
                    alphaPrev = alpha
                    beta = (cc*alphaPrev*0.5d0)**2
                    alpha = 1.0d0/(dd - beta/alphaPrev)
                    gamma = alpha*beta/alphaPrev
                end if
                call cheb_combine(blk, alpha, gamma)
            end if
            call prof_toc(proj_prof, PROF_SWEEP, t0)
            t0 = prof_tic()
            call exchange_scalar_halos(c, phi, blk, ifaceRow=.true.)
            ! Re-mirror the outlet phi ghosts EVERY iteration: the exchange's
            ! tangential extension can write physical halos, so do not rely on
            ! them staying zero.
            if (anyOutlet) call apply_scalar_bc(blk, bc, phi, phiMode)
            call prof_toc(proj_prof, PROF_PHI_EXCHANGE, t0)
            t0 = prof_tic()
            call jacobi_apply(ps, blk, dt_gamma, ibm, outLow, outHigh, refd)
            call prof_toc(proj_prof, PROF_APPLY, t0)
            t0 = prof_tic()
            call apply_bc(blk, bc)
            call prof_toc(proj_prof, PROF_PROJ_BC, t0)
            t0 = prof_tic()
            if (iIter == ps%nIter) then
                ! Last iteration: the full shell, all four variables, with the
                ! cross-level transfers -- this is what the next substage reads.
                call exchange_halos(c, blk, [VAR_U, VAR_V, VAR_W, VAR_P])
            else if (c%nPeers == 0) then
                ! Between iterations only the divergence stencil reads the
                ! velocity halo, so refresh just q(nb+1) per dim, normal
                ! component (see comm.f90 sync_divergence_halos). Single-rank
                ! only: with peers the same planes would have to come over MPI,
                ! and the message is the whole point of the saving, so that
                ! needs the entry list partitioned first.
                call sync_divergence_halos(c, blk)
            else
                call exchange_halos(c, blk, [VAR_U, VAR_V, VAR_W], interp=.false.)
            end if
            call prof_toc(proj_prof, PROF_PROJ_VEL_EXCHANGE, t0)
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

    ! Chebyshev increment buffer (same bounds as phi), allocated on first use.
    subroutine allocate_delta(blk)
        type(block_set_type), intent(in) :: blk
        if (.not. allocated(delta)) then
            allocate(delta(lbound(phi,1):ubound(phi,1), lbound(phi,2):ubound(phi,2), &
                           lbound(phi,3):ubound(phi,3), blk%nBlocks))
            delta = 0.0d0
#ifdef USE_OPENMP_OFFLOAD
            !$omp target enter data map(to: delta)
#endif
        end if
    end subroutine allocate_delta

    ! Chebyshev combine: phi := alpha*phi + gamma*delta (phi currently holds the
    ! preconditioned residual z), then delta := phi (the new increment delta_k),
    ! over interior cells. After this phi == delta_k and jacobi_apply adds it
    ! exactly as it would the Jacobi increment.
    subroutine cheb_combine(blk, alpha, gamma)
        type(block_set_type), intent(in) :: blk
        real(C_DOUBLE), intent(in) :: alpha, gamma
        integer(C_INT) :: i, j, k, b, nx, ny, nz
        real(C_DOUBLE) :: s
        nx = blk%nb(1); ny = blk%nb(2); nz = blk%nb(3)
#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do collapse(4) &
        !$omp& map(to: alpha, gamma, nx, ny, nz) map(tofrom: phi, delta) private(i,j,k,b,s)
#endif
        do b = 1_C_INT, blk%nBlocks
        do k = 1_C_INT, nz
            do j = 1_C_INT, ny
                do i = 1_C_INT, nx
                    s = alpha*phi(i,j,k,b) + gamma*delta(i,j,k,b)
                    phi(i,j,k,b) = s
                    delta(i,j,k,b) = s
                end do
            end do
        end do
        end do
#ifdef USE_OPENMP_OFFLOAD
        !$omp end target teams distribute parallel do
#endif
    end subroutine cheb_combine

    ! Jacobi step 1: phi(i,j,k) = -omega * div / denom for every interior cell,
    ! from the frozen velocity. denom is the projection diagonal (sum of the
    ! cell's non-pinned face metrics); pinned faces (walls, closed faces) leave
    ! the diagonal exactly as in the red-black scheme.
    ! omega: the damping factor. Plain Jacobi passes ps%omega (the increment IS
    ! phi = -omega*div/denom); Chebyshev passes 1.0 so phi holds the pure
    ! diagonal-preconditioned residual z = -div/denom (the damping then comes
    ! from the Chebyshev coefficients).
    subroutine jacobi_compute_phi(blk, ibm, omega, outLow, outHigh, refd)
        type(block_set_type), intent(inout) :: blk
        type(ibm_type), intent(in) :: ibm
        real(C_DOUBLE), intent(in) :: omega
        logical(C_BOOL), intent(in) :: outLow(3), outHigh(3), refd(3)

        real(C_DOUBLE) :: denom, div
        real(C_DOUBLE) :: mu_u_i, mu_u_ip, mu_v_j, mu_v_jp, mu_w_k, mu_w_kp
        integer(C_INT) :: i, ip, j, jp, k, kp, b, nBlocks, nx, ny, nz

        nx = blk%nb(1); ny = blk%nb(2); nz = blk%nb(3)
        nBlocks = blk%nBlocks

#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do collapse(4) &
        !$omp& map(to: omega, nx, ny, nz, outLow(1:3), outHigh(1:3), refd(1:3), &
        !$omp& blk%physLow, blk%physHigh, &
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
                    ! and the regular 1/h otherwise. An OUTLET physical face counts
                    ! 2*d1f (face_grad_denom): its phi ghost is MIRRORED, so the
                    ! flux sensitivity per unit phi_i doubles -- the matching
                    ! correction in jacobi_apply uses the regular d1f against the
                    ! mirrored ghost. Keep the pair consistent: SPD lives there.
                    denom = (face_grad_denom(blk%physLow(1,b), i == 1_C_INT, blk%d1x(i,VAR_U,b), outLow(1), refd(1))*mu_u_i &
                           + face_grad_denom(blk%physHigh(1,b), i == nx, blk%d1x(ip,VAR_U,b), outHigh(1), refd(1))*mu_u_ip)*blk%d1x(i,VAR_P,b) &
                          + (face_grad_denom(blk%physLow(2,b), j == 1_C_INT, blk%d1y(j,VAR_V,b), outLow(2), refd(2))*mu_v_j &
                           + face_grad_denom(blk%physHigh(2,b), j == ny, blk%d1y(jp,VAR_V,b), outHigh(2), refd(2))*mu_v_jp)*blk%d1y(j,VAR_P,b) &
                          + (face_grad_denom(blk%physLow(3,b), k == 1_C_INT, blk%d1z(k,VAR_W,b), outLow(3), refd(3))*mu_w_k &
                           + face_grad_denom(blk%physHigh(3,b), k == nz, blk%d1z(kp,VAR_W,b), outHigh(3), refd(3))*mu_w_kp)*blk%d1z(k,VAR_P,b)

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
    subroutine jacobi_apply(ps, blk, dt_gamma, ibm, outLow, outHigh, refd)
        type(pressure_solver_type), intent(in) :: ps
        type(block_set_type), intent(inout) :: blk
        real(C_DOUBLE), intent(in) :: dt_gamma
        type(ibm_type), intent(in) :: ibm
        logical(C_BOOL), intent(in) :: outLow(3), outHigh(3), refd(3)

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
        ! An OUTLET physical face (either side) is also corrected, with the
        ! REGULAR d1f against the MIRRORED phi ghost (face_grad_corr): e.g.
        ! q(nx+1) += (phi(nx) - (-phi(nx)))*d1f*mu = 2*phi(nx)*d1f*mu -- the
        ! half-cell Dirichlet gradient whose 2*d1f sensitivity the denominator
        ! already counts. The predictor never writes the outlet face, so this
        ! correction (+ the initial value) is its entire evolution: the standard
        ! do-nothing Dirichlet-pressure outlet.
#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do collapse(4) &
        !$omp& map(to: nx, ny, nz, outLow(1:3), outHigh(1:3), refd(1:3), &
        !$omp& blk%physLow, blk%physHigh, blk%d1x, blk%d1y, blk%d1z, ibm%mu) &
        !$omp& map(tofrom: blk%q, phi) private(i,ip,j,jp,k,kp,b,cf)
#endif
        do b = 1_C_INT, nBlocks
        do k = 1_C_INT, nz
            do j = 1_C_INT, ny
                do i = 1_C_INT, nx
                    ip = i + 1; jp = j + 1; kp = k + 1

                    cf = face_grad_corr(blk%physLow(1,b), i == 1_C_INT, blk%d1x(i,VAR_U,b), outLow(1), refd(1))
                    if (cf /= 0.0d0) blk%q(i,j,k,VAR_U,b) = blk%q(i,j,k,VAR_U,b) &
                        + (phi(i-1,j,k,b) - phi(i,j,k,b))*cf*ibm%mu(i,j,k,VAR_U,b)
                    cf = face_grad_corr(blk%physLow(2,b), j == 1_C_INT, blk%d1y(j,VAR_V,b), outLow(2), refd(2))
                    if (cf /= 0.0d0) blk%q(i,j,k,VAR_V,b) = blk%q(i,j,k,VAR_V,b) &
                        + (phi(i,j-1,k,b) - phi(i,j,k,b))*cf*ibm%mu(i,j,k,VAR_V,b)
                    cf = face_grad_corr(blk%physLow(3,b), k == 1_C_INT, blk%d1z(k,VAR_W,b), outLow(3), refd(3))
                    if (cf /= 0.0d0) blk%q(i,j,k,VAR_W,b) = blk%q(i,j,k,VAR_W,b) &
                        + (phi(i,j,k-1,b) - phi(i,j,k,b))*cf*ibm%mu(i,j,k,VAR_W,b)

                    ! High faces: the owned 2:1 interface ones, plus an outlet
                    ! physical face (FACE_PHYS + declared outlet).
                    if (i == nx .and. (is_interface(blk%physHigh(1,b)) .or. &
                        (outHigh(1) .and. blk%physHigh(1,b) == FACE_PHYS))) &
                        blk%q(ip,j,k,VAR_U,b) = blk%q(ip,j,k,VAR_U,b) &
                            + (phi(i,j,k,b) - phi(ip,j,k,b)) &
                              *face_grad_corr(blk%physHigh(1,b), .true., blk%d1x(ip,VAR_U,b), outHigh(1), refd(1))*ibm%mu(ip,j,k,VAR_U,b)
                    if (j == ny .and. (is_interface(blk%physHigh(2,b)) .or. &
                        (outHigh(2) .and. blk%physHigh(2,b) == FACE_PHYS))) &
                        blk%q(i,jp,k,VAR_V,b) = blk%q(i,jp,k,VAR_V,b) &
                            + (phi(i,j,k,b) - phi(i,jp,k,b)) &
                              *face_grad_corr(blk%physHigh(2,b), .true., blk%d1y(jp,VAR_V,b), outHigh(2), refd(2))*ibm%mu(i,jp,k,VAR_V,b)
                    if (k == nz .and. (is_interface(blk%physHigh(3,b)) .or. &
                        (outHigh(3) .and. blk%physHigh(3,b) == FACE_PHYS))) &
                        blk%q(i,j,kp,VAR_W,b) = blk%q(i,j,kp,VAR_W,b) &
                            + (phi(i,j,k,b) - phi(i,j,kp,b)) &
                              *face_grad_corr(blk%physHigh(3,b), .true., blk%d1z(kp,VAR_W,b), outHigh(3), refd(3))*ibm%mu(i,j,kp,VAR_W,b)
                end do
            end do
        end do
        end do
#ifdef USE_OPENMP_OFFLOAD
        !$omp end target teams distribute parallel do
#endif
    end subroutine jacobi_apply

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
    ! refined = the face's normal direction halves per level ([blocks]
    ! refine_dims; always true in xyz mode). A 2:1 face normal to an
    ! UNREFINED direction shares its node line on both sides (the
    ! xz-quadtree y face): the normal gradient is the UNIFORM stencil d1f
    ! -- only the tangential resolution (i.e. the phi ghost content:
    ! 4-sample restrict / injection) differs. Denominator and correction
    ! then both count the SAME d1f: the SPD pair.
    pure real(C_DOUBLE) function face_grad(fk, atBnd, d1f, refined)
!$omp declare target
        integer(C_INT), intent(in) :: fk
        logical, intent(in) :: atBnd
        real(C_DOUBLE), intent(in) :: d1f
        logical(C_BOOL), intent(in) :: refined

        if (.not. atBnd) then
            face_grad = d1f
        else if (fk == FACE_PHYS .or. fk == FACE_CLOSED) then
            face_grad = 0.0d0
        else if (.not. refined) then
            face_grad = d1f
        else if (fk == FACE_COARSE) then
            face_grad = (2.0d0/3.0d0)*d1f
        else if (fk == FACE_FINE) then
            face_grad = (4.0d0/3.0d0)*d1f
        else
            face_grad = d1f
        end if
    end function face_grad

    ! Outlet-aware metrics: the Dirichlet-pressure outlet is the ONE face kind
    ! where the Jacobi DENOMINATOR and the velocity-face CORRECTION need
    ! different values, because the correction reads the MIRRORED phi ghost
    ! (phi_ghost = -phi_i, i.e. phi held at 0 ON the face):
    !   denominator:  2*d1f  (d/dphi_i of the corrected face flux),
    !   correction:     d1f  (applied to phi_i - phi_ghost = 2*phi_i).
    ! The product is the same half-cell Dirichlet gradient on both sides of
    ! the operator -- change one without the other and the projection loses
    ! SPD at exactly this face row. Every other face kind defers to face_grad
    ! (whose branches are validated and locked).
    pure real(C_DOUBLE) function face_grad_denom(fk, atBnd, d1f, outlet, refined)
!$omp declare target
        integer(C_INT), intent(in) :: fk
        logical, intent(in) :: atBnd
        real(C_DOUBLE), intent(in) :: d1f
        logical(C_BOOL), intent(in) :: outlet, refined

        if (atBnd .and. outlet .and. fk == FACE_PHYS) then
            face_grad_denom = 2.0d0*d1f
        else
            face_grad_denom = face_grad(fk, atBnd, d1f, refined)
        end if
    end function face_grad_denom

    pure real(C_DOUBLE) function face_grad_corr(fk, atBnd, d1f, outlet, refined)
!$omp declare target
        integer(C_INT), intent(in) :: fk
        logical, intent(in) :: atBnd
        real(C_DOUBLE), intent(in) :: d1f
        logical(C_BOOL), intent(in) :: outlet, refined

        if (atBnd .and. outlet .and. fk == FACE_PHYS) then
            face_grad_corr = d1f
        else
            face_grad_corr = face_grad(fk, atBnd, d1f, refined)
        end if
    end function face_grad_corr

    ! -------------------------------------------------------------------------
    ! Red-black SOR projection ([pressure] solver = redblack): the original
    ! coloured Gauss-Seidel scheme, restored as a selectable single-level
    ! alternative to damped Jacobi. In-place over-relaxed sweeps; the same-level
    ! halo exchange refreshes the velocity each colour (pressure only at the
    ! end). It reuses the SPD diagonal metric face_grad_denom for BOTH the
    ! divergence denominator AND the in-place face corrections: applying phi in
    ! place (no ghost read) means the correction metric equals the diagonal
    ! metric (d1f interior, 0 pinned, 2*d1f for the Dirichlet-p outlet), so the
    ! outlet is handled exactly as in the Jacobi path and, on outlet-free grids,
    ! the sweep is identical to the original red-black scheme.
    subroutine redblack_projection(ps, blk, dt_gamma, ibm, bc, c)
        type(pressure_solver_type), intent(in) :: ps
        type(block_set_type), intent(inout) :: blk
        real(C_DOUBLE), intent(in) :: dt_gamma
        type(ibm_type), intent(in) :: ibm
        type(boundary_type), intent(in) :: bc
        type(comm_type), intent(inout) :: c

        integer(C_INT) :: iIter, color, dir
        logical(C_BOOL) :: outLow(3), outHigh(3), refd(3)
        real(C_DOUBLE) :: t0

        t0 = prof_tic()
        refd = blk%refMask == 1_C_INT
        do dir = 1, 3
            outLow(dir)  = bc%facePatchType(boundary_face_id(int(dir), 0)) == PATCH_OUTLET
            outHigh(dir) = bc%facePatchType(boundary_face_id(int(dir), 1)) == PATCH_OUTLET
        end do
        call prof_toc(proj_prof, PROF_PROJ_SETUP, t0)

        ! The coupled sweep updates pressure AND velocity in one kernel, so
        ! proj_timing's `apply` and `phi_exchange` stay zero here: the red-black
        ! path spends everything in `sweep` and the velocity exchange. That
        ! asymmetry against the Jacobi path is real, not missing instrumentation.
        do iIter = 1_C_INT, ps%nIter
            do color = 1_C_INT, 0_C_INT, -1_C_INT
                t0 = prof_tic()
                call redblack_sweep(ps, blk, dt_gamma, ibm, color, outLow, outHigh, refd)
                call prof_toc(proj_prof, PROF_SWEEP, t0)
                t0 = prof_tic()
                call apply_bc(blk, bc)
                call prof_toc(proj_prof, PROF_PROJ_BC, t0)
                t0 = prof_tic()
                if (iIter == ps%nIter .and. color == 0_C_INT) then
                    call exchange_halos(c, blk, [VAR_U, VAR_V, VAR_W, VAR_P])
                else
                    call exchange_halos(c, blk, [VAR_U, VAR_V, VAR_W], interp=.false.)
                end if
                call prof_toc(proj_prof, PROF_PROJ_VEL_EXCHANGE, t0)
            end do
        end do
    end subroutine redblack_projection

    ! One red-black sweep of `color`. Each block sweeps from its halo layer (0)
    ! redundantly with the owning neighbour except on physical boundaries -- the
    ! rank-level red-black trick one level down -- so results are independent of
    ! the block/rank layout; parity is anchored to the global index space via
    ! the block origin. `sor` = ps%omega (the config "sor"; over-relaxation > 1
    ! is allowed here, unlike the Jacobi damping).
    subroutine redblack_sweep(ps, blk, dt_gamma, ibm, color, outLow, outHigh, refd)
        type(pressure_solver_type), intent(in) :: ps
        type(block_set_type), intent(inout) :: blk
        real(C_DOUBLE), intent(in) :: dt_gamma
        type(ibm_type), intent(in) :: ibm
        integer(C_INT), intent(in) :: color
        logical(C_BOOL), intent(in) :: outLow(3), outHigh(3), refd(3)

        real(C_DOUBLE) :: phi, denom, idt, sor, div
        real(C_DOUBLE) :: gLo1, gHi1, gLo2, gHi2, gLo3, gHi3
        real(C_DOUBLE) :: mu_u_i, mu_u_ip, mu_v_j, mu_v_jp, mu_w_k, mu_w_kp
        integer(C_INT) :: i, ip, j, jp, k, kp, b, nBlocks, nLowerHaloDirections, iColor, nColorX
        integer(C_INT) :: iLo, jLo, kLo, colorOffset
        integer(C_INT) :: hi(1:3)

        hi = blk%nb(1:3)
        sor = ps%omega
        nColorX = (hi(1) + 2_C_INT)/2_C_INT
        nBlocks = blk%nBlocks
        idt = 1.0_C_DOUBLE/dt_gamma

#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do collapse(4) &
        !$omp& map(to: color, nColorX, sor, idt, dt_gamma, hi(1:3), &
        !$omp& outLow(1:3), outHigh(1:3), refd(1:3), &
        !$omp& blk%origin, blk%physLow, blk%physHigh, &
        !$omp& blk%d1x, blk%d1y, blk%d1z, ibm%mu) &
        !$omp& map(tofrom: blk%q) &
        !$omp& private(i,ip,j,jp,k,kp,b,iColor,iLo,jLo,kLo,colorOffset, &
        !$omp& phi,denom,div,nLowerHaloDirections, gLo1,gHi1,gLo2,gHi2,gLo3,gHi3, &
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

                    ip = i + 1; jp = j + 1; kp = k + 1

                    mu_u_i  = ibm%mu(i,j,k,VAR_U,b);  mu_u_ip = ibm%mu(ip,j,k,VAR_U,b)
                    mu_v_j  = ibm%mu(i,j,k,VAR_V,b);  mu_v_jp = ibm%mu(i,jp,k,VAR_V,b)
                    mu_w_k  = ibm%mu(i,j,k,VAR_W,b);  mu_w_kp = ibm%mu(i,j,kp,VAR_W,b)

                    ! SPD face metric (0 pinned, 2*d1f outlet, d1f otherwise); the
                    ! in-place correction reuses the same value as the diagonal.
                    gLo1 = face_grad_denom(blk%physLow(1,b),  i == 1_C_INT, blk%d1x(i,VAR_U,b),  outLow(1),  refd(1))
                    gHi1 = face_grad_denom(blk%physHigh(1,b), i == hi(1),   blk%d1x(ip,VAR_U,b), outHigh(1), refd(1))
                    gLo2 = face_grad_denom(blk%physLow(2,b),  j == 1_C_INT, blk%d1y(j,VAR_V,b),  outLow(2),  refd(2))
                    gHi2 = face_grad_denom(blk%physHigh(2,b), j == hi(2),   blk%d1y(jp,VAR_V,b), outHigh(2), refd(2))
                    gLo3 = face_grad_denom(blk%physLow(3,b),  k == 1_C_INT, blk%d1z(k,VAR_W,b),  outLow(3),  refd(3))
                    gHi3 = face_grad_denom(blk%physHigh(3,b), k == hi(3),   blk%d1z(kp,VAR_W,b), outHigh(3), refd(3))

                    denom = (gLo1*mu_u_i + gHi1*mu_u_ip)*blk%d1x(i,VAR_P,b) &
                          + (gLo2*mu_v_j + gHi2*mu_v_jp)*blk%d1y(j,VAR_P,b) &
                          + (gLo3*mu_w_k + gHi3*mu_w_kp)*blk%d1z(k,VAR_P,b)

                    div = (blk%q(ip,j,k,VAR_U,b)-blk%q(i,j,k,VAR_U,b))*blk%d1x(i,VAR_P,b) &
                        + (blk%q(i,jp,k,VAR_V,b)-blk%q(i,j,k,VAR_V,b))*blk%d1y(j,VAR_P,b) &
                        + (blk%q(i,j,kp,VAR_W,b)-blk%q(i,j,k,VAR_W,b))*blk%d1z(k,VAR_P,b)

                    phi = -sor*div/denom

                    blk%q(i,j,k,VAR_P,b) = blk%q(i,j,k,VAR_P,b) + phi*idt

                    ! Symmetric correction: every non-pinned face of the cell is
                    ! moved by phi*metric*mu; the neighbour of the other colour
                    ! moves the shared face from its side, so over both colours a
                    ! face sees (phi_below - phi_above)*d1f -- the same operator as
                    ! the Jacobi apply, here in place.
                    blk%q(i,j,k,VAR_U,b)  = blk%q(i,j,k,VAR_U,b)  - phi*gLo1*mu_u_i
                    blk%q(ip,j,k,VAR_U,b) = blk%q(ip,j,k,VAR_U,b) + phi*gHi1*mu_u_ip
                    blk%q(i,j,k,VAR_V,b)  = blk%q(i,j,k,VAR_V,b)  - phi*gLo2*mu_v_j
                    blk%q(i,jp,k,VAR_V,b) = blk%q(i,jp,k,VAR_V,b) + phi*gHi2*mu_v_jp
                    blk%q(i,j,k,VAR_W,b)  = blk%q(i,j,k,VAR_W,b)  - phi*gLo3*mu_w_k
                    blk%q(i,j,kp,VAR_W,b) = blk%q(i,j,kp,VAR_W,b) + phi*gHi3*mu_w_kp
                END DO
            END DO
        END DO
        END DO
#ifdef USE_OPENMP_OFFLOAD
        !$omp end target teams distribute parallel do
#endif
    end subroutine redblack_sweep

end module pressure_solver
