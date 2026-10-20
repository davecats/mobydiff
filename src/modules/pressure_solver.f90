module pressure_solver
    use, intrinsic :: iso_c_binding
    use :: init, only: dns_type, VAR_U, VAR_V, VAR_W, VAR_P
    use :: blocks, only: block_set_type, FACE_PHYS, FACE_CLOSED, FACE_COARSE, FACE_FINE
    use :: ibmm, only: ibm_type, ibm_body_blocks
    use :: boundary, only: boundary_type, apply_bc, boundary_face_id, PATCH_OUTLET
    use :: comm, only: comm_type, exchange_halos, exchange_scalar_halos, sync_divergence_halos, &
        comm_allreduce_max_int
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
    ! Its ghost plane behind a Dirichlet-pressure OUTLET face holds phi_face
    ! = 0 (the increment of a HELD pressure) for the whole run: zeroed at
    ! allocation and written by nothing -- the proof is in pressure_projection.
    real(C_DOUBLE), allocatable :: phi(:,:,:,:)

    ! Chebyshev search increment delta_k = alpha_k z_k + gamma_k delta_{k-1}
    ! (the actual per-iteration update applied to p and the velocity faces, so
    ! the apply step is identical to Jacobi's -- it just adds delta instead of
    ! phi). Allocated only when Chebyshev is enabled.
    real(C_DOUBLE), allocatable :: delta(:,:,:,:)

    ! RECIPROCAL of the projection diagonal, one value per interior cell.
    !
    ! The diagonal is a function of the face kinds, the grid metrics and ibm%mu
    ! ONLY -- and update_ibm_mu runs once per RK substage, before the projection
    ! -- so it is identical across all niter iterations and is formed at most
    ! once per substage instead of once per iteration (and, since 2026-09-17,
    ! only on the blocks whose mu changes: see below). That alone removed six
    ! branchy face_grad calls, six mu loads and ~11 fp64 operations from
    ! the per-iteration kernel (measured: -26 % on the sweep).
    !
    ! It is stored INVERTED so the iteration kernel multiplies instead of
    ! divides. That is a deliberate NUMERICS change, not a refactor: a*(1/b)
    ! differs from a/b by up to 1 ulp, so the projection is NOT bit-identical to
    ! code predating 2026-09-02, and "bit-exact vs the previous binary" is no
    ! longer available as a gate for it. What it buys, measured with ncu: the
    ! fp64 division was 50.5 % of jacobi_compute_phi (562.94 -> 278.85 us), and
    ! removing it moves the kernel from SM-bound to DRAM-bound at 89 % of peak,
    ! i.e. essentially optimal. fp64 divide is a long software sequence on a GPU
    ! whose fp64 rate is a fraction of its fp32 rate.
    !
    ! Why a last-bit change here is safe: phi is an iterative INCREMENT and this
    ! is its preconditioner diagonal, so the perturbation is many orders below
    ! the projection's own residual after niter iterations. The operator stays
    ! SPD -- the diagonal is strictly positive, only its final bit moves -- and
    ! the correction metric (the same face_grad values, cfLow/cfHigh) is not
    ! involved. A degenerate all-faces-pinned cell produced
    ! Inf/NaN before (division by an exactly zero diagonal) and still does, now
    ! via 1/0; the failure mode is unchanged in kind.
    real(C_DOUBLE), allocatable :: rdenom(:,:,:,:)
    ! The blocks compute_rdenom still has to visit -- see the metric tables.
    ! Unallocated until the first projection, which visits every block; narrowed
    ! afterwards to the blocks holding IBM coefficients. Empty means nothing can
    ! change any more, which is the whole of a body-free case.
    integer, allocatable :: rdenomBlocks(:)
    integer, save :: nRdenomBlocks = 0
    logical, save :: dnsHasTerminal = .false.

    ! STATIC face metric tables, precomputed per (face-normal index, direction,
    ! block) and shared by the two volume kernels of the projection.
    !
    ! Every metric these kernels need at a face is face_grad (face kind,
    ! at-boundary, d1f, outlet, refined) -- a function of the block and the
    ! index along the face NORMAL
    ! only, never of the two tangential indices the collapse(4) loops also span.
    ! Evaluating them inside a kernel therefore recomputed the same branchy
    ! expression nb^2 times per value and, worse, kept
    ! physLow/physHigh/d1x/d1y/d1z/outLow/outHigh/refd live across the whole loop
    ! body. That is what made both kernels register-limited, and cost real time:
    ! jacobi_apply's face kernel went 88 -> 64 registers, 29.6 -> 46.0 %
    ! occupancy, 42.5 -> 60.5 % of peak DRAM and -29.8 % on the clock with its
    ! traffic unchanged to the last digit (measured, ncu, both binaries on one
    ! node: overheadTest/results_apply_registers_2026-09-14.md).
    !
    ! Hoisting is bit-exact by construction -- identical expression, identical
    ! inputs, evaluated once instead of nb^2 times.
    !
    !   cfLow/cfHigh   jacobi_apply's velocity-face CORRECTION metric: the
    !                  diagonal's own face_grad values (ONE metric, see
    !                  face_grad), except that cfLow(idx,d,b) carries the
    !                  interface zeroing (2:1 interface faces are corrected by
    !                  interface_correct, not there) and cfHigh(d,b) is zero
    !                  unless the block's high face in d is an OUTLET physical
    !                  face -- the only high face that kernel owns, so one
    !                  value per block suffices.
    !   dnLow/dnHigh   compute_rdenom's DIAGONAL metric, both faces, every index
    !                  (the high face enters the diagonal at every cell, not
    !                  only at the block boundary -- hence the full index range
    !                  here where cfHigh needs one value).
    !   d1P            the cell's own divergence metric d1?(idx,VAR_P,b), the
    !                  last thing that kept all three d1 arrays live in
    !                  compute_rdenom. It MUST stay a separate factor: the
    !                  diagonal is (dnLow*mu + dnHigh*mu)*d1P, and folding d1P
    !                  into dnLow/dnHigh would distribute the multiply and move
    !                  the last bits.
    !
    ! All of them are STATIC: face kinds come from the leaf table and the metrics
    ! from the node lines, and neither changes during a run. So they are formed
    ! ONCE on the host and mapped once.
    !
    ! rdenom follows ibm%mu, which update_ibm_mu rewrites every substage -- but
    ! only where there is a body to rewrite it for. mu, the penalization factor
    ! of dt*coef (ibm.f90 penal_incr_factor), is EXACTLY 1.0 wherever coef is
    ! zero, whatever dt does, so a block holding no
    ! coefficient has a dt-independent rdenom and needs computing once, like the
    ! tables above. The first projection visits every block; afterwards only the
    ! body blocks (rdenomBlocks, from ibm_body_blocks), which on a body-free
    ! case is none at all. Skipping is bit-exact by construction: same inputs,
    ! same expression, same result -- the recomputation was writing back the
    ! values already there.
    real(C_DOUBLE), allocatable :: cfLow(:,:,:)   ! (1:maxval(nb), 3, nBlocks)
    real(C_DOUBLE), allocatable :: cfHigh(:,:)    ! (3, nBlocks)
    real(C_DOUBLE), allocatable :: dnLow(:,:,:)   ! (1:maxval(nb), 3, nBlocks)
    real(C_DOUBLE), allocatable :: dnHigh(:,:,:)  ! (1:maxval(nb), 3, nBlocks)
    real(C_DOUBLE), allocatable :: d1P(:,:,:)     ! (1:maxval(nb), 3, nBlocks)

    ! The 2:1 interface faces of THIS rank, listed once (init_interface_faces):
    ! per face its block, the direction of its normal, its index along that
    ! direction (1 = the block's low face, nb+1 = its high halo face) and the
    ! correction metric face_grad(kind, .true., d1f, refd(d)) -- a function of
    ! (block, direction, side) alone, the cfLow argument once more. The list
    ! is static like the tables (face kinds from the leaf table, metrics from
    ! the node lines) and is what makes interface_correct ONE launch over the
    ! faces that exist, instead of three launches over every block testing the
    ! face kind at every point (2026-10-21). nIfaces = 0 on a rank that holds no
    ! interface face, and the kernel does not launch there.
    integer(C_INT), allocatable :: ifBlock(:), ifDir(:), ifIdx(:)
    real(C_DOUBLE), allocatable :: ifCoef(:)
    integer(C_INT) :: nIfaces = 0_C_INT
    logical, save :: ifaceListBuilt = .false.

    ! Does this case have any 2:1 interface face? Face kinds come from the leaf
    ! table and never change during a run, so the answer is computed once and
    ! cached. It must be a GLOBAL answer: the red-black path drives a collective
    ! phi exchange from it, and ranks that disagreed would post mismatched
    ! messages.
    logical, save :: ifaceKnown = .false.
    logical, save :: ifaceAny = .false.

contains

    subroutine init_pressure_solver(ps, dns, bc, has_terminal)
        type(pressure_solver_type), intent(inout) :: ps
        type(dns_type), intent(in) :: dns
        type(boundary_type), intent(in) :: bc
        logical, intent(in), optional :: has_terminal

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

        ! Red-black SOR: 2:1 refinement is supported since R1 (the sweep stores
        ! its increment, a per-colour cross-level exchange transmits it and
        ! interface_correct applies it -- the same operator the Jacobi path
        ! uses). Chebyshev stays exclusive: it needs a stationary linear
        ! operator. Like the original scheme it needs even global sizes in
        ! periodic directions for a consistent colouring.
        if (ps%method == PRESSURE_REDBLACK) then
            if (ps%cheb) error stop &
                "[pressure] solver = redblack is incompatible with accel = chebyshev"
            do dir = 1, 3
                if (bc%isPeriodic(dir) .and. mod(dns%globalSize(dir), 2_C_INT) /= 0_C_INT) &
                    error stop "red-black pressure solver requires even global sizes in periodic directions"
            end do
        end if
    end subroutine init_pressure_solver

    subroutine pressure_projection(ps, blk, dt_gamma, ibm, bc, c)
        type(pressure_solver_type), intent(in) :: ps
        type(block_set_type), intent(inout) :: blk
        real(C_DOUBLE), intent(in) :: dt_gamma
        type(ibm_type), intent(in) :: ibm
        type(boundary_type), intent(in) :: bc
        type(comm_type), intent(inout) :: c

        real(C_DOUBLE) :: t0
        integer(C_INT) :: iIter, dir
        real(C_DOUBLE) :: omega
        real(C_DOUBLE) :: dd, cc, alpha, alphaPrev, beta, gamma
        logical(C_BOOL) :: outLow(3), outHigh(3), refd(3)
        logical :: hasIface

        if (ps%method == PRESSURE_REDBLACK) then
            call redblack_projection(ps, blk, dt_gamma, ibm, bc, c)
            return
        end if

        t0 = prof_tic()
        call allocate_phi(blk)
        call allocate_rdenom(blk)
        if (ps%cheb) call allocate_delta(blk)

        ! Per-direction refinement flags ([blocks] refine_dims): a 2:1 face
        ! normal to an unrefined direction uses the uniform gradient metric
        ! (face_grad), all true in xyz mode.
        refd = blk%refMask == 1_C_INT
        hasIface = interfaces_present(blk, c)
        if (hasIface) call init_interface_faces(blk, refd)

        ! Dirichlet-pressure outlet faces (declared [boundary] _patch = outlet).
        ! The flags are per DOMAIN face; inside the kernels they act only where
        ! the block's own face is FACE_PHYS (i.e. the block face IS that domain
        ! face). The projection owns and corrects the outlet face flux against
        ! phi_face = 0, the increment of a HELD outlet pressure -- this is what
        ! keeps the Poisson problem solvable under net prescribed inflow.
        ! That face value is DATA, not a boundary row: phi's ghost plane behind
        ! an outlet face (index 0 or nb+1 along the outlet normal, tangential
        ! indices 1..nb) is zero from allocate_phi and nothing writes it
        ! afterwards --
        !   jacobi_compute_phi and redblack_sweep write interior cells (the
        !   sweep also index 0, but on OPEN faces only: iLo = 1 at a physical
        !   face);
        !   the scalar exchange writes the plane 0 / nb+1 of a direction only
        !   through an entry with off(d) = -1 / +1 (comm.f90 entry_boxes),
        !   which needs a leaf across that face -- an outlet is a non-periodic
        !   domain face and has none; an entry with off(d) = 0 reaches those
        !   indices only through the tangential extension, i.e. at HALO
        !   indices of its own direction: the edge and corner ghosts, which no
        !   kernel reads (jacobi_apply and interface_correct read phi at
        !   tangential indices 1..nb only).
        ! So the correction at an outlet face reads a zero ghost through the
        ! same 2*d1f metric as the diagonal (face_grad): one metric, SPD by
        ! construction. Until 2026-10-21 the ghost was instead re-MIRRORED
        ! (-phi_i) after every exchange and corrected with d1f -- the same bits,
        ! (phi_i - (-phi_i))*d1f == (phi_i - 0)*(2*d1f) in IEEE, at 36 surface
        ! launches per step.
        call outlet_flags(bc, outLow, outHigh)

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
        ! The diagonal is constant across the iterations (see rdenom), so it is
        ! formed once here rather than inside the loop -- and on a body-free rank
        ! it is constant across SUBSTAGES too, so it is formed once per run.
        call init_face_metrics(blk, outLow, outHigh, refd)
        if (.not. allocated(rdenomBlocks)) then
            ! First projection: every block, then narrow to the body ones.
            allocate(rdenomBlocks(blk%nBlocks))
            rdenomBlocks = [(dir, dir = 1, int(blk%nBlocks))]
            nRdenomBlocks = int(blk%nBlocks)
            dnsHasTerminal = c%has_terminal
            call enter_rdenom_blocks()
            call compute_rdenom(blk, ibm)
            call narrow_rdenom_blocks(ibm)
        else if (nRdenomBlocks > 0) then
            call compute_rdenom(blk, ibm)
        end if
        call prof_toc(proj_prof, PROF_PROJ_SETUP, t0)

        do iIter = 1_C_INT, ps%nIter
            t0 = prof_tic()
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
            end if
            ! The Chebyshev combine rides the same kernel (one pass over phi
            ! and delta instead of two, one launch instead of two per
            ! iteration; 2026-10-08).
            call jacobi_compute_phi(blk, omega, ps%cheb, alpha, gamma)
            call prof_toc(proj_prof, PROF_SWEEP, t0)
            t0 = prof_tic()
            call exchange_scalar_halos(c, phi, ifaceRow=.true.)
            call prof_toc(proj_prof, PROF_PHI_EXCHANGE, t0)
            t0 = prof_tic()
            call jacobi_apply(blk, dt_gamma, ibm)
            if (hasIface) call interface_correct(blk, ibm)
            call prof_toc(proj_prof, PROF_APPLY, t0)
            if (iIter == ps%nIter) then
                call projection_finish(blk, bc, c)
            else
                ! Between iterations only the divergence stencil reads the
                ! velocity halo, so refresh just q(nb+1) per dim, normal
                ! component (see comm.f90 sync_divergence_halos): the same-rank
                ! prefix as one kernel, the off-rank prefix as one message per
                ! peer, a sixth of the copy-only round's bytes at the same
                ! three kernel launches.
                t0 = prof_tic()
                call sync_divergence_halos(c, blk)
                call prof_toc(proj_prof, PROF_PROJ_VEL_EXCHANGE, t0)
            end if
        end do
    end subroutine pressure_projection

    ! The Dirichlet-pressure outlet faces (declared [boundary] _patch = outlet),
    ! per DOMAIN face and direction; inside the kernels they act only where
    ! the block's own face is FACE_PHYS. Both projections start from these.
    subroutine outlet_flags(bc, outLow, outHigh)
        type(boundary_type), intent(in) :: bc
        logical(C_BOOL), intent(out) :: outLow(3), outHigh(3)
        integer :: dir

        do dir = 1, 3
            outLow(dir)  = bc%facePatchType(boundary_face_id(dir, 0)) == PATCH_OUTLET
            outHigh(dir) = bc%facePatchType(boundary_face_id(dir, 1)) == PATCH_OUTLET
        end do
    end subroutine outlet_flags

    ! The end of EITHER projection: the physical-face ghosts ONCE, after the
    ! last correction and before the full exchange whose tangential extension
    ! copies them into neighbouring blocks' corner halos (numerics review F5,
    ! 2026-09-27) -- nothing between iterations reads them (the divergence
    ! stencil reads only the face-normal velocity at q(nb+1), pinned Dirichlet
    ! faces are masked out of the correction, phi has its own ghosts, and a
    ! Neumann NORMAL velocity, the one write that was not idempotent, is a
    ! config error since F4) -- then the full shell, all four variables with
    ! the cross-level transfers: what the next substage reads.
    subroutine projection_finish(blk, bc, c)
        type(block_set_type), intent(inout) :: blk
        type(boundary_type), intent(in) :: bc
        type(comm_type), intent(inout) :: c
        real(C_DOUBLE) :: t0

        t0 = prof_tic()
        call apply_bc(blk, bc)
        call prof_toc(proj_prof, PROF_PROJ_BC, t0)
        t0 = prof_tic()
        call exchange_halos(c, blk, [VAR_U, VAR_V, VAR_W, VAR_P])
        call prof_toc(proj_prof, PROF_PROJ_VEL_EXCHANGE, t0)
    end subroutine projection_finish

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

    ! rdenomBlocks lives on the device for the whole run (compute_rdenom indexes
    ! through it every substage on a body case), so it is mapped once here and
    ! remapped by narrow_rdenom_blocks when the list shrinks.
    subroutine enter_rdenom_blocks()
#ifdef USE_OPENMP_OFFLOAD
        !$omp target enter data map(to: rdenomBlocks)
#endif
    end subroutine enter_rdenom_blocks

    ! After the first projection has filled every block, keep only the blocks
    ! whose rdenom can still change -- those holding IBM coefficients. On a
    ! body-free rank that is none, and compute_rdenom never runs again.
    subroutine narrow_rdenom_blocks(ibm)
        type(ibm_type), intent(in) :: ibm

        integer, allocatable :: body(:)
        integer :: nBody

        call ibm_body_blocks(ibm, body, nBody)
#ifdef USE_OPENMP_OFFLOAD
        !$omp target exit data map(delete: rdenomBlocks)
#endif
        deallocate(rdenomBlocks)
        allocate(rdenomBlocks(max(1, nBody)))
        rdenomBlocks = 0
        if (nBody > 0) rdenomBlocks(1:nBody) = body(1:nBody)
        nRdenomBlocks = nBody
        call enter_rdenom_blocks()
        deallocate(body)
        ! Printed like the trip force's block list, and for the same reason: the
        ! fraction is the whole value of the change and a silent 100% would look
        ! exactly like a silent 0%.
        if (dnsHasTerminal) print '(a,i0,a,i0)', &
            " projection: rdenom recomputed on body blocks ", nBody, "/", size(rdenom,4)
    end subroutine narrow_rdenom_blocks

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

    ! Reciprocal-diagonal buffer, same bounds as phi so the two index identically.
    subroutine allocate_rdenom(blk)
        type(block_set_type), intent(in) :: blk
        if (.not. allocated(rdenom)) then
            allocate(rdenom(lbound(phi,1):ubound(phi,1), lbound(phi,2):ubound(phi,2), &
                            lbound(phi,3):ubound(phi,3), blk%nBlocks))
            rdenom = 0.0d0
#ifdef USE_OPENMP_OFFLOAD
            !$omp target enter data map(to: rdenom)
#endif
        end if
    end subroutine allocate_rdenom

    ! Fill the static face metric tables (see their declaration) and map them.
    ! Called from pressure_projection once the outlet/refinement flags are known;
    ! it returns immediately on every later call because nothing it reads can
    ! change during a run.
    subroutine init_face_metrics(blk, outLow, outHigh, refd)
        type(block_set_type), intent(in) :: blk
        logical(C_BOOL), intent(in) :: outLow(3), outHigh(3), refd(3)

        integer(C_INT) :: b, d, idx, n

        if (allocated(cfLow)) return
        allocate(cfLow(maxval(blk%nb),3,blk%nBlocks), cfHigh(3,blk%nBlocks), &
                 dnLow(maxval(blk%nb),3,blk%nBlocks), dnHigh(maxval(blk%nb),3,blk%nBlocks), &
                 d1P(maxval(blk%nb),3,blk%nBlocks))
        cfLow = 0.0d0
        cfHigh = 0.0d0
        dnLow = 0.0d0
        dnHigh = 0.0d0
        d1P = 0.0d0

        do b = 1_C_INT, blk%nBlocks
            do d = 1_C_INT, 3_C_INT
                n = blk%nb(d)
                do idx = 1_C_INT, n
                    dnLow(idx,d,b) = face_grad(blk%physLow(d,b), idx == 1_C_INT, &
                        d1_normal(blk, d, idx, b), outLow(d), refd(d))
                    dnHigh(idx,d,b) = face_grad(blk%physHigh(d,b), idx == n, &
                        d1_normal(blk, d, idx + 1_C_INT, b), outHigh(d), refd(d))
                    d1P(idx,d,b) = d1_divergence(blk, d, idx, b)
                end do
                ! The correction metric is the diagonal's, except: interface
                ! faces belong to interface_correct (zero rather than a second
                ! predicate, exactly as the kernel used to do), and the high
                ! face is the neighbour's low face (filled by the velocity
                ! exchange) unless it is an outlet physical face.
                cfLow(1:n,d,b) = dnLow(1:n,d,b)
                if (is_interface(blk%physLow(d,b))) cfLow(1,d,b) = 0.0d0
                if (outHigh(d) .and. blk%physHigh(d,b) == FACE_PHYS) cfHigh(d,b) = dnHigh(n,d,b)
            end do
        end do

#ifdef USE_OPENMP_OFFLOAD
        !$omp target enter data map(to: cfLow, cfHigh, dnLow, dnHigh, d1P)
#endif
    end subroutine init_face_metrics

    ! Fill the 2:1 interface face list (see its declaration) and map it. Like
    ! init_face_metrics it runs once: nothing it reads changes during a run.
    ! The order is (block, direction, low face, high face); it does not matter,
    ! every entry writes its own face and reads phi only. An interface face is
    ! never a physical face, so face_grad's outlet flag is simply false here.
    subroutine init_interface_faces(blk, refd)
        type(block_set_type), intent(in) :: blk
        logical(C_BOOL), intent(in) :: refd(3)

        integer(C_INT) :: b, d, n, e

        if (ifaceListBuilt) return
        ifaceListBuilt = .true.

        n = 0_C_INT
        do b = 1_C_INT, blk%nBlocks
            do d = 1_C_INT, 3_C_INT
                if (is_interface(blk%physLow(d,b)))  n = n + 1_C_INT
                if (is_interface(blk%physHigh(d,b))) n = n + 1_C_INT
            end do
        end do
        nIfaces = n
        ! At least one element so the device map is well-formed on a rank
        ! without interface faces (the kernel never reads it there).
        allocate(ifBlock(max(1_C_INT, n)), ifDir(max(1_C_INT, n)), &
                 ifIdx(max(1_C_INT, n)), ifCoef(max(1_C_INT, n)))
        ifBlock = 0_C_INT; ifDir = 0_C_INT; ifIdx = 0_C_INT; ifCoef = 0.0d0

        e = 0_C_INT
        do b = 1_C_INT, blk%nBlocks
            do d = 1_C_INT, 3_C_INT
                if (is_interface(blk%physLow(d,b))) then
                    e = e + 1_C_INT
                    ifBlock(e) = b; ifDir(e) = d; ifIdx(e) = 1_C_INT
                    ifCoef(e) = face_grad(blk%physLow(d,b), .true., &
                        d1_normal(blk, d, 1_C_INT, b), .false._C_BOOL, refd(d))
                end if
                if (is_interface(blk%physHigh(d,b))) then
                    e = e + 1_C_INT
                    ifBlock(e) = b; ifDir(e) = d; ifIdx(e) = blk%nb(d) + 1_C_INT
                    ifCoef(e) = face_grad(blk%physHigh(d,b), .true., &
                        d1_normal(blk, d, blk%nb(d) + 1_C_INT, b), .false._C_BOOL, refd(d))
                end if
            end do
        end do

#ifdef USE_OPENMP_OFFLOAD
        !$omp target enter data map(to: ifBlock, ifDir, ifIdx, ifCoef)
#endif
    end subroutine init_interface_faces

    ! The face metric d1f of direction d at face index idx of block b, for the
    ! velocity component NORMAL to that face, and the cell's own divergence
    ! metric in direction d. Host-side helpers for init_face_metrics only.
    real(C_DOUBLE) function d1_normal(blk, d, idx, b)
        type(block_set_type), intent(in) :: blk
        integer(C_INT), intent(in) :: d, idx, b

        if (d == 1_C_INT) then
            d1_normal = blk%d1x(idx,VAR_U,b)
        else if (d == 2_C_INT) then
            d1_normal = blk%d1y(idx,VAR_V,b)
        else
            d1_normal = blk%d1z(idx,VAR_W,b)
        end if
    end function d1_normal

    real(C_DOUBLE) function d1_divergence(blk, d, idx, b)
        type(block_set_type), intent(in) :: blk
        integer(C_INT), intent(in) :: d, idx, b

        if (d == 1_C_INT) then
            d1_divergence = blk%d1x(idx,VAR_P,b)
        else if (d == 2_C_INT) then
            d1_divergence = blk%d1y(idx,VAR_P,b)
        else
            d1_divergence = blk%d1z(idx,VAR_P,b)
        end if
    end function d1_divergence

    ! The RECIPROCAL projection diagonal: each face's pressure-gradient metric,
    ! summed, then inverted.
    ! face_grad returns 0 for a pinned wall face, the coarse-fine gradient 1/d
    ! for a 2:1 interface face (the symmetric composite stencil), and the regular
    ! 1/h otherwise. An OUTLET physical face counts 2*d1f: phi is held at 0 ON
    ! the face, so the flux sensitivity per unit phi_i is the half-cell
    ! gradient's -- and jacobi_apply corrects that face with the same 2*d1f
    ! against the zero ghost (one metric: SPD lives there).
    !
    ! Called once per substage on the blocks listed in rdenomBlocks (every
    ! block on the first call, the body blocks afterwards -- see the rdenom
    ! declaration), and stored inverted so the iteration kernel multiplies
    ! rather than divides.

    subroutine compute_rdenom(blk, ibm)
        type(block_set_type), intent(in) :: blk
        type(ibm_type), intent(in) :: ibm

        real(C_DOUBLE) :: mu_u_i, mu_u_ip, mu_v_j, mu_v_jp, mu_w_k, mu_w_kp
        integer(C_INT) :: i, ip, j, jp, k, kp, b, bb, nTodo, nx, ny, nz

        nx = blk%nb(1); ny = blk%nb(2); nz = blk%nb(3)
        nTodo = int(nRdenomBlocks, C_INT)
        if (nTodo <= 0_C_INT) return

        ! dnLow/dnHigh/d1P are the face_grad and divergence metrics
        ! precomputed per (normal index, direction, block) -- see their
        ! declaration for why that is bit-exact and why d1P stays a separate
        ! factor. Only ibm%mu, which changes every substage, is still read here.
#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do collapse(4) &
        !$omp& map(to: nx, ny, nz, nTodo, dnLow, dnHigh, d1P, ibm%mu, rdenomBlocks) &
        !$omp& map(tofrom: rdenom) &
        !$omp& private(i,ip,j,jp,k,kp,b,bb, &
        !$omp& mu_u_i,mu_u_ip,mu_v_j,mu_v_jp,mu_w_k,mu_w_kp)
#endif
        do bb = 1_C_INT, nTodo
        do k = 1_C_INT, nz
            do j = 1_C_INT, ny
                do i = 1_C_INT, nx
                    b = rdenomBlocks(bb)
                    ip = i + 1; jp = j + 1; kp = k + 1

                    mu_u_i  = ibm%mu(i,j,k,VAR_U,b);  mu_u_ip = ibm%mu(ip,j,k,VAR_U,b)
                    mu_v_j  = ibm%mu(i,j,k,VAR_V,b);  mu_v_jp = ibm%mu(i,jp,k,VAR_V,b)
                    mu_w_k  = ibm%mu(i,j,k,VAR_W,b);  mu_w_kp = ibm%mu(i,j,kp,VAR_W,b)

                    rdenom(i,j,k,b) = 1.0d0/( &
                            (dnLow(i,1,b)*mu_u_i + dnHigh(i,1,b)*mu_u_ip)*d1P(i,1,b) &
                          + (dnLow(j,2,b)*mu_v_j + dnHigh(j,2,b)*mu_v_jp)*d1P(j,2,b) &
                          + (dnLow(k,3,b)*mu_w_k + dnHigh(k,3,b)*mu_w_kp)*d1P(k,3,b))
                end do
            end do
        end do
        end do
#ifdef USE_OPENMP_OFFLOAD
        !$omp end target teams distribute parallel do
#endif
    end subroutine compute_rdenom

    ! Jacobi step 1: phi(i,j,k) = -omega * div / denom for every interior cell,
    ! from the frozen velocity. denom is the projection diagonal (sum of the
    ! cell's non-pinned face metrics); pinned faces (walls, closed faces) leave
    ! the diagonal exactly as in the red-black scheme.
    ! omega: the damping factor. Plain Jacobi passes ps%omega (the increment IS
    ! phi = -omega*div/denom); Chebyshev passes 1.0 so phi holds the pure
    ! diagonal-preconditioned residual z = -div/denom (the damping then comes
    ! from the Chebyshev coefficients).
    ! Chebyshev (cheb = .true.): the preconditioned residual z computed here is
    ! combined in place, phi := alpha*z + gamma*delta, and delta := phi (the
    ! new increment delta_k), so jacobi_apply adds phi exactly as it would the
    ! Jacobi increment. One pass over phi and delta; the combine was its own
    ! kernel (cheb_combine) until 2026-10-08.
    subroutine jacobi_compute_phi(blk, omega, cheb, alpha, gamma)
        type(block_set_type), intent(inout) :: blk
        real(C_DOUBLE), intent(in) :: omega, alpha, gamma
        logical, intent(in) :: cheb

        real(C_DOUBLE) :: div, z, s
        integer(C_INT) :: i, ip, j, jp, k, kp, b, nBlocks, nx, ny, nz

        nx = blk%nb(1); ny = blk%nb(2); nz = blk%nb(3)
        nBlocks = blk%nBlocks

        ! d1P is the same divergence metric d1?(idx,VAR_P,b), gathered into one
        ! table by init_face_metrics -- one array base here instead of three.
        ! Two loops, one per mode, because delta exists only under Chebyshev
        ! (an unallocated array cannot sit in a map clause); the divergence and
        ! z lines are the same text in both.
        if (cheb) then
#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do collapse(4) &
        !$omp& map(to: omega, alpha, gamma, nx, ny, nz, d1P) &
        !$omp& map(tofrom: phi, delta, blk%q, rdenom) &
        !$omp& private(i,ip,j,jp,k,kp,b,div,z,s)
#endif
        do b = 1_C_INT, nBlocks
        do k = 1_C_INT, nz
            do j = 1_C_INT, ny
                do i = 1_C_INT, nx
                    ip = i + 1; jp = j + 1; kp = k + 1
                    div = (blk%q(ip,j,k,VAR_U,b)-blk%q(i,j,k,VAR_U,b))*d1P(i,1,b) &
                        + (blk%q(i,jp,k,VAR_V,b)-blk%q(i,j,k,VAR_V,b))*d1P(j,2,b) &
                        + (blk%q(i,j,kp,VAR_W,b)-blk%q(i,j,k,VAR_W,b))*d1P(k,3,b)
                    z = -omega*div*rdenom(i,j,k,b)
                    s = alpha*z + gamma*delta(i,j,k,b)
                    phi(i,j,k,b) = s
                    delta(i,j,k,b) = s
                end do
            end do
        end do
        end do
#ifdef USE_OPENMP_OFFLOAD
        !$omp end target teams distribute parallel do
#endif
        else
#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do collapse(4) &
        !$omp& map(to: omega, nx, ny, nz, d1P) &
        !$omp& map(tofrom: phi, blk%q, rdenom) &
        !$omp& private(i,ip,j,jp,k,kp,b,div)
#endif
        do b = 1_C_INT, nBlocks
        do k = 1_C_INT, nz
            do j = 1_C_INT, ny
                do i = 1_C_INT, nx
                    ip = i + 1; jp = j + 1; kp = k + 1

                    div = (blk%q(ip,j,k,VAR_U,b)-blk%q(i,j,k,VAR_U,b))*d1P(i,1,b) &
                        + (blk%q(i,jp,k,VAR_V,b)-blk%q(i,j,k,VAR_V,b))*d1P(j,2,b) &
                        + (blk%q(i,j,kp,VAR_W,b)-blk%q(i,j,k,VAR_W,b))*d1P(k,3,b)

                    ! Multiply by the stored RECIPROCAL: see the rdenom
                    ! declaration for why this is not bit-identical to a divide
                    ! and why that is the right trade here.
                    phi(i,j,k,b) = -omega*div*rdenom(i,j,k,b)
                end do
            end do
        end do
        end do
#ifdef USE_OPENMP_OFFLOAD
        !$omp end target teams distribute parallel do
#endif
        end if
    end subroutine jacobi_compute_phi

    ! Jacobi step 3: apply the frozen increment. Pressure gets p += phi/dt_gamma;
    ! each velocity face is corrected once from the two cells it separates,
    !   q_face += (phi_below - phi_above) * d1(face) * mu(face),
    ! the exact sum of the two adjacent red-black corrections, but read from the
    ! frozen phi (with the halo layer filled by the scalar exchange) so there is
    ! no in-place race and no colouring. Each cell corrects its own LOW faces
    ! (1..nb) and, on a declared outlet, the high face nb+1 against the
    ! zero phi ghost (cfHigh); every other high halo face is the
    ! neighbour's low face, filled by the velocity exchange. Pinned faces have
    ! cfLow = 0 and are left untouched.
    subroutine jacobi_apply(blk, dt_gamma, ibm)
        type(block_set_type), intent(inout) :: blk
        real(C_DOUBLE), intent(in) :: dt_gamma
        type(ibm_type), intent(in) :: ibm

        real(C_DOUBLE) :: idt, cf
        integer(C_INT) :: i, ip, j, jp, k, kp, b, nBlocks, nx, ny, nz

        nx = blk%nb(1); ny = blk%nb(2); nz = blk%nb(3)
        nBlocks = blk%nBlocks
        idt = 1.0_C_DOUBLE/dt_gamma

        ! ONE kernel for the pressure update AND the velocity-face corrections
        ! (two until 2026-10-08: the pressure update streamed phi and p on its
        ! own, 10 % of the device time of a Jacobi step for one extra read of
        ! phi and 18 launches per step; the cell's own phi is loaded here
        ! anyway). The two write disjoint fields and read phi only, so the
        ! fusion changes no operand and no order -- bit-exact by construction.
        !
        ! Velocity-face corrections. Each cell corrects its LOW faces (1..nb):
        ! q_face += (phi_below - phi_self) * grad(face) * mu, with grad(face) from
        ! face_grad (0 if pinned, ifGrad if a 2:1 interface, 1/h otherwise). The
        ! HIGH face is normally the neighbour's low face (filled by the exchange),
        ! so it is corrected here ONLY when it is itself a 2:1 interface face --
        ! then this block owns it, reconstructing it from the halo phi (the coarse
        ! halo = avg of the fine increments / the fine halo = the coarse increment,
        ! both from the scalar exchange), which keeps the two sides conservative.
        ! An OUTLET physical face (either side) is also corrected, with the
        ! half-cell Dirichlet metric 2*d1f against the ZERO phi ghost (the
        ! increment of the held outlet pressure; see pressure_projection):
        ! q(nx+1) += (phi(nx) - 0)*2*d1f*mu -- the same 2*d1f the denominator
        ! counts, so the operator is SPD at that row by construction. The
        ! predictor advances the outlet face with that same half-cell gradient
        ! of the stored pressure (step.f90 predict_outlet_faces), so phi is a
        ! true increment there too and p += phi/dt_gamma keeps the boundary
        ! cell tied to the outlet value. cfLow/cfHigh hold face_grad for every
        ! face this kernel owns, precomputed per (block, normal index) -- see
        ! their declaration for why that is both bit-exact and the point of the
        ! kernel's shape.
#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do collapse(4) &
        !$omp& map(to: idt, nx, ny, nz, cfLow, cfHigh, ibm%mu) &
        !$omp& map(tofrom: blk%q, phi) private(i,ip,j,jp,k,kp,b,cf)
#endif
        do b = 1_C_INT, nBlocks
        do k = 1_C_INT, nz
            do j = 1_C_INT, ny
                do i = 1_C_INT, nx
                    blk%q(i,j,k,VAR_P,b) = blk%q(i,j,k,VAR_P,b) + phi(i,j,k,b)*idt
                    cf = cfLow(i,1,b)
                    if (cf /= 0.0d0) blk%q(i,j,k,VAR_U,b) = blk%q(i,j,k,VAR_U,b) &
                        + (phi(i-1,j,k,b) - phi(i,j,k,b))*cf*ibm%mu(i,j,k,VAR_U,b)
                    cf = cfLow(j,2,b)
                    if (cf /= 0.0d0) blk%q(i,j,k,VAR_V,b) = blk%q(i,j,k,VAR_V,b) &
                        + (phi(i,j-1,k,b) - phi(i,j,k,b))*cf*ibm%mu(i,j,k,VAR_V,b)
                    cf = cfLow(k,3,b)
                    if (cf /= 0.0d0) blk%q(i,j,k,VAR_W,b) = blk%q(i,j,k,VAR_W,b) &
                        + (phi(i,j,k-1,b) - phi(i,j,k,b))*cf*ibm%mu(i,j,k,VAR_W,b)

                    ! High faces: only an OUTLET physical face (FACE_PHYS +
                    ! declared outlet) is owned here, and cfHigh is zero unless
                    ! it is one; the 2:1 interface ones moved to
                    ! interface_correct. Reading cfHigh under the plane
                    ! predicate keeps it off the common path.
                    if (i == nx) then
                        cf = cfHigh(1,b); ip = i + 1
                        if (cf /= 0.0d0) blk%q(ip,j,k,VAR_U,b) = blk%q(ip,j,k,VAR_U,b) &
                            + (phi(i,j,k,b) - phi(ip,j,k,b))*cf*ibm%mu(ip,j,k,VAR_U,b)
                    end if
                    if (j == ny) then
                        cf = cfHigh(2,b); jp = j + 1
                        if (cf /= 0.0d0) blk%q(i,jp,k,VAR_V,b) = blk%q(i,jp,k,VAR_V,b) &
                            + (phi(i,j,k,b) - phi(i,jp,k,b))*cf*ibm%mu(i,jp,k,VAR_V,b)
                    end if
                    if (k == nz) then
                        cf = cfHigh(3,b); kp = k + 1
                        if (cf /= 0.0d0) blk%q(i,j,kp,VAR_W,b) = blk%q(i,j,kp,VAR_W,b) &
                            + (phi(i,j,k,b) - phi(i,j,kp,b))*cf*ibm%mu(i,j,kp,VAR_W,b)
                    end if
                end do
            end do
        end do
        end do
#ifdef USE_OPENMP_OFFLOAD
        !$omp end target teams distribute parallel do
#endif
    end subroutine jacobi_apply

    ! See the ifaceAny declaration: global, cached, computed once.
    logical function interfaces_present(blk, c)
        type(block_set_type), intent(in) :: blk
        type(comm_type), intent(in) :: c

        integer(C_INT) :: flag(1)
        integer :: b, d

        if (.not. ifaceKnown) then
            flag(1) = 0_C_INT
            scan: do b = 1, int(blk%nBlocks)
                do d = 1, 3
                    if (is_interface(blk%physLow(d,b)) .or. &
                        is_interface(blk%physHigh(d,b))) then
                        flag(1) = 1_C_INT
                        exit scan
                    end if
                end do
            end do scan
            call comm_allreduce_max_int(c, flag)
            ifaceAny = flag(1) == 1_C_INT
            ifaceKnown = .true.
        end if
        interfaces_present = ifaceAny
    end function interfaces_present

    ! The 2:1 interface-face correction, shared by BOTH smoothers:
    !
    !   q_face += (phi_neighbour - phi_self) * face_grad * mu
    !
    ! for the interface faces only -- the low face at index 1 and the high halo
    ! face at index nb+1. The formula and its operand order are exactly the ones
    ! jacobi_apply used to carry inline, which is what makes the Jacobi path
    ! bit-exact across the extraction; face_grad comes from the ifCoef table,
    ! formed once on the host from the same arguments (the cfLow argument: the
    ! same expression on the same inputs, evaluated once instead of per point).
    !
    ! Why it is a separate pass rather than part of a volume kernel: the faces
    ! it touches are planes, so this is nb^2 work against the sweep's nb^3, and
    ! red-black needs it applied AFTER a phi exchange that cannot happen inside
    ! its sweep.
    !
    ! ONE launch over the face list (three per call until 2026-10-21, one per
    ! direction, each over every block with the face kind tested per point):
    ! the collapsed loop runs over (face, t2, t1) with t1/t2 spanning the
    ! largest tangential extents of the three orientations, and a face whose
    ! own extents are smaller skips the excess -- nothing on cubic blocks.
    ! Low and high faces are one statement: with f the face index, the
    ! difference is phi(f-1) - phi(f) on either side (phi(0) - phi(1) at the
    ! low face, phi(nb) - phi(nb+1) at the high one).
    !
    ! An interface face is never an outlet face, so the plain face_grad is
    ! the metric here (the outlet-aware pair lives in the tables of the
    ! volume kernels).
    !
    ! Interior cells only (1..nb in the tangential indices). Tangential-halo
    ! copies of an interface face are repaired by the same-level velocity
    ! exchange that follows, so the edge/corner phi ghosts (plain injection,
    ! not the ifaceRow restrict) are never read here.
    subroutine interface_correct(blk, ibm)
        type(block_set_type), intent(inout) :: blk
        type(ibm_type), intent(in) :: ibm

        integer(C_INT) :: e, t1, t2, n1, n2, nF, nx, ny, nz
        integer(C_INT) :: i, j, k, io, jo, ko, b, d, f, var
        real(C_DOUBLE) :: cf
        logical :: inFace

        nF = nIfaces
        if (nF <= 0_C_INT) return
        nx = blk%nb(1); ny = blk%nb(2); nz = blk%nb(3)
        ! Tangential extents: x faces (ny, nz), y faces (nx, nz), z faces
        ! (nx, ny); t1 and t2 run over the largest of each.
        n1 = max(nx, ny); n2 = max(ny, nz)

#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do collapse(3) &
        !$omp& map(to: nF, n1, n2, nx, ny, nz, ifBlock, ifDir, ifIdx, ifCoef, ibm%mu) &
        !$omp& map(tofrom: blk%q, phi) private(t1,t2,i,j,k,io,jo,ko,b,d,f,var,cf,inFace)
#endif
        do e = 1_C_INT, nF
        do t2 = 1_C_INT, n2
            do t1 = 1_C_INT, n1
                b = ifBlock(e); d = ifDir(e); f = ifIdx(e); cf = ifCoef(e)
                ! (i,j,k) the face point, (io,jo,ko) the cell below it along d.
                if (d == 1_C_INT) then
                    inFace = t1 <= ny .and. t2 <= nz
                    i = f; j = t1; k = t2; io = f - 1_C_INT; jo = t1; ko = t2
                    var = VAR_U
                else if (d == 2_C_INT) then
                    inFace = t1 <= nx .and. t2 <= nz
                    i = t1; j = f; k = t2; io = t1; jo = f - 1_C_INT; ko = t2
                    var = VAR_V
                else
                    inFace = t1 <= nx .and. t2 <= ny
                    i = t1; j = t2; k = f; io = t1; jo = t2; ko = f - 1_C_INT
                    var = VAR_W
                end if
                if (inFace) blk%q(i,j,k,var,b) = blk%q(i,j,k,var,b) &
                    + (phi(io,jo,ko,b) - phi(i,j,k,b))*cf*ibm%mu(i,j,k,var,b)
            end do
        end do
        end do
#ifdef USE_OPENMP_OFFLOAD
        !$omp end target teams distribute parallel do
#endif
    end subroutine interface_correct

    ! A 2:1 coarse-fine interface face (either orientation).
    pure logical function is_interface(fk)
!$omp declare target
        integer(C_INT), intent(in) :: fk

        is_interface = fk == FACE_COARSE .or. fk == FACE_FINE
    end function is_interface

    ! Pressure-gradient metric for one cell face: THE metric of the projection
    ! operator -- the diagonal (compute_rdenom, the red-black sweep), the
    ! Jacobi velocity-face correction (cfLow/cfHigh) and the red-black in-place
    ! correction all use this one value, which is what makes the operator SPD
    ! at every face row. atBnd marks a face on the
    ! block boundary (where the face kind fk applies); interior faces always use
    ! the regular metric d1f. A pinned wall/closed face contributes 0 -- except
    ! a declared Dirichlet-pressure OUTLET (outlet: the domain face's flag,
    ! acting where the block face is FACE_PHYS), whose metric is the half-cell
    ! gradient 2*d1f against phi_face = 0 (the Jacobi path reads that zero from
    ! phi's ghost, the red-black path reads no ghost: pressure_projection). A 2:1
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
    pure real(C_DOUBLE) function face_grad(fk, atBnd, d1f, outlet, refined)
!$omp declare target
        integer(C_INT), intent(in) :: fk
        logical, intent(in) :: atBnd
        real(C_DOUBLE), intent(in) :: d1f
        logical(C_BOOL), intent(in) :: outlet, refined

        if (.not. atBnd) then
            face_grad = d1f
        else if (outlet .and. fk == FACE_PHYS) then
            face_grad = 2.0d0*d1f
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

    ! -------------------------------------------------------------------------
    ! Red-black SOR projection ([pressure] solver = redblack): the original
    ! coloured Gauss-Seidel scheme, restored as a selectable single-level
    ! alternative to damped Jacobi. In-place over-relaxed sweeps; the same-level
    ! halo exchange refreshes the velocity each colour (pressure only at the
    ! end). It uses the one metric face_grad for BOTH the divergence
    ! denominator AND the in-place face corrections (d1f interior, 0 pinned,
    ! 2*d1f for the Dirichlet-p outlet; the Jacobi tables hold the same
    ! values), so the outlet is handled exactly as in the Jacobi path and, on
    ! outlet-free grids, the sweep is identical to the original red-black
    ! scheme.
    subroutine redblack_projection(ps, blk, dt_gamma, ibm, bc, c)
        type(pressure_solver_type), intent(in) :: ps
        type(block_set_type), intent(inout) :: blk
        real(C_DOUBLE), intent(in) :: dt_gamma
        type(ibm_type), intent(in) :: ibm
        type(boundary_type), intent(in) :: bc
        type(comm_type), intent(inout) :: c

        integer(C_INT) :: iIter, color, dir
        logical(C_BOOL) :: outLow(3), outHigh(3), refd(3)
        logical :: hasIface
        real(C_DOUBLE) :: t0

        t0 = prof_tic()
        refd = blk%refMask == 1_C_INT
        hasIface = interfaces_present(blk, c)
        ! phi exists on this path ONLY to cross a level jump; a single-level
        ! red-black run allocates and maps nothing.
        if (hasIface) then
            call allocate_phi(blk)
            call init_interface_faces(blk, refd)
        end if
        call outlet_flags(bc, outLow, outHigh)
        call prof_toc(proj_prof, PROF_PROJ_SETUP, t0)

        ! The coupled sweep updates pressure AND velocity in one kernel, so on a
        ! SINGLE-LEVEL grid proj_timing's `apply` and `phi_exchange` stay zero
        ! here and the red-black path spends everything in `sweep` and the
        ! velocity exchange. That asymmetry against the Jacobi path is real, not
        ! missing instrumentation. With a 2:1 interface the two buckets fill:
        ! `phi_exchange` is the per-colour cross-level transmission of the
        ! increment and `apply` is interface_correct.
        do iIter = 1_C_INT, ps%nIter
            do color = 1_C_INT, 0_C_INT, -1_C_INT
                t0 = prof_tic()
                ! Zero per COLOUR, not per iteration: the patch must see this
                ! colour's increment alone (see redblack_sweep).
                if (hasIface) call zero_phi(blk)
                call redblack_sweep(ps, blk, dt_gamma, ibm, color, outLow, outHigh, refd, hasIface)
                call prof_toc(proj_prof, PROF_SWEEP, t0)
                if (hasIface) then
                    t0 = prof_tic()
                    call exchange_scalar_halos(c, phi, ifaceRow=.true.)
                    call prof_toc(proj_prof, PROF_PHI_EXCHANGE, t0)
                    t0 = prof_tic()
                    call interface_correct(blk, ibm)
                    call prof_toc(proj_prof, PROF_APPLY, t0)
                end if
                if (iIter == ps%nIter .and. color == 0_C_INT) then
                    call projection_finish(blk, bc, c)
                else
                    t0 = prof_tic()
                    call exchange_halos(c, blk, [VAR_U, VAR_V, VAR_W], cross_level=.false.)
                    call prof_toc(proj_prof, PROF_PROJ_VEL_EXCHANGE, t0)
                end if
            end do
        end do
    end subroutine redblack_projection

    ! Zero phi over the whole array, halos included: the sweep writes only its
    ! own colour's interior cells, and the interface patch must not read a stale
    ! value anywhere else.
    subroutine zero_phi(blk)
        type(block_set_type), intent(in) :: blk
        integer(C_INT) :: i, j, k, b

#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do collapse(4) map(tofrom: phi) &
        !$omp& private(i,j,k,b)
#endif
        do b = 1_C_INT, blk%nBlocks
        do k = lbound(phi,3), ubound(phi,3)
            do j = lbound(phi,2), ubound(phi,2)
                do i = lbound(phi,1), ubound(phi,1)
                    phi(i,j,k,b) = 0.0d0
                end do
            end do
        end do
        end do
#ifdef USE_OPENMP_OFFLOAD
        !$omp end target teams distribute parallel do
#endif
    end subroutine zero_phi

    ! One red-black sweep of `color`. Each block sweeps from its halo layer (0)
    ! redundantly with the owning neighbour except on physical boundaries -- the
    ! rank-level red-black trick one level down -- so results are independent of
    ! the block/rank layout; parity is anchored to the global index space via
    ! the block origin. `sor` = ps%omega (the config "sor"; over-relaxation > 1
    ! is allowed here, unlike the Jacobi damping).
    ! storePhi: with a 2:1 interface present the increment has to survive the
    ! sweep so interface_correct can transmit it across the level jump, so it is
    ! written to the module phi array as well as applied in place. Interface
    ! faces are then kept in the DENOMINATOR but NOT corrected here -- the patch
    ! owns them, which keeps the interface formula in exactly one place and
    ! avoids stray writes into tangential halo copies of an interface face.
    subroutine redblack_sweep(ps, blk, dt_gamma, ibm, color, outLow, outHigh, refd, storePhi)
        type(pressure_solver_type), intent(in) :: ps
        type(block_set_type), intent(inout) :: blk
        real(C_DOUBLE), intent(in) :: dt_gamma
        type(ibm_type), intent(in) :: ibm
        integer(C_INT), intent(in) :: color
        logical(C_BOOL), intent(in) :: outLow(3), outHigh(3), refd(3)
        logical, intent(in) :: storePhi

        real(C_DOUBLE) :: phiC, denom, idt, sor, div
        real(C_DOUBLE) :: gLo1, gHi1, gLo2, gHi2, gLo3, gHi3
        real(C_DOUBLE) :: cLo1, cHi1, cLo2, cHi2, cLo3, cHi3
        logical :: store
        real(C_DOUBLE) :: mu_u_i, mu_u_ip, mu_v_j, mu_v_jp, mu_w_k, mu_w_kp
        integer(C_INT) :: i, ip, j, jp, k, kp, b, nBlocks, nLowerHaloDirections, iColor, nColorX
        integer(C_INT) :: iLo, jLo, kLo, colorOffset
        integer(C_INT) :: hi(1:3)

        hi = blk%nb(1:3)
        sor = ps%omega
        nColorX = (hi(1) + 2_C_INT)/2_C_INT
        nBlocks = blk%nBlocks
        idt = 1.0_C_DOUBLE/dt_gamma
        store = storePhi

#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do collapse(4) &
        !$omp& map(to: color, nColorX, sor, idt, dt_gamma, hi(1:3), store, &
        !$omp& outLow(1:3), outHigh(1:3), refd(1:3), &
        !$omp& blk%origin, blk%physLow, blk%physHigh, &
        !$omp& blk%d1x, blk%d1y, blk%d1z, ibm%mu) &
        !$omp& map(tofrom: blk%q, phi) &
        !$omp& private(i,ip,j,jp,k,kp,b,iColor,iLo,jLo,kLo,colorOffset, &
        !$omp& phiC,denom,div,nLowerHaloDirections, gLo1,gHi1,gLo2,gHi2,gLo3,gHi3, &
        !$omp& cLo1,cHi1,cLo2,cHi2,cLo3,cHi3, &
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
                    gLo1 = face_grad(blk%physLow(1,b),  i == 1_C_INT, blk%d1x(i,VAR_U,b),  outLow(1),  refd(1))
                    gHi1 = face_grad(blk%physHigh(1,b), i == hi(1),   blk%d1x(ip,VAR_U,b), outHigh(1), refd(1))
                    gLo2 = face_grad(blk%physLow(2,b),  j == 1_C_INT, blk%d1y(j,VAR_V,b),  outLow(2),  refd(2))
                    gHi2 = face_grad(blk%physHigh(2,b), j == hi(2),   blk%d1y(jp,VAR_V,b), outHigh(2), refd(2))
                    gLo3 = face_grad(blk%physLow(3,b),  k == 1_C_INT, blk%d1z(k,VAR_W,b),  outLow(3),  refd(3))
                    gHi3 = face_grad(blk%physHigh(3,b), k == hi(3),   blk%d1z(kp,VAR_W,b), outHigh(3), refd(3))

                    denom = (gLo1*mu_u_i + gHi1*mu_u_ip)*blk%d1x(i,VAR_P,b) &
                          + (gLo2*mu_v_j + gHi2*mu_v_jp)*blk%d1y(j,VAR_P,b) &
                          + (gLo3*mu_w_k + gHi3*mu_w_kp)*blk%d1z(k,VAR_P,b)

                    div = (blk%q(ip,j,k,VAR_U,b)-blk%q(i,j,k,VAR_U,b))*blk%d1x(i,VAR_P,b) &
                        + (blk%q(i,jp,k,VAR_V,b)-blk%q(i,j,k,VAR_V,b))*blk%d1y(j,VAR_P,b) &
                        + (blk%q(i,j,kp,VAR_W,b)-blk%q(i,j,k,VAR_W,b))*blk%d1z(k,VAR_P,b)

                    phiC = -sor*div/denom

                    blk%q(i,j,k,VAR_P,b) = blk%q(i,j,k,VAR_P,b) + phiC*idt

                    ! Interface faces stay in the denominator above but are not
                    ! corrected in place: interface_correct applies them from the
                    ! exchanged phi. Every other face keeps the diagonal metric,
                    ! so with no interface present cX == gX and the arithmetic is
                    ! untouched.
                    cLo1 = merge(0.0d0, gLo1, i == 1_C_INT     .and. is_interface(blk%physLow(1,b)))
                    cHi1 = merge(0.0d0, gHi1, i == hi(1)       .and. is_interface(blk%physHigh(1,b)))
                    cLo2 = merge(0.0d0, gLo2, j == 1_C_INT     .and. is_interface(blk%physLow(2,b)))
                    cHi2 = merge(0.0d0, gHi2, j == hi(2)       .and. is_interface(blk%physHigh(2,b)))
                    cLo3 = merge(0.0d0, gLo3, k == 1_C_INT     .and. is_interface(blk%physLow(3,b)))
                    cHi3 = merge(0.0d0, gHi3, k == hi(3)       .and. is_interface(blk%physHigh(3,b)))

                    ! Symmetric correction: every non-pinned face of the cell is
                    ! moved by phi*metric*mu; the neighbour of the other colour
                    ! moves the shared face from its side, so over both colours a
                    ! face sees (phi_below - phi_above)*d1f -- the same operator as
                    ! the Jacobi apply, here in place.
                    blk%q(i,j,k,VAR_U,b)  = blk%q(i,j,k,VAR_U,b)  - phiC*cLo1*mu_u_i
                    blk%q(ip,j,k,VAR_U,b) = blk%q(ip,j,k,VAR_U,b) + phiC*cHi1*mu_u_ip
                    blk%q(i,j,k,VAR_V,b)  = blk%q(i,j,k,VAR_V,b)  - phiC*cLo2*mu_v_j
                    blk%q(i,jp,k,VAR_V,b) = blk%q(i,jp,k,VAR_V,b) + phiC*cHi2*mu_v_jp
                    blk%q(i,j,k,VAR_W,b)  = blk%q(i,j,k,VAR_W,b)  - phiC*cLo3*mu_w_k
                    blk%q(i,j,kp,VAR_W,b) = blk%q(i,j,kp,VAR_W,b) + phiC*cHi3*mu_w_kp

                    ! Only this colour's cells hold a non-zero increment (phi was
                    ! zeroed before the sweep), so each interface face receives
                    ! the red half in one pass and the black half in the other --
                    ! the two sum to the full two-sided correction, no double
                    ! counting.
                    if (store) phi(i,j,k,b) = phiC
                END DO
            END DO
        END DO
        END DO
#ifdef USE_OPENMP_OFFLOAD
        !$omp end target teams distribute parallel do
#endif
    end subroutine redblack_sweep

end module pressure_solver
