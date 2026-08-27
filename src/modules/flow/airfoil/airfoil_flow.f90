module airfoil_flow
    ! Quasi-2D immersed-boundary airfoil/cylinder case ([case] name = airfoil,
    ! docs/next_session_airfoil.md, phases A1/A2). The angle of attack turns
    ! into a Dirichlet freestream on the inlet faces (x_min, y_min, y_max) and
    ! a Dirichlet-pressure outlet at x_max, composed purely of PATCH TYPES
    ! (resolve_face_bcs derives the per-variable rows); geometry comes from
    ! the standard file-based IBM path. Runtime statistics are C_L(t)/C_D(t)
    ! from the CONTROL-VOLUME momentum budget over the box `[case.airfoil]
    ! cv_box` (see cv_forces). The penalization integral F = int coef*u dV
    ! that this case used to report is GONE: it is exact bookkeeping only
    ! while the solid interior is present, and production runs remove the
    ! buried core ([blocks] remove_solid, the default), which puts its share
    ! of the pressure-dominated loading outside the coef bookkeeping. Recover
    ! it from history if ever needed; do not reinstate it as a fallback.
    use, intrinsic :: iso_c_binding
    use :: flow_case_base, only: case_type
    use :: generic_flow, only: set_generic_defaults
    use :: init, only: dns_type, grid_type, VAR_U, VAR_V, VAR_W, VAR_P
    use :: blocks, only: block_set_type
    use :: boundary, only: boundary_type, boundary_face_id, PATCH_INLET, PATCH_OUTLET
    use :: pressure_solver, only: pressure_solver_type
    use :: ibmm, only: ibm_type
    use :: turbulence, only: turb_type, TURB_NONE
    use :: comm, only: comm_type, comm_allreduce_sum, comm_allreduce_max
    use :: case_config_helpers, only: next_config_entry, to_lower, clean_config_string
    implicit none

    private

    character(len=*), parameter :: AIRFOIL_CASE_NAME = "airfoil"
    ! per-block CV accumulators: momC momL fluxC fluxL
    integer, parameter :: NCV = 4
    ! per-block accumulators of the p_inf pass: p*dA and dA on the upstream border
    integer, parameter :: NPI = 2

    type, extends(case_type), public :: airfoil_case_type
        real(C_DOUBLE) :: aoa = 0.0d0            ! angle of attack [deg]
        real(C_DOUBLE) :: u_inf = 1.0d0
        real(C_DOUBLE) :: chord = 1.0d0
        ! [case.airfoil] span = z (default) | y: the periodic (extrusion)
        ! direction. span = y is the xz-quadtree orientation
        ! (docs/next_session_refine2d.md): chord along x, LIFT along z,
        ! freestream (U cos a, 0, U sin a), inlets x_min/z_min/z_max,
        ! outlet x_max — [blocks] refine_dims = xz then never refines the
        ! span. liftDim/spanDim are derived at config read.
        integer :: spanDim = 3
        integer :: liftDim = 2
        integer :: force_sample_interval = 10
        character(len=256) :: runtime_file = "forces.txt"
        logical :: header_written = .false.
        ! [case.airfoil] cv_box = c0 c1 l0 l1: the control volume for the
        ! runtime force budget, in physical coordinates along the CHORD (x)
        ! and LIFT axes; the span is the full periodic extent. REQUIRED
        ! unless force_sample_interval = 0 (there is no sane default: the
        ! case does not know where the immersed body sits, that comes from
        ! the IBM file).
        ! Borders are SNAPPED at setup to the node line of the coarsest
        ! level they cross, so every border coincides with a cell face at
        ! every level it crosses and the box stays closed.
        real(C_DOUBLE) :: cv_box(4) = 0.0d0
        logical :: cv_set = .false.
        ! d/dt of the momentum inside the box, differenced between the last
        ! two SAMPLES (hence over force_sample_interval steps).
        real(C_DOUBLE) :: mom_prev(2) = 0.0d0
        real(C_DOUBLE) :: t_prev = 0.0d0
        logical :: have_prev = .false.
        ! [case.airfoil] steady_tol: stop the run once the flow is steady.
        ! The measure is the unsteady term of the budget expressed in the
        ! SAME units as the reported coefficients, |2 dmom/dt| / qref, so
        ! the threshold is read like a C_L/C_D increment (1e-4 = the
        ! unsteady term no longer moves the fourth digit). 0 disables it.
        ! It is only testable once d/dt exists, i.e. from the second sample.
        real(C_DOUBLE) :: steady_tol = 0.0d0
        ! Consecutive qualifying samples required. An oscillating flow's
        ! dmom/dt passes through ZERO twice per shedding period, so a
        ! single-sample test would stop a perfectly unsteady run at a
        ! turning point; requiring a run of samples rejects that.
        integer :: steady_samples = 3
        integer :: steady_hits = 0
    contains
        procedure :: read_config => airfoil_read_config
        procedure :: apply_defaults => airfoil_apply_defaults
        procedure :: setup_after_grid => airfoil_setup_after_grid
        procedure :: initialise_fields => airfoil_initialise_fields
        procedure :: after_step => airfoil_after_step
        procedure :: finalize => airfoil_finalize
    end type airfoil_case_type

    public :: create_airfoil_case, AIRFOIL_CASE_NAME

contains

    subroutine create_airfoil_case(flow)
        class(case_type), allocatable, intent(out) :: flow

        allocate(airfoil_case_type :: flow)
        flow%name = AIRFOIL_CASE_NAME
    end subroutine create_airfoil_case

    ! Freestream direction cosines from the angle of attack: ux along the
    ! chord (x), ul along the LIFT direction (y for span = z, z for
    ! span = y).
    subroutine freestream(this, ux, ul)
        class(airfoil_case_type), intent(in) :: this
        real(C_DOUBLE), intent(out) :: ux, ul
        real(C_DOUBLE) :: a

        a = this%aoa*(4.0d0*atan(1.0d0))/180.0d0
        ux = this%u_inf*cos(a)
        ul = this%u_inf*sin(a)
    end subroutine freestream

    subroutine airfoil_apply_defaults(this, dns, g, bc, c, ps)
        class(airfoil_case_type), intent(inout) :: this
        type(dns_type), intent(inout) :: dns
        type(grid_type), intent(inout) :: g
        type(boundary_type), intent(inout) :: bc
        type(comm_type), intent(inout) :: c
        type(pressure_solver_type), intent(inout) :: ps

        real(C_DOUBLE) :: ux, ul, vel(3)
        integer :: face_id, dir, side

        call set_generic_defaults(dns, g, bc, c, ps)
        this%name = AIRFOIL_CASE_NAME
        call freestream(this, ux, ul)
        vel = 0.0d0
        vel(1) = ux
        vel(this%liftDim) = ul

        dns%ibm_enabled = .true.
        dns%forcing = 0.0d0             ! no volume forcing in airfoil runs
        dns%initial_velocity = vel

        ! Freestream composition in ONE vocabulary -- patch types (the ini is
        ! parsed AFTER apply_defaults, so explicit [boundary] keys still win):
        ! x_min and both LIFT-direction faces are Dirichlet inlets carrying
        ! (U cos a) e_x + (U sin a) e_lift, x_max the Dirichlet-pressure
        ! outlet; the SPAN direction is periodic (quasi-2D, n_span = nb).
        ! Large |aoa| needs a taller domain -- the lift-direction faces are
        ! far-field Dirichlet, the standard penalization freestream.
        bc%isPeriodic = .false.
        bc%isPeriodic(this%spanDim) = .true.
        do dir = 1, 3
            if (dir == this%spanDim) cycle
            do side = 0, 1
                face_id = boundary_face_id(dir, side)
                if (dir == 1 .and. side == 1) then
                    bc%facePatchType(face_id) = PATCH_OUTLET
                else
                    bc%facePatchType(face_id) = PATCH_INLET
                    bc%faceBcDefaultValue(VAR_U,face_id) = vel(1)
                    bc%faceBcDefaultValue(VAR_V,face_id) = vel(2)
                    bc%faceBcDefaultValue(VAR_W,face_id) = vel(3)
                end if
            end do
        end do
    end subroutine airfoil_apply_defaults

    subroutine airfoil_setup_after_grid(this, blk, dns, g, bc, c)
        class(airfoil_case_type), intent(inout) :: this
        type(block_set_type), intent(inout) :: blk
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        type(boundary_type), intent(in) :: bc
        type(comm_type), intent(in) :: c

        integer :: ib

        ! The control-volume budget is the ONLY force statistic, so the box
        ! is required whenever forces are sampled -- there is no sane
        ! default (the case does not know where the immersed body sits).
        ! force_sample_interval = 0 is the explicit opt-out.
        if (.not. this%cv_set) then
            if (this%force_sample_interval <= 0) return
            error stop "[case.airfoil] cv_box is required for runtime forces &
                &(set it to c0 c1 l0 l1 around the body, or set &
                &force_sample_interval = 0 to disable forces)"
        end if
        do ib = 1, 4
            call snap_border(this, blk, c, ib)
        end do
        if (this%cv_box(1) >= this%cv_box(2) .or. this%cv_box(3) >= this%cv_box(4)) &
            error stop "[case.airfoil] cv_box: empty control volume after snapping"
        if (c%has_terminal) write(*,'(A,4(1X,ES14.6))') &
            " airfoil CV box snapped to (c0 c1 l0 l1):", this%cv_box
    end subroutine airfoil_setup_after_grid

    ! Snap one CV border onto a cell face of the COARSEST level it crosses.
    ! A coarse node is a node at every finer level (the level lines subdivide
    ! by midpoints), so the snapped border is a face everywhere it runs and
    ! the box closes exactly -- without this a border laid on a fine-only
    ! face would simply have no face to integrate where it crosses coarser
    ! blocks. Reduced with max() alone: first the minimum level via its
    ! negative, then the face coordinate (every block at that level yields
    ! the same candidate, so the max returns it).
    subroutine snap_border(this, blk, c, ib)
        class(airfoil_case_type), intent(inout) :: this
        type(block_set_type), intent(in) :: blk
        type(comm_type), intent(in) :: c
        integer, intent(in) :: ib

        real(C_DOUBLE) :: r(1), coord, cand, dbest, d, lo, hi, tlo, thi, t0, t1, cf
        integer :: b, i, n, minlev, var, pass
        logical :: chord, spanY, cross

        coord = this%cv_box(ib)
        chord = (ib <= 2)
        spanY = (this%spanDim == 2)
        if (chord) then
            var = VAR_U; t0 = this%cv_box(3); t1 = this%cv_box(4)
        else
            var = merge(VAR_W, VAR_V, spanY); t0 = this%cv_box(1); t1 = this%cv_box(2)
        end if
        n = merge(int(blk%nb(1)), merge(int(blk%nb(3)), int(blk%nb(2)), spanY), chord)

        do pass = 1, 2
            r(1) = -1.0d30
            do b = 1, blk%nBlocks
                if (pass == 2 .and. blk%level(b) /= minlev) cycle
                ! does the border run through this block, inside the box?
                if (chord) then
                    lo = blk%x(1,var,b); hi = blk%x(n+1,var,b)
                    tlo = merge(blk%z(1,var,b), blk%y(1,var,b), spanY)
                    thi = merge(blk%z(int(blk%nb(3)),var,b), blk%y(int(blk%nb(2)),var,b), spanY)
                else
                    lo = merge(blk%z(1,var,b), blk%y(1,var,b), spanY)
                    hi = merge(blk%z(n+1,var,b), blk%y(n+1,var,b), spanY)
                    tlo = blk%x(1,var,b); thi = blk%x(int(blk%nb(1)),var,b)
                end if
                if (coord < lo .or. coord > hi) cycle
                if (thi < t0 .or. tlo > t1) cycle
                if (pass == 1) then
                    r(1) = max(r(1), -real(blk%level(b), C_DOUBLE))
                else
                    dbest = 1.0d30; cand = coord
                    do i = 1, n + 1
                        if (chord) then
                            cf = blk%x(i,var,b)
                        else
                            cf = merge(blk%z(i,var,b), blk%y(i,var,b), spanY)
                        end if
                        d = abs(cf - coord)
                        if (d < dbest) then
                            dbest = d; cand = cf
                        end if
                    end do
                    r(1) = max(r(1), cand)
                end if
            end do
            call comm_allreduce_max(c, r)
            if (pass == 1) then
                if (r(1) < -1.0d29) error stop "[case.airfoil] cv_box: a border crosses no block"
                minlev = nint(-r(1))
            else
                this%cv_box(ib) = r(1)
            end if
        end do
    end subroutine snap_border

    ! Uniform freestream everywhere (impulsive start; the IBM damps the body
    ! interior within the first steps).
    subroutine airfoil_initialise_fields(this, blk, dns, g, bc, c)
        class(airfoil_case_type), intent(inout) :: this
        type(block_set_type), intent(inout) :: blk
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        type(boundary_type), intent(in) :: bc
        type(comm_type), intent(in) :: c

        real(C_DOUBLE) :: ux, ul, vel(3)

        call freestream(this, ux, ul)
        vel = 0.0d0
        vel(1) = ux
        vel(this%liftDim) = ul
        blk%q(:,:,:,VAR_U,:) = vel(1)
        blk%q(:,:,:,VAR_V,:) = vel(2)
        blk%q(:,:,:,VAR_W,:) = vel(3)
    end subroutine airfoil_initialise_fields

    ! Sample C_L(t)/C_D(t) on the END-OF-STEP field -- after the last
    ! projection exchange, so every halo the budget reads is current.
    subroutine airfoil_after_step(this, blk, dns, g, c, ibm, turb)
        class(airfoil_case_type), intent(inout) :: this
        type(block_set_type), intent(inout) :: blk
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        type(comm_type), intent(in) :: c
        type(ibm_type), intent(in) :: ibm
        type(turb_type), intent(in) :: turb

        if (this%force_sample_interval <= 0) return
        if (modulo(int(dns%step_current), this%force_sample_interval) /= 0) return
        call cv_forces(this, blk, dns, c, turb)
    end subroutine airfoil_after_step

    ! CONTROL-VOLUME momentum budget -- the PRODUCTION force statistic.
    !
    !   F = - d/dt int_V u dV - oint_S [ u (u.n) + (p - p_inf) n - tau.n ] dS
    !   tau = (nu + nu_t) (grad u + grad u^T)
    !
    ! over the box [c0,c1] x [l0,l1] (chord x lift) spanning the full
    ! periodic extent, so the two span faces cancel and only the four
    ! lateral borders are integrated. It never touches the body, so it is
    ! INDIFFERENT to [blocks] remove_solid -- which is why it replaced the
    ! penalization integral outright rather than being an alternative to it.
    !
    ! FLUX-EXACT: each border is snapped at setup to a face of the coarsest
    ! level it crosses, and only INTERIOR faces (i = 1..nb, the west face of
    ! cell i) are integrated. Every physical face is then counted exactly
    ! once -- by the block on its east side, at that block's own level -- so
    ! a border crossing a 2:1 interface is neither double counted nor
    ! missed, and a uniform field integrates to exactly zero. Halos supply
    ! the i-1/i+1 neighbours, so the velocity interpolations and gradients
    ! are central everywhere -- the block-edge one-sided fallbacks of the
    ! offline tools/cv_forces.py are not needed. Only the PRESSURE is
    ! one-sided on a face that is a block's low edge, because that halo
    ! carries the blended 2:1 ghost rather than the neighbour cell value
    ! (see the branch below; it costs ~5e-5 in C_D on the Re 40 cylinder).
    !
    ! p_inf (the area-weighted mean pressure of the upstream border) is
    ! subtracted PER FACE, in a first reduction pass over that border alone.
    ! For a closed box oint n dS = 0, so a constant cancels analytically and
    ! it would be tempting to subtract it from the assembled signed areas
    ! instead -- that is catastrophic in floating point: the stored pn
    ! carries the run's accumulated projection drift (O(1e3) on the Re 40
    ! cylinder against a force of O(0.2)), so the per-face terms must be
    ! O(p - p_inf) for the sum to have any significant digits left.
    !
    ! The unsteady term is differenced between the last two samples, so it
    ! costs one extra momentum reduction and no extra state: at the first
    ! sample after a (re)start it is unknown and the budget is reported
    ! without it (flagged by dmdt = 0 in the file).
    subroutine cv_forces(this, blk, dns, c, turb)
        class(airfoil_case_type), intent(inout) :: this
        type(block_set_type), intent(inout) :: blk
        type(dns_type), intent(in) :: dns
        type(comm_type), intent(in) :: c
        type(turb_type), intent(in) :: turb

        real(C_DOUBLE), allocatable :: fb(:,:)
        real(C_DOUBLE) :: acc(NCV), p_inf, fc, fl, dmdt(2), a, cl, cd, qref
        real(C_DOUBLE) :: c0, c1, l0, l1, nu, tol, coord, sgn
        real(C_DOUBLE) :: lc, cc, dA, un, pf, ut, nue, dundn, dundt, dutdn
        real(C_DOUBLE) :: tnn, tnt, fn, ft, dl2, wcp, wcm, vol
        real(C_DOUBLE) :: sC, sL, fC_, fL_, unsteady
        ! plain local copy: mapping a component of the polymorphic `this`
        ! into a target region is not portable
        real(C_DOUBLE) :: box(4)
        integer(C_INT) :: i, j, k, b, m, ii, nbx, nby, nbz, nBlocks, ib, varL
        integer(C_INT) :: dj, dk
        logical :: spanY, hasNut, have_dmdt

        nbx = blk%nb(1); nby = blk%nb(2); nbz = blk%nb(3)
        nBlocks = blk%nBlocks
        box = this%cv_box
        c0 = box(1); c1 = box(2)
        l0 = box(3); l1 = box(4)
        nu = 1.0d0/dns%re
        spanY = (this%spanDim == 2)
        if (spanY) then
            varL = VAR_W; dj = 0_C_INT; dk = 1_C_INT
        else
            varL = VAR_V; dj = 1_C_INT; dk = 0_C_INT
        end if
        hasNut = (turb%model /= TURB_NONE)
        ! a face matches a border when it is within this fraction of the
        ! finest spacing -- borders are snapped, so the match is exact up to
        ! round-off in the node-line arithmetic
        tol = 1.0d-9*max(1.0d0, abs(c1 - c0))

        ! ---- pass 1: the reference pressure, from the upstream border ----
        p_inf = cv_pinf(blk, c, box, tol, spanY, nbx, nby, nbz)

        allocate(fb(NCV, nBlocks))

#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute &
        !$omp& map(to: nbx, nby, nbz, c0, c1, l0, l1, nu, tol, varL, dj, dk, spanY, hasNut, &
        !$omp&        p_inf, box, blk%x, blk%y, blk%z, blk%d1x, blk%d1y, blk%d1z, &
        !$omp&        blk%q, turb%nut) &
        !$omp& map(from: fb(1:NCV,1:nBlocks)) &
        !$omp& private(i,j,k,m,ii,ib,sC,sL,fC_,fL_,coord,sgn,lc,cc,dA, &
        !$omp&         un,pf,ut,nue,dundn,dundt,dutdn,tnn,tnt,fn,ft,dl2,wcp,wcm,vol)
#endif
        do b = 1_C_INT, nBlocks
            sC = 0.0d0; sL = 0.0d0
            fC_ = 0.0d0; fL_ = 0.0d0

            ! ---- momentum inside the box (each component on its own grid) ----
            do k = 1_C_INT, nbz
                do j = 1_C_INT, nby
                    do i = 1_C_INT, nbx
                        lc = merge(blk%z(k,VAR_U,b), blk%y(j,VAR_U,b), spanY)
                        if (blk%x(i,VAR_U,b) > c0 .and. blk%x(i,VAR_U,b) < c1 .and. &
                            lc > l0 .and. lc < l1) then
                            vol = 1.0d0/(blk%d1x(i,VAR_U,b)*blk%d1y(j,VAR_U,b)*blk%d1z(k,VAR_U,b))
                            sC = sC + blk%q(i,j,k,VAR_U,b)*vol
                        end if
                        lc = merge(blk%z(k,varL,b), blk%y(j,varL,b), spanY)
                        if (blk%x(i,varL,b) > c0 .and. blk%x(i,varL,b) < c1 .and. &
                            lc > l0 .and. lc < l1) then
                            vol = 1.0d0/(blk%d1x(i,varL,b)*blk%d1y(j,varL,b)*blk%d1z(k,varL,b))
                            sL = sL + blk%q(i,j,k,varL,b)*vol
                        end if
                    end do
                end do
            end do

            ! ---- the four lateral borders ----
            do ib = 1_C_INT, 4_C_INT
                coord = box(ib)
                sgn = merge(-1.0d0, 1.0d0, ib == 1_C_INT .or. ib == 3_C_INT)
                if (ib <= 2_C_INT) then
                    ! CHORD-normal border: u faces
                    ii = -1_C_INT
                    do i = 1_C_INT, nbx
                        if (abs(blk%x(i,VAR_U,b) - coord) < tol) ii = i
                    end do
                    if (ii > 0_C_INT) then
                        do k = 1_C_INT, nbz
                            do j = 1_C_INT, nby
                                lc = merge(blk%z(k,VAR_U,b), blk%y(j,VAR_U,b), spanY)
                                if (lc <= l0 .or. lc >= l1) cycle
                                dA = 1.0d0/(blk%d1y(j,VAR_U,b)*blk%d1z(k,VAR_U,b))
                                un = blk%q(ii,j,k,VAR_U,b)
                                ! Interpolate from the OWNING side only. The
                                ! halo pressure is not the neighbour's cell
                                ! value at a 2:1 interface -- it is the blended
                                ! ghost (2 p_C + p_f)/3 the projection needs --
                                ! so a central average across the face reads a
                                ! modified quantity. Second-order one-sided
                                ! from the two interior cells east of the face
                                ! whenever the halo would be involved.
                                if (ii >= 2_C_INT) then
                                    pf = 0.5d0*(blk%q(ii-1,j,k,VAR_P,b) + blk%q(ii,j,k,VAR_P,b))
                                else
                                    pf = 1.5d0*blk%q(ii,j,k,VAR_P,b) - 0.5d0*blk%q(ii+1,j,k,VAR_P,b)
                                end if
                                ut = 0.25d0*(blk%q(ii-1,j,k,varL,b) + blk%q(ii,j,k,varL,b) &
                                           + blk%q(ii-1,j+dj,k+dk,varL,b) + blk%q(ii,j+dj,k+dk,varL,b))
                                nue = nu
                                if (hasNut) nue = nue &
                                    + 0.5d0*(turb%nut(ii-1,j,k,b) + turb%nut(ii,j,k,b))
                                dundn = (blk%q(ii+1,j,k,VAR_U,b) - blk%q(ii-1,j,k,VAR_U,b)) &
                                      / (blk%x(ii+1,VAR_U,b) - blk%x(ii-1,VAR_U,b))
                                dl2 = merge(blk%z(k+dk,VAR_U,b) - blk%z(k-dk,VAR_U,b), &
                                            blk%y(j+dj,VAR_U,b) - blk%y(j-dj,VAR_U,b), spanY)
                                dundt = (blk%q(ii,j+dj,k+dk,VAR_U,b) &
                                       - blk%q(ii,j-dj,k-dk,VAR_U,b))/dl2
                                wcp = 0.5d0*(blk%q(ii,j,k,varL,b) + blk%q(ii,j+dj,k+dk,varL,b))
                                wcm = 0.5d0*(blk%q(ii-1,j,k,varL,b) + blk%q(ii-1,j+dj,k+dk,varL,b))
                                dutdn = (wcp - wcm)*blk%d1x(ii,VAR_U,b)
                                tnn = 2.0d0*nue*dundn
                                tnt = nue*(dundt + dutdn)
                                fn = (un*un + (pf - p_inf) - tnn)*dA
                                ft = (un*ut - tnt)*dA
                                fC_ = fC_ - sgn*fn
                                fL_ = fL_ - sgn*ft
                            end do
                        end do
                    end if
                else
                    ! LIFT-normal border: varL faces; the normal index is k
                    ! when the span is y, j otherwise
                    ii = -1_C_INT
                    do m = 1_C_INT, merge(nbz, nby, spanY)
                        lc = merge(blk%z(m,varL,b), blk%y(m,varL,b), spanY)
                        if (abs(lc - coord) < tol) ii = m
                    end do
                    if (ii > 0_C_INT) then
                        do m = 1_C_INT, merge(nby, nbz, spanY)
                            do i = 1_C_INT, nbx
                                j = merge(m, ii, spanY)
                                k = merge(ii, m, spanY)
                                cc = blk%x(i,varL,b)
                                if (cc <= c0 .or. cc >= c1) cycle
                                dA = merge(1.0d0/(blk%d1x(i,varL,b)*blk%d1y(j,varL,b)), &
                                           1.0d0/(blk%d1x(i,varL,b)*blk%d1z(k,varL,b)), spanY)
                                un = blk%q(i,j,k,varL,b)
                                if (ii >= 2_C_INT) then
                                    pf = 0.5d0*(blk%q(i,j-dj,k-dk,VAR_P,b) + blk%q(i,j,k,VAR_P,b))
                                else
                                    pf = 1.5d0*blk%q(i,j,k,VAR_P,b) &
                                       - 0.5d0*blk%q(i,j+dj,k+dk,VAR_P,b)
                                end if
                                ut = 0.25d0*(blk%q(i,j-dj,k-dk,VAR_U,b) + blk%q(i,j,k,VAR_U,b) &
                                           + blk%q(i+1,j-dj,k-dk,VAR_U,b) + blk%q(i+1,j,k,VAR_U,b))
                                nue = nu
                                if (hasNut) nue = nue &
                                    + 0.5d0*(turb%nut(i,j-dj,k-dk,b) + turb%nut(i,j,k,b))
                                dl2 = merge(blk%z(k+dk,varL,b) - blk%z(k-dk,varL,b), &
                                            blk%y(j+dj,varL,b) - blk%y(j-dj,varL,b), spanY)
                                dundn = (blk%q(i,j+dj,k+dk,varL,b) &
                                       - blk%q(i,j-dj,k-dk,varL,b))/dl2
                                dundt = (blk%q(i+1,j,k,varL,b) - blk%q(i-1,j,k,varL,b)) &
                                      / (blk%x(i+1,varL,b) - blk%x(i-1,varL,b))
                                wcp = 0.5d0*(blk%q(i,j,k,VAR_U,b) + blk%q(i+1,j,k,VAR_U,b))
                                wcm = 0.5d0*(blk%q(i,j-dj,k-dk,VAR_U,b) &
                                           + blk%q(i+1,j-dj,k-dk,VAR_U,b))
                                dutdn = (wcp - wcm)*merge(blk%d1z(k,varL,b), blk%d1y(j,varL,b), spanY)
                                tnn = 2.0d0*nue*dundn
                                tnt = nue*(dundt + dutdn)
                                fn = (un*un + (pf - p_inf) - tnn)*dA
                                ft = (un*ut - tnt)*dA
                                fL_ = fL_ - sgn*fn
                                fC_ = fC_ - sgn*ft
                            end do
                        end do
                    end if
                end if
            end do

            fb(1,b) = sC;  fb(2,b) = sL
            fb(3,b) = fC_; fb(4,b) = fL_
        end do
#ifdef USE_OPENMP_OFFLOAD
        !$omp end target teams distribute
#endif

        call cv_reduce(blk, c, fb, NCV, acc)

        fc = acc(3)
        fl = acc(4)

        dmdt = 0.0d0
        have_dmdt = this%have_prev .and. dns%t_current > this%t_prev
        if (have_dmdt) then
            dmdt(1) = (acc(1) - this%mom_prev(1))/(dns%t_current - this%t_prev)
            dmdt(2) = (acc(2) - this%mom_prev(2))/(dns%t_current - this%t_prev)
        end if
        this%mom_prev(1) = acc(1); this%mom_prev(2) = acc(2)
        this%t_prev = dns%t_current
        this%have_prev = .true.

        fc = fc - dmdt(1)
        fl = fl - dmdt(2)

        ! Drag along the freestream, lift normal to it (rho = 1):
        ! C = 2 F.e / (U_inf^2 * chord * L_span).
        a = this%aoa*(4.0d0*atan(1.0d0))/180.0d0
        qref = this%u_inf**2*this%chord*dns%leng(this%spanDim)
        cd = 2.0d0*( fc*cos(a) + fl*sin(a))/qref
        cl = 2.0d0*(-fc*sin(a) + fl*cos(a))/qref

        ! The steady measure travels with the coefficients: it is what
        ! `steady_tol` is compared against, so a run's trace carries the
        ! evidence for where its own tolerance should sit (and, afterwards,
        ! why it stopped when it did). Zero until d/dt exists.
        unsteady = 0.0d0
        if (have_dmdt) unsteady = 2.0d0*max(abs(dmdt(1)), abs(dmdt(2)))/qref

        if (c%has_terminal) call append_forces(this, dns, cl, cd, unsteady)

        ! Steady-state stop. The measure is the unsteady term in the same
        ! units as the reported coefficients, so `steady_tol` reads like a
        ! C_L/C_D increment. dmdt comes out of an exact allreduce, so every
        ! rank decides identically -- this must stay OUTSIDE has_terminal or
        ! the ranks part company at the loop exit.
        if (this%steady_tol > 0.0d0 .and. have_dmdt) then
            if (unsteady < this%steady_tol) then
                this%steady_hits = this%steady_hits + 1
            else
                this%steady_hits = 0
            end if
            if (this%steady_hits >= this%steady_samples) then
                this%stop_requested = .true.
                if (c%has_terminal) write(*,'(A,ES12.4,A,I0,A,ES12.4,A)') &
                    " airfoil: steady state reached (|2 dmom/dt|/qref =", unsteady, &
                    " over ", this%steady_samples, " samples, tol ", this%steady_tol, &
                    "); writing the final field and stopping"
            end if
        end if
    end subroutine cv_forces

    ! Pass 1 of the CV budget: p_inf, the area-weighted mean face pressure of
    ! the UPSTREAM (chord-low) border. It has to be subtracted PER FACE (see
    ! cv_forces), which means it must be known before the flux sum -- hence
    ! its own reduction. Same face set, same ownership rule and the same
    ! pressure interpolation as the flux pass, and the same deterministic
    ! reduction: p_inf enters the reported force, so it must be rank-count
    ! independent too.
    function cv_pinf(blk, c, box, tol, spanY, nbx, nby, nbz) result(p_inf)
        type(block_set_type), intent(inout) :: blk
        type(comm_type), intent(in) :: c
        real(C_DOUBLE), intent(in) :: box(4), tol
        logical, intent(in) :: spanY
        integer(C_INT), intent(in) :: nbx, nby, nbz
        real(C_DOUBLE) :: p_inf

        real(C_DOUBLE), allocatable :: pb(:,:)
        real(C_DOUBLE) :: acc(NPI), c0, l0, l1, pA_, aW_, lc, dA, pf
        integer(C_INT) :: i, j, k, b, ii, nBlocks

        c0 = box(1); l0 = box(3); l1 = box(4)
        nBlocks = blk%nBlocks
        allocate(pb(NPI, nBlocks))

#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute &
        !$omp& map(to: nbx, nby, nbz, c0, l0, l1, tol, spanY, &
        !$omp&        blk%x, blk%y, blk%z, blk%d1y, blk%d1z, blk%q) &
        !$omp& map(from: pb(1:NPI,1:nBlocks)) &
        !$omp& private(i,j,k,ii,pA_,aW_,lc,dA,pf)
#endif
        do b = 1_C_INT, nBlocks
            pA_ = 0.0d0; aW_ = 0.0d0
            ii = -1_C_INT
            do i = 1_C_INT, nbx
                if (abs(blk%x(i,VAR_U,b) - c0) < tol) ii = i
            end do
            if (ii > 0_C_INT) then
                do k = 1_C_INT, nbz
                    do j = 1_C_INT, nby
                        lc = merge(blk%z(k,VAR_U,b), blk%y(j,VAR_U,b), spanY)
                        if (lc <= l0 .or. lc >= l1) cycle
                        dA = 1.0d0/(blk%d1y(j,VAR_U,b)*blk%d1z(k,VAR_U,b))
                        if (ii >= 2_C_INT) then
                            pf = 0.5d0*(blk%q(ii-1,j,k,VAR_P,b) + blk%q(ii,j,k,VAR_P,b))
                        else
                            pf = 1.5d0*blk%q(ii,j,k,VAR_P,b) - 0.5d0*blk%q(ii+1,j,k,VAR_P,b)
                        end if
                        pA_ = pA_ + pf*dA
                        aW_ = aW_ + dA
                    end do
                end do
            end if
            pb(1,b) = pA_; pb(2,b) = aW_
        end do
#ifdef USE_OPENMP_OFFLOAD
        !$omp end target teams distribute
#endif

        call cv_reduce(blk, c, pb, NPI, acc)
        p_inf = 0.0d0
        if (acc(2) > 0.0d0) p_inf = acc(1)/acc(2)
    end function cv_pinf

    ! Deterministic reduction of per-block partials: scatter into the global
    ! block table (one contributor per entry, so the allreduce is exact) and
    ! sum in global-id order on every rank -- the result is independent of
    ! the rank count by construction.
    subroutine cv_reduce(blk, c, fb, n, acc)
        type(block_set_type), intent(in) :: blk
        type(comm_type), intent(in) :: c
        real(C_DOUBLE), intent(in) :: fb(:,:)
        integer, intent(in) :: n
        real(C_DOUBLE), intent(out) :: acc(n)

        real(C_DOUBLE), allocatable :: fbg(:)
        integer :: b

        allocate(fbg(n*blk%nBlocksGlobal))
        fbg = 0.0d0
        do b = 1, int(blk%nBlocks)
            fbg(n*blk%globalId(b)+1:n*blk%globalId(b)+n) = fb(:,b)
        end do
        call comm_allreduce_sum(c, fbg)
        acc = 0.0d0
        do b = 1, int(blk%nBlocksGlobal)
            acc = acc + fbg(n*(b-1)+1:n*(b-1)+n)
        end do
    end subroutine cv_reduce

    subroutine append_forces(this, dns, cl, cd, unsteady)
        class(airfoil_case_type), intent(inout) :: this
        type(dns_type), intent(in) :: dns
        real(C_DOUBLE), intent(in) :: cl, cd
        ! |2 dmom/dt|/qref, the budget's own unsteady term: the quantity
        ! `steady_tol` thresholds, logged so convergence is auditable.
        real(C_DOUBLE), intent(in) :: unsteady

        integer :: unit, stat
        character(len=*), parameter :: header = "iteration time cl cd dmomdt"

        if (.not. this%header_written) then
            write(*,'(A)') header
            open(newunit=unit, file=trim(this%runtime_file), status="replace", &
                action="write", iostat=stat)
            if (stat == 0) then
                write(unit,'(A)') header
                close(unit)
            else
                print *, "warning: could not open runtime data file: ", trim(this%runtime_file)
            end if
            this%header_written = .true.
        end if

        write(*,'(I10,4(1X,ES16.8))') int(dns%step_current), dns%t_current, cl, cd, unsteady
        open(newunit=unit, file=trim(this%runtime_file), status="old", &
            position="append", action="write", iostat=stat)
        if (stat == 0) then
            write(unit,'(I10,4(1X,ES16.8))') int(dns%step_current), dns%t_current, cl, cd, unsteady
            close(unit)
        end if
    end subroutine append_forces

    subroutine airfoil_finalize(this, dns, g, c)
        class(airfoil_case_type), intent(inout) :: this
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        type(comm_type), intent(in) :: c
    end subroutine airfoil_finalize

    subroutine airfoil_read_config(this, input_file, has_terminal)
        class(airfoil_case_type), intent(inout) :: this
        character(len=*), intent(in) :: input_file
        logical, intent(in), optional :: has_terminal

        integer :: unit, stat, line_no
        character(len=512) :: key, value
        character(len=64) :: section
        logical :: exists, ok

        section = ""
        line_no = 0
        inquire(file=trim(input_file), exist=exists)
        if (.not. exists) return

        open(newunit=unit, file=trim(input_file), status="old", action="read", iostat=stat)
        if (stat /= 0) return

        do
            call next_config_entry(unit, section, key, value, line_no, ok)
            if (.not. ok) exit
            if (trim(section) /= "case.airfoil") cycle
            call apply_airfoil_case_value(this, to_lower(key), value, line_no, has_terminal)
        end do

        close(unit)
    end subroutine airfoil_read_config

    subroutine apply_airfoil_case_value(this, key, value, line_no, has_terminal)
        type(airfoil_case_type), intent(inout) :: this
        character(len=*), intent(in) :: key, value
        integer, intent(in) :: line_no
        logical, intent(in), optional :: has_terminal

        integer :: int_value, stat
        real(C_DOUBLE) :: real_value
        logical :: terminal

        terminal = .true.
        if (present(has_terminal)) terminal = has_terminal

        select case (trim(key))
        case ("aoa")
            read(value, *, iostat=stat) real_value
            if (stat == 0) this%aoa = real_value
        case ("u_inf")
            read(value, *, iostat=stat) real_value
            if (stat == 0) this%u_inf = real_value
        case ("chord")
            read(value, *, iostat=stat) real_value
            if (stat == 0) this%chord = real_value
        case ("force_sample_interval")
            read(value, *, iostat=stat) int_value
            if (stat == 0) this%force_sample_interval = int_value
        case ("cv_box")
            ! c0 c1 l0 l1: the control volume for the runtime force budget.
            block
                real(C_DOUBLE) :: v(4)
                read(value, *, iostat=stat) v
                if (stat == 0) then
                    this%cv_box = v
                    this%cv_set = .true.
                else
                    print *, "[case.airfoil] cv_box needs 4 reals (c0 c1 l0 l1) at line", line_no
                    error stop "invalid [case.airfoil] cv_box"
                end if
            end block
        case ("steady_tol")
            ! stop once |2 dmom/dt|/qref (the unsteady term in coefficient
            ! units) stays below this for steady_samples consecutive samples
            read(value, *, iostat=stat) real_value
            if (stat == 0) this%steady_tol = real_value
        case ("steady_samples")
            read(value, *, iostat=stat) int_value
            if (stat == 0) then
                if (int_value < 1) then
                    print *, "[case.airfoil] steady_samples must be >= 1 at line", line_no
                    error stop "invalid [case.airfoil] steady_samples"
                end if
                this%steady_samples = int_value
            end if
        case ("runtime_file")
            this%runtime_file = clean_config_string(value)
        case ("span")
            select case (trim(to_lower(clean_config_string(value))))
            case ("z")
                this%spanDim = 3
                this%liftDim = 2
            case ("y")
                this%spanDim = 2
                this%liftDim = 3
            case default
                print *, "[case.airfoil] span must be z or y at line", line_no
                error stop "invalid [case.airfoil] span"
            end select
        case default
            if (terminal) print *, "warning: unknown airfoil case key on input line", line_no, ": ", trim(key)
        end select
    end subroutine apply_airfoil_case_value

end module airfoil_flow
