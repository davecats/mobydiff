!--------------------------!
!                          !
!     Passive scalars      !
!                          !
!--------------------------!
!
! Advected/diffused cell-centred scalars s_1..s_N (docs/history/next_session_scalar.md).
!
! DESIGN (increment S0/S1): a scalar is an EXTRA VARIABLE of blk%q, stored at
! the PRESSURE point, so
!   * the metric arrays never grow -- every scalar stencil reads the VAR_P
!     column of blk%x/d1x/...;
!   * the 2:1-interface halo exchange serves scalars for free: comm.f90 treats
!     any variable index > 3 as cell-centred (8-cell restrict average, injected
!     prolong, the (2 coarse + fine)/3 ghost blend), i.e. the pressure
!     treatment, and scalars ride the SAME MPI message as u,v,w,p;
!   * output/restart is one file, one collective write.
! `[scalar] count = 0` is bit-exact with a scalar-free build BY CONSTRUCTION:
! dns%nVar == NVAR reproduces every allocation shape, no kernel here is called,
! the boundary tables have no scalar column and no extra exchange is issued.
!
! LANDMINE: the scalar index DIFFERS between the arrays -- q slot VAR_S0+is,
! qs/oldrhs slot SCR_S0+is (pressure has no scratch plane).
!
! Per-scalar state lives in scalar_type as flat ALLOCATABLE arrays mapped to
! the device once (enter_scalar_data), the turb_type/sst_type/bodyforce_type
! pattern. There is deliberately NO fixed bound on the scalar count.
module scalar
    use, intrinsic :: iso_c_binding
    use :: init, only: dns_type, VAR_U, VAR_V, VAR_W, VAR_P, &
        VAR_S0, SCR_S0, SOLID_FACE_THRESHOLD
    use :: blocks, only: block_set_type, FACE_CLOSED
    use :: boundary, only: boundary_type, boundary_face_id, apply_bc, set_scalar_bc_rows, BC_DIRICHLET, &
        BC_NEUMANN, PATCH_GENERIC, PATCH_WALL, PATCH_INLET, PATCH_OUTLET
    use :: comm, only: comm_type, exchange_halos
    use :: turbulence, only: turb_type, turbulence_is_enabled
    use :: ibmm, only: penal_incr_factor, penal_state_factor, penal_state_minus_incr
    use :: scalar_types, only: SC_CONV_ADV, SC_CONV_DIV, SC_CONV_SKEW, SC_IBM_ADIABATIC, &
        SC_IBM_DIRICHLET, SC_INIT_LINEAR_Y, SC_SRC_VELOCITY, SRC_DIR_NAME, &
        scalar_conjugate_enabled, scalar_type, scalars_enabled
    use :: scalar_wallfn, only: jayatilleke_p, plain_face_diffusivity, thermal_yplus
    use :: scalar_conjugate, only: band_face_total, band_material, check_conjugate_refinement, conjugate_face_diffusivity, &
        conjugate_face_local_k, conjugate_face_total, conjugate_tangential, init_scalar_conjugate
    use :: scalar_config, only: deallocate_config, validate_conjugate_config
    implicit none
    private
    public :: destroy_scalar, enter_scalar_data, exit_scalar_data, init_scalar, init_scalar_fields
    public :: scalar_finish, scalar_sync, scalar_transport

contains

    ! Resolve the patch-derived boundary rows and build the per-block metric
    ! tables. Runs after the config/restart metadata is final (the patch types
    ! and the block set must both exist).
    subroutine init_scalar(sc, blk, bc, wall_function, has_terminal)
        type(scalar_type), intent(inout) :: sc
        type(block_set_type), intent(in) :: blk
        type(boundary_type), intent(inout) :: bc
        ! [rans] wall_treatment = wall_function: allocate the wall-cell y+
        ! field and report the thermal wall function's per-scalar constants.
        logical, intent(in) :: wall_function
        logical, intent(in) :: has_terminal

        integer :: is, dir, side, face_id, nx, ny, nz, i, b

        if (.not. scalars_enabled(sc)) return

        ! Patch-derived defaults, set-if-unset (the A0 face concept): a wall
        ! is adiabatic, an inlet holds the scalar's inlet value, an outlet
        ! and a generic patch are zero-gradient. BC_OUTFLOW is a
        ! face-staggered concept and does not apply to a cell-centred scalar.
        do dir = 1, 3
            do side = 0, 1
                face_id = boundary_face_id(dir, side)
                select case (bc%facePatchType(face_id))
                case (PATCH_INLET)
                    do is = 1, int(sc%n)
                        call resolve_scalar_row(sc, is, face_id, BC_DIRICHLET, &
                            sc%inlet(is), has_terminal, strict=.true.)
                    end do
                case (PATCH_OUTLET, PATCH_GENERIC)
                    do is = 1, int(sc%n)
                        call resolve_scalar_row(sc, is, face_id, BC_NEUMANN, 0.0d0, &
                            has_terminal, strict=.true.)
                    end do
                case (PATCH_WALL)
                    ! A wall is the ONE face kind where both scalar types are
                    ! physical: adiabatic (Neumann 0, the default) and
                    ! isothermal / fixed-concentration (Dirichlet) are equally
                    ! valid, so an explicit _type key is HONOURED here rather
                    ! than treated as a contradiction (Section 5 of
                    ! docs/history/next_session_scalar.md: "unless the ini gives a
                    ! type/value"). Inlet and outlet stay strict: an inlet must
                    ! impose a value, an outlet must be zero-gradient.
                    do is = 1, int(sc%n)
                        call resolve_scalar_row(sc, is, face_id, BC_NEUMANN, 0.0d0, &
                            has_terminal, strict=.false.)
                    end do
                end select
            end do
        end do
        ! The scalar rows join the boundary tables as their own columns, so
        ! apply_bc serves them exactly as it serves u, v, w, p (init order:
        ! init_boundary_faces, this, then enter_boundary_data maps them).
        do is = 1, int(sc%n)
            call set_scalar_bc_rows(bc, blk, int(VAR_S0) + is, sc%bcType(is,:), sc%bcValue(is,:))
        end do

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))
        if (allocated(sc%invDx)) deallocate(sc%invDx, sc%invDy, sc%invDz)
        allocate(sc%invDx(0:nx+1,blk%nBlocks), sc%invDy(0:ny+1,blk%nBlocks), &
                 sc%invDz(0:nz+1,blk%nBlocks))
        if (allocated(sc%cdx)) deallocate(sc%cdx, sc%cdy, sc%cdz)
        allocate(sc%cdx(0:nx+1,blk%nBlocks), sc%cdy(0:ny+1,blk%nBlocks), &
                 sc%cdz(0:nz+1,blk%nBlocks))
        do b = 1, int(blk%nBlocks)
            do i = 0, nx+1
                sc%invDx(i,b) = inv_delta(blk%x(i,VAR_P,b) - blk%x(i-1,VAR_P,b))
                sc%cdx(i,b) = inv_delta(blk%x(i+1,VAR_P,b) - blk%x(i-1,VAR_P,b))
            end do
            do i = 0, ny+1
                sc%invDy(i,b) = inv_delta(blk%y(i,VAR_P,b) - blk%y(i-1,VAR_P,b))
                sc%cdy(i,b) = inv_delta(blk%y(i+1,VAR_P,b) - blk%y(i-1,VAR_P,b))
            end do
            do i = 0, nz+1
                sc%invDz(i,b) = inv_delta(blk%z(i,VAR_P,b) - blk%z(i-1,VAR_P,b))
                sc%cdz(i,b) = inv_delta(blk%z(i+1,VAR_P,b) - blk%z(i-1,VAR_P,b))
            end do
        end do

        if (allocated(sc%nutNone)) deallocate(sc%nutNone)
        allocate(sc%nutNone(0:0,0:0,0:0,1))
        sc%nutNone = 0.0d0

        ! Conjugate signed distance: full size (ghost-inclusive) only with a
        ! conjugate scalar, a 1-cell dummy otherwise -- the nutNone idiom, so
        ! the device maps stay uniform and every access is mode-guarded. It
        ! is allocated HERE (before enter_scalar_data maps it) but FILLED
        ! later by init_scalar_conjugate, which needs the IBM coefficients.
        if (allocated(sc%phi)) deallocate(sc%phi)
        if (allocated(sc%vfrac)) deallocate(sc%vfrac)
        if (scalar_conjugate_enabled(sc)) then
            allocate(sc%phi(0:nx+1,0:ny+1,0:nz+1,blk%nBlocks))
            allocate(sc%vfrac(0:nx+1,0:ny+1,0:nz+1,blk%nBlocks))
        else
            allocate(sc%phi(0:0,0:0,0:0,1))
            allocate(sc%vfrac(0:0,0:0,0:0,1))
        end if
        sc%phi = 0.0d0
        sc%vfrac = 1.0d0

        ! Thermal wall function (S5a): the per-scalar constants are pure
        ! functions of (Pr, Pr_t) and always cheap, so they are computed
        ! unconditionally; the y+ FIELD is full size only under wall
        ! functions (a 1-cell dummy otherwise, the nutNone idiom).
        if (allocated(sc%wfP)) deallocate(sc%wfP, sc%wfYpt)
        allocate(sc%wfP(int(sc%n)), sc%wfYpt(int(sc%n)))
        do is = 1, int(sc%n)
            sc%wfP(is) = jayatilleke_p(sc%pr(is)/sc%prt(is))
            sc%wfYpt(is) = thermal_yplus(sc%pr(is), sc%prt(is), sc%wfP(is))
        end do
        if (allocated(sc%wfYplus)) deallocate(sc%wfYplus)
        if (wall_function) then
            allocate(sc%wfYplus(0:nx+1,0:ny+1,0:nz+1,blk%nBlocks))
        else
            allocate(sc%wfYplus(0:0,0:0,0:0,1))
        end if
        sc%wfYplus = 0.0d0

        if (allocated(sc%varList)) deallocate(sc%varList)
        allocate(sc%varList(int(sc%n)))
        do is = 1, int(sc%n)
            sc%varList(is) = VAR_S0 + int(is, C_INT)
        end do

        if (has_terminal) then
            print *, "passive scalars:", sc%n
            do is = 1, int(sc%n)
                print '(a,i0,a,a,a,f8.4,a,es10.3,a)', "   s", is, " '", trim(sc%name(is)), &
                    "'  pr =", sc%pr(is), "  source =", sc%source(is), &
                    merge(" * u_" // SRC_DIR_NAME(sc%srcDir(is)) // " (Kasagi)", &
                          "               ", sc%srcType(is) == SC_SRC_VELOCITY)
                if (wall_function) print '(a,f8.4,a,f8.3)', &
                    "      thermal wall function: P =", sc%wfP(is), "  y+_T =", sc%wfYpt(is)
            end do
        end if
    end subroutine init_scalar

    ! Set-if-unset with the A0 contradiction rule: on a strict face an
    ! explicit _type key that disagrees with the declared patch is a hard
    ! config error; on a permissive one (a wall) it simply wins.
    subroutine resolve_scalar_row(sc, is, face_id, want, value, terminal, strict)
        type(scalar_type), intent(inout) :: sc
        integer, intent(in) :: is, face_id
        integer(C_INT), intent(in) :: want
        real(C_DOUBLE), intent(in) :: value
        logical, intent(in) :: terminal
        logical, intent(in) :: strict

        if (sc%bcTypeSet(is,face_id)) then
            if (strict .and. sc%bcType(is,face_id) /= want) then
                if (terminal) print '(a,i0,a,i0)', &
                    " error: explicit [scalar.N] _type key contradicts the declared " // &
                    "patch type on face ", face_id, ", scalar ", is
                error stop "[scalar] BC type contradicts the declared patch type"
            end if
        else
            sc%bcType(is,face_id) = want
        end if
        if (.not. sc%bcValueSet(is,face_id)) sc%bcValue(is,face_id) = value
    end subroutine resolve_scalar_row

    real(C_DOUBLE) function inv_delta(d) result(inv)
        real(C_DOUBLE), intent(in) :: d

        inv = 1.0d0/sign(max(abs(d), 1.0d-300), d)
    end function inv_delta

    ! Initial condition in the q scalar slots (host side, before the device
    ! map). Always called, including on restart: read_field then overwrites
    ! the scalars whose dataset the file carries, so an absent dataset (a
    ! newly added scalar, an older restart) keeps these values -- the RANS
    ! named-scalar restart precedent.
    subroutine init_scalar_fields(sc, dns, blk)
        type(scalar_type), intent(in) :: sc
        type(dns_type), intent(in) :: dns
        type(block_set_type), intent(inout) :: blk

        integer :: is, i, j, k, b, nx, ny, nz
        real(C_DOUBLE) :: s0, s1, y

        if (.not. scalars_enabled(sc)) return

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))

        do is = 1, int(sc%n)
            blk%q(:,:,:,VAR_S0+is,:) = sc%initValue(is)
            if (sc%initProfile(is) /= SC_INIT_LINEAR_Y) cycle
            ! linear_y: the straight line between the y faces' Dirichlet
            ! values (the conduction/heated-channel convenience IC).
            s0 = sc%bcValue(is, boundary_face_id(2,0))
            s1 = sc%bcValue(is, boundary_face_id(2,1))
            do b = 1, int(blk%nBlocks)
                do k = 0, nz+1
                    do j = 0, ny+1
                        y = blk%y(j,VAR_P,b)
                        do i = 0, nx+1
                            blk%q(i,j,k,VAR_S0+is,b) = s0 + (s1 - s0)*y/dns%leng(2)
                        end do
                    end do
                end do
            end do
        end do
    end subroutine init_scalar_fields

    subroutine enter_scalar_data(sc)
        type(scalar_type), intent(inout) :: sc

        if (.not. scalars_enabled(sc)) return

#ifdef USE_OPENMP_OFFLOAD
        !$omp target enter data map(to: sc)
        !$omp target enter data map(to: sc%pr, sc%prt, sc%prtModel, sc%source, &
        !$omp& sc%srcType, sc%srcDir, sc%plainList, sc%conjList, &
        !$omp& sc%ibmValue, sc%ibmMode, &
        !$omp& sc%invDx, sc%invDy, sc%invDz, sc%nutNone, &
        !$omp& sc%cdx, sc%cdy, sc%cdz, &
        !$omp& sc%wfP, sc%wfYpt, sc%wfYplus, sc%phi, sc%vfrac, sc%solidK, sc%solidC, &
        !$omp& sc%solidSource, sc%contactR, sc%tangCorr, &
        !$omp& sc%solidDepth, sc%outerK, sc%outerC)
#endif
    end subroutine enter_scalar_data

    subroutine exit_scalar_data(sc)
        type(scalar_type), intent(inout) :: sc

        if (.not. scalars_enabled(sc)) return

#ifdef USE_OPENMP_OFFLOAD
        !$omp target exit data map(delete: sc%pr, sc%prt, sc%prtModel, sc%source, &
        !$omp& sc%srcType, sc%srcDir, sc%plainList, sc%conjList, &
        !$omp& sc%ibmValue, sc%ibmMode, &
        !$omp& sc%invDx, sc%invDy, sc%invDz, sc%nutNone, &
        !$omp& sc%cdx, sc%cdy, sc%cdz, &
        !$omp& sc%wfP, sc%wfYpt, sc%wfYplus, sc%phi, sc%vfrac, sc%solidK, sc%solidC, &
        !$omp& sc%solidSource, sc%contactR, sc%tangCorr, &
        !$omp& sc%solidDepth, sc%outerK, sc%outerC)
        !$omp target exit data map(delete: sc)
#endif
    end subroutine exit_scalar_data

    subroutine destroy_scalar(sc)
        type(scalar_type), intent(inout) :: sc

        call deallocate_config(sc)
        if (allocated(sc%varList)) deallocate(sc%varList)
        if (allocated(sc%nutNone)) deallocate(sc%nutNone)
        if (allocated(sc%wfP)) deallocate(sc%wfP)
        if (allocated(sc%wfYpt)) deallocate(sc%wfYpt)
        if (allocated(sc%wfYplus)) deallocate(sc%wfYplus)
        if (allocated(sc%invDx)) deallocate(sc%invDx)
        if (allocated(sc%invDy)) deallocate(sc%invDy)
        if (allocated(sc%invDz)) deallocate(sc%invDz)
        if (allocated(sc%cdx)) deallocate(sc%cdx)
        if (allocated(sc%cdy)) deallocate(sc%cdy)
        if (allocated(sc%cdz)) deallocate(sc%cdz)
        if (allocated(sc%phi)) deallocate(sc%phi)
        if (allocated(sc%vfrac)) deallocate(sc%vfrac)
        sc%n = 0_C_INT
        sc%nConjugate = 0_C_INT
        sc%nPlain = 0_C_INT
    end subroutine destroy_scalar

    ! One RK3 substage of every scalar: TWO fused kernels, one per
    ! immersed-body family (collapse(4) over b,k,j,i with an inner loop over
    ! that family's scalar LIST, so the GPU sees one launch per family that
    ! is present, whatever the scalar count). scalar_transport_plain serves
    ! the dirichlet and adiabatic modes -- and every scalar of a body-free
    ! run -- scalar_transport_conjugate the conjugate mode. They were ONE
    ! kernel until 2026-10-20 (docs/next_session_code_simplification.md item
    ! 1): the per-scalar mode flags were re-evaluated per cell and the ~100
    ! private variables of the dormant branches cost occupancy on both
    ! families. Every per-cell expression moved VERBATIM, so each kernel is
    ! bit-exact with the fused one on its own scalars (gated at nofma: the
    ! 9-case suite, conjugate C1-C3 and validation/scalar/mixed{,les}.ini,
    ! which carry one scalar of every mode). Both read the START-of-substage
    ! q and write only their own qs/oldrhs scalar slots, so their order is
    ! immaterial.
    !
    ! Convection: 2nd-order central in DIVERGENCE form on the p cell's own six
    ! face velocities -- no velocity interpolation is needed on the staggered
    ! mesh, and the flux form telescopes, so sum(s dV) changes only by the
    ! boundary flux.
    !
    ! [scalar] convection subtracts f*s*(div u)|stencil built from the SAME
    ! face velocities, f = 0 / 1/2 / 1 for divergence / skew / advective
    ! (docs/history/next_session_scalar.md Section 2). NOTE what each buys; the
    ! scalar and the momentum kernel used to share one key, which was
    ! misleading, because f differs between them: the FULL
    ! subtraction is the ADVECTIVE form u.grad s, which preserves a UNIFORM
    ! scalar exactly for ANY advecting field (div u never enters), at the cost
    ! of the exact global conservation the divergence form has. Subtracting
    ! HALF would instead be the skew-symmetric form (the momentum kernel's
    ! choice, energy-neutral in sum s^2 for any advecting field) but leaves
    ! 1/2 s div u on a uniform field. The projection makes div u small either
    ! way; the trade is uniform-preservation vs. sum s^2 neutrality.
    !
    ! There is deliberately NO upwind option (module header of rans.f90 and
    ! docs/history/next_session_scalar.md Section 2): central/skew only until the
    ! shared TVD increment lands.
    !
    ! Diffusion: face-flux form on the module's own inverse centre-to-centre
    ! distances, with the face diffusivity
    !
    !   D_face = 1/(Re Pr) + 1/2 (nu_t,L + nu_t,R)/Pr_t(face)
    !
    ! reading turb%nut -- the ONE blended eddy viscosity that LES (WALE /
    ! Smagorinsky), RANS (SST) and IDDES all write, so a single code path
    ! serves all three models (increment S2). Pr_t is the per-scalar constant
    ! or the Kays-Crawford correlation, selected by [scalar.N] prt_model.
    !
    ! Immersed body (increment S3), per scalar via [scalar.N] ibm_wall:
    !   * dirichlet -- volume penalization toward [scalar.N] ibm_value with
    !     ibm.f90's integrating factors at x = dt_gamma coef_p/Pr, formed
    !     inline (unconditionally stable, so no dt restriction; no mu array
    !     per scalar). The stored coefficient carries the 1/Re scaling, so
    !     coef_p/Pr is exactly the scalar's own 1/(Re Pr) -- ONE cell-centred
    !     coefficient array serves every scalar, each with its own Pr. Inside
    !     the body coef_p = SOLID/Re, so the factors -> 0 and the cell holds
    !     ibm_value. Diffusive fluxes are NOT masked: the solid cell holds
    !     the wall value and delivers the flux. First-order/staircase, like
    !     the velocity penalization.
    !   * adiabatic -- no penalization; instead BOTH the convective and the
    !     diffusive flux are masked on each of the six p-cell faces whose
    !     staggered velocity coefficient is solid (the solw/sole/... test
    !     rans.f90 uses for k). The mask is symmetric across a face, so the
    !     flux form still telescopes and the fluid conserves sum(s dV)
    !     exactly.
    !   * conjugate (increment C1, docs/history/next_session_conjugate.md) -- the
    !     solid is a REAL unknown. No penalization either; instead
    !       - the face diffusivity becomes conjugate_face_diffusivity, which
    !         is kappa*dm within one material and the distance-weighted
    !         harmonic mean across a cut face. Every UNCUT FLUID face keeps
    !         today's expression verbatim, which is why the feature is
    !         untouched away from the body;
    !       - CONVECTION is hard-masked on every face whose staggered
    !         velocity node is solid AND on every cut face (u.n = 0 on the
    !         interface exactly, so this is consistent to second order).
    !         Relying on the penalised velocity being "small" is not
    !         acceptable here: the solid now carries a genuine temperature
    !         field that must not be advected. The skew term's divergence is
    !         then built from the SAME masked face velocities, so a uniform
    !         scalar is still preserved exactly;
    !       - the flux divergence is divided by the local capacity C
    !         (1 in the fluid, C_s = (rho c)_s/(rho c)_f in the solid), and
    !         the solid takes its own volumetric source;
    !       - nu_t enters NEITHER the solid NOR a cut face.
    !     The interface term is written as a FACE FLUX, never a cell source,
    !     which is what keeps sum(C s dV) conserved to round-off whatever the
    !     scheme's accuracy -- and what lets the C2 tangential correction and
    !     the C4 wedge model drop in without touching anything around them.
    ! With no immersed body (useIbm off) none of the branches is taken and
    ! the arithmetic is the S2 kernel's, byte for byte.
    !
    ! The nut ghost cells follow the momentum SGS convention: the halo
    ! exchange fills them across block and rank boundaries, and they stay at
    ! zero on a physical face (a resolved wall has nu_t -> 0 there anyway).
    !
    ! Under [rans] wall_treatment = wall_function (increment S5a) the face
    ! diffusivity at WALL CELLS comes from the thermal wall function instead
    ! (Kader/Jayatilleke, above) -- the momentum wall function's wall-cell
    ! nu_t is deliberately NOT reused there, because it carries the momentum
    ! sublayer's resistance, not the thermal one. Every other cell keeps
    ! nu_t/Pr_t, and with wall functions off the branch is not entered.
    subroutine scalar_transport(sc, blk, dns, turb, coef, dt_alpha, dt_beta, dt_gamma)
        type(scalar_type), intent(in) :: sc
        type(block_set_type), intent(inout) :: blk
        type(dns_type), intent(in) :: dns
        type(turb_type), intent(in) :: turb
        ! ibm%coef: with scalars on it always carries the cell-centred VAR_P
        ! column (init_ibm), so there is no dummy-array case here -- the
        ! body is switched off by dns%ibm_enabled alone.
        real(C_DOUBLE), intent(in) :: coef(0:,0:,0:,1:,1:)
        real(C_DOUBLE), intent(in) :: dt_alpha, dt_beta, dt_gamma

        if (.not. scalars_enabled(sc)) return

        ! With no turbulence model turb%nut does not exist; hand the kernels
        ! the 1-cell dummy and switch the eddy term off, which makes the
        ! arithmetic byte-identical to the molecular-only (S1) kernel.
        if (turbulence_is_enabled(turb) .and. allocated(turb%nut)) then
            if (sc%nPlain > 0_C_INT) call scalar_transport_plain(sc, blk, dns, &
                dt_alpha, dt_beta, dt_gamma, turb%nut, .true., coef)
            if (sc%nConjugate > 0_C_INT) call scalar_transport_conjugate(sc, blk, dns, &
                dt_alpha, dt_beta, dt_gamma, turb%nut, .true., coef)
        else
            if (sc%nPlain > 0_C_INT) call scalar_transport_plain(sc, blk, dns, &
                dt_alpha, dt_beta, dt_gamma, sc%nutNone, .false., coef)
            if (sc%nConjugate > 0_C_INT) call scalar_transport_conjugate(sc, blk, dns, &
                dt_alpha, dt_beta, dt_gamma, sc%nutNone, .false., coef)
        end if
    end subroutine scalar_transport

    ! The dirichlet / adiabatic family: every scalar of sc%plainList. With no
    ! immersed body (useIbm off) neither branch is taken and the arithmetic
    ! is the S2 kernel's, byte for byte.
    subroutine scalar_transport_plain(sc, blk, dns, dt_alpha, dt_beta, dt_gamma, &
            nut, useNut, coef)
        type(scalar_type), intent(in) :: sc
        type(block_set_type), intent(inout) :: blk
        type(dns_type), intent(in) :: dns
        real(C_DOUBLE), intent(in) :: dt_alpha, dt_beta, dt_gamma
        real(C_DOUBLE), intent(in) :: nut(0:,0:,0:,1:)
        logical, intent(in) :: useNut
        real(C_DOUBLE), intent(in) :: coef(0:,0:,0:,1:,1:)

        integer :: i, j, k, b, il, is, nx, ny, nz, nBlocks, nList, var, scr
        real(C_DOUBLE) :: ire, re, uw, ue, vs, vn, wb, wt, divu
        real(C_DOUBLE) :: s0, conv, diff, rhs, srcVal, dm, fw, fe, ss, xpen
        ! The 7-point nu_t and wall-function y+ stencils, loaded once per
        ! cell (guarded: each array is a 1-cell dummy when its model is
        ! off) and handed to plain_face_diffusivity pairwise.
        real(C_DOUBLE) :: nutc, nutw, nute, nuts, nutn, nutb, nutt
        real(C_DOUBLE) :: ypc, ypw, ype, yps, ypn, ypb, ypt
        real(C_DOUBLE) :: dxw, dxe, dys, dyn, dzb, dzt
        real(C_DOUBLE) :: convFac
        logical :: nondiv, useIbm, adiab, wallfn
        logical :: clw, cle, cls, cln, clb, clt
        logical :: solw, sole, sols, soln, solb, solt
        logical :: mw, me, ms, mn, mb, mt

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))
        nBlocks = int(blk%nBlocks)
        nList = int(sc%nPlain)
        re = dns%re
        ire = 1.0d0/re
        ! f = 0 (divergence) skips the subtraction entirely, so the default
        ! path keeps its arithmetic byte for byte -- the branch is not entered
        ! rather than multiplying by a zero factor.
        nondiv = sc%convMode /= SC_CONV_DIV
        convFac = 0.0d0
        if (sc%convMode == SC_CONV_SKEW) convFac = 0.5d0
        if (sc%convMode == SC_CONV_ADV)  convFac = 1.0d0
        useIbm = logical(dns%ibm_enabled)
        ! S5a: the thermal wall function replaces the face eddy diffusivity
        ! AT WALL CELLS. Off (every resolved-wall and every non-RANS run) the
        ! branch is not entered and the S2/S3 arithmetic is reproduced byte
        ! for byte -- sc%wfYplus is then a 1-cell dummy, mapped but unread.
        wallfn = dns%rans_wall_treatment == 1_C_INT

        !$omp target teams distribute parallel do collapse(4) &
        !$omp& map(to: ire, re, dt_alpha, dt_beta, dt_gamma, nondiv, convFac, useNut, useIbm, &
        !$omp& nList, nx, ny, nz, sc%plainList, &
        !$omp& blk%q, blk%d1x, blk%d1y, blk%d1z, blk%physLow, blk%physHigh, nut, coef, &
        !$omp& sc%pr, sc%prt, sc%prtModel, sc%source, sc%srcType, sc%srcDir, &
        !$omp& sc%invDx, sc%invDy, sc%invDz, &
        !$omp& sc%ibmMode, sc%ibmValue, sc%wfP, sc%wfYpt, sc%wfYplus, wallfn) &
        !$omp& map(tofrom: blk%qs, blk%oldrhs) &
        !$omp& private(i,j,k,b,il,is,var,scr,uw,ue,vs,vn,wb,wt,divu, &
        !$omp& s0,conv,diff,rhs,srcVal,dm,fw,fe,ss,xpen,adiab, &
        !$omp& clw,cle,cls,cln,clb,clt,solw,sole,sols,soln,solb,solt, &
        !$omp& mw,me,ms,mn,mb,mt, &
        !$omp& nutc,nutw,nute,nuts,nutn,nutb,nutt,ypc,ypw,ype,yps,ypn,ypb,ypt, &
        !$omp& dxw,dxe,dys,dyn,dzb,dzt)
        do b = 1, nBlocks
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    ! The p cell's six face velocities: exactly the six
                    ! staggered components, no interpolation.
                    uw = blk%q(i,  j,  k,  VAR_U,b)
                    ue = blk%q(i+1,j,  k,  VAR_U,b)
                    vs = blk%q(i,  j,  k,  VAR_V,b)
                    vn = blk%q(i,  j+1,k,  VAR_V,b)
                    wb = blk%q(i,  j,  k,  VAR_W,b)
                    wt = blk%q(i,  j,  k+1,VAR_W,b)
                    divu = (ue - uw)*blk%d1x(i,VAR_P,b) &
                         + (vn - vs)*blk%d1y(j,VAR_P,b) &
                         + (wt - wb)*blk%d1z(k,VAR_P,b)

                    ! A FACE_CLOSED face borders a block removed inside the
                    ! immersed body: its halo is held at zero, so the
                    ! DIFFUSIVE flux must be masked (the convective one
                    ! already vanishes -- the face velocity is pinned to 0).
                    ! FACE_PHYS faces are NOT masked: their ghosts carry the
                    ! boundary condition.
                    clw = i == 1  .and. blk%physLow(1,b)  == FACE_CLOSED
                    cle = i == nx .and. blk%physHigh(1,b) == FACE_CLOSED
                    cls = j == 1  .and. blk%physLow(2,b)  == FACE_CLOSED
                    cln = j == ny .and. blk%physHigh(2,b) == FACE_CLOSED
                    clb = k == 1  .and. blk%physLow(3,b)  == FACE_CLOSED
                    clt = k == nz .and. blk%physHigh(3,b) == FACE_CLOSED

                    ! Solid staggered faces of this p cell (adiabatic body
                    ! mode only; scalar-independent, so hoisted out of the
                    ! scalar loop). rans.f90's test for k, verbatim.
                    solw = .false.; sole = .false.; sols = .false.
                    soln = .false.; solb = .false.; solt = .false.
                    if (useIbm) then
                        solw = abs(coef(i,  j,  k,  VAR_U,b)) > SOLID_FACE_THRESHOLD
                        sole = abs(coef(i+1,j,  k,  VAR_U,b)) > SOLID_FACE_THRESHOLD
                        sols = abs(coef(i,  j,  k,  VAR_V,b)) > SOLID_FACE_THRESHOLD
                        soln = abs(coef(i,  j+1,k,  VAR_V,b)) > SOLID_FACE_THRESHOLD
                        solb = abs(coef(i,  j,  k,  VAR_W,b)) > SOLID_FACE_THRESHOLD
                        solt = abs(coef(i,  j,  k+1,VAR_W,b)) > SOLID_FACE_THRESHOLD
                    end if

                    ! The nu_t stencil (and under wall functions the y+ one):
                    ! scalar-independent, loaded once per cell and guarded,
                    ! since each array is a 1-cell dummy when its model is
                    ! off. plain_face_diffusivity forms the 1/2 (L + R) face
                    ! average the momentum SGS correction uses, from the pairs
                    ! (Pr_t may be per-scalar, so the face value cannot be).
                    nutc = 0.0d0; nutw = 0.0d0; nute = 0.0d0; nuts = 0.0d0
                    nutn = 0.0d0; nutb = 0.0d0; nutt = 0.0d0
                    ypc = 0.0d0; ypw = 0.0d0; ype = 0.0d0; yps = 0.0d0
                    ypn = 0.0d0; ypb = 0.0d0; ypt = 0.0d0
                    if (useNut) then
                        nutc = nut(i,j,k,b)
                        nutw = nut(i-1,j,k,b); nute = nut(i+1,j,k,b)
                        nuts = nut(i,j-1,k,b); nutn = nut(i,j+1,k,b)
                        nutb = nut(i,j,k-1,b); nutt = nut(i,j,k+1,b)
                        if (wallfn) then
                            ypc = sc%wfYplus(i,j,k,b)
                            ypw = sc%wfYplus(i-1,j,k,b); ype = sc%wfYplus(i+1,j,k,b)
                            yps = sc%wfYplus(i,j-1,k,b); ypn = sc%wfYplus(i,j+1,k,b)
                            ypb = sc%wfYplus(i,j,k-1,b); ypt = sc%wfYplus(i,j,k+1,b)
                        end if
                    end if

                    do il = 1, nList
                        is = int(sc%plainList(il))
                        var = VAR_S0 + is
                        scr = SCR_S0 + is
                        s0 = blk%q(i,j,k,var,b)

                        ! An adiabatic body seals this cell's solid faces:
                        ! convective AND diffusive flux masked. Every other
                        ! configuration leaves adiab .false., so the six
                        ! masks below collapse to the FACE_CLOSED flags and
                        ! the S2 expressions are reproduced exactly.
                        adiab = useIbm .and. sc%ibmMode(is) == SC_IBM_ADIABATIC
                        mw = clw .or. (adiab .and. solw)
                        me = cle .or. (adiab .and. sole)
                        ms = cls .or. (adiab .and. sols)
                        mn = cln .or. (adiab .and. soln)
                        mb = clb .or. (adiab .and. solb)
                        mt = clt .or. (adiab .and. solt)

                        if (adiab) then
                            conv = (merge(0.0d0, ue*0.5d0*(s0 + blk%q(i+1,j,k,var,b)), me) &
                                  - merge(0.0d0, uw*0.5d0*(blk%q(i-1,j,k,var,b) + s0), mw)) &
                                    *blk%d1x(i,VAR_P,b) &
                                 + (merge(0.0d0, vn*0.5d0*(s0 + blk%q(i,j+1,k,var,b)), mn) &
                                  - merge(0.0d0, vs*0.5d0*(blk%q(i,j-1,k,var,b) + s0), ms)) &
                                    *blk%d1y(j,VAR_P,b) &
                                 + (merge(0.0d0, wt*0.5d0*(s0 + blk%q(i,j,k+1,var,b)), mt) &
                                  - merge(0.0d0, wb*0.5d0*(blk%q(i,j,k-1,var,b) + s0), mb)) &
                                    *blk%d1z(k,VAR_P,b)
                        else
                            conv = (ue*0.5d0*(s0 + blk%q(i+1,j,k,var,b)) &
                                  - uw*0.5d0*(blk%q(i-1,j,k,var,b) + s0))*blk%d1x(i,VAR_P,b) &
                                 + (vn*0.5d0*(s0 + blk%q(i,j+1,k,var,b)) &
                                  - vs*0.5d0*(blk%q(i,j-1,k,var,b) + s0))*blk%d1y(j,VAR_P,b) &
                                 + (wt*0.5d0*(s0 + blk%q(i,j,k+1,var,b)) &
                                  - wb*0.5d0*(blk%q(i,j,k-1,var,b) + s0))*blk%d1z(k,VAR_P,b)
                        end if
                        ! The skew/advective term subtracts the UNMASKED
                        ! divergence here, the adiabatic mode included (its
                        ! shipped, gated arithmetic; the conjugate kernel is
                        ! the one that builds it from the masked velocities).
                        if (nondiv) conv = conv - convFac*s0*divu

                        ! Molecular + eddy diffusivity on each face, ONE
                        ! definition (plain_face_diffusivity, shared with the
                        ! statistics kernels). With useNut off every face
                        ! keeps dm exactly, so the S1 arithmetic is reproduced
                        ! bit-for-bit.
                        dm = ire/sc%pr(is)
                        dxw = plain_face_diffusivity(dm, useNut, wallfn, nutw, nutc, ypw, ypc, &
                            sc%pr(is), sc%prt(is), sc%prtModel(is), re, sc%wfP(is), sc%wfYpt(is))
                        dxe = plain_face_diffusivity(dm, useNut, wallfn, nutc, nute, ypc, ype, &
                            sc%pr(is), sc%prt(is), sc%prtModel(is), re, sc%wfP(is), sc%wfYpt(is))
                        dys = plain_face_diffusivity(dm, useNut, wallfn, nuts, nutc, yps, ypc, &
                            sc%pr(is), sc%prt(is), sc%prtModel(is), re, sc%wfP(is), sc%wfYpt(is))
                        dyn = plain_face_diffusivity(dm, useNut, wallfn, nutc, nutn, ypc, ypn, &
                            sc%pr(is), sc%prt(is), sc%prtModel(is), re, sc%wfP(is), sc%wfYpt(is))
                        dzb = plain_face_diffusivity(dm, useNut, wallfn, nutb, nutc, ypb, ypc, &
                            sc%pr(is), sc%prt(is), sc%prtModel(is), re, sc%wfP(is), sc%wfYpt(is))
                        dzt = plain_face_diffusivity(dm, useNut, wallfn, nutc, nutt, ypc, ypt, &
                            sc%pr(is), sc%prt(is), sc%prtModel(is), re, sc%wfP(is), sc%wfYpt(is))

                        fw = merge(0.0d0, dxw*(s0 - blk%q(i-1,j,k,var,b)) &
                            *sc%invDx(i,b), mw)
                        fe = merge(0.0d0, dxe*(blk%q(i+1,j,k,var,b) - s0) &
                            *sc%invDx(i+1,b), me)
                        diff = (fe - fw)*blk%d1x(i,VAR_P,b)
                        fw = merge(0.0d0, dys*(s0 - blk%q(i,j-1,k,var,b)) &
                            *sc%invDy(j,b), ms)
                        fe = merge(0.0d0, dyn*(blk%q(i,j+1,k,var,b) - s0) &
                            *sc%invDy(j+1,b), mn)
                        diff = diff + (fe - fw)*blk%d1y(j,VAR_P,b)
                        fw = merge(0.0d0, dzb*(s0 - blk%q(i,j,k-1,var,b)) &
                            *sc%invDz(k,b), mb)
                        fe = merge(0.0d0, dzt*(blk%q(i,j,k+1,var,b) - s0) &
                            *sc%invDz(k+1,b), mt)
                        diff = diff + (fe - fw)*blk%d1z(k,VAR_P,b)

                        ! The volumetric source. UNIFORM is the S0 constant
                        ! and reproduces the old arithmetic operand for
                        ! operand, so that path stays bit-exact. VELOCITY is
                        ! Kasagi's fully-developed channel source, beta*u_x,
                        ! built from the p cell's OWN two staggered faces --
                        ! which the convection term above already loaded, so
                        ! it costs one add and one multiply. u is the
                        ! INSTANTANEOUS velocity, which is the whole point:
                        ! the term comes from u.grad(beta x), not from a mean.
                        ! source_dir picks the pair; the x branch is the
                        ! default and is the original expression operand for
                        ! operand, so an ini without the key is bit-exact.
                        srcVal = sc%source(is)
                        if (sc%srcType(is) == SC_SRC_VELOCITY) then
                            if (sc%srcDir(is) == 2_C_INT) then
                                srcVal = sc%source(is)*0.5d0*(vs + vn)
                            else if (sc%srcDir(is) == 3_C_INT) then
                                srcVal = sc%source(is)*0.5d0*(wb + wt)
                            else
                                srcVal = sc%source(is)*0.5d0*(uw + ue)
                            end if
                        end if

                        rhs = -conv + diff + srcVal
                        ss = s0 + dt_alpha*rhs + dt_beta*blk%oldrhs(i,j,k,scr,b)
                        ! Volume penalization toward the body value (dirichlet
                        ! mode), with the momentum predictor's integrating
                        ! factors (ibm.f90 penal_*): the substage solves
                        ! ds/dt = R - lambda (s - s_body), lambda = coef_p/Pr,
                        ! so the start-of-substage value s0 takes the STATE
                        ! factor, the increment the INCREMENT factor, and the
                        ! body value the remainder 1 - state (EXACTLY 1 in a
                        ! solid cell, which therefore holds s_body). A fluid
                        ! cell (coef_p = 0, factors exactly 1) is skipped: the
                        ! divides are paid in the body's cells only.
                        ! oldrhs keeps the UNpenalized rhs, exactly as the
                        ! momentum predictor does.
                        if (useIbm .and. sc%ibmMode(is) == SC_IBM_DIRICHLET) then
                            xpen = dt_gamma*coef(i,j,k,VAR_P,b)
                            if (xpen /= 0.0d0) then
                                xpen = xpen/sc%pr(is)
                                ss = ss*penal_incr_factor(xpen) + penal_state_minus_incr(xpen)*s0 &
                                     + (1.0d0 - penal_state_factor(xpen))*sc%ibmValue(is)
                            end if
                        end if
                        blk%qs(i,j,k,scr,b) = ss
                        blk%oldrhs(i,j,k,scr,b) = rhs
                    end do
                end do
            end do
        end do
        end do
        !$omp end target teams distribute parallel do
    end subroutine scalar_transport_plain

    ! The conjugate family: every scalar of sc%conjList. Three things the
    ! config guarantees are built in rather than tested per cell
    ! (validate_conjugate_config): there IS an immersed body, [scalar]
    ! convection is advective or skew (the divergence form is refused, so the
    ! div(u) subtraction always runs), and wall functions are off. A
    ! conjugate scalar is never adiabatic, so the diffusive masks are the
    ! FACE_CLOSED flags alone and the convective masks add the solid and the
    ! cut faces.
    subroutine scalar_transport_conjugate(sc, blk, dns, dt_alpha, dt_beta, dt_gamma, &
            nut, useNut, coef)
        type(scalar_type), intent(in) :: sc
        type(block_set_type), intent(inout) :: blk
        type(dns_type), intent(in) :: dns
        real(C_DOUBLE), intent(in) :: dt_alpha, dt_beta, dt_gamma
        real(C_DOUBLE), intent(in) :: nut(0:,0:,0:,1:)
        logical, intent(in) :: useNut
        real(C_DOUBLE), intent(in) :: coef(0:,0:,0:,1:,1:)

        integer :: i, j, k, b, il, is, nx, ny, nz, nBlocks, nList, var, scr
        real(C_DOUBLE) :: ire, re, uw, ue, vs, vn, wb, wt, divu, divuse
        real(C_DOUBLE) :: s0, conv, diff, rhs, srcVal, dm, fw, fe, ss
        real(C_DOUBLE) :: dxw, dxe, dys, dyn, dzb, dzt
        ! Conjugate interface (C1): the six neighbour signed distances and
        ! this cell's, plus the per-scalar solid properties.
        real(C_DOUBLE) :: phc, phw, phe, phs, phn, phb, pht, ks, rc
        ! The solid's second material band (F2): its depth and outer kappa,
        ! and this cell's band capacity and band source. Without a band the
        ! last two are the scalar's solid_rhocp / solid_source, so the rhs
        ! expression below is unchanged operand for operand.
        real(C_DOUBLE) :: dsh, ko, csb, ssb
        ! Tangential correction (C2): the face-centred gradients of T and phi
        ! in (face-normal, tangential 1, tangential 2) form, and the six face
        ! corrections themselves. Zero unless the face is CUT and the scalar
        ! asked for the correction, so the C1 flux line is untouched
        ! elsewhere.
        real(C_DOUBLE) :: gtd, gt1, gt2, gpd, gp1, gp2
        real(C_DOUBLE) :: crw, cre, crs, crn, crb, crt
        real(C_DOUBLE) :: ntw, nte, nts, ntn, ntb, ntt
        real(C_DOUBLE) :: convFac
        logical :: solc, tang, banded, nutOn
        logical :: cutw, cute, cuts, cutn, cutb, cutt
        logical :: clw, cle, cls, cln, clb, clt
        logical :: solw, sole, sols, soln, solb, solt
        logical :: cmw, cme, cms, cmn, cmb, cmt

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))
        nBlocks = int(blk%nBlocks)
        nList = int(sc%nConjugate)
        re = dns%re
        ire = 1.0d0/re
        ! The divergence form (f = 0) is refused for a conjugate scalar, so
        ! the subtraction below is unconditional.
        convFac = 0.0d0
        if (sc%convMode == SC_CONV_SKEW) convFac = 0.5d0
        if (sc%convMode == SC_CONV_ADV)  convFac = 1.0d0

        !$omp target teams distribute parallel do collapse(4) &
        !$omp& map(to: ire, re, dt_alpha, dt_beta, dt_gamma, convFac, useNut, &
        !$omp& nList, nx, ny, nz, sc%conjList, &
        !$omp& sc%cdx, sc%cdy, sc%cdz, sc%tangCorr, &
        !$omp& blk%q, blk%d1x, blk%d1y, blk%d1z, blk%physLow, blk%physHigh, nut, coef, &
        !$omp& sc%pr, sc%prt, sc%prtModel, sc%source, sc%srcType, sc%srcDir, &
        !$omp& sc%invDx, sc%invDy, sc%invDz, &
        !$omp& sc%phi, sc%vfrac, sc%solidK, sc%solidC, sc%solidSource, sc%contactR, &
        !$omp& sc%solidDepth, sc%outerK, sc%outerC) &
        !$omp& map(tofrom: blk%qs, blk%oldrhs) &
        !$omp& private(i,j,k,b,il,is,var,scr,uw,ue,vs,vn,wb,wt,divu,divuse, &
        !$omp& s0,conv,diff,rhs,srcVal,dm,fw,fe,ss,solc,ks,rc,tang, &
        !$omp& banded,dsh,ko,csb,ssb,nutOn, &
        !$omp& gtd,gt1,gt2,gpd,gp1,gp2,crw,cre,crs,crn,crb,crt, &
        !$omp& phc,phw,phe,phs,phn,phb,pht, &
        !$omp& cutw,cute,cuts,cutn,cutb,cutt, &
        !$omp& clw,cle,cls,cln,clb,clt,solw,sole,sols,soln,solb,solt, &
        !$omp& cmw,cme,cms,cmn,cmb,cmt, &
        !$omp& ntw,nte,nts,ntn,ntb,ntt,dxw,dxe,dys,dyn,dzb,dzt)
        do b = 1, nBlocks
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    ! The p cell's six face velocities: exactly the six
                    ! staggered components, no interpolation.
                    uw = blk%q(i,  j,  k,  VAR_U,b)
                    ue = blk%q(i+1,j,  k,  VAR_U,b)
                    vs = blk%q(i,  j,  k,  VAR_V,b)
                    vn = blk%q(i,  j+1,k,  VAR_V,b)
                    wb = blk%q(i,  j,  k,  VAR_W,b)
                    wt = blk%q(i,  j,  k+1,VAR_W,b)
                    divu = (ue - uw)*blk%d1x(i,VAR_P,b) &
                         + (vn - vs)*blk%d1y(j,VAR_P,b) &
                         + (wt - wb)*blk%d1z(k,VAR_P,b)

                    ! A FACE_CLOSED face borders a block removed inside the
                    ! immersed body: its halo is held at zero, so the
                    ! DIFFUSIVE flux must be masked. FACE_PHYS faces are NOT
                    ! masked: their ghosts carry the boundary condition.
                    clw = i == 1  .and. blk%physLow(1,b)  == FACE_CLOSED
                    cle = i == nx .and. blk%physHigh(1,b) == FACE_CLOSED
                    cls = j == 1  .and. blk%physLow(2,b)  == FACE_CLOSED
                    cln = j == ny .and. blk%physHigh(2,b) == FACE_CLOSED
                    clb = k == 1  .and. blk%physLow(3,b)  == FACE_CLOSED
                    clt = k == nz .and. blk%physHigh(3,b) == FACE_CLOSED

                    ! Solid staggered faces of this p cell -- nothing is
                    ! advected across them. rans.f90's test for k, verbatim.
                    solw = abs(coef(i,  j,  k,  VAR_U,b)) > SOLID_FACE_THRESHOLD
                    sole = abs(coef(i+1,j,  k,  VAR_U,b)) > SOLID_FACE_THRESHOLD
                    sols = abs(coef(i,  j,  k,  VAR_V,b)) > SOLID_FACE_THRESHOLD
                    soln = abs(coef(i,  j+1,k,  VAR_V,b)) > SOLID_FACE_THRESHOLD
                    solb = abs(coef(i,  j,  k,  VAR_W,b)) > SOLID_FACE_THRESHOLD
                    solt = abs(coef(i,  j,  k+1,VAR_W,b)) > SOLID_FACE_THRESHOLD

                    ! Conjugate interface geometry (C1): the six neighbour
                    ! signed distances and the cut flags. Pure geometry --
                    ! scalar-independent, so hoisted out of the scalar loop
                    ! (the level-set weight itself needs no normal at all,
                    ! the obliquity lemma; see conjugate_face_diffusivity).
                    phc = sc%phi(i,  j,  k,  b)
                    phw = sc%phi(i-1,j,  k,  b)
                    phe = sc%phi(i+1,j,  k,  b)
                    phs = sc%phi(i,  j-1,k,  b)
                    phn = sc%phi(i,  j+1,k,  b)
                    phb = sc%phi(i,  j,  k-1,b)
                    pht = sc%phi(i,  j,  k+1,b)
                    solc = phc < 0.0d0
                    cutw = (phw < 0.0d0) .neqv. solc
                    cute = (phe < 0.0d0) .neqv. solc
                    cuts = (phs < 0.0d0) .neqv. solc
                    cutn = (phn < 0.0d0) .neqv. solc
                    cutb = (phb < 0.0d0) .neqv. solc
                    cutt = (pht < 0.0d0) .neqv. solc

                    ! The CONVECTIVE masks: the solid conducts (diffusion is
                    ! masked at FACE_CLOSED faces only) while nothing at all
                    ! is advected across a solid or a cut face.
                    cmw = clw .or. solw .or. cutw
                    cme = cle .or. sole .or. cute
                    cms = cls .or. sols .or. cuts
                    cmn = cln .or. soln .or. cutn
                    cmb = clb .or. solb .or. cutb
                    cmt = clt .or. solt .or. cutt

                    ! Face eddy viscosities: the same 1/2 (L + R) average the
                    ! momentum SGS correction uses for its face-normal terms.
                    ! Scalar-independent, so they are hoisted out of the
                    ! scalar loop; Pr_t is not (it may be per-scalar).
                    ntw = 0.0d0; nte = 0.0d0; nts = 0.0d0
                    ntn = 0.0d0; ntb = 0.0d0; ntt = 0.0d0
                    if (useNut) then
                        ntw = 0.5d0*(nut(i-1,j,k,b) + nut(i,j,k,b))
                        nte = 0.5d0*(nut(i,j,k,b) + nut(i+1,j,k,b))
                        nts = 0.5d0*(nut(i,j-1,k,b) + nut(i,j,k,b))
                        ntn = 0.5d0*(nut(i,j,k,b) + nut(i,j+1,k,b))
                        ntb = 0.5d0*(nut(i,j,k-1,b) + nut(i,j,k,b))
                        ntt = 0.5d0*(nut(i,j,k,b) + nut(i,j,k+1,b))
                    end if

                    do il = 1, nList
                        is = int(sc%conjList(il))
                        var = VAR_S0 + is
                        scr = SCR_S0 + is
                        s0 = blk%q(i,j,k,var,b)
                        tang = sc%tangCorr(is) /= 0_C_INT
                        banded = sc%solidDepth(is) > 0.0d0

                        conv = (merge(0.0d0, ue*0.5d0*(s0 + blk%q(i+1,j,k,var,b)), cme) &
                              - merge(0.0d0, uw*0.5d0*(blk%q(i-1,j,k,var,b) + s0), cmw)) &
                                *blk%d1x(i,VAR_P,b) &
                             + (merge(0.0d0, vn*0.5d0*(s0 + blk%q(i,j+1,k,var,b)), cmn) &
                              - merge(0.0d0, vs*0.5d0*(blk%q(i,j-1,k,var,b) + s0), cms)) &
                                *blk%d1y(j,VAR_P,b) &
                             + (merge(0.0d0, wt*0.5d0*(s0 + blk%q(i,j,k+1,var,b)), cmt) &
                              - merge(0.0d0, wb*0.5d0*(blk%q(i,j,k-1,var,b) + s0), cmb)) &
                                *blk%d1z(k,VAR_P,b)
                        ! The skew/advective term must subtract the divergence
                        ! of the SAME face velocities the convective flux
                        ! used, or a sealed solid cell would be "advected" by
                        ! the penalised velocity's residual divergence. With
                        ! the masked form a uniform scalar stays uniform
                        ! exactly, in the solid as in the fluid.
                        divuse = (merge(0.0d0, ue, cme) - merge(0.0d0, uw, cmw)) &
                                    *blk%d1x(i,VAR_P,b) &
                               + (merge(0.0d0, vn, cmn) - merge(0.0d0, vs, cms)) &
                                    *blk%d1y(j,VAR_P,b) &
                               + (merge(0.0d0, wt, cmt) - merge(0.0d0, wb, cmb)) &
                                    *blk%d1z(k,VAR_P,b)
                        conv = conv - convFac*s0*divuse

                        ! THE conjugate scheme: one face coefficient, ONE
                        ! definition (conjugate_face_total, shared with the
                        ! statistics kernels): kappa*dm within one material
                        ! (kappa = 1 in the fluid, so a fluid-fluid face is
                        ! the plain kernel's line unchanged), across a cut
                        ! face the distance-weighted harmonic mean built on
                        ! the level-set fraction of the arm, the band's
                        ! coefficient (F2) for a banded scalar; and the eddy
                        ! part on a FLUID cell's UNCUT faces only.
                        dm = ire/sc%pr(is)
                        ks = sc%solidK(is)
                        rc = sc%contactR(is)
                        nutOn = useNut .and. .not. solc
                        ! The band branch (F2) is entered only by a scalar
                        ! that asked for one, so the single-material calls
                        ! stay byte for byte what they were.
                        if (banded) then
                            dsh = sc%solidDepth(is)
                            ko = sc%outerK(is)
                            dxw = band_face_total(dm, phw, phc, ks, ko, dsh, rc, sc%invDx(i,b), &
                                nutOn .and. .not. cutw, ntw, sc%pr(is), sc%prt(is), sc%prtModel(is), re)
                            dxe = band_face_total(dm, phc, phe, ks, ko, dsh, rc, sc%invDx(i+1,b), &
                                nutOn .and. .not. cute, nte, sc%pr(is), sc%prt(is), sc%prtModel(is), re)
                            dys = band_face_total(dm, phs, phc, ks, ko, dsh, rc, sc%invDy(j,b), &
                                nutOn .and. .not. cuts, nts, sc%pr(is), sc%prt(is), sc%prtModel(is), re)
                            dyn = band_face_total(dm, phc, phn, ks, ko, dsh, rc, sc%invDy(j+1,b), &
                                nutOn .and. .not. cutn, ntn, sc%pr(is), sc%prt(is), sc%prtModel(is), re)
                            dzb = band_face_total(dm, phb, phc, ks, ko, dsh, rc, sc%invDz(k,b), &
                                nutOn .and. .not. cutb, ntb, sc%pr(is), sc%prt(is), sc%prtModel(is), re)
                            dzt = band_face_total(dm, phc, pht, ks, ko, dsh, rc, sc%invDz(k+1,b), &
                                nutOn .and. .not. cutt, ntt, sc%pr(is), sc%prt(is), sc%prtModel(is), re)
                        else
                            dxw = conjugate_face_total(dm, phw, phc, ks, rc, sc%invDx(i,b), &
                                nutOn .and. .not. cutw, ntw, sc%pr(is), sc%prt(is), sc%prtModel(is), re)
                            dxe = conjugate_face_total(dm, phc, phe, ks, rc, sc%invDx(i+1,b), &
                                nutOn .and. .not. cute, nte, sc%pr(is), sc%prt(is), sc%prtModel(is), re)
                            dys = conjugate_face_total(dm, phs, phc, ks, rc, sc%invDy(j,b), &
                                nutOn .and. .not. cuts, nts, sc%pr(is), sc%prt(is), sc%prtModel(is), re)
                            dyn = conjugate_face_total(dm, phc, phn, ks, rc, sc%invDy(j+1,b), &
                                nutOn .and. .not. cutn, ntn, sc%pr(is), sc%prt(is), sc%prtModel(is), re)
                            dzb = conjugate_face_total(dm, phb, phc, ks, rc, sc%invDz(k,b), &
                                nutOn .and. .not. cutb, ntb, sc%pr(is), sc%prt(is), sc%prtModel(is), re)
                            dzt = conjugate_face_total(dm, phc, pht, ks, rc, sc%invDz(k+1,b), &
                                nutOn .and. .not. cutt, ntt, sc%pr(is), sc%prt(is), sc%prtModel(is), re)
                        end if
                        ! The cell's capacity and source: the solid's own pair,
                        ! except that the OUTER band (F2) carries its own
                        ! capacity and NO source -- it stands in for whatever
                        ! is beyond the wall, and a solid_source that fired
                        ! there would have nowhere to send its heat (an
                        ! insulated band would simply run away).
                        csb = sc%solidC(is)
                        ssb = sc%solidSource(is)
                        if (banded) then
                            if (band_material(phc, dsh) == 2) then
                                csb = sc%outerC(is)
                                ssb = 0.0d0
                            end if
                        end if
                        crw = 0.0d0; cre = 0.0d0; crs = 0.0d0
                        crn = 0.0d0; crb = 0.0d0; crt = 0.0d0
                        ! THE C2 TANGENTIAL CORRECTION (off by default;
                        ! [scalar.N] tangential_correction). The exact
                        ! cut-face flux is
                        !   F = k_face (T_R - T_L)/h_d + s_t (k_loc - k_face)
                        ! and C1 keeps only the first term. The second is
                        ! IDENTICALLY ZERO away from the interface --
                        ! k_loc == k_face there -- so it is evaluated at
                        ! CUT FACES ONLY, and the flux line below stays
                        ! byte-for-byte C1's everywhere else.
                        !
                        ! Both gradients come from the central differences
                        ! the stencil already reaches: the face-normal
                        ! component is the arm difference, each tangential
                        ! component the mean of the two cells' own central
                        ! differences. Written in the fixed low -> high
                        ! cell order, so the two cells sharing a face
                        ! compute the SAME number to the last bit and the
                        ! correction stays a conservative face flux.
                        !
                        ! The tangential terms reach the EDGE ghosts of
                        ! the 3^3 neighbourhood, which the 26-neighbour
                        ! exchange fills -- 1 == 4 ranks is exact, gated.
                        ! A cut face may not sit on a 2:1 block face
                        ! (check_conjugate_refinement), but a cut face one
                        ! cell inside a refined block can still read a
                        ! halo value that came from a level transfer;
                        ! refine_body's one-block buffer is what keeps
                        ! that away from the surface in practice.
                        if (tang) then
                            if (cutw .and. .not. clw) then
                                gtd = (s0 - blk%q(i-1,j,k,var,b))*sc%invDx(i,b)
                                gt1 = 0.5d0*((blk%q(i-1,j+1,k,var,b) - blk%q(i-1,j-1,k,var,b)) &
                                           + (blk%q(i,  j+1,k,var,b) - blk%q(i,  j-1,k,var,b))) &
                                        *sc%cdy(j,b)
                                gt2 = 0.5d0*((blk%q(i-1,j,k+1,var,b) - blk%q(i-1,j,k-1,var,b)) &
                                           + (blk%q(i,  j,k+1,var,b) - blk%q(i,  j,k-1,var,b))) &
                                        *sc%cdz(k,b)
                                gpd = (phc - phw)*sc%invDx(i,b)
                                gp1 = 0.5d0*((sc%phi(i-1,j+1,k,b) - sc%phi(i-1,j-1,k,b)) &
                                           + (sc%phi(i,  j+1,k,b) - sc%phi(i,  j-1,k,b))) &
                                        *sc%cdy(j,b)
                                gp2 = 0.5d0*((sc%phi(i-1,j,k+1,b) - sc%phi(i-1,j,k-1,b)) &
                                           + (sc%phi(i,  j,k+1,b) - sc%phi(i,  j,k-1,b))) &
                                        *sc%cdz(k,b)
                                crw = conjugate_tangential(gtd, gt1, gt2, gpd, gp1, gp2, &
                                        dxw/dm, phw, phc, ks, sc%invDx(i,b)) &
                                    *(conjugate_face_local_k(dm, phw, phc, ks) - dxw)
                            end if
                            if (cute .and. .not. cle) then
                                gtd = (blk%q(i+1,j,k,var,b) - s0)*sc%invDx(i+1,b)
                                gt1 = 0.5d0*((blk%q(i,  j+1,k,var,b) - blk%q(i,  j-1,k,var,b)) &
                                           + (blk%q(i+1,j+1,k,var,b) - blk%q(i+1,j-1,k,var,b))) &
                                        *sc%cdy(j,b)
                                gt2 = 0.5d0*((blk%q(i,  j,k+1,var,b) - blk%q(i,  j,k-1,var,b)) &
                                           + (blk%q(i+1,j,k+1,var,b) - blk%q(i+1,j,k-1,var,b))) &
                                        *sc%cdz(k,b)
                                gpd = (phe - phc)*sc%invDx(i+1,b)
                                gp1 = 0.5d0*((sc%phi(i,  j+1,k,b) - sc%phi(i,  j-1,k,b)) &
                                           + (sc%phi(i+1,j+1,k,b) - sc%phi(i+1,j-1,k,b))) &
                                        *sc%cdy(j,b)
                                gp2 = 0.5d0*((sc%phi(i,  j,k+1,b) - sc%phi(i,  j,k-1,b)) &
                                           + (sc%phi(i+1,j,k+1,b) - sc%phi(i+1,j,k-1,b))) &
                                        *sc%cdz(k,b)
                                cre = conjugate_tangential(gtd, gt1, gt2, gpd, gp1, gp2, &
                                        dxe/dm, phc, phe, ks, sc%invDx(i+1,b)) &
                                    *(conjugate_face_local_k(dm, phc, phe, ks) - dxe)
                            end if
                            if (cuts .and. .not. cls) then
                                gtd = (s0 - blk%q(i,j-1,k,var,b))*sc%invDy(j,b)
                                gt1 = 0.5d0*((blk%q(i+1,j-1,k,var,b) - blk%q(i-1,j-1,k,var,b)) &
                                           + (blk%q(i+1,j,  k,var,b) - blk%q(i-1,j,  k,var,b))) &
                                        *sc%cdx(i,b)
                                gt2 = 0.5d0*((blk%q(i,j-1,k+1,var,b) - blk%q(i,j-1,k-1,var,b)) &
                                           + (blk%q(i,j,  k+1,var,b) - blk%q(i,j,  k-1,var,b))) &
                                        *sc%cdz(k,b)
                                gpd = (phc - phs)*sc%invDy(j,b)
                                gp1 = 0.5d0*((sc%phi(i+1,j-1,k,b) - sc%phi(i-1,j-1,k,b)) &
                                           + (sc%phi(i+1,j,  k,b) - sc%phi(i-1,j,  k,b))) &
                                        *sc%cdx(i,b)
                                gp2 = 0.5d0*((sc%phi(i,j-1,k+1,b) - sc%phi(i,j-1,k-1,b)) &
                                           + (sc%phi(i,j,  k+1,b) - sc%phi(i,j,  k-1,b))) &
                                        *sc%cdz(k,b)
                                crs = conjugate_tangential(gtd, gt1, gt2, gpd, gp1, gp2, &
                                        dys/dm, phs, phc, ks, sc%invDy(j,b)) &
                                    *(conjugate_face_local_k(dm, phs, phc, ks) - dys)
                            end if
                            if (cutn .and. .not. cln) then
                                gtd = (blk%q(i,j+1,k,var,b) - s0)*sc%invDy(j+1,b)
                                gt1 = 0.5d0*((blk%q(i+1,j,  k,var,b) - blk%q(i-1,j,  k,var,b)) &
                                           + (blk%q(i+1,j+1,k,var,b) - blk%q(i-1,j+1,k,var,b))) &
                                        *sc%cdx(i,b)
                                gt2 = 0.5d0*((blk%q(i,j,  k+1,var,b) - blk%q(i,j,  k-1,var,b)) &
                                           + (blk%q(i,j+1,k+1,var,b) - blk%q(i,j+1,k-1,var,b))) &
                                        *sc%cdz(k,b)
                                gpd = (phn - phc)*sc%invDy(j+1,b)
                                gp1 = 0.5d0*((sc%phi(i+1,j,  k,b) - sc%phi(i-1,j,  k,b)) &
                                           + (sc%phi(i+1,j+1,k,b) - sc%phi(i-1,j+1,k,b))) &
                                        *sc%cdx(i,b)
                                gp2 = 0.5d0*((sc%phi(i,j,  k+1,b) - sc%phi(i,j,  k-1,b)) &
                                           + (sc%phi(i,j+1,k+1,b) - sc%phi(i,j+1,k-1,b))) &
                                        *sc%cdz(k,b)
                                crn = conjugate_tangential(gtd, gt1, gt2, gpd, gp1, gp2, &
                                        dyn/dm, phc, phn, ks, sc%invDy(j+1,b)) &
                                    *(conjugate_face_local_k(dm, phc, phn, ks) - dyn)
                            end if
                            if (cutb .and. .not. clb) then
                                gtd = (s0 - blk%q(i,j,k-1,var,b))*sc%invDz(k,b)
                                gt1 = 0.5d0*((blk%q(i+1,j,k-1,var,b) - blk%q(i-1,j,k-1,var,b)) &
                                           + (blk%q(i+1,j,k,  var,b) - blk%q(i-1,j,k,  var,b))) &
                                        *sc%cdx(i,b)
                                gt2 = 0.5d0*((blk%q(i,j+1,k-1,var,b) - blk%q(i,j-1,k-1,var,b)) &
                                           + (blk%q(i,j+1,k,  var,b) - blk%q(i,j-1,k,  var,b))) &
                                        *sc%cdy(j,b)
                                gpd = (phc - phb)*sc%invDz(k,b)
                                gp1 = 0.5d0*((sc%phi(i+1,j,k-1,b) - sc%phi(i-1,j,k-1,b)) &
                                           + (sc%phi(i+1,j,k,  b) - sc%phi(i-1,j,k,  b))) &
                                        *sc%cdx(i,b)
                                gp2 = 0.5d0*((sc%phi(i,j+1,k-1,b) - sc%phi(i,j-1,k-1,b)) &
                                           + (sc%phi(i,j+1,k,  b) - sc%phi(i,j-1,k,  b))) &
                                        *sc%cdy(j,b)
                                crb = conjugate_tangential(gtd, gt1, gt2, gpd, gp1, gp2, &
                                        dzb/dm, phb, phc, ks, sc%invDz(k,b)) &
                                    *(conjugate_face_local_k(dm, phb, phc, ks) - dzb)
                            end if
                            if (cutt .and. .not. clt) then
                                gtd = (blk%q(i,j,k+1,var,b) - s0)*sc%invDz(k+1,b)
                                gt1 = 0.5d0*((blk%q(i+1,j,k,  var,b) - blk%q(i-1,j,k,  var,b)) &
                                           + (blk%q(i+1,j,k+1,var,b) - blk%q(i-1,j,k+1,var,b))) &
                                        *sc%cdx(i,b)
                                gt2 = 0.5d0*((blk%q(i,j+1,k,  var,b) - blk%q(i,j-1,k,  var,b)) &
                                           + (blk%q(i,j+1,k+1,var,b) - blk%q(i,j-1,k+1,var,b))) &
                                        *sc%cdy(j,b)
                                gpd = (pht - phc)*sc%invDz(k+1,b)
                                gp1 = 0.5d0*((sc%phi(i+1,j,k,  b) - sc%phi(i-1,j,k,  b)) &
                                           + (sc%phi(i+1,j,k+1,b) - sc%phi(i-1,j,k+1,b))) &
                                        *sc%cdx(i,b)
                                gp2 = 0.5d0*((sc%phi(i,j+1,k,  b) - sc%phi(i,j-1,k,  b)) &
                                           + (sc%phi(i,j+1,k+1,b) - sc%phi(i,j-1,k+1,b))) &
                                        *sc%cdy(j,b)
                                crt = conjugate_tangential(gtd, gt1, gt2, gpd, gp1, gp2, &
                                        dzt/dm, phc, pht, ks, sc%invDz(k+1,b)) &
                                    *(conjugate_face_local_k(dm, phc, pht, ks) - dzt)
                            end if
                        end if

                        fw = merge(0.0d0, dxw*(s0 - blk%q(i-1,j,k,var,b)) &
                            *sc%invDx(i,b), clw)
                        fe = merge(0.0d0, dxe*(blk%q(i+1,j,k,var,b) - s0) &
                            *sc%invDx(i+1,b), cle)
                        diff = (fe - fw)*blk%d1x(i,VAR_P,b)
                        fw = merge(0.0d0, dys*(s0 - blk%q(i,j-1,k,var,b)) &
                            *sc%invDy(j,b), cls)
                        fe = merge(0.0d0, dyn*(blk%q(i,j+1,k,var,b) - s0) &
                            *sc%invDy(j+1,b), cln)
                        diff = diff + (fe - fw)*blk%d1y(j,VAR_P,b)
                        fw = merge(0.0d0, dzb*(s0 - blk%q(i,j,k-1,var,b)) &
                            *sc%invDz(k,b), clb)
                        fe = merge(0.0d0, dzt*(blk%q(i,j,k+1,var,b) - s0) &
                            *sc%invDz(k+1,b), clt)
                        diff = diff + (fe - fw)*blk%d1z(k,VAR_P,b)
                        ! The C2 correction enters as its own flux divergence,
                        ! ADDED to the finished C1 term rather than folded into
                        ! it: the fused expression above then stays byte for
                        ! byte what it was, so the correction is inert by
                        ! CONSTRUCTION when off (the bodyforce lesson -- not a
                        ! "+ 0.0" argument). It is still a face flux, so
                        ! conservation is untouched: the neighbour subtracts
                        ! the same number.
                        if (tang) diff = diff + (cre - crw)*blk%d1x(i,VAR_P,b) &
                                              + (crn - crs)*blk%d1y(j,VAR_P,b) &
                                              + (crt - crb)*blk%d1z(k,VAR_P,b)

                        ! The volumetric source (see the plain kernel).
                        srcVal = sc%source(is)
                        if (sc%srcType(is) == SC_SRC_VELOCITY) then
                            if (sc%srcDir(is) == 2_C_INT) then
                                srcVal = sc%source(is)*0.5d0*(vs + vn)
                            else if (sc%srcDir(is) == 3_C_INT) then
                                srcVal = sc%source(is)*0.5d0*(wb + wt)
                            else
                                srcVal = sc%source(is)*0.5d0*(uw + ue)
                            end if
                        end if

                        ! Divide the flux divergence by the cell's own
                        ! volumetric capacity, and let the solid carry its
                        ! own volumetric source (Joule heating, ...).
                        !
                        ! C3: the capacity is FLUID-FRACTION WEIGHTED,
                        !   C_cell = f + (1 - f) C_s,
                        ! because a cut cell physically contains both
                        ! materials. It cannot change any steady state (it
                        ! divides an rhs that vanishes there -- C1 gated
                        ! exactly that), and it is what makes the TRANSIENT
                        ! second order: with the pointwise value a cut cell
                        ! carries one material's whole heat capacity, an
                        ! O(1) error on an O(h) band, i.e. first order.
                        ! The volumetric source is weighted with the SAME
                        ! fraction, as the power density f C_f S_f +
                        ! (1-f) C_s S_s it is (C_f = 1); at f = 1 and f = 0
                        ! it reduces to the pointwise C1 statement.
                        rhs = (-conv + diff &
                               + sc%vfrac(i,j,k,b)*srcVal &
                               + (1.0d0 - sc%vfrac(i,j,k,b))*csb*ssb) &
                            /(sc%vfrac(i,j,k,b) &
                              + (1.0d0 - sc%vfrac(i,j,k,b))*csb)
                        ! No penalization: the solid cells are unknowns, not
                        ! boundary values. oldrhs keeps the rhs, as the
                        ! momentum predictor does.
                        ss = s0 + dt_alpha*rhs + dt_beta*blk%oldrhs(i,j,k,scr,b)
                        blk%qs(i,j,k,scr,b) = ss
                        blk%oldrhs(i,j,k,scr,b) = rhs
                    end do
                end do
            end do
        end do
        end do
        !$omp end target teams distribute parallel do
    end subroutine scalar_transport_conjugate

    ! Substage tail: qs -> q for the scalar slots, then the ghost/halo sync.
    ! Runs AFTER the pressure projection returns -- the scalar never enters
    ! the projection loop.
    subroutine scalar_finish(sc, blk, bc, c)
        type(scalar_type), intent(in) :: sc
        type(block_set_type), intent(inout) :: blk
        type(boundary_type), intent(in) :: bc
        type(comm_type), intent(inout) :: c

        integer :: i, j, k, b, is, nx, ny, nz, nBlocks, nScal

        if (.not. scalars_enabled(sc)) return

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))
        nBlocks = int(blk%nBlocks)
        nScal = int(sc%n)

        !$omp target teams distribute parallel do collapse(4) &
        !$omp& map(to: nScal, blk%qs) map(tofrom: blk%q) &
        !$omp& private(i,j,k,b,is)
        do b = 1, nBlocks
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    do is = 1, nScal
                        blk%q(i,j,k,VAR_S0+is,b) = blk%qs(i,j,k,SCR_S0+is,b)
                    end do
                end do
            end do
        end do
        end do
        !$omp end target teams distribute parallel do

        call scalar_sync(sc, blk, bc, c)
    end subroutine scalar_finish

    ! Physical ghosts + one batched halo exchange over all scalar variables
    ! (they ride the same boundary rows, the same entries and the same
    ! message as u,v,w,p).
    subroutine scalar_sync(sc, blk, bc, c)
        type(scalar_type), intent(in) :: sc
        type(block_set_type), intent(inout) :: blk
        type(boundary_type), intent(in) :: bc
        type(comm_type), intent(inout) :: c

        if (.not. scalars_enabled(sc)) return

        call apply_bc(blk, bc, vars=sc%varList)
        call exchange_halos(c, blk, sc%varList)
    end subroutine scalar_sync

end module scalar
