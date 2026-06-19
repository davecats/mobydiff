module pressure_solver
    use, intrinsic :: iso_c_binding
    use :: init, only: dns_type, VAR_U, VAR_V, VAR_W, VAR_P
    use :: blocks, only: block_set_type, FACE_PHYS, FACE_CLOSED, FACE_FINE, FACE_COARSE
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

    ! Pressure at the start of the current projection. The interface velocity is
    ! reconstructed from the in-projection pressure CHANGE (p - pStart): the
    ! predictor already applied the start-pressure interface gradient, so the
    ! change completes it to the full new-pressure gradient without double count.
    real(C_DOUBLE), allocatable :: pStart(:,:,:,:)

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

    subroutine pressure_projection(ps, blk, dt_gamma, ibm, bc, c)
        type(pressure_solver_type), intent(in) :: ps
        type(block_set_type), intent(inout) :: blk
        real(C_DOUBLE), intent(in) :: dt_gamma
        type(ibm_type), intent(in) :: ibm
        type(boundary_type), intent(in) :: bc
        type(comm_type), intent(inout) :: c

        integer(C_INT) :: iIter, color

        ! 2:1 refined grids (nLevels > 1) use the composite projection: one
        ! coupled SPD system whose interface rows are relaxed in situ. The owner
        ! of each shared face (the block ABOVE it, holding it as its interior low
        ! face) reconstructs the face from the pressures in the sweep
        ! (q = qs - (p_above - p_below)*ifGrad); the cross-level coupling is
        ! carried by the per-colour pressure + velocity exchange. The single-
        ! level path below is the plain red-black SOR (no interface faces, so the
        ! sweep's interface reconstruction is inert and the two paths coincide).
        if (blk%nLevels > 1_C_INT) then
            ! The predictor velocity halos (including the cross-interface transfer)
            ! are already current from the main-loop exchange right before this
            ! call; no velocity has changed since, so composite_qs_setup can read
            ! the slaved interface faces directly.
            call snapshot_start_pressure(blk)
            call composite_qs_setup(blk)                          ! freeze the slaved interface qs
            do iIter = 1_C_INT, ps%nIter
                do color = 1_C_INT, 0_C_INT, -1_C_INT
                    call redblack_sweep(ps, blk, dt_gamma, ibm, color)
                    call apply_bc(blk, bc)
                    ! Refresh the cross-level pressure (RESTRICT + inject) and
                    ! reconcile the interface velocity in one exchange. Between
                    ! sweeps the sweep reads a neighbour pressure only across a 2:1
                    ! interface, so pressure is transferred on the interface entries
                    ! only; the final exchange ships it everywhere so the next
                    ! substage's momentum sees current same-level halo pressures.
                    if (iIter == ps%nIter .and. color == 0_C_INT) then
                        ! Final exchange, two passes (E3 + two-phase). The momentum
                        ! predictor of the next substage reads these halos, so the
                        ! coarse->fine velocity is prolonged by linear interpolation
                        ! (O(h^2)). The relaxation above used injection, so its
                        ! contraction is unaffected. Pass A is a plain injection
                        ! exchange that refreshes every coarse halo from the
                        ! post-sweep interiors (same-level copies + restricts);
                        ! pass B then does the linear PROLONG, whose 2-point stencil
                        ! reads those now-current coarse halos (exact second order to
                        ! the patch edges, no clamp, generic across orientations and
                        ! wedged coarse blocks). Pass B re-runs copies/restricts
                        ! idempotently, so the order within it is immaterial.
                        call exchange_halos(c, blk, [VAR_U, VAR_V, VAR_W, VAR_P])
                        call exchange_halos(c, blk, [VAR_U, VAR_V, VAR_W], linear_prolong=.true.)
                    else
                        ! Per-colour: velocity injection (the relaxation contracts
                        ! only with injection halos), but the PROLONG pressure halo
                        ! is smoothed tangentially so the owned-face reconstruction
                        ! reads a non-staircase coarse pressure (fixes the fine-owns
                        ! consistency defect; see comm%pLinProlong).
                        call exchange_halos(c, blk, [VAR_U, VAR_V, VAR_W, VAR_P], &
                            p_interface_only=.true., p_linear_prolong=.true.)
                    end if
                end do
            end do
            return
        end if

        do iIter = 1_C_INT, ps%nIter
            do color = 1_C_INT, 0_C_INT, -1_C_INT
                call redblack_sweep(ps, blk, dt_gamma, ibm, color)
                call apply_bc(blk, bc)
                if (iIter == ps%nIter .and. color == 0_C_INT) then
                    call exchange_halos(c, blk, [VAR_U, VAR_V, VAR_W, VAR_P])
                else
                    call exchange_halos(c, blk, [VAR_U, VAR_V, VAR_W], interp=.false.)
                end if
            end do
        end do
    end subroutine pressure_projection

    ! Allocate the device-resident start-of-projection pressure store once.
    subroutine allocate_pstart(blk)
        type(block_set_type), intent(in) :: blk
        if (.not. allocated(pStart)) then
            ! Match blk%q's bounds exactly: a section like blk%q(:,:,:,VAR_P,:)
            ! would re-base every dimension to 1, but the kernels index pStart
            ! from 0 (the halo layer), so allocate with the real 0-based bounds.
            allocate(pStart(lbound(blk%q,1):ubound(blk%q,1), &
                            lbound(blk%q,2):ubound(blk%q,2), &
                            lbound(blk%q,3):ubound(blk%q,3), &
                            lbound(blk%q,5):ubound(blk%q,5)))
#ifdef USE_OPENMP_OFFLOAD
            !$omp target enter data map(alloc: pStart)
#endif
        end if
    end subroutine allocate_pstart

    ! Snapshot the start-of-projection pressure into pStart (device kernel, no
    ! host round-trip). Used by the interface reconstruction in the sweep.
    subroutine snapshot_start_pressure(blk)
        type(block_set_type), intent(in) :: blk
        integer(C_INT) :: i, j, k, b, nx, ny, nz
        call allocate_pstart(blk)
        nx = blk%nb(1); ny = blk%nb(2); nz = blk%nb(3)
#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do collapse(4) &
        !$omp& map(to: blk%q) map(tofrom: pStart) private(i,j,k,b)
#endif
        do b = 1_C_INT, blk%nBlocks
            do k = 0_C_INT, nz+2
                do j = 0_C_INT, ny+2
                    do i = 0_C_INT, nx+2
                        pStart(i,j,k,b) = blk%q(i,j,k,VAR_P,b)
                    end do
                end do
            end do
        end do
#ifdef USE_OPENMP_OFFLOAD
        !$omp end target teams distribute parallel do
#endif
    end subroutine snapshot_start_pressure

    ! Freeze the slaved interface predictor velocity. The owner of a shared face
    ! is the block above it (its interior low face); the block below holds the
    ! face as a high-side halo whose value is the transferred (RESTRICT/PROLONG)
    ! owner velocity, just written by exchange_halos. Freeze it into qs there so
    ! the sweep reconstruction q = qs - (p_above - p_below)*ifGrad starts from
    ! the owner velocity. The owner keeps its own predictor qs at the low face.
    subroutine composite_qs_setup(blk)
        type(block_set_type), intent(inout) :: blk
        integer(C_INT) :: b, i, j, k, nx, ny, nz
        nx = blk%nb(1); ny = blk%nb(2); nz = blk%nb(3)
#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do collapse(4) &
        !$omp& map(to: blk%physHigh) map(tofrom: blk%qs, blk%q) private(i,j,k,b)
#endif
        do b = 1_C_INT, blk%nBlocks
            do k = 1_C_INT, nz
                do j = 1_C_INT, ny
                    do i = 1_C_INT, nx
                        if (i == nx .and. is_interface(blk%physHigh(1,b))) &
                            blk%qs(nx+1,j,k,VAR_U,b) = blk%q(nx+1,j,k,VAR_U,b)
                        if (j == ny .and. is_interface(blk%physHigh(2,b))) &
                            blk%qs(i,ny+1,k,VAR_V,b) = blk%q(i,ny+1,k,VAR_V,b)
                        if (k == nz .and. is_interface(blk%physHigh(3,b))) &
                            blk%qs(i,j,nz+1,VAR_W,b) = blk%q(i,j,nz+1,VAR_W,b)
                    end do
                end do
            end do
        end do
#ifdef USE_OPENMP_OFFLOAD
        !$omp end target teams distribute parallel do
#endif
    end subroutine composite_qs_setup

    subroutine redblack_sweep(ps, blk, dt_gamma, ibm, color)
        type(pressure_solver_type), intent(in) :: ps
        type(block_set_type), intent(inout) :: blk
        real(C_DOUBLE), intent(in) :: dt_gamma
        type(ibm_type), intent(in) :: ibm
        integer(C_INT), intent(in) :: color

        real(C_DOUBLE) :: phi, denom, div, idt, sor
        real(C_DOUBLE) :: gradU_i, gradU_ip, gradV_j, gradV_jp, gradW_k, gradW_kp
        real(C_DOUBLE) :: mu_u_i, mu_u_ip, mu_v_j, mu_v_jp, mu_w_k, mu_w_kp
        ! True when the cell sits on a 2:1 coarse-fine interface at the low/high
        ! face in x/y/z (enables the pressure-based face reconstruction below).
        logical :: ifLoX, ifHiX, ifLoY, ifHiY, ifLoZ, ifHiZ
        integer(C_INT) :: i, ip, j, jp, k, kp, b, nBlocks, nLowerHaloDirections
        integer(C_INT) :: iColor, nColorX, iLo, jLo, kLo, colorOffset
        integer(C_INT) :: hi(1:3)   ! block interior size (cells per dimension)

        ! Each block sweeps from 0 (its halo layer, redundantly with the
        ! neighbour that owns those cells) except on physical boundaries,
        ! exactly the rank-level scheme one level down. Red-black parity is
        ! anchored to the global index space through the block origin.
        hi = blk%nb(1:3)
        sor = ps%sor
        nColorX = (hi(1) + 2_C_INT)/2_C_INT
        nBlocks = blk%nBlocks
        idt = 1.0_C_DOUBLE/dt_gamma
        ! pStart is device-resident (enter data, once); the map(to:) below then
        ! finds it present and does NOT copy it per sweep. It is read only at
        ! interface faces (none on single-level grids), but must be mapped.
        call allocate_pstart(blk)

#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do collapse(4) &
        !$omp& map(to: color, nColorX, sor, idt, dt_gamma, hi(1:3), &
        !$omp& blk%origin, blk%physLow, blk%physHigh, blk%ifGrad, blk%qs, pStart, &
        !$omp& blk%d1x, blk%d1y, blk%d1z, ibm%mu) &
        !$omp& map(tofrom: blk%q) &
        !$omp& private(i,ip,j,jp,k,kp,b,iColor,iLo,jLo,kLo,colorOffset, &
        !$omp& phi,denom,div,nLowerHaloDirections, &
        !$omp& gradU_i,gradU_ip,gradV_j,gradV_jp,gradW_k,gradW_kp, &
        !$omp& ifLoX,ifHiX,ifLoY,ifHiY,ifLoZ,ifHiZ, &
        !$omp& mu_u_i,mu_u_ip,mu_v_j,mu_v_jp,mu_w_k,mu_w_kp)
#endif
        DO b = 1_C_INT, nBlocks
        DO k = 0_C_INT, hi(3)
            DO j = 0_C_INT, hi(2)
                DO iColor = 0_C_INT, nColorX-1_C_INT
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

                    mu_u_i  = ibm%mu(i,j,k,VAR_U,b)
                    mu_u_ip = ibm%mu(ip,j,k,VAR_U,b)
                    mu_v_j  = ibm%mu(i,j,k,VAR_V,b)
                    mu_v_jp = ibm%mu(i,jp,k,VAR_V,b)
                    mu_w_k  = ibm%mu(i,j,k,VAR_W,b)
                    mu_w_kp = ibm%mu(i,j,kp,VAR_W,b)

                    ! A cell face is a 2:1 interface when its block-boundary side
                    ! carries FACE_FINE/FACE_COARSE. There the cell-centre spacing
                    ! d1 is replaced by the adjoint interface gradient ifGrad
                    ! (1/gap over the true coarse-fine cell-centre distance), and
                    ! the face is reconstructed from the pressures rather than
                    ! incrementally corrected. Inert on single-level grids.
                    ifLoX = i == 1_C_INT     .and. is_interface(blk%physLow(1,b))
                    ifHiX = i == hi(1)       .and. is_interface(blk%physHigh(1,b))
                    ifLoY = j == 1_C_INT     .and. is_interface(blk%physLow(2,b))
                    ifHiY = j == hi(2)       .and. is_interface(blk%physHigh(2,b))
                    ifLoZ = k == 1_C_INT     .and. is_interface(blk%physLow(3,b))
                    ifHiZ = k == hi(3)       .and. is_interface(blk%physHigh(3,b))
                    gradU_i  = merge(blk%ifGrad(1,b), blk%d1x(i,VAR_U,b),  ifLoX)
                    gradU_ip = merge(blk%ifGrad(2,b), blk%d1x(ip,VAR_U,b), ifHiX)
                    gradV_j  = merge(blk%ifGrad(3,b), blk%d1y(j,VAR_V,b),  ifLoY)
                    gradV_jp = merge(blk%ifGrad(4,b), blk%d1y(jp,VAR_V,b), ifHiY)
                    gradW_k  = merge(blk%ifGrad(5,b), blk%d1z(k,VAR_W,b),  ifLoZ)
                    gradW_kp = merge(blk%ifGrad(6,b), blk%d1z(kp,VAR_W,b), ifHiZ)

                    ! The diagonal drops only faces pinned to zero flux forever
                    ! (walls, closed faces). A 2:1 interface flux is a live
                    ! unknown and stays in the diagonal with its ifGrad coefficient.
                    denom = (merge(0.0d0, mu_u_i*gradU_i, &
                                      face_pinned(blk%physLow(1,b)) .and. i == 1_C_INT) &
                           + merge(0.0d0, mu_u_ip*gradU_ip, &
                                      face_pinned(blk%physHigh(1,b)) .and. i == hi(1)))*blk%d1x(i,VAR_P,b) &
                          + (merge(0.0d0, mu_v_j*gradV_j, &
                                      face_pinned(blk%physLow(2,b)) .and. j == 1_C_INT) &
                           + merge(0.0d0, mu_v_jp*gradV_jp, &
                                      face_pinned(blk%physHigh(2,b)) .and. j == hi(2)))*blk%d1y(j,VAR_P,b) &
                          + (merge(0.0d0, mu_w_k*gradW_k, &
                                      face_pinned(blk%physLow(3,b)) .and. k == 1_C_INT) &
                           + merge(0.0d0, mu_w_kp*gradW_kp, &
                                      face_pinned(blk%physHigh(3,b)) .and. k == hi(3)))*blk%d1z(k,VAR_P,b)

                    div = (blk%q(ip,j,k,VAR_U,b)-blk%q(i,j,k,VAR_U,b))*blk%d1x(i,VAR_P,b) &
                        + (blk%q(i,jp,k,VAR_V,b)-blk%q(i,j,k,VAR_V,b))*blk%d1y(j,VAR_P,b) &
                        + (blk%q(i,j,kp,VAR_W,b)-blk%q(i,j,k,VAR_W,b))*blk%d1z(k,VAR_P,b)

                    phi = -sor*div/denom
                    blk%q(i,j,k,VAR_P,b) = blk%q(i,j,k,VAR_P,b) + phi*idt

                    ! Reconstruct each interface face this cell owns/holds from
                    ! the UPDATED pressures: q = qs - dt_gamma*(p_above-p_below)*
                    ! ifGrad, using the in-projection pressure change p - pStart.
                    if (ifLoX) blk%q(i,j,k,VAR_U,b) = blk%qs(i,j,k,VAR_U,b) - dt_gamma*mu_u_i*blk%ifGrad(1,b) &
                        *((blk%q(i,j,k,VAR_P,b)-pStart(i,j,k,b)) - (blk%q(i-1,j,k,VAR_P,b)-pStart(i-1,j,k,b)))
                    if (ifHiX) blk%q(ip,j,k,VAR_U,b) = blk%qs(ip,j,k,VAR_U,b) - dt_gamma*mu_u_ip*blk%ifGrad(2,b) &
                        *((blk%q(ip,j,k,VAR_P,b)-pStart(ip,j,k,b)) - (blk%q(i,j,k,VAR_P,b)-pStart(i,j,k,b)))
                    if (ifLoY) blk%q(i,j,k,VAR_V,b) = blk%qs(i,j,k,VAR_V,b) - dt_gamma*mu_v_j*blk%ifGrad(3,b) &
                        *((blk%q(i,j,k,VAR_P,b)-pStart(i,j,k,b)) - (blk%q(i,j-1,k,VAR_P,b)-pStart(i,j-1,k,b)))
                    if (ifHiY) blk%q(i,jp,k,VAR_V,b) = blk%qs(i,jp,k,VAR_V,b) - dt_gamma*mu_v_jp*blk%ifGrad(4,b) &
                        *((blk%q(i,jp,k,VAR_P,b)-pStart(i,jp,k,b)) - (blk%q(i,j,k,VAR_P,b)-pStart(i,j,k,b)))
                    if (ifLoZ) blk%q(i,j,k,VAR_W,b) = blk%qs(i,j,k,VAR_W,b) - dt_gamma*mu_w_k*blk%ifGrad(5,b) &
                        *((blk%q(i,j,k,VAR_P,b)-pStart(i,j,k,b)) - (blk%q(i,j,k-1,VAR_P,b)-pStart(i,j,k-1,b)))
                    if (ifHiZ) blk%q(i,j,kp,VAR_W,b) = blk%qs(i,j,kp,VAR_W,b) - dt_gamma*mu_w_kp*blk%ifGrad(6,b) &
                        *((blk%q(i,j,kp,VAR_P,b)-pStart(i,j,kp,b)) - (blk%q(i,j,k,VAR_P,b)-pStart(i,j,k,b)))

                    ! Every other face is corrected to drive div -> 0 (both this
                    ! block's and the neighbour's copy of a shared face relax;
                    ! the exchange reconciles them). Pinned and reconstructed
                    ! interface faces are skipped.
                    blk%q(i,j,k,VAR_U,b) = blk%q(i,j,k,VAR_U,b) &
                        - merge(0.0d0, phi*blk%d1x(i,VAR_U,b)*mu_u_i, &
                                (face_pinned(blk%physLow(1,b)) .and. i == 1_C_INT) .or. ifLoX)
                    blk%q(ip,j,k,VAR_U,b) = blk%q(ip,j,k,VAR_U,b) &
                        + merge(0.0d0, phi*blk%d1x(ip,VAR_U,b)*mu_u_ip, &
                                (face_pinned(blk%physHigh(1,b)) .and. i == hi(1)) .or. ifHiX)
                    blk%q(i,j,k,VAR_V,b) = blk%q(i,j,k,VAR_V,b) &
                        - merge(0.0d0, phi*blk%d1y(j,VAR_V,b)*mu_v_j, &
                                (face_pinned(blk%physLow(2,b)) .and. j == 1_C_INT) .or. ifLoY)
                    blk%q(i,jp,k,VAR_V,b) = blk%q(i,jp,k,VAR_V,b) &
                        + merge(0.0d0, phi*blk%d1y(jp,VAR_V,b)*mu_v_jp, &
                                (face_pinned(blk%physHigh(2,b)) .and. j == hi(2)) .or. ifHiY)
                    blk%q(i,j,k,VAR_W,b) = blk%q(i,j,k,VAR_W,b) &
                        - merge(0.0d0, phi*blk%d1z(k,VAR_W,b)*mu_w_k, &
                                (face_pinned(blk%physLow(3,b)) .and. k == 1_C_INT) .or. ifLoZ)
                    blk%q(i,j,kp,VAR_W,b) = blk%q(i,j,kp,VAR_W,b) &
                        + merge(0.0d0, phi*blk%d1z(kp,VAR_W,b)*mu_w_kp, &
                                (face_pinned(blk%physHigh(3,b)) .and. k == hi(3)) .or. ifHiZ)
                END DO
            END DO
        END DO
        END DO
#ifdef USE_OPENMP_OFFLOAD
        !$omp end target teams distribute parallel do
#endif
    end subroutine redblack_sweep

    pure logical function is_interface(fk)
!$omp declare target
        integer(C_INT), intent(in) :: fk
        is_interface = fk == FACE_FINE .or. fk == FACE_COARSE
    end function is_interface

    ! Pinned faces carry zero flux forever (physical walls, closed faces against
    ! removed blocks): they leave both the diagonal and the corrections.
    pure logical function face_pinned(fk)
!$omp declare target
        integer(C_INT), intent(in) :: fk
        face_pinned = fk == FACE_PHYS .or. fk == FACE_CLOSED
    end function face_pinned

end module pressure_solver
