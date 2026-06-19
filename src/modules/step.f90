!--------------------------!
!                          !
!       Time-stepper       !
!          module          !
!                          !
!--------------------------!
!
! authors: Dr.-Ing. Davide Gatti
!          B.Sc. Ahmet Cumhur
!
! date:    28.04.26
!

module step
    use, intrinsic :: iso_c_binding
    use :: init, only: dns_type, VAR_U, VAR_V, VAR_W, VAR_P, &
        CFL_COURANT, CFL_PECLET, NCFL
    use :: blocks, only: block_set_type, FACE_PHYS, FACE_CLOSED, FACE_COARSE, FACE_FINE
    use :: ibmm, only: ibm_type
    use :: comm, only: comm_type, comm_allreduce_max
    use :: les_model, only: les_type, les_is_enabled, les_profile_type, &
        les_wall_seconds, add_les_profile, LES_PROF_SGS
    implicit none

    real(C_DOUBLE), parameter :: rk_alpha(3) = [64.0d0/120.0d0,  50.0d0/120.0d0,  90.0d0/120.0d0]
    real(C_DOUBLE), parameter :: rk_beta(3)  = [ 0.0d0,         -34.0d0/120.0d0, -50.0d0/120.0d0]
    real(C_DOUBLE), parameter :: rk_gamma(3) = [64.0d0/120.0d0,  16.0d0/120.0d0,  40.0d0/120.0d0]

contains

    subroutine precompute_peclet_rate(dns, blk, c)
        type(dns_type), intent(inout) :: dns
        type(block_set_type), intent(in) :: blk
        type(comm_type), intent(in) :: c

        integer :: i, b, nx, ny, nz
        real(C_DOUBLE) :: ire, local_rate(1)

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))
        ire = 1.0d0/dns%re

        local_rate = 0.0d0
        do b = 1, int(blk%nBlocks)
            do i = 0, nx+1
                local_rate(1) = max(local_rate(1), ire*blk%d1x(i,VAR_P,b)**2)
            end do
            do i = 0, ny+1
                local_rate(1) = max(local_rate(1), ire*blk%d1y(i,VAR_P,b)**2)
            end do
            do i = 0, nz+1
                local_rate(1) = max(local_rate(1), ire*blk%d1z(i,VAR_P,b)**2)
            end do
        end do

        call comm_allreduce_max(c, local_rate)
        dns%peclet_rate = local_rate(1)
    end subroutine precompute_peclet_rate

    subroutine momentum(blk, dns, dt_alpha, dt_beta, dt_gamma, ibm, les, les_prof)
        type(block_set_type), intent(inout) :: blk
        type(dns_type),   intent(in)    :: dns
        real(C_DOUBLE),   intent(in)    :: dt_alpha, dt_beta, dt_gamma
        type(ibm_type),   intent(in)    :: ibm
        type(les_type),   intent(in), optional :: les
        type(les_profile_type), intent(inout), optional :: les_prof

        integer :: i,j,k,b,ip,im,kp,km,jp,jm
        integer :: nx, ny, nz, nBlocks
        ! First predicted face index per component: 2 if the low face is pinned
        ! by a wall / physical BC (skip it), else 1 (a halo or interface face).
        integer :: uStartX, vStartY, wStartZ

        real(C_DOUBLE) :: diff_ux,diff_uy,diff_uz
        real(C_DOUBLE) :: diff_vx,diff_vy,diff_vz
        real(C_DOUBLE) :: diff_wx,diff_wy,diff_wz

        real(C_DOUBLE) :: uu_p,uu_m,uv_p,uv_m,uw_p,uw_m
        real(C_DOUBLE) :: vu_p,vu_m,vv_p,vv_m,vw_p,vw_m
        real(C_DOUBLE) :: wu_p,wu_m,ww_p,ww_m,wv_p,wv_m

        real(C_DOUBLE) :: dpx,dpy,dpz,rhsu,rhsv,rhsw
        real(C_DOUBLE) :: mu_u, mu_v, mu_w
        real(C_DOUBLE) :: ire
        real(C_DOUBLE) :: forcing(1:3)
        logical :: use_les

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))
        nBlocks = int(blk%nBlocks)
        ire = 1.0d0/dns%re
        forcing = dns%forcing
        use_les = .false.
        if (present(les)) use_les = les_is_enabled(les) .and. allocated(les%nut)

        ! Predictor for all staggered velocity components. The face on a
        ! physical lower boundary is held by apply_bc: each block skips it via
        ! its own physLow mask.
        !$omp target teams distribute parallel do collapse(4) &
        !$omp& map(to: dt_alpha, dt_beta, dt_gamma, &
        !$omp& ire, forcing(1:3), &
        !$omp& blk%physLow, blk%physHigh, blk%ifGrad, blk%d1x, blk%d1y, blk%d1z, &
        !$omp& blk%lapXm, blk%lapX0, blk%lapXp, blk%lapYm, blk%lapY0, blk%lapYp, &
        !$omp& blk%lapZm, blk%lapZ0, blk%lapZp, blk%q, ibm%mu) &
        !$omp& map(tofrom: blk%qs, blk%oldrhs) &
        !$omp& private(i,j,k,b,ip,im,jp,jm,kp,km,uStartX,vStartY,wStartZ, &
        !$omp& uu_p,uu_m,uv_p,uv_m,uw_p,uw_m, &
        !$omp& vu_p,vu_m,vv_p,vv_m,vw_p,vw_m,wu_p,wu_m,ww_p,ww_m,wv_p,wv_m, &
        !$omp& diff_ux,diff_uy,diff_uz,diff_vx,diff_vy,diff_vz,diff_wx,diff_wy,diff_wz, &
        !$omp& dpx,dpy,dpz,rhsu,rhsv,rhsw,mu_u,mu_v,mu_w)
        do b = 1, nBlocks
        do k = 1, nz+1
            do j = 1, ny+1
                do i = 1, nx+1
                    ip = i+1
                    im = i-1
                    jp = j+1
                    jm = j-1
                    kp = k+1
                    km = k-1
                    ! Only physical walls and closed faces are pinned (skipped).
                    ! The interior low face of a 2:1 interface is the owner's
                    ! shared face and IS predicted; the high face is a halo
                    ! supplied by the exchange, so every block ends its loop at nb.
                    uStartX = merge(2, 1, blk%physLow(1,b) == FACE_PHYS .or. blk%physLow(1,b) == FACE_CLOSED)
                    vStartY = merge(2, 1, blk%physLow(2,b) == FACE_PHYS .or. blk%physLow(2,b) == FACE_CLOSED)
                    wStartZ = merge(2, 1, blk%physLow(3,b) == FACE_PHYS .or. blk%physLow(3,b) == FACE_CLOSED)

                    if (i >= uStartX .and. i <= nx .and. j <= ny .and. k <= nz) then
                        uu_p = (blk%q(i,j,k,VAR_U,b) + blk%q(ip,j,k,VAR_U,b))**2
                        uu_m = (blk%q(im,j,k,VAR_U,b) + blk%q(i,j,k,VAR_U,b))**2
                        ! Convective-flux reflux at a fine-owns 2:1 interface (i=1,
                        ! coarse below): the stored coarse-side normal velocity u(0)
                        ! is injected at the wrong position (the coarse face one fine
                        ! cell deeper). Its correct value at the fine halo position is
                        ! 0.5*(u(0)+u(1)); use it in the interface convective flux so
                        ! the owned-face momentum matches the coarse side. The stored
                        ! halo (read by the divergence) is untouched -> mass exact.
                        if (i == 1 .and. blk%physLow(1,b) == FACE_COARSE) &
                            uu_m = (0.5d0*blk%q(im,j,k,VAR_U,b) + 1.5d0*blk%q(i,j,k,VAR_U,b))**2

                        uv_p = (blk%q(i,j,k,VAR_U,b) + blk%q(i,jp,k,VAR_U,b)) &
                             * (blk%q(im,jp,k,VAR_V,b) + blk%q(i,jp,k,VAR_V,b))
                        uv_m = (blk%q(i,jm,k,VAR_U,b) + blk%q(i,j,k,VAR_U,b)) &
                             * (blk%q(im,j,k,VAR_V,b) + blk%q(i,j,k,VAR_V,b))

                        uw_p = (blk%q(i,j,k,VAR_U,b) + blk%q(i,j,kp,VAR_U,b)) &
                             * (blk%q(im,j,kp,VAR_W,b) + blk%q(i,j,kp,VAR_W,b))
                        uw_m = (blk%q(i,j,km,VAR_U,b) + blk%q(i,j,k,VAR_U,b)) &
                             * (blk%q(im,j,k,VAR_W,b) + blk%q(i,j,k,VAR_W,b))

                        diff_ux = blk%lapXm(i,VAR_U,b)*blk%q(im,j,k,VAR_U,b) &
                                + blk%lapX0(i,VAR_U,b)*blk%q(i,j,k,VAR_U,b) &
                                + blk%lapXp(i,VAR_U,b)*blk%q(ip,j,k,VAR_U,b)
                        diff_uy = blk%lapYm(j,VAR_U,b)*blk%q(i,jm,k,VAR_U,b) &
                                + blk%lapY0(j,VAR_U,b)*blk%q(i,j,k,VAR_U,b) &
                                + blk%lapYp(j,VAR_U,b)*blk%q(i,jp,k,VAR_U,b)
                        diff_uz = blk%lapZm(k,VAR_U,b)*blk%q(i,j,km,VAR_U,b) &
                                + blk%lapZ0(k,VAR_U,b)*blk%q(i,j,k,VAR_U,b) &
                                + blk%lapZp(k,VAR_U,b)*blk%q(i,j,kp,VAR_U,b)

                        ! At the owner's low interface face the adjoint gradient
                        ! ifGrad (over the coarse-fine gap) replaces d1, matching
                        ! the projection so predictor and correction stay consistent.
                        dpx = (blk%q(i,j,k,VAR_P,b)-blk%q(im,j,k,VAR_P,b)) &
                            * merge(blk%ifGrad(1,b), blk%d1x(i,VAR_U,b), i == 1 .and. &
                              (blk%physLow(1,b) == FACE_FINE .or. blk%physLow(1,b) == FACE_COARSE))

                        rhsu = ( &
                            -0.25d0*( (uu_p-uu_m)*blk%d1x(i,VAR_U,b) &
                                     +(uv_p-uv_m)*blk%d1y(j,VAR_U,b) &
                                     +(uw_p-uw_m)*blk%d1z(k,VAR_U,b)) &
                            + forcing(VAR_U) &
                            + ire*(diff_ux + diff_uy + diff_uz) )

                        blk%qs(i,j,k,VAR_U,b) = blk%q(i,j,k,VAR_U,b) + dt_alpha*rhsu &
                            + dt_beta*blk%oldrhs(i,j,k,VAR_U,b) - dt_gamma*dpx

                        mu_u = ibm%mu(i,j,k,VAR_U,b)
                        blk%qs(i,j,k,VAR_U,b) = blk%qs(i,j,k,VAR_U,b)*mu_u

                        blk%oldrhs(i,j,k,VAR_U,b) = rhsu
                    end if

                    if (j >= vStartY .and. j <= ny .and. i <= nx .and. k <= nz) then
                        vu_p = (blk%q(i,j,k,VAR_V,b) + blk%q(ip,j,k,VAR_V,b)) &
                             * (blk%q(ip,jm,k,VAR_U,b) + blk%q(ip,j,k,VAR_U,b))
                        vu_m = (blk%q(im,j,k,VAR_V,b) + blk%q(i,j,k,VAR_V,b)) &
                             * (blk%q(i,jm,k,VAR_U,b) + blk%q(i,j,k,VAR_U,b))

                        vv_p = (blk%q(i,j,k,VAR_V,b) + blk%q(i,jp,k,VAR_V,b))**2
                        vv_m = (blk%q(i,jm,k,VAR_V,b) + blk%q(i,j,k,VAR_V,b))**2
                        if (j == 1 .and. blk%physLow(2,b) == FACE_COARSE) &   ! fine-owns reflux (see uu_m)
                            vv_m = (0.5d0*blk%q(i,jm,k,VAR_V,b) + 1.5d0*blk%q(i,j,k,VAR_V,b))**2

                        vw_p = (blk%q(i,j,k,VAR_V,b) + blk%q(i,j,kp,VAR_V,b)) &
                             * (blk%q(i,jm,kp,VAR_W,b) + blk%q(i,j,kp,VAR_W,b))
                        vw_m = (blk%q(i,j,km,VAR_V,b) + blk%q(i,j,k,VAR_V,b)) &
                             * (blk%q(i,jm,k,VAR_W,b) + blk%q(i,j,k,VAR_W,b))

                        diff_vx = blk%lapXm(i,VAR_V,b)*blk%q(im,j,k,VAR_V,b) &
                                + blk%lapX0(i,VAR_V,b)*blk%q(i,j,k,VAR_V,b) &
                                + blk%lapXp(i,VAR_V,b)*blk%q(ip,j,k,VAR_V,b)
                        diff_vy = blk%lapYm(j,VAR_V,b)*blk%q(i,jm,k,VAR_V,b) &
                                + blk%lapY0(j,VAR_V,b)*blk%q(i,j,k,VAR_V,b) &
                                + blk%lapYp(j,VAR_V,b)*blk%q(i,jp,k,VAR_V,b)
                        diff_vz = blk%lapZm(k,VAR_V,b)*blk%q(i,j,km,VAR_V,b) &
                                + blk%lapZ0(k,VAR_V,b)*blk%q(i,j,k,VAR_V,b) &
                                + blk%lapZp(k,VAR_V,b)*blk%q(i,j,kp,VAR_V,b)

                        dpy = (blk%q(i,j,k,VAR_P,b)-blk%q(i,jm,k,VAR_P,b)) &
                            * merge(blk%ifGrad(3,b), blk%d1y(j,VAR_V,b), j == 1 .and. &
                              (blk%physLow(2,b) == FACE_FINE .or. blk%physLow(2,b) == FACE_COARSE))

                        rhsv = ( &
                            -0.25d0*((vu_p-vu_m)*blk%d1x(i,VAR_V,b) &
                                    +(vv_p-vv_m)*blk%d1y(j,VAR_V,b) &
                                    +(vw_p-vw_m)*blk%d1z(k,VAR_V,b)) &
                            + forcing(VAR_V) &
                            + ire*(diff_vx + diff_vy + diff_vz) )

                        blk%qs(i,j,k,VAR_V,b) = blk%q(i,j,k,VAR_V,b) + dt_alpha*rhsv &
                            + dt_beta*blk%oldrhs(i,j,k,VAR_V,b) - dt_gamma*dpy

                        mu_v = ibm%mu(i,j,k,VAR_V,b)
                        blk%qs(i,j,k,VAR_V,b) = blk%qs(i,j,k,VAR_V,b)*mu_v

                        blk%oldrhs(i,j,k,VAR_V,b) = rhsv
                    end if

                    if (k >= wStartZ .and. k <= nz .and. i <= nx .and. j <= ny) then
                        wu_p = (blk%q(i,j,k,VAR_W,b) + blk%q(ip,j,k,VAR_W,b)) &
                             * (blk%q(ip,j,km,VAR_U,b) + blk%q(ip,j,k,VAR_U,b))
                        wu_m = (blk%q(im,j,k,VAR_W,b) + blk%q(i,j,k,VAR_W,b)) &
                             * (blk%q(i,j,km,VAR_U,b) + blk%q(i,j,k,VAR_U,b))

                        ww_p = (blk%q(i,j,k,VAR_W,b) + blk%q(i,j,kp,VAR_W,b))**2
                        ww_m = (blk%q(i,j,km,VAR_W,b) + blk%q(i,j,k,VAR_W,b))**2
                        if (k == 1 .and. blk%physLow(3,b) == FACE_COARSE) &   ! fine-owns reflux (see uu_m)
                            ww_m = (0.5d0*blk%q(i,j,km,VAR_W,b) + 1.5d0*blk%q(i,j,k,VAR_W,b))**2

                        wv_p = (blk%q(i,j,k,VAR_W,b) + blk%q(i,jp,k,VAR_W,b)) &
                             * (blk%q(i,jp,km,VAR_V,b) + blk%q(i,jp,k,VAR_V,b))
                        wv_m = (blk%q(i,jm,k,VAR_W,b) + blk%q(i,j,k,VAR_W,b)) &
                             * (blk%q(i,j,km,VAR_V,b) + blk%q(i,j,k,VAR_V,b))

                        diff_wx = blk%lapXm(i,VAR_W,b)*blk%q(im,j,k,VAR_W,b) &
                                + blk%lapX0(i,VAR_W,b)*blk%q(i,j,k,VAR_W,b) &
                                + blk%lapXp(i,VAR_W,b)*blk%q(ip,j,k,VAR_W,b)
                        diff_wy = blk%lapYm(j,VAR_W,b)*blk%q(i,jm,k,VAR_W,b) &
                                + blk%lapY0(j,VAR_W,b)*blk%q(i,j,k,VAR_W,b) &
                                + blk%lapYp(j,VAR_W,b)*blk%q(i,jp,k,VAR_W,b)
                        diff_wz = blk%lapZm(k,VAR_W,b)*blk%q(i,j,km,VAR_W,b) &
                                + blk%lapZ0(k,VAR_W,b)*blk%q(i,j,k,VAR_W,b) &
                                + blk%lapZp(k,VAR_W,b)*blk%q(i,j,kp,VAR_W,b)

                        dpz = (blk%q(i,j,k,VAR_P,b)-blk%q(i,j,km,VAR_P,b)) &
                            * merge(blk%ifGrad(5,b), blk%d1z(k,VAR_W,b), k == 1 .and. &
                              (blk%physLow(3,b) == FACE_FINE .or. blk%physLow(3,b) == FACE_COARSE))

                        rhsw = ( &
                            -0.25d0*( (wu_p-wu_m)*blk%d1x(i,VAR_W,b) &
                                     +(wv_p-wv_m)*blk%d1y(j,VAR_W,b) &
                                     +(ww_p-ww_m)*blk%d1z(k,VAR_W,b)) &
                            + forcing(VAR_W) &
                            + ire*(diff_wx + diff_wy + diff_wz) )

                        blk%qs(i,j,k,VAR_W,b) = blk%q(i,j,k,VAR_W,b) + dt_alpha*rhsw &
                            + dt_beta*blk%oldrhs(i,j,k,VAR_W,b) - dt_gamma*dpz

                        mu_w = ibm%mu(i,j,k,VAR_W,b)
                        blk%qs(i,j,k,VAR_W,b) = blk%qs(i,j,k,VAR_W,b)*mu_w

                        blk%oldrhs(i,j,k,VAR_W,b) = rhsw
                    end if
                end do
            end do
        end do
        end do
        !$omp end target teams distribute parallel do

        if (present(les)) then
            if (use_les) then
                if (present(les_prof)) then
                    call add_les_momentum_correction(blk, dns, dt_alpha, ibm, les, les_prof)
                else
                    call add_les_momentum_correction(blk, dns, dt_alpha, ibm, les)
                end if
            end if
        end if

        ! Commit the predictor: copy qs -> q on every interior face (1..nb per
        ! component). Interface high faces (nb+1) are halos refreshed by the
        ! exchange, so they are not committed here.
        !$omp target teams distribute parallel do collapse(4) &
        !$omp& map(to: blk%qs) &
        !$omp& map(tofrom: blk%q) &
        !$omp& private(i,j,k,b)
        do b = 1, nBlocks
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    blk%q(i,j,k,VAR_U,b) = blk%qs(i,j,k,VAR_U,b)
                    blk%q(i,j,k,VAR_V,b) = blk%qs(i,j,k,VAR_V,b)
                    blk%q(i,j,k,VAR_W,b) = blk%qs(i,j,k,VAR_W,b)
                end do
            end do
        end do
        end do
        !$omp end target teams distribute parallel do

    end subroutine momentum

    ! Truncation-error probe (MOBY_TRUNC). With blk%q set to the EXACT TGV field
    ! plus the scheme's interface halos, print the discrete u-momentum operator
    ! terms along a y-z RMS at each x-index of the first fine-owns block
    ! (physLow(1)==FACE_COARSE): convection C, pressure P=-dp/dx, viscous D, and
    ! the inviscid balance C+P. For the TGV the exact convection cancels the exact
    ! pressure gradient, so C+P -> 0 with the grid in a consistent scheme; an O(1)
    ! C+P at the interface (i=1) that does not shrink with h is the order-0 defect,
    ! and the C-vs-P split says whether it is convective or pressure.
    subroutine truncation_probe(blk, dns)
        type(block_set_type), intent(inout) :: blk
        type(dns_type), intent(in) :: dns
        integer :: i,j,k,b,ip,im,jp,jm,kp,km, nx,ny,nz, bsel, npl, kk, kind
        real(C_DOUBLE) :: ire, uu_p,uu_m,uv_p,uv_m,uw_p,uw_m, conv, pterm, dterm
        real(C_DOUBLE) :: dux,duy,duz, sC,sP,sCP,sD, grad
        logical :: iface

        ire = 1.0d0/dns%re
        nx = int(blk%nb(1)); ny = int(blk%nb(2)); nz = int(blk%nb(3))
        npl = ny*nz
#ifdef USE_OPENMP_OFFLOAD
        !$omp target update from(blk%q)
#endif
      do kk = 1, 2
        kind = merge(FACE_COARSE, FACE_FINE, kk == 1)   ! fine-owns, then coarse-owns
        bsel = -1
        do b = 1, int(blk%nBlocks)
            if (blk%physLow(1,b) == kind) then; bsel = b; exit; end if
        end do
        if (bsel < 0) cycle
        b = bsel
        print '(a,a,3i5,a,es12.4)', merge("TRUNC FINE-OWNS  ", "TRUNC COARSE-OWNS", kk == 1), &
            " block origin", int(blk%origin(:,b)), "  x_interface=", blk%x(1,VAR_U,b)
        print '(a)', "   i        x        |conv|        |pres|       |conv+pres|     |visc|"
        do i = 1, min(nx, 8)
            ip=i+1; im=i-1
            iface = (i == 1)
            grad = merge(blk%ifGrad(1,b), blk%d1x(i,VAR_U,b), iface)
            sC=0; sP=0; sCP=0; sD=0
            do k = 1, nz
                kp=k+1; km=k-1
                do j = 1, ny
                    jp=j+1; jm=j-1
                    uu_p = (blk%q(i,j,k,VAR_U,b)+blk%q(ip,j,k,VAR_U,b))**2
                    uu_m = (blk%q(im,j,k,VAR_U,b)+blk%q(i,j,k,VAR_U,b))**2
                    if (iface) uu_m = (0.5d0*blk%q(im,j,k,VAR_U,b)+1.5d0*blk%q(i,j,k,VAR_U,b))**2
                    uv_p = (blk%q(i,j,k,VAR_U,b)+blk%q(i,jp,k,VAR_U,b)) &
                         * (blk%q(im,jp,k,VAR_V,b)+blk%q(i,jp,k,VAR_V,b))
                    uv_m = (blk%q(i,jm,k,VAR_U,b)+blk%q(i,j,k,VAR_U,b)) &
                         * (blk%q(im,j,k,VAR_V,b)+blk%q(i,j,k,VAR_V,b))
                    uw_p = (blk%q(i,j,k,VAR_U,b)+blk%q(i,j,kp,VAR_U,b)) &
                         * (blk%q(im,j,kp,VAR_W,b)+blk%q(i,j,kp,VAR_W,b))
                    uw_m = (blk%q(i,j,km,VAR_U,b)+blk%q(i,j,k,VAR_U,b)) &
                         * (blk%q(im,j,k,VAR_W,b)+blk%q(i,j,k,VAR_W,b))
                    conv = -0.25d0*((uu_p-uu_m)*blk%d1x(i,VAR_U,b) &
                                   +(uv_p-uv_m)*blk%d1y(j,VAR_U,b) &
                                   +(uw_p-uw_m)*blk%d1z(k,VAR_U,b))
                    pterm = -(blk%q(i,j,k,VAR_P,b)-blk%q(im,j,k,VAR_P,b))*grad
                    dux = blk%lapXm(i,VAR_U,b)*blk%q(im,j,k,VAR_U,b) &
                        + blk%lapX0(i,VAR_U,b)*blk%q(i,j,k,VAR_U,b) &
                        + blk%lapXp(i,VAR_U,b)*blk%q(ip,j,k,VAR_U,b)
                    duy = blk%lapYm(j,VAR_U,b)*blk%q(i,jm,k,VAR_U,b) &
                        + blk%lapY0(j,VAR_U,b)*blk%q(i,j,k,VAR_U,b) &
                        + blk%lapYp(j,VAR_U,b)*blk%q(i,jp,k,VAR_U,b)
                    duz = blk%lapZm(k,VAR_U,b)*blk%q(i,j,km,VAR_U,b) &
                        + blk%lapZ0(k,VAR_U,b)*blk%q(i,j,k,VAR_U,b) &
                        + blk%lapZp(k,VAR_U,b)*blk%q(i,j,kp,VAR_U,b)
                    dterm = ire*(dux+duy+duz)
                    sC = sC + conv*conv; sP = sP + pterm*pterm
                    sCP = sCP + (conv+pterm)**2; sD = sD + dterm*dterm
                end do
            end do
            print '(i5,f10.4,4es14.4)', i, blk%x(i,VAR_U,b), &
                sqrt(sC/npl), sqrt(sP/npl), sqrt(sCP/npl), sqrt(sD/npl)
        end do
      end do
    end subroutine truncation_probe


    subroutine add_les_momentum_correction(blk, dns, dt_alpha, ibm, les, les_prof)
        type(block_set_type), intent(inout) :: blk
        type(dns_type), intent(in) :: dns
        real(C_DOUBLE), intent(in) :: dt_alpha
        type(ibm_type), intent(in) :: ibm
        type(les_type), intent(in) :: les
        type(les_profile_type), intent(inout), optional :: les_prof

        integer :: i, j, k, b, ip, im, jp, jm, kp, km
        integer :: nx, ny, nz, nBlocks
        ! First predicted face index per component: 2 if the low face is pinned
        ! by a wall / physical BC (skip it), else 1 (a halo or interface face).
        integer :: uStartX, vStartY, wStartZ
        real(C_DOUBLE) :: sgs_u, sgs_v, sgs_w
        real(C_DOUBLE) :: tau_xp, tau_xm, tau_yp, tau_ym, tau_zp, tau_zm
        real(C_DOUBLE) :: nut_edge, nut0, nut1, wx, wy, mu_u, mu_v, mu_w
        real(C_DOUBLE) :: profile_start

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))
        nBlocks = int(blk%nBlocks)

        if (present(les_prof)) profile_start = les_wall_seconds()

        !$omp target teams distribute parallel do collapse(4) &
        !$omp& map(to: dt_alpha, &
        !$omp& blk%physLow, blk%d1x, blk%d1y, blk%d1z, blk%q, ibm%mu, les%nut, &
        !$omp& les%u_from_p_x, les%v_from_p_y, les%w_from_p_z, &
        !$omp& les%inv_dx, les%inv_dy, les%inv_dz) &
        !$omp& map(tofrom: blk%qs, blk%oldrhs) &
        !$omp& private(i,j,k,b,ip,im,jp,jm,kp,km,uStartX,vStartY,wStartZ, &
        !$omp& sgs_u,sgs_v,sgs_w,tau_xp,tau_xm, &
        !$omp& tau_yp,tau_ym,tau_zp,tau_zm,nut_edge,nut0,nut1,wx,wy,mu_u,mu_v,mu_w)
        do b = 1, nBlocks
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    ip = i + 1
                    im = i - 1
                    jp = j + 1
                    jm = j - 1
                    kp = k + 1
                    km = k - 1
                    ! Only physical walls and closed faces are pinned. Both
                    ! sides of a 2:1 interface predict the shared face (the
                    ! restriction overwrites the coarse copy, Phase 3c).
                    uStartX = merge(2, 1, blk%physLow(1,b) == FACE_PHYS .or. blk%physLow(1,b) == FACE_CLOSED)
                    vStartY = merge(2, 1, blk%physLow(2,b) == FACE_PHYS .or. blk%physLow(2,b) == FACE_CLOSED)
                    wStartZ = merge(2, 1, blk%physLow(3,b) == FACE_PHYS .or. blk%physLow(3,b) == FACE_CLOSED)

                    if (i >= uStartX) then
                        tau_xp = 2.0d0*les%nut(i,j,k,b) &
                               * (blk%q(ip,j,k,VAR_U,b) - blk%q(i,j,k,VAR_U,b))*blk%d1x(i,VAR_P,b)
                        tau_xm = 2.0d0*les%nut(im,j,k,b) &
                               * (blk%q(i,j,k,VAR_U,b) - blk%q(im,j,k,VAR_U,b))*blk%d1x(im,VAR_P,b)

                        wx = les%u_from_p_x(i,b)
                        wy = les%v_from_p_y(jp,b)
                        nut0 = (1.0d0 - wx)*les%nut(im,j,k,b) + wx*les%nut(i,j,k,b)
                        nut1 = (1.0d0 - wx)*les%nut(im,jp,k,b) + wx*les%nut(i,jp,k,b)
                        nut_edge = (1.0d0 - wy)*nut0 + wy*nut1
                        tau_yp = nut_edge*( &
                            (blk%q(i,jp,k,VAR_U,b) - blk%q(i,j,k,VAR_U,b))*les%inv_dy(jp,VAR_U,b) &
                          + (blk%q(i,jp,k,VAR_V,b) - blk%q(im,jp,k,VAR_V,b))*les%inv_dx(i,VAR_V,b) )

                        wx = les%u_from_p_x(i,b)
                        wy = les%v_from_p_y(j,b)
                        nut0 = (1.0d0 - wx)*les%nut(im,jm,k,b) + wx*les%nut(i,jm,k,b)
                        nut1 = (1.0d0 - wx)*les%nut(im,j,k,b) + wx*les%nut(i,j,k,b)
                        nut_edge = (1.0d0 - wy)*nut0 + wy*nut1
                        tau_ym = nut_edge*( &
                            (blk%q(i,j,k,VAR_U,b) - blk%q(i,jm,k,VAR_U,b))*les%inv_dy(j,VAR_U,b) &
                          + (blk%q(i,j,k,VAR_V,b) - blk%q(im,j,k,VAR_V,b))*les%inv_dx(i,VAR_V,b) )

                        wx = les%u_from_p_x(i,b)
                        wy = les%w_from_p_z(kp,b)
                        nut0 = (1.0d0 - wx)*les%nut(im,j,k,b) + wx*les%nut(i,j,k,b)
                        nut1 = (1.0d0 - wx)*les%nut(im,j,kp,b) + wx*les%nut(i,j,kp,b)
                        nut_edge = (1.0d0 - wy)*nut0 + wy*nut1
                        tau_zp = nut_edge*( &
                            (blk%q(i,j,kp,VAR_U,b) - blk%q(i,j,k,VAR_U,b))*les%inv_dz(kp,VAR_U,b) &
                          + (blk%q(i,j,kp,VAR_W,b) - blk%q(im,j,kp,VAR_W,b))*les%inv_dx(i,VAR_W,b) )

                        wx = les%u_from_p_x(i,b)
                        wy = les%w_from_p_z(k,b)
                        nut0 = (1.0d0 - wx)*les%nut(im,j,km,b) + wx*les%nut(i,j,km,b)
                        nut1 = (1.0d0 - wx)*les%nut(im,j,k,b) + wx*les%nut(i,j,k,b)
                        nut_edge = (1.0d0 - wy)*nut0 + wy*nut1
                        tau_zm = nut_edge*( &
                            (blk%q(i,j,k,VAR_U,b) - blk%q(i,j,km,VAR_U,b))*les%inv_dz(k,VAR_U,b) &
                          + (blk%q(i,j,k,VAR_W,b) - blk%q(im,j,k,VAR_W,b))*les%inv_dx(i,VAR_W,b) )

                        sgs_u = (tau_xp - tau_xm)*blk%d1x(i,VAR_U,b) &
                              + (tau_yp - tau_ym)*blk%d1y(j,VAR_U,b) &
                              + (tau_zp - tau_zm)*blk%d1z(k,VAR_U,b)

                        mu_u = ibm%mu(i,j,k,VAR_U,b)
                        blk%qs(i,j,k,VAR_U,b) = blk%qs(i,j,k,VAR_U,b) + dt_alpha*sgs_u*mu_u
                        blk%oldrhs(i,j,k,VAR_U,b) = blk%oldrhs(i,j,k,VAR_U,b) + sgs_u
                    end if

                    if (j >= vStartY) then
                        wx = les%u_from_p_x(ip,b)
                        wy = les%v_from_p_y(j,b)
                        nut0 = (1.0d0 - wx)*les%nut(i,jm,k,b) + wx*les%nut(ip,jm,k,b)
                        nut1 = (1.0d0 - wx)*les%nut(i,j,k,b) + wx*les%nut(ip,j,k,b)
                        nut_edge = (1.0d0 - wy)*nut0 + wy*nut1
                        tau_xp = nut_edge*( &
                            (blk%q(ip,j,k,VAR_U,b) - blk%q(ip,jm,k,VAR_U,b))*les%inv_dy(j,VAR_U,b) &
                          + (blk%q(ip,j,k,VAR_V,b) - blk%q(i,j,k,VAR_V,b))*les%inv_dx(ip,VAR_V,b) )

                        wx = les%u_from_p_x(i,b)
                        wy = les%v_from_p_y(j,b)
                        nut0 = (1.0d0 - wx)*les%nut(im,jm,k,b) + wx*les%nut(i,jm,k,b)
                        nut1 = (1.0d0 - wx)*les%nut(im,j,k,b) + wx*les%nut(i,j,k,b)
                        nut_edge = (1.0d0 - wy)*nut0 + wy*nut1
                        tau_xm = nut_edge*( &
                            (blk%q(i,j,k,VAR_U,b) - blk%q(i,jm,k,VAR_U,b))*les%inv_dy(j,VAR_U,b) &
                          + (blk%q(i,j,k,VAR_V,b) - blk%q(im,j,k,VAR_V,b))*les%inv_dx(i,VAR_V,b) )

                        tau_yp = 2.0d0*les%nut(i,j,k,b) &
                               * (blk%q(i,jp,k,VAR_V,b) - blk%q(i,j,k,VAR_V,b))*blk%d1y(j,VAR_P,b)
                        tau_ym = 2.0d0*les%nut(i,jm,k,b) &
                               * (blk%q(i,j,k,VAR_V,b) - blk%q(i,jm,k,VAR_V,b))*blk%d1y(jm,VAR_P,b)

                        wx = les%v_from_p_y(j,b)
                        wy = les%w_from_p_z(kp,b)
                        nut0 = (1.0d0 - wx)*les%nut(i,jm,k,b) + wx*les%nut(i,j,k,b)
                        nut1 = (1.0d0 - wx)*les%nut(i,jm,kp,b) + wx*les%nut(i,j,kp,b)
                        nut_edge = (1.0d0 - wy)*nut0 + wy*nut1
                        tau_zp = nut_edge*( &
                            (blk%q(i,j,kp,VAR_V,b) - blk%q(i,j,k,VAR_V,b))*les%inv_dz(kp,VAR_V,b) &
                          + (blk%q(i,j,kp,VAR_W,b) - blk%q(i,jm,kp,VAR_W,b))*les%inv_dy(j,VAR_W,b) )

                        wx = les%v_from_p_y(j,b)
                        wy = les%w_from_p_z(k,b)
                        nut0 = (1.0d0 - wx)*les%nut(i,jm,km,b) + wx*les%nut(i,j,km,b)
                        nut1 = (1.0d0 - wx)*les%nut(i,jm,k,b) + wx*les%nut(i,j,k,b)
                        nut_edge = (1.0d0 - wy)*nut0 + wy*nut1
                        tau_zm = nut_edge*( &
                            (blk%q(i,j,k,VAR_V,b) - blk%q(i,j,km,VAR_V,b))*les%inv_dz(k,VAR_V,b) &
                          + (blk%q(i,j,k,VAR_W,b) - blk%q(i,jm,k,VAR_W,b))*les%inv_dy(j,VAR_W,b) )

                        sgs_v = (tau_xp - tau_xm)*blk%d1x(i,VAR_V,b) &
                              + (tau_yp - tau_ym)*blk%d1y(j,VAR_V,b) &
                              + (tau_zp - tau_zm)*blk%d1z(k,VAR_V,b)

                        mu_v = ibm%mu(i,j,k,VAR_V,b)
                        blk%qs(i,j,k,VAR_V,b) = blk%qs(i,j,k,VAR_V,b) + dt_alpha*sgs_v*mu_v
                        blk%oldrhs(i,j,k,VAR_V,b) = blk%oldrhs(i,j,k,VAR_V,b) + sgs_v
                    end if

                    if (k >= wStartZ) then
                        wx = les%u_from_p_x(ip,b)
                        wy = les%w_from_p_z(k,b)
                        nut0 = (1.0d0 - wx)*les%nut(i,j,km,b) + wx*les%nut(ip,j,km,b)
                        nut1 = (1.0d0 - wx)*les%nut(i,j,k,b) + wx*les%nut(ip,j,k,b)
                        nut_edge = (1.0d0 - wy)*nut0 + wy*nut1
                        tau_xp = nut_edge*( &
                            (blk%q(ip,j,k,VAR_U,b) - blk%q(ip,j,km,VAR_U,b))*les%inv_dz(k,VAR_U,b) &
                          + (blk%q(ip,j,k,VAR_W,b) - blk%q(i,j,k,VAR_W,b))*les%inv_dx(ip,VAR_W,b) )

                        wx = les%u_from_p_x(i,b)
                        wy = les%w_from_p_z(k,b)
                        nut0 = (1.0d0 - wx)*les%nut(im,j,km,b) + wx*les%nut(i,j,km,b)
                        nut1 = (1.0d0 - wx)*les%nut(im,j,k,b) + wx*les%nut(i,j,k,b)
                        nut_edge = (1.0d0 - wy)*nut0 + wy*nut1
                        tau_xm = nut_edge*( &
                            (blk%q(i,j,k,VAR_U,b) - blk%q(i,j,km,VAR_U,b))*les%inv_dz(k,VAR_U,b) &
                          + (blk%q(i,j,k,VAR_W,b) - blk%q(im,j,k,VAR_W,b))*les%inv_dx(i,VAR_W,b) )

                        wx = les%v_from_p_y(jp,b)
                        wy = les%w_from_p_z(k,b)
                        nut0 = (1.0d0 - wx)*les%nut(i,j,km,b) + wx*les%nut(i,jp,km,b)
                        nut1 = (1.0d0 - wx)*les%nut(i,j,k,b) + wx*les%nut(i,jp,k,b)
                        nut_edge = (1.0d0 - wy)*nut0 + wy*nut1
                        tau_yp = nut_edge*( &
                            (blk%q(i,jp,k,VAR_V,b) - blk%q(i,jp,km,VAR_V,b))*les%inv_dz(k,VAR_V,b) &
                          + (blk%q(i,jp,k,VAR_W,b) - blk%q(i,j,k,VAR_W,b))*les%inv_dy(jp,VAR_W,b) )

                        wx = les%v_from_p_y(j,b)
                        wy = les%w_from_p_z(k,b)
                        nut0 = (1.0d0 - wx)*les%nut(i,jm,km,b) + wx*les%nut(i,j,km,b)
                        nut1 = (1.0d0 - wx)*les%nut(i,jm,k,b) + wx*les%nut(i,j,k,b)
                        nut_edge = (1.0d0 - wy)*nut0 + wy*nut1
                        tau_ym = nut_edge*( &
                            (blk%q(i,j,k,VAR_V,b) - blk%q(i,j,km,VAR_V,b))*les%inv_dz(k,VAR_V,b) &
                          + (blk%q(i,j,k,VAR_W,b) - blk%q(i,jm,k,VAR_W,b))*les%inv_dy(j,VAR_W,b) )

                        tau_zp = 2.0d0*les%nut(i,j,k,b) &
                               * (blk%q(i,j,kp,VAR_W,b) - blk%q(i,j,k,VAR_W,b))*blk%d1z(k,VAR_P,b)
                        tau_zm = 2.0d0*les%nut(i,j,km,b) &
                               * (blk%q(i,j,k,VAR_W,b) - blk%q(i,j,km,VAR_W,b))*blk%d1z(km,VAR_P,b)

                        sgs_w = (tau_xp - tau_xm)*blk%d1x(i,VAR_W,b) &
                              + (tau_yp - tau_ym)*blk%d1y(j,VAR_W,b) &
                              + (tau_zp - tau_zm)*blk%d1z(k,VAR_W,b)

                        mu_w = ibm%mu(i,j,k,VAR_W,b)
                        blk%qs(i,j,k,VAR_W,b) = blk%qs(i,j,k,VAR_W,b) + dt_alpha*sgs_w*mu_w
                        blk%oldrhs(i,j,k,VAR_W,b) = blk%oldrhs(i,j,k,VAR_W,b) + sgs_w
                    end if
                end do
            end do
        end do
        end do
        !$omp end target teams distribute parallel do

        if (present(les_prof)) call add_les_profile(les_prof, LES_PROF_SGS, &
            les_wall_seconds() - profile_start)
    end subroutine add_les_momentum_correction

    subroutine get_timestep_rates(blk, dns, rates, les)
        type(block_set_type), intent(inout) :: blk
        type(dns_type),   intent(in)    :: dns
        real(C_DOUBLE), intent(out) :: rates(1:NCFL)
        type(les_type), intent(in), optional :: les

        integer :: i,j,k,b
        integer :: nx, ny, nz, nBlocks
        real(C_DOUBLE) :: cfl_rate, peclet_rate, ire, nu_eff
        logical :: use_les

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))
        nBlocks = int(blk%nBlocks)
        cfl_rate = 0.0d0
        use_les = .false.
        if (present(les)) use_les = les_is_enabled(les) .and. allocated(les%nut)

        !$omp target teams distribute parallel do collapse(4) reduction(max:cfl_rate) &
        !$omp& map(to: blk%d1x, blk%d1y, blk%d1z, blk%q) &
        !$omp& private(i,j,k,b)
        do b = 1, nBlocks
        do k = 0, nz+1
            do j = 0, ny+1
                do i = 0, nx+1
                    cfl_rate = max(cfl_rate, abs(blk%q(i,j,k,VAR_U,b)*blk%d1x(i,VAR_U,b)))
                    cfl_rate = max(cfl_rate, abs(blk%q(i,j,k,VAR_V,b)*blk%d1y(j,VAR_V,b)))
                    cfl_rate = max(cfl_rate, abs(blk%q(i,j,k,VAR_W,b)*blk%d1z(k,VAR_W,b)))
                end do
            end do
        end do
        end do
        !$omp end target teams distribute parallel do

        rates(CFL_COURANT) = cfl_rate
        rates(CFL_PECLET) = dns%peclet_rate
        if (.not. use_les) return

        peclet_rate = dns%peclet_rate
        ire = 1.0d0/dns%re

        !$omp target teams distribute parallel do collapse(4) reduction(max:peclet_rate) &
        !$omp& map(to: blk%d1x, blk%d1y, blk%d1z, les%nut) &
        !$omp& private(i,j,k,b,nu_eff)
        do b = 1, nBlocks
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    nu_eff = ire + max(0.0d0, les%nut(i,j,k,b))
                    peclet_rate = max(peclet_rate, nu_eff*blk%d1x(i,VAR_P,b)**2)
                    peclet_rate = max(peclet_rate, nu_eff*blk%d1y(j,VAR_P,b)**2)
                    peclet_rate = max(peclet_rate, nu_eff*blk%d1z(k,VAR_P,b)**2)
                end do
            end do
        end do
        end do
        !$omp end target teams distribute parallel do

        rates(CFL_PECLET) = peclet_rate
    end subroutine get_timestep_rates

    logical function run_should_continue(dns, loop_steps)
        type(dns_type), intent(in) :: dns
        integer(C_INT), intent(in) :: loop_steps

        real(C_DOUBLE) :: time_tol

        run_should_continue = .true.
        if (dns%nsteps > 0_C_INT .and. loop_steps >= dns%nsteps) run_should_continue = .false.
        if (dns%t_final > 0.0d0) then
            time_tol = max(1.0d-12, 100.0d0*epsilon(1.0d0)*max(1.0d0, abs(dns%t_final)))
            if (dns%t_current >= dns%t_final - time_tol) run_should_continue = .false.
        end if
    end function run_should_continue

    subroutine trim_dt_for_final_time(dns)
        type(dns_type), intent(inout) :: dns

        real(C_DOUBLE) :: remaining

        if (dns%t_final <= 0.0d0) return

        remaining = dns%t_final - dns%t_current
        dns%dt = min(dns%dt, max(0.0d0, remaining))
    end subroutine trim_dt_for_final_time

    subroutine update_timestep_limits(blk, dns, c, les)
        type(block_set_type), intent(inout) :: blk
        type(dns_type), intent(inout) :: dns
        type(comm_type), intent(in) :: c
        type(les_type), intent(in), optional :: les

        real(C_DOUBLE) :: rates(1:NCFL), next_dt
        logical :: have_limit

        if (dns%cflmax <= 0.0d0 .and. dns%pecletmax <= 0.0d0) return

        if (present(les)) then
            call get_timestep_rates(blk, dns, rates, les)
        else
            call get_timestep_rates(blk, dns, rates)
        end if
        call comm_allreduce_max(c, rates)

        ! An active limit caps dt at min(dtmax, the limits) -- seed dtmax so dt
        ! can also grow toward it, not just shrink from the current value.
        next_dt = dns%dtmax
        have_limit = .false.
        if (dns%cflmax > 0.0d0 .and. rates(CFL_COURANT) > 0.0d0) then
            next_dt = min(next_dt, dns%cflmax/rates(CFL_COURANT))
            have_limit = .true.
        end if
        if (dns%pecletmax > 0.0d0 .and. rates(CFL_PECLET) > 0.0d0) then
            next_dt = min(next_dt, dns%pecletmax/rates(CFL_PECLET))
            have_limit = .true.
        end if

        if (have_limit) dns%dt = next_dt
        dns%cfl(CFL_COURANT) = dns%dt*rates(CFL_COURANT)
        dns%cfl(CFL_PECLET) = dns%dt*rates(CFL_PECLET)
    end subroutine update_timestep_limits

end module step
