!--------------------------!
!                          !
!   Geometry-agnostic      !
!   wall distance          !
!                          !
!--------------------------!
!
! IDDES phase T1b (docs/next_session_iddes.md): unsigned Euclidean distance
! to the surface of ANY analytic immersed body, computed from its
! inside/outside indicator alone (isInBody in production, passed as a
! procedure argument so gates can drive the same machinery with a second
! geometry). Host-only init code — nothing here runs in the time loop.
!
! Method:
!  1. Surface point cloud: every finest-refinement-level grid-line segment
!     between adjacent cell centres whose endpoints straddle the indicator
!     is bisected to a surface point (the finest level so refine_body
!     leaves see the surface at their own h). Deterministic two-pass scan
!     (count, prefix, fill) over z-planes, so the cloud is independent of
!     the OpenMP thread count and of the rank layout (every rank builds
!     the same global cloud).
!  2. Nearest cloud point: a kd-tree over the cloud (median split, bbox
!     pruning). A uniform-bin ring search was considered instead but
!     degenerates to O((d/bin)^3) bin visits for cells far from the
!     surface; the tree is O(log n) everywhere and needs no bin tuning on
!     stretched grids. Periodic directions are handled by querying the
!     +-L images of the point and folding the result to the minimum image.
!  3. Polish: the raw cloud distance overestimates by O(s^2/d) (s = local
!     sample spacing) — worst at the near-wall cells the RANS model cares
!     about. The surface is re-sampled on a shrinking 3x3x3 lattice around
!     the current nearest point (bisecting the straddling lattice
!     segments), halving the spacing until the distance error bound drops
!     below the requested tolerance ([rans] dwall_tol).
!
! In periodic directions the indicator must be length-periodic — the same
! assumption the IBM coefficient machinery already makes when it evaluates
! isInBody at halo coordinates beyond the domain.

module walldist
    use, intrinsic :: iso_c_binding
    use :: init, only: dns_type, cell_center_at
    use :: ibmm, only: ibm_type
    implicit none

    private
    public :: body_indicator_i, walldist_type
    public :: build_walldist, destroy_walldist, walldist_distance

    ! The indicator signature is exactly isInBody's (ibm.f90).
    abstract interface
        logical function body_indicator_i(xIN, ibm, dns)
            import :: C_DOUBLE, ibm_type, dns_type
            real(C_DOUBLE), intent(in) :: xIN(1:3)
            type(ibm_type), intent(in) :: ibm
            type(dns_type), intent(in) :: dns
        end function body_indicator_i
    end interface

    integer, parameter :: LEAF_SIZE = 12         ! kd-tree leaf bucket
    integer, parameter :: MAX_POLISH_LEVELS = 80
    integer, parameter :: MAX_BISECT_ITER = 200
    ! Bisection cannot resolve the surface position below round-off of the
    ! domain coordinates; this floor keeps tolerances meaningful.
    real(C_DOUBLE), parameter :: BISECT_FLOOR = 1.0d-13
    real(C_DOUBLE), parameter :: CLOUD_BISECT_TOL = 1.0d-12

    type :: walldist_type
        integer :: nCloud = 0
        ! Surface samples (kd-tree order) and the local finest-grid spacing
        ! that generated each one (the polish start scale).
        real(C_DOUBLE), allocatable :: pts(:,:)     ! (3, nCloud)
        real(C_DOUBLE), allocatable :: scale(:)     ! (nCloud)
        ! kd-tree nodes: point range [lo, hi], children (0 = leaf), split
        ! dimension and per-node bounding box.
        integer :: nNodes = 0
        integer, allocatable :: nodeLo(:), nodeHi(:)
        integer, allocatable :: nodeLeft(:), nodeRight(:)
        integer, allocatable :: nodeDim(:)
        real(C_DOUBLE), allocatable :: bbMin(:,:), bbMax(:,:)  ! (3, nNodes)
        ! Domain topology for the minimum-image metric.
        real(C_DOUBLE) :: leng(3) = 0.0d0
        real(C_DOUBLE) :: x0(3) = 0.0d0
        logical :: periodic(3) = .false.
    end type walldist_type

contains

    ! Build the surface cloud + kd-tree from the finest-level node lines
    ! (blk%line? column nLevels) and the indicator. n = finest-level cell
    ! counts per direction.
    subroutine build_walldist(w, f, ibm, dns, lineX, lineY, lineZ, n, leng, periodic, has_terminal)
        type(walldist_type), intent(inout) :: w
        procedure(body_indicator_i) :: f
        type(ibm_type), intent(in) :: ibm
        type(dns_type), intent(in) :: dns
        real(C_DOUBLE), intent(in) :: lineX(0:), lineY(0:), lineZ(0:)
        integer, intent(in) :: n(3)
        real(C_DOUBLE), intent(in) :: leng(3)
        logical, intent(in) :: periodic(3)
        logical, intent(in) :: has_terminal

        real(C_DOUBLE), allocatable :: xs(:), ys(:), zs(:)   ! cell centres 0:n+1
        real(C_DOUBLE), allocatable :: wx(:), wy(:), wz(:)   ! local cell widths
        integer, allocatable :: cnt(:), offs(:)
        logical(C_BOOL) :: perb(3)
        integer :: lo(3), hi(3), i, kz, total

        call destroy_walldist(w)
        w%leng = leng
        w%periodic = periodic
        w%x0 = [lineX(0), lineY(0), lineZ(0)]
        perb = logical(periodic, C_BOOL)

        ! Sample range per direction: interior centres for periodic
        ! directions (the ghost row is the wrap image of an interior one),
        ! ghost-inclusive for walls (bodies clipped by the domain boundary).
        do i = 1, 3
            lo(i) = merge(1, 0, periodic(i))
            hi(i) = merge(n(i), n(i) + 1, periodic(i))
        end do

        allocate(xs(0:n(1)+1), ys(0:n(2)+1), zs(0:n(3)+1))
        do i = 0, n(1)+1
            xs(i) = cell_center_at(lineX, n(1), leng(1), i, perb(1))
        end do
        do i = 0, n(2)+1
            ys(i) = cell_center_at(lineY, n(2), leng(2), i, perb(2))
        end do
        do i = 0, n(3)+1
            zs(i) = cell_center_at(lineZ, n(3), leng(3), i, perb(3))
        end do
        allocate(wx(0:n(1)), wy(0:n(2)), wz(0:n(3)))
        wx = xs(1:n(1)+1) - xs(0:n(1))
        wy = ys(1:n(2)+1) - ys(0:n(2))
        wz = zs(1:n(3)+1) - zs(0:n(3))

        ! Pass 1: count crossings per z-plane; pass 2: fill each plane's
        ! points at its prefix offset (fixed in-plane order, so the cloud
        ! is deterministic under any thread count).
        allocate(cnt(lo(3):hi(3)), offs(lo(3):hi(3)))
        !$omp parallel do schedule(dynamic)
        do kz = lo(3), hi(3)
            call scan_plane(f, ibm, dns, kz, xs, ys, zs, wx, wy, wz, &
                lo, hi, n, .false., 0, w, cnt(kz))
        end do
        !$omp end parallel do

        total = 0
        do kz = lo(3), hi(3)
            offs(kz) = total
            total = total + cnt(kz)
        end do
        w%nCloud = total
        if (has_terminal) print '(A,I0,A)', " wall-distance surface cloud: ", total, " points"
        if (total == 0) then
            if (has_terminal) print *, "warning: the body indicator produced no surface", &
                " crossings; wall distance left at its no-wall value"
            return
        end if

        allocate(w%pts(3, total), w%scale(total))
        !$omp parallel do schedule(dynamic)
        do kz = lo(3), hi(3)
            call scan_plane(f, ibm, dns, kz, xs, ys, zs, wx, wy, wz, &
                lo, hi, n, .true., offs(kz), w, cnt(kz))
        end do
        !$omp end parallel do

        call kd_build(w)
    end subroutine build_walldist

    subroutine destroy_walldist(w)
        type(walldist_type), intent(inout) :: w

        if (allocated(w%pts)) deallocate(w%pts)
        if (allocated(w%scale)) deallocate(w%scale)
        if (allocated(w%nodeLo)) deallocate(w%nodeLo, w%nodeHi, w%nodeLeft, &
            w%nodeRight, w%nodeDim, w%bbMin, w%bbMax)
        w%nCloud = 0
        w%nNodes = 0
    end subroutine destroy_walldist

    ! Scan one z-plane: in-plane x/y segments at zs(kz), plus the z segments
    ! between planes kz and kz+1 (when kz <= n(3)). do_fill = .false. only
    ! counts; .true. writes points starting at slot offs+1.
    subroutine scan_plane(f, ibm, dns, kz, xs, ys, zs, wx, wy, wz, lo, hi, n, do_fill, offs, w, cnt)
        procedure(body_indicator_i) :: f
        type(ibm_type), intent(in) :: ibm
        type(dns_type), intent(in) :: dns
        integer, intent(in) :: kz, lo(3), hi(3), n(3), offs
        real(C_DOUBLE), intent(in) :: xs(0:), ys(0:), zs(0:)
        real(C_DOUBLE), intent(in) :: wx(0:), wy(0:), wz(0:)
        logical, intent(in) :: do_fill
        type(walldist_type), intent(inout) :: w
        integer, intent(out) :: cnt

        logical, allocatable :: e0(:,:), e1(:,:)
        real(C_DOUBLE) :: xA(3), xB(3), scal
        integer :: i, j

        allocate(e0(lo(1):n(1)+1, lo(2):n(2)+1))
        do j = lo(2), n(2)+1
            do i = lo(1), n(1)+1
                e0(i,j) = f([xs(i), ys(j), zs(kz)], ibm, dns)
            end do
        end do
        if (kz <= n(3)) then
            allocate(e1(lo(1):n(1)+1, lo(2):n(2)+1))
            do j = lo(2), n(2)+1
                do i = lo(1), n(1)+1
                    e1(i,j) = f([xs(i), ys(j), zs(kz+1)], ibm, dns)
                end do
            end do
        end if

        cnt = 0
        ! x segments
        do j = lo(2), hi(2)
            do i = lo(1), n(1)
                if (e0(i,j) .neqv. e0(i+1,j)) then
                    cnt = cnt + 1
                    if (do_fill) then
                        xA = [xs(i),   ys(j), zs(kz)]
                        xB = [xs(i+1), ys(j), zs(kz)]
                        scal = max(wx(i), wy(clampi(j, lo(2), n(2))), &
                                          wz(clampi(kz, lo(3), n(3))))
                        call store_point(w, f, ibm, dns, xA, xB, scal, offs + cnt)
                    end if
                end if
            end do
        end do
        ! y segments
        do j = lo(2), n(2)
            do i = lo(1), hi(1)
                if (e0(i,j) .neqv. e0(i,j+1)) then
                    cnt = cnt + 1
                    if (do_fill) then
                        xA = [xs(i), ys(j),   zs(kz)]
                        xB = [xs(i), ys(j+1), zs(kz)]
                        scal = max(wx(clampi(i, lo(1), n(1))), wy(j), &
                                   wz(clampi(kz, lo(3), n(3))))
                        call store_point(w, f, ibm, dns, xA, xB, scal, offs + cnt)
                    end if
                end if
            end do
        end do
        ! z segments toward plane kz+1
        if (kz <= n(3)) then
            do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                    if (e0(i,j) .neqv. e1(i,j)) then
                        cnt = cnt + 1
                        if (do_fill) then
                            xA = [xs(i), ys(j), zs(kz)]
                            xB = [xs(i), ys(j), zs(kz+1)]
                            scal = max(wx(clampi(i, lo(1), n(1))), &
                                       wy(clampi(j, lo(2), n(2))), wz(kz))
                            call store_point(w, f, ibm, dns, xA, xB, scal, offs + cnt)
                        end if
                    end if
                end do
            end do
        end if
    end subroutine scan_plane

    integer function clampi(v, vlo, vhi) result(r)
        integer, intent(in) :: v, vlo, vhi
        r = min(max(v, vlo), vhi)
    end function clampi

    ! Bisect the straddling segment to a surface point and store it, folded
    ! into the primary domain along periodic directions (so a single +-L
    ! image query covers the minimum-image metric).
    subroutine store_point(w, f, ibm, dns, xA, xB, scal, slot)
        type(walldist_type), intent(inout) :: w
        procedure(body_indicator_i) :: f
        type(ibm_type), intent(in) :: ibm
        type(dns_type), intent(in) :: dns
        real(C_DOUBLE), intent(in) :: xA(3), xB(3), scal
        integer, intent(in) :: slot

        real(C_DOUBLE) :: p(3)
        integer :: d

        call bisect_indicator(f, ibm, dns, xA, xB, CLOUD_BISECT_TOL, p)
        do d = 1, 3
            if (w%periodic(d)) p(d) = w%x0(d) + modulo(p(d) - w%x0(d), w%leng(d))
        end do
        w%pts(:, slot) = p
        w%scale(slot) = scal
    end subroutine store_point

    ! Generic host-side bisection on the indicator (the device-resident
    ! bisection in ibm.f90 is hardwired to isInBody; procedure arguments do
    ! not survive OpenMP target offload, hence this host twin).
    subroutine bisect_indicator(f, ibm, dns, xA_in, xB_in, tol, xS)
        procedure(body_indicator_i) :: f
        type(ibm_type), intent(in) :: ibm
        type(dns_type), intent(in) :: dns
        real(C_DOUBLE), intent(in) :: xA_in(3), xB_in(3), tol
        real(C_DOUBLE), intent(out) :: xS(3)

        real(C_DOUBLE) :: xA(3), xB(3), xM(3)
        logical :: sideA
        integer :: it

        xA = xA_in
        xB = xB_in
        sideA = f(xA, ibm, dns)
        do it = 1, MAX_BISECT_ITER
            if (sqrt(sum((xB - xA)**2)) <= tol) exit
            xM = 0.5d0*(xA + xB)
            if (f(xM, ibm, dns) .eqv. sideA) then
                xA = xM
            else
                xB = xM
            end if
        end do
        xS = 0.5d0*(xA + xB)
    end subroutine bisect_indicator

    !========================
    ! kd-tree
    !========================

    subroutine kd_build(w)
        type(walldist_type), intent(inout) :: w

        integer :: maxNodes, root

        maxNodes = 4*(w%nCloud/LEAF_SIZE + 2)
        allocate(w%nodeLo(maxNodes), w%nodeHi(maxNodes), w%nodeLeft(maxNodes), &
            w%nodeRight(maxNodes), w%nodeDim(maxNodes))
        allocate(w%bbMin(3, maxNodes), w%bbMax(3, maxNodes))
        w%nNodes = 0
        root = kd_build_node(w, 1, w%nCloud)
    end subroutine kd_build

    recursive integer function kd_build_node(w, plo, phi) result(id)
        type(walldist_type), intent(inout) :: w
        integer, intent(in) :: plo, phi

        integer :: d, mid, splitDim
        real(C_DOUBLE) :: bmin(3), bmax(3)

        w%nNodes = w%nNodes + 1
        id = w%nNodes
        w%nodeLo(id) = plo
        w%nodeHi(id) = phi
        do d = 1, 3
            bmin(d) = minval(w%pts(d, plo:phi))
            bmax(d) = maxval(w%pts(d, plo:phi))
        end do
        w%bbMin(:, id) = bmin
        w%bbMax(:, id) = bmax

        if (phi - plo + 1 <= LEAF_SIZE) then
            w%nodeLeft(id) = 0
            w%nodeRight(id) = 0
            w%nodeDim(id) = 0
            return
        end if

        splitDim = maxloc(bmax - bmin, dim=1)
        mid = (plo + phi)/2
        call kd_select(w, plo, phi, mid, splitDim)
        w%nodeDim(id) = splitDim
        w%nodeLeft(id) = kd_build_node(w, plo, mid)
        w%nodeRight(id) = kd_build_node(w, mid + 1, phi)
    end function kd_build_node

    ! Quickselect: after the call pts(dim, plo:k) <= pts(dim, k+1:phi).
    subroutine kd_select(w, plo_in, phi_in, k, dim)
        type(walldist_type), intent(inout) :: w
        integer, intent(in) :: plo_in, phi_in, k, dim

        integer :: plo, phi, i, j
        real(C_DOUBLE) :: pivot

        plo = plo_in
        phi = phi_in
        do while (plo < phi)
            ! median-of-three pivot (deterministic)
            pivot = median3(w%pts(dim, plo), w%pts(dim, (plo+phi)/2), w%pts(dim, phi))
            i = plo
            j = phi
            do
                do while (w%pts(dim, i) < pivot)
                    i = i + 1
                end do
                do while (w%pts(dim, j) > pivot)
                    j = j - 1
                end do
                if (i >= j) exit
                call swap_point(w, i, j)
                i = i + 1
                j = j - 1
            end do
            if (k <= j) then
                phi = j
            else
                plo = j + 1
            end if
        end do
    end subroutine kd_select

    real(C_DOUBLE) function median3(a, b, c) result(m)
        real(C_DOUBLE), intent(in) :: a, b, c
        m = max(min(a, b), min(max(a, b), c))
    end function median3

    subroutine swap_point(w, a, b)
        type(walldist_type), intent(inout) :: w
        integer, intent(in) :: a, b

        real(C_DOUBLE) :: t(3), ts

        t = w%pts(:, a); w%pts(:, a) = w%pts(:, b); w%pts(:, b) = t
        ts = w%scale(a); w%scale(a) = w%scale(b); w%scale(b) = ts
    end subroutine swap_point

    real(C_DOUBLE) function bbox_dist2(w, id, x) result(d2)
        type(walldist_type), intent(in) :: w
        integer, intent(in) :: id
        real(C_DOUBLE), intent(in) :: x(3)

        real(C_DOUBLE) :: dd(3)

        dd = max(0.0d0, w%bbMin(:, id) - x, x - w%bbMax(:, id))
        d2 = sum(dd*dd)
    end function bbox_dist2

    recursive subroutine kd_query(w, id, x, best2, bestIdx)
        type(walldist_type), intent(in) :: w
        integer, intent(in) :: id
        real(C_DOUBLE), intent(in) :: x(3)
        real(C_DOUBLE), intent(inout) :: best2
        integer, intent(inout) :: bestIdx

        integer :: p, near, far
        real(C_DOUBLE) :: d2

        if (bbox_dist2(w, id, x) >= best2) return

        if (w%nodeLeft(id) == 0) then
            do p = w%nodeLo(id), w%nodeHi(id)
                d2 = sum((w%pts(:, p) - x)**2)
                if (d2 < best2) then
                    best2 = d2
                    bestIdx = p
                end if
            end do
            return
        end if

        if (x(w%nodeDim(id)) <= w%bbMax(w%nodeDim(id), w%nodeLeft(id))) then
            near = w%nodeLeft(id)
            far = w%nodeRight(id)
        else
            near = w%nodeRight(id)
            far = w%nodeLeft(id)
        end if
        call kd_query(w, near, x, best2, bestIdx)
        call kd_query(w, far, x, best2, bestIdx)
    end subroutine kd_query

    !========================
    ! distance query + polish
    !========================

    ! Unsigned distance from x to the indicator's surface, to absolute
    ! tolerance ~tol: nearest cloud point (minimum image over periodic
    ! directions), then local surface re-sampling around it.
    real(C_DOUBLE) function walldist_distance(w, f, ibm, dns, x, tol) result(d)
        type(walldist_type), intent(in) :: w
        procedure(body_indicator_i) :: f
        type(ibm_type), intent(in) :: ibm
        type(dns_type), intent(in) :: dns
        real(C_DOUBLE), intent(in) :: x(3), tol

        integer :: kx, ky, kz, m(3), bestIdx, dd
        real(C_DOUBLE) :: best2, xq(3), p(3), dx

        if (w%nCloud == 0) error stop "walldist_distance: empty surface cloud"

        do dd = 1, 3
            m(dd) = merge(1, 0, w%periodic(dd))
        end do

        best2 = huge(1.0d0)
        bestIdx = 0
        do kz = -m(3), m(3)
            do ky = -m(2), m(2)
                do kx = -m(1), m(1)
                    xq = x + [kx, ky, kz]*w%leng
                    if (bbox_dist2(w, 1, xq) >= best2) cycle
                    call kd_query(w, 1, xq, best2, bestIdx)
                end do
            end do
        end do

        ! Fold the winner to the minimum image in x's frame.
        p = w%pts(:, bestIdx)
        do dd = 1, 3
            if (w%periodic(dd)) then
                dx = p(dd) - x(dd)
                p(dd) = x(dd) + dx - w%leng(dd)*nint(dx/w%leng(dd))
            end if
        end do

        d = polish_distance(f, ibm, dns, x, p, w%scale(bestIdx), tol)
    end function walldist_distance

    ! Shrinking-neighbourhood polish: re-sample the surface on a 3x3x3
    ! lattice of spacing sigma around the current nearest point p (bisecting
    ! the straddling lattice segments), keep the closest crossing, halve
    ! sigma. The residual tangential offset after a level is O(sigma), so
    ! the distance error is bounded by min(2 sigma, 2 sigma^2/d); the loop
    ! exits when that bound reaches tol. Cells far from the surface satisfy
    ! the bound immediately and skip the polish entirely.
    real(C_DOUBLE) function polish_distance(f, ibm, dns, x, p_in, sigma0, tol) result(d)
        procedure(body_indicator_i) :: f
        type(ibm_type), intent(in) :: ibm
        type(dns_type), intent(in) :: dns
        real(C_DOUBLE), intent(in) :: x(3), p_in(3), sigma0, tol

        logical :: s(-1:1, -1:1, -1:1)
        real(C_DOUBLE) :: p(3), pNew(3), q(3), xA(3), xB(3)
        real(C_DOUBLE) :: sigma, tb, dNew, dq
        integer :: lvl, a, b, c

        p = p_in
        d = sqrt(sum((x - p)**2))
        sigma = sigma0
        tb = max(0.25d0*tol, BISECT_FLOOR)

        do lvl = 1, MAX_POLISH_LEVELS
            if (min(2.0d0*sigma, 2.0d0*sigma*sigma/max(d, tiny(1.0d0))) <= tol) exit
            if (sigma <= BISECT_FLOOR) exit

            do c = -1, 1
                do b = -1, 1
                    do a = -1, 1
                        s(a,b,c) = f(p + sigma*[a, b, c], ibm, dns)
                    end do
                end do
            end do

            dNew = d
            pNew = p
            do c = -1, 1
                do b = -1, 1
                    do a = -1, 0
                        ! x-, y- and z-oriented lattice segments
                        if (s(a,b,c) .neqv. s(a+1,b,c)) then
                            xA = p + sigma*[a, b, c]
                            xB = p + sigma*[a+1, b, c]
                            call bisect_indicator(f, ibm, dns, xA, xB, tb, q)
                            dq = sqrt(sum((x - q)**2))
                            if (dq < dNew) then
                                dNew = dq
                                pNew = q
                            end if
                        end if
                        if (s(b,a,c) .neqv. s(b,a+1,c)) then
                            xA = p + sigma*[b, a, c]
                            xB = p + sigma*[b, a+1, c]
                            call bisect_indicator(f, ibm, dns, xA, xB, tb, q)
                            dq = sqrt(sum((x - q)**2))
                            if (dq < dNew) then
                                dNew = dq
                                pNew = q
                            end if
                        end if
                        if (s(b,c,a) .neqv. s(b,c,a+1)) then
                            xA = p + sigma*[b, c, a]
                            xB = p + sigma*[b, c, a+1]
                            call bisect_indicator(f, ibm, dns, xA, xB, tb, q)
                            dq = sqrt(sum((x - q)**2))
                            if (dq < dNew) then
                                dNew = dq
                                pNew = q
                            end if
                        end if
                    end do
                end do
            end do

            d = dNew
            p = pNew
            sigma = 0.5d0*sigma
        end do
    end function polish_distance

end module walldist
