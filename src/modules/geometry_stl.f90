!--------------------------!
!                          !
!   STL geometry           !
!   (moby_prepare, P1)     !
!                          !
!--------------------------!
!
! Prepare/solve split P1 (docs/prepare_solve_strategy.md): an inside/outside
! indicator for watertight binary-STL bodies with exactly the analytic
! isInBody signature (body_indicator_i, ibm.f90), so STL geometry flows
! through the SAME host machinery as analytic bodies -- classify_* masks,
! the graded coefficient bisection (set_ibm_coeff_host) and the walldist
! wall distance. Host-only preprocessing code; nothing here runs in the
! solver's time loop, and the loaded body lives in module state because the
! shared indicator signature has no slot for it (one body per run, like the
! analytic isInBody).
!
! Method (find-inside.cpl of the AMPHIBIOUS CPL solver is the template):
!  - triangles from one or more binary STL files, float32 vertices widened
!    to float64 exactly (the quantized coordinates ARE the as-built
!    geometry -- the mobygeom convention);
!  - a BVH over the triangles (median split on the centroid along the
!    largest extent, in-place index partition, small leaves);
!  - point classification by ray-parity casting along a fixed direction
!    table with majority vote (3 rays, escalating to 11 on disagreement).
!    Parity is orientation-independent, so no normal fixing is needed --
!    but the mesh must be watertight; a ray through a triangle edge or in
!    a triangle's plane is DEGENERATE and retried with a deterministically
!    perturbed direction (never a random one: the classification must be
!    reproducible across runs, ranks and thread counts);
!  - periodic directions by minimum-image testing: the point's -L/0/+L
!    images (pruned by the mesh bounding box) are each tested against the
!    as-stored triangles, so bodies straddling a periodic boundary
!    classify correctly.
!
! Everything is read-only after stl_geometry_load, so the indicator is
! safe to call from OpenMP host-parallel loops.

module geometry_stl
    use, intrinsic :: iso_c_binding
    use :: init, only: dns_type, VAR_P
    use :: ibmm, only: ibm_type
    use :: blocks, only: block_set_type
    implicit none

    private
    public :: stl_geometry_load, stl_geometry_destroy, stl_geometry_loaded
    public :: stl_is_in_body, stl_fill_dwall, stl_cull_box

    integer, parameter :: BVH_LEAF = 4        ! max triangles per leaf
    integer, parameter :: MAX_ATTEMPTS = 32   ! degenerate-ray retries
    integer, parameter :: NDIRS = 11          ! fixed ray-direction table
    integer, parameter :: STACK_SIZE = 128    ! BVH traversal (median split
                                              ! keeps the tree balanced)
    ! Ray-parallel cutoff, relative to the triangle's area scale (|det| is
    ! the (e1, e2, d) volume and grows with the triangle, so an absolute
    ! cutoff misbehaves on domain-sized slab faces).
    real(C_DOUBLE), parameter :: DET_EPS = 1.0d-10
    real(C_DOUBLE), parameter :: HIT_EPS = 1.0d-12   ! min ray parameter t
    ! Barycentric edge zone flagged degenerate (a crossing through a shared
    ! edge double-counts). Kept TIGHT: barycentric units scale with the
    ! triangle, and near-surface query points hit at tiny t, where the
    ! retry perturbation moves the hit point by only ~t*amplitude.
    real(C_DOUBLE), parameter :: EDGE_TOL = 1.0d-7

    ! The loaded body: v0 and the two edge vectors per triangle (the
    ! Moeller-Trumbore form), the BVH flat arrays, and the domain topology
    ! for the periodic images.
    integer :: nTri = 0
    real(C_DOUBLE), allocatable :: triV0(:,:), triE1(:,:), triE2(:,:)  ! (3, nTri)
    real(C_DOUBLE), allocatable :: triArea2(:)   ! |e1 x e2|, the det scale
    integer :: nNodes = 0
    integer, allocatable :: bvhLeft(:), bvhRight(:), bvhFirst(:), bvhCount(:)
    real(C_DOUBLE), allocatable :: bvhMin(:,:), bvhMax(:,:)             ! (3, nodes)
    integer, allocatable :: triIdx(:)
    real(C_DOUBLE) :: bbLo(3) = 0.0d0, bbHi(3) = 0.0d0, bbPad(3) = 0.0d0
    real(C_DOUBLE) :: leng(3) = 0.0d0
    logical :: isPeriodic(3) = .false.
    ! Distance queries image a periodic dim ONLY when the mesh is narrower
    ! than the cell there: an STL spanning the full cell (the padded wall
    ! slabs) is already its own periodic continuation, and its overhanging
    ! skin is INTERIOR to the periodic union of the solid -- imaging it
    ! would report the distance to a fictitious wall (mobygeom's no-image
    ! convention for exactly these files). Membership (parity, OR over
    ! images) is correct either way and always images.
    logical :: imageDim(3) = .false.
    real(C_DOUBLE) :: rayDir(3, NDIRS)

contains

    logical function stl_geometry_loaded()
        stl_geometry_loaded = nTri > 0
    end function stl_geometry_loaded

    ! Load one or more binary STL files ([ibm] stl_file, one path per
    ! occurrence) and build the BVH. The optional transform is
    ! v*scale + translate on the float64-widened vertices (mobygeom's
    ! convention and operation order, so transformed geometries agree
    ! bitwise). lengIn/periodicIn give the domain topology for
    ! minimum-image queries.
    subroutine stl_geometry_load(files, scale, translate, lengIn, periodicIn, &
            has_terminal)
        character(len=*), intent(in) :: files(:)
        real(C_DOUBLE), intent(in) :: scale, translate(3)
        real(C_DOUBLE), intent(in) :: lengIn(3)
        logical, intent(in) :: periodicIn(3)
        logical, intent(in) :: has_terminal

        integer :: f, n, total, offset, pos

        call stl_geometry_destroy()
        leng = lengIn
        isPeriodic = periodicIn

        total = 0
        do f = 1, size(files)
            call stl_count_triangles(trim(files(f)), n)
            total = total + n
        end do
        if (total == 0) error stop "stl_file: no triangles found"

        allocate(triV0(3, total), triE1(3, total), triE2(3, total), triArea2(total))

        offset = 0
        do f = 1, size(files)
            call stl_read_triangles(trim(files(f)), scale, translate, offset, &
                has_terminal)
        end do
        ! Drop exactly-degenerate triangles (zero area): they poison both
        ! the parity det scale and the distance query's edge divisions.
        nTri = 0
        do pos = 1, total
            n = pos   ! reuse: source slot
            triArea2(n) = norm2(cross3(triE1(:, n), triE2(:, n)))
            if (triArea2(n) > 1.0d-30) then
                nTri = nTri + 1
                triV0(:, nTri) = triV0(:, n)
                triE1(:, nTri) = triE1(:, n)
                triE2(:, nTri) = triE2(:, n)
                triArea2(nTri) = triArea2(n)
            end if
        end do
        if (nTri < total .and. has_terminal) &
            print '(A,I0,A)', " STL geometry: dropped ", total - nTri, " degenerate triangles"
        if (nTri == 0) error stop "stl_file: only degenerate triangles"

        call init_ray_directions()
        call build_bvh()

        bbLo = bvhMin(:, 1)   ! root node = the whole soup
        bbHi = bvhMax(:, 1)
        ! Points exactly on the bounding box must not be culled.
        bbPad = 1.0d-9*max(bbHi - bbLo, 1.0d0)
        imageDim = isPeriodic .and. (bbHi - bbLo < leng)

        if (has_terminal) then
            print '(A,I0,A)', " STL geometry: ", nTri, " triangles"
            print '(A,3(ES13.5),A,3(ES13.5))', "   bbox lo:", bbLo, "   hi:", bbHi
        end if
    end subroutine stl_geometry_load

    subroutine stl_geometry_destroy()
        if (allocated(triV0)) deallocate(triV0, triE1, triE2, triArea2)
        if (allocated(bvhLeft)) deallocate(bvhLeft, bvhRight, bvhFirst, bvhCount)
        if (allocated(bvhMin)) deallocate(bvhMin, bvhMax)
        if (allocated(triIdx)) deallocate(triIdx)
        nTri = 0
        nNodes = 0
    end subroutine stl_geometry_destroy

    ! Binary STL: 80-byte header, uint32 triangle count, then 50 bytes per
    ! triangle (float32 normal + 3 float32 vertices + uint16 attribute).
    ! ASCII STL ("solid" header, size not matching binary) is parsed by
    ! keyword; its decimal vertices go straight to float64 -- the same
    ! rounding trimesh/numpy give mobygeom, so the two see one geometry.
    logical function stl_is_ascii(file_name) result(is_ascii)
        character(len=*), intent(in) :: file_name

        integer :: unit, ios
        integer(C_INT32_T) :: count32
        integer(C_INT64_T) :: file_size
        character(len=5) :: head

        open(newunit=unit, file=trim(file_name), access="stream", &
            form="unformatted", status="old", action="read", iostat=ios)
        if (ios /= 0) then
            print *, "error: cannot open STL file: ", trim(file_name)
            error stop
        end if
        inquire(unit=unit, size=file_size)
        head = ""
        count32 = 0
        if (file_size >= 84_C_INT64_T) then
            read(unit, pos=1) head
            read(unit, pos=81) count32
        end if
        close(unit)

        if (file_size == 84_C_INT64_T + 50_C_INT64_T*int(count32, C_INT64_T)) then
            is_ascii = .false.
        else if (head == "solid") then
            is_ascii = .true.
        else
            print *, "error: corrupt binary STL (size mismatch): ", trim(file_name)
            error stop
        end if
    end function stl_is_ascii

    subroutine stl_count_triangles(file_name, n)
        character(len=*), intent(in) :: file_name
        integer, intent(out) :: n

        integer :: unit, ios
        integer(C_INT32_T) :: count32
        character(len=64) :: line

        if (stl_is_ascii(file_name)) then
            n = 0
            open(newunit=unit, file=trim(file_name), status="old", action="read")
            do
                read(unit, '(A)', iostat=ios) line
                if (ios /= 0) exit
                line = adjustl(line)
                if (line(1:6) == "facet ") n = n + 1
            end do
            close(unit)
        else
            open(newunit=unit, file=trim(file_name), access="stream", &
                form="unformatted", status="old", action="read")
            read(unit, pos=81) count32
            close(unit)
            n = int(count32)
        end if
    end subroutine stl_count_triangles

    subroutine stl_read_triangles(file_name, scale, translate, offset, has_terminal)
        character(len=*), intent(in) :: file_name
        real(C_DOUBLE), intent(in) :: scale, translate(3)
        integer, intent(inout) :: offset
        logical, intent(in) :: has_terminal

        integer :: unit, i, n, ios, nv
        integer(C_INT32_T) :: count32
        real(C_FLOAT) :: rec(12)
        integer(C_INT16_T) :: attr
        real(C_DOUBLE) :: v(3,3)
        character(len=256) :: line

        n = 0
        if (stl_is_ascii(file_name)) then
            open(newunit=unit, file=trim(file_name), status="old", action="read")
            nv = 0
            do
                read(unit, '(A)', iostat=ios) line
                if (ios /= 0) exit
                line = adjustl(line)
                if (line(1:7) /= "vertex ") cycle
                nv = nv + 1
                read(line(8:), *, iostat=ios) v(:, nv)
                if (ios /= 0) then
                    print *, "error: bad ASCII STL vertex in: ", trim(file_name)
                    error stop
                end if
                if (nv == 3) then
                    n = n + 1
                    call store_triangle(offset + n, v, scale, translate)
                    nv = 0
                end if
            end do
            close(unit)
            if (nv /= 0) error stop "ASCII STL: vertex count not a multiple of 3"
        else
            open(newunit=unit, file=trim(file_name), access="stream", &
                form="unformatted", status="old", action="read")
            read(unit, pos=81) count32
            n = int(count32)
            do i = 1, n
                read(unit) rec, attr
                ! rec(1:3) is the stored normal -- parity casting never
                ! needs it. float32 -> float64 widening is exact.
                v(:,1) = real(rec(4:6), C_DOUBLE)
                v(:,2) = real(rec(7:9), C_DOUBLE)
                v(:,3) = real(rec(10:12), C_DOUBLE)
                call store_triangle(offset + i, v, scale, translate)
            end do
            close(unit)
        end if
        if (has_terminal) print '(3A,I0,A)', " STL file ", trim(file_name), ": ", n, " triangles"
        offset = offset + n
    end subroutine stl_read_triangles

    ! The transform runs in float64 on the widened/parsed vertices --
    ! mobygeom's convention and operation order (v*scale + translate).
    subroutine store_triangle(t, v, scale, translate)
        integer, intent(in) :: t
        real(C_DOUBLE), intent(in) :: v(3,3), scale, translate(3)

        real(C_DOUBLE) :: w(3,3)

        w(:,1) = v(:,1)*scale + translate
        w(:,2) = v(:,2)*scale + translate
        w(:,3) = v(:,3)*scale + translate
        triV0(:, t) = w(:,1)
        triE1(:, t) = w(:,2) - w(:,1)
        triE2(:, t) = w(:,3) - w(:,1)
    end subroutine store_triangle

    ! Conservative solid-possible box for classification culling: a point
    ! outside it cannot be inside the body through ANY image the indicator
    ! tests (in imaged periodic dims the box widens by +-L, the hull of
    ! the image intervals).
    subroutine stl_cull_box(lo, hi)
        real(C_DOUBLE), intent(out) :: lo(3), hi(3)

        lo = bbLo - bbPad - merge(leng, [0.0d0, 0.0d0, 0.0d0], imageDim)
        hi = bbHi + bbPad + merge(leng, [0.0d0, 0.0d0, 0.0d0], imageDim)
    end subroutine stl_cull_box

    ! Fixed ray directions (find-inside.cpl's table), normalized. The
    ! first three decide unanimous points; disagreement escalates to the
    ! full table with majority vote.
    subroutine init_ray_directions()
        real(C_DOUBLE) :: d(3, NDIRS)
        integer :: i

        d(:, 1) = [ 0.12345d0,  0.54321d0,  0.98765d0]
        d(:, 2) = [-0.54321d0,  0.98765d0, -0.12345d0]
        d(:, 3) = [-0.98765d0, -0.12345d0,  0.54321d0]
        d(:, 4) = [ 0.31415d0, -0.92653d0,  0.58979d0]
        d(:, 5) = [-0.89793d0,  0.23846d0, -0.26433d0]
        d(:, 6) = [ 0.83279d0,  0.50288d0, -0.41971d0]
        d(:, 7) = [-0.69314d0, -0.71828d0,  0.31415d0]
        d(:, 8) = [ 0.14159d0,  0.26535d0, -0.89793d0]
        d(:, 9) = [-0.23846d0, -0.26433d0,  0.83279d0]
        d(:,10) = [ 0.50288d0, -0.41971d0,  0.69314d0]
        d(:,11) = [-0.71828d0,  0.31415d0, -0.14159d0]
        do i = 1, NDIRS
            rayDir(:, i) = d(:, i)/norm2(d(:, i))
        end do
    end subroutine init_ray_directions

    !========================
    ! BVH
    !========================

    subroutine build_bvh()
        integer :: maxNodes, t

        maxNodes = 4*nTri + 2
        allocate(bvhLeft(maxNodes), bvhRight(maxNodes))
        allocate(bvhFirst(maxNodes), bvhCount(maxNodes))
        allocate(bvhMin(3, maxNodes), bvhMax(3, maxNodes))
        allocate(triIdx(nTri))
        do t = 1, nTri
            triIdx(t) = t
        end do

        nNodes = 1
        bvhFirst(1) = 1
        bvhCount(1) = nTri
        bvhLeft(1) = 0
        bvhRight(1) = 0
        call update_node_bounds(1)
        call subdivide_node(1)
    end subroutine build_bvh

    subroutine update_node_bounds(node)
        integer, intent(in) :: node

        integer :: i, t

        bvhMin(:, node) = huge(1.0d0)
        bvhMax(:, node) = -huge(1.0d0)
        do i = bvhFirst(node), bvhFirst(node) + bvhCount(node) - 1
            t = triIdx(i)
            bvhMin(:, node) = min(bvhMin(:, node), triV0(:, t), &
                triV0(:, t) + triE1(:, t), triV0(:, t) + triE2(:, t))
            bvhMax(:, node) = max(bvhMax(:, node), triV0(:, t), &
                triV0(:, t) + triE1(:, t), triV0(:, t) + triE2(:, t))
        end do
    end subroutine update_node_bounds

    recursive subroutine subdivide_node(node)
        integer, intent(in) :: node

        integer :: axis, i, j, tmp, leftChild, rightChild, leftCount
        real(C_DOUBLE) :: cmin(3), cmax(3), extent(3), split, c

        if (bvhCount(node) <= BVH_LEAF) return

        ! Split on the largest centroid extent at its midpoint.
        cmin = huge(1.0d0)
        cmax = -huge(1.0d0)
        do i = bvhFirst(node), bvhFirst(node) + bvhCount(node) - 1
            associate (t => triIdx(i))
                cmin = min(cmin, triV0(:, t) + (triE1(:, t) + triE2(:, t))/3.0d0)
                cmax = max(cmax, triV0(:, t) + (triE1(:, t) + triE2(:, t))/3.0d0)
            end associate
        end do
        extent = cmax - cmin
        axis = maxloc(extent, dim=1)
        split = 0.5d0*(cmin(axis) + cmax(axis))

        i = bvhFirst(node)
        j = bvhFirst(node) + bvhCount(node) - 1
        do while (i <= j)
            c = triV0(axis, triIdx(i)) &
                + (triE1(axis, triIdx(i)) + triE2(axis, triIdx(i)))/3.0d0
            if (c < split) then
                i = i + 1
            else
                tmp = triIdx(i)
                triIdx(i) = triIdx(j)
                triIdx(j) = tmp
                j = j - 1
            end if
        end do
        leftCount = i - bvhFirst(node)
        ! Degenerate split (all centroids on one side): halve by count.
        if (leftCount == 0 .or. leftCount == bvhCount(node)) &
            leftCount = bvhCount(node)/2

        leftChild = nNodes + 1
        rightChild = nNodes + 2
        nNodes = nNodes + 2
        bvhFirst(leftChild) = bvhFirst(node)
        bvhCount(leftChild) = leftCount
        bvhLeft(leftChild) = 0
        bvhRight(leftChild) = 0
        bvhFirst(rightChild) = bvhFirst(node) + leftCount
        bvhCount(rightChild) = bvhCount(node) - leftCount
        bvhLeft(rightChild) = 0
        bvhRight(rightChild) = 0
        bvhLeft(node) = leftChild
        bvhRight(node) = rightChild
        bvhCount(node) = 0

        call update_node_bounds(leftChild)
        call update_node_bounds(rightChild)
        call subdivide_node(leftChild)
        call subdivide_node(rightChild)
    end subroutine subdivide_node

    !========================
    ! Ray-parity inside test
    !========================

    ! The shared indicator signature: ibm and dns ride along unused (the
    ! geometry is module state).
    logical function stl_is_in_body(xIN, ibm, dns) result(inside)
        real(C_DOUBLE), intent(in) :: xIN(1:3)
        type(ibm_type), intent(in) :: ibm
        type(dns_type), intent(in) :: dns

        integer :: kx, ky, kz, k(3)
        real(C_DOUBLE) :: xs(3)

        inside = .false.
        ! Minimum-image sweep: a body straddling a periodic boundary is
        ! seen by the neighbouring image of the query point. The bbox cull
        ! keeps this at one parity test for almost every point.
        do kz = -merge(1, 0, isPeriodic(3)), merge(1, 0, isPeriodic(3))
            do ky = -merge(1, 0, isPeriodic(2)), merge(1, 0, isPeriodic(2))
                do kx = -merge(1, 0, isPeriodic(1)), merge(1, 0, isPeriodic(1))
                    k = [kx, ky, kz]
                    xs = xIN + real(k, C_DOUBLE)*leng
                    if (any(xs < bbLo - bbPad) .or. any(xs > bbHi + bbPad)) cycle
                    if (parity_inside(xs)) then
                        inside = .true.
                        return
                    end if
                end do
            end do
        end do
    end function stl_is_in_body

    ! Majority vote over the direction table; each degenerate cast retries
    ! with a deterministic perturbation.
    logical function parity_inside(x) result(inside)
        real(C_DOUBLE), intent(in) :: x(3)

        integer :: i, attempts, vote, votesIn, votesOut
        real(C_DOUBLE) :: d(3)

        votesIn = 0
        votesOut = 0
        do i = 1, NDIRS
            d = rayDir(:, i)
            attempts = 0
            vote = cast_ray(x, d)
            do while (vote < 0)
                attempts = attempts + 1
                if (attempts > MAX_ATTEMPTS) then
                    print *, "error: STL ray casting degenerate at point", x, &
                        " -- is the mesh watertight? (mobygeom stress-stl-watertightness)"
                    error stop
                end if
                ! LARGE deterministic rotation: parity is direction-
                ! independent, and a near-surface origin hits at tiny t
                ! where a small perturbation cannot leave the edge zone.
                d = rayDir(:, i) + 0.37d0*real(attempts, C_DOUBLE) &
                    *rayDir(:, 1 + mod(i + attempts - 1, NDIRS))
                d = d/norm2(d)
                vote = cast_ray(x, d)
            end do
            if (vote == 1) then
                votesIn = votesIn + 1
            else
                votesOut = votesOut + 1
            end if
            ! Unanimous after three rays: decided.
            if (i == 3 .and. (votesIn == 3 .or. votesOut == 3)) exit
        end do
        inside = votesIn > votesOut
    end function parity_inside

    ! One parity cast: 1 = odd crossings (inside), 0 = even, -1 = the ray
    ! grazed an edge or a near-parallel triangle (degenerate, retry).
    integer function cast_ray(o, d) result(vote)
        real(C_DOUBLE), intent(in) :: o(3), d(3)

        integer :: stack(STACK_SIZE), sp, node, i, hit, crossings

        crossings = 0
        sp = 1
        stack(1) = 1
        do while (sp > 0)
            node = stack(sp)
            sp = sp - 1
            if (.not. ray_hits_box(o, d, node)) cycle
            if (bvhCount(node) > 0) then
                do i = bvhFirst(node), bvhFirst(node) + bvhCount(node) - 1
                    hit = ray_triangle(o, d, triIdx(i))
                    if (hit < 0) then
                        vote = -1
                        return
                    end if
                    crossings = crossings + hit
                end do
            else
                if (sp + 2 > STACK_SIZE) error stop "STL BVH traversal stack overflow"
                stack(sp + 1) = bvhLeft(node)
                stack(sp + 2) = bvhRight(node)
                sp = sp + 2
            end if
        end do
        vote = mod(crossings, 2)
    end function cast_ray

    logical function ray_hits_box(o, d, node) result(hits)
        real(C_DOUBLE), intent(in) :: o(3), d(3)
        integer, intent(in) :: node

        real(C_DOUBLE) :: invD, t1, t2, tmin, tmax
        integer :: a

        tmin = -huge(1.0d0)
        tmax = huge(1.0d0)
        do a = 1, 3
            invD = 1.0d0/sign(max(abs(d(a)), 1.0d-300), d(a))
            t1 = (bvhMin(a, node) - o(a))*invD
            t2 = (bvhMax(a, node) - o(a))*invD
            tmin = max(tmin, min(t1, t2))
            tmax = min(tmax, max(t1, t2))
        end do
        hits = tmax >= tmin .and. tmax >= 0.0d0
    end function ray_hits_box

    ! Moeller-Trumbore. 1 = clean crossing (t > 0), 0 = miss or behind,
    ! -1 = degenerate (near-parallel or a hit in the EDGE_TOL edge zone).
    integer function ray_triangle(o, d, t) result(hit)
        real(C_DOUBLE), intent(in) :: o(3), d(3)
        integer, intent(in) :: t

        real(C_DOUBLE) :: pvec(3), qvec(3), tvec(3), det, invDet, u, v, tt

        pvec = cross3(d, triE2(:, t))
        det = dot_product(triE1(:, t), pvec)
        if (abs(det) < DET_EPS*triArea2(t)) then
            hit = -1
            return
        end if
        invDet = 1.0d0/det
        tvec = o - triV0(:, t)
        u = dot_product(tvec, pvec)*invDet
        if (u < 0.0d0 .or. u > 1.0d0) then
            hit = 0
            return
        end if
        qvec = cross3(tvec, triE1(:, t))
        v = dot_product(d, qvec)*invDet
        if (v < 0.0d0 .or. u + v > 1.0d0) then
            hit = 0
            return
        end if
        tt = dot_product(triE2(:, t), qvec)*invDet
        if (tt <= HIT_EPS) then
            hit = 0
            return
        end if
        if (u < EDGE_TOL .or. u > 1.0d0 - EDGE_TOL .or. &
            v < EDGE_TOL .or. v > 1.0d0 - EDGE_TOL .or. &
            u + v > 1.0d0 - EDGE_TOL) then
            hit = -1
            return
        end if
        hit = 1
    end function ray_triangle

    pure function cross3(a, b) result(c)
        real(C_DOUBLE), intent(in) :: a(3), b(3)
        real(C_DOUBLE) :: c(3)

        c(1) = a(2)*b(3) - a(3)*b(2)
        c(2) = a(3)*b(1) - a(1)*b(3)
        c(3) = a(1)*b(2) - a(2)*b(1)
    end function cross3

    !========================
    ! Exact wall distance
    !========================

    ! Per-leaf cell-centred dwall tiles: the exact (BVH nearest
    ! point-triangle) unsigned distance to the STL surface -- the same
    ! query mobygeom's dwall_blocks uses (igl), so the two agree to
    ! round-off. The indicator-driven walldist machinery is NOT used for
    ! STL bodies: its surface-cloud + polish bisections cost millions of
    ! near-surface parity casts, and the exact query is both faster and
    ! more accurate here.
    subroutine stl_fill_dwall(dwall, blk)
        real(C_DOUBLE), intent(inout) :: dwall(0:,0:,0:,1:)
        type(block_set_type), intent(in) :: blk

        integer :: i, j, k, b, nx, ny, nz
        real(C_DOUBLE) :: xA(3)

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))

        !$omp parallel do collapse(2) private(i,j,k,b,xA)
        do b = 1, int(blk%nBlocks)
        do k = 0, nz+1
            do j = 0, ny+1
                do i = 0, nx+1
                    xA(1) = blk%x(i, VAR_P, b)
                    xA(2) = blk%y(j, VAR_P, b)
                    xA(3) = blk%z(k, VAR_P, b)
                    dwall(i,j,k,b) = stl_distance(xA)
                end do
            end do
        end do
        end do
        !$omp end parallel do
    end subroutine stl_fill_dwall

    ! Unsigned distance to the surface, minimum over the periodic images
    ! of the query point (pruned by the running best against each image's
    ! root bounding box).
    real(C_DOUBLE) function stl_distance(xIN) result(d)
        real(C_DOUBLE), intent(in) :: xIN(3)

        integer :: kx, ky, kz
        real(C_DOUBLE) :: xs(3), best2

        best2 = huge(1.0d0)
        do kz = -merge(1, 0, imageDim(3)), merge(1, 0, imageDim(3))
            do ky = -merge(1, 0, imageDim(2)), merge(1, 0, imageDim(2))
                do kx = -merge(1, 0, imageDim(1)), merge(1, 0, imageDim(1))
                    xs = xIN + real([kx, ky, kz], C_DOUBLE)*leng
                    if (box_dist2(xs, 1) >= best2) cycle
                    call nearest_tri_dist2(xs, best2)
                end do
            end do
        end do
        d = sqrt(best2)
    end function stl_distance

    real(C_DOUBLE) function box_dist2(x, node) result(d2)
        real(C_DOUBLE), intent(in) :: x(3)
        integer, intent(in) :: node

        real(C_DOUBLE) :: dd(3)

        dd = max(bvhMin(:, node) - x, x - bvhMax(:, node), 0.0d0)
        d2 = dot_product(dd, dd)
    end function box_dist2

    subroutine nearest_tri_dist2(x, best2)
        real(C_DOUBLE), intent(in) :: x(3)
        real(C_DOUBLE), intent(inout) :: best2

        integer :: stack(STACK_SIZE), sp, node, i
        real(C_DOUBLE) :: d2

        sp = 1
        stack(1) = 1
        do while (sp > 0)
            node = stack(sp)
            sp = sp - 1
            if (box_dist2(x, node) >= best2) cycle
            if (bvhCount(node) > 0) then
                do i = bvhFirst(node), bvhFirst(node) + bvhCount(node) - 1
                    d2 = point_tri_dist2(x, triIdx(i))
                    if (d2 < best2) best2 = d2
                end do
            else
                if (sp + 2 > STACK_SIZE) error stop "STL BVH traversal stack overflow"
                ! Visit the nearer child first for tighter pruning.
                if (box_dist2(x, bvhLeft(node)) <= box_dist2(x, bvhRight(node))) then
                    stack(sp + 1) = bvhRight(node)
                    stack(sp + 2) = bvhLeft(node)
                else
                    stack(sp + 1) = bvhLeft(node)
                    stack(sp + 2) = bvhRight(node)
                end if
                sp = sp + 2
            end if
        end do
    end subroutine nearest_tri_dist2

    ! Squared distance from a point to triangle t (Eberly's region
    ! decomposition on the parametrisation v0 + s e1 + t e2).
    real(C_DOUBLE) function point_tri_dist2(p, t) result(dist2)
        real(C_DOUBLE), intent(in) :: p(3)
        integer, intent(in) :: t

        real(C_DOUBLE) :: dv(3), a, b, c, d, e, det, s, tt, inv, tmp0, tmp1, numer, denom

        dv = triV0(:, t) - p
        a = dot_product(triE1(:, t), triE1(:, t))
        b = dot_product(triE1(:, t), triE2(:, t))
        c = dot_product(triE2(:, t), triE2(:, t))
        d = dot_product(triE1(:, t), dv)
        e = dot_product(triE2(:, t), dv)
        det = max(a*c - b*b, 0.0d0)
        s = b*e - c*d
        tt = b*d - a*e

        if (s + tt <= det) then
            if (s < 0.0d0) then
                if (tt < 0.0d0) then                        ! region 4
                    if (d < 0.0d0) then
                        tt = 0.0d0
                        s = min(max(-d/a, 0.0d0), 1.0d0)
                    else
                        s = 0.0d0
                        tt = min(max(-e/c, 0.0d0), 1.0d0)
                    end if
                else                                        ! region 3
                    s = 0.0d0
                    tt = min(max(-e/c, 0.0d0), 1.0d0)
                end if
            else if (tt < 0.0d0) then                       ! region 5
                tt = 0.0d0
                s = min(max(-d/a, 0.0d0), 1.0d0)
            else                                            ! region 0
                inv = 1.0d0/max(det, tiny(1.0d0))
                s = s*inv
                tt = tt*inv
            end if
        else
            if (s < 0.0d0) then                             ! region 2
                tmp0 = b + d
                tmp1 = c + e
                if (tmp1 > tmp0) then
                    numer = tmp1 - tmp0
                    denom = a - 2.0d0*b + c
                    s = min(max(numer/max(denom, tiny(1.0d0)), 0.0d0), 1.0d0)
                    tt = 1.0d0 - s
                else
                    s = 0.0d0
                    tt = min(max(-e/c, 0.0d0), 1.0d0)
                end if
            else if (tt < 0.0d0) then                       ! region 6
                tmp0 = b + e
                tmp1 = a + d
                if (tmp1 > tmp0) then
                    numer = tmp1 - tmp0
                    denom = a - 2.0d0*b + c
                    tt = min(max(numer/max(denom, tiny(1.0d0)), 0.0d0), 1.0d0)
                    s = 1.0d0 - tt
                else
                    tt = 0.0d0
                    s = min(max(-d/a, 0.0d0), 1.0d0)
                end if
            else                                            ! region 1
                numer = (c + e) - (b + d)
                denom = a - 2.0d0*b + c
                s = min(max(numer/max(denom, tiny(1.0d0)), 0.0d0), 1.0d0)
                tt = 1.0d0 - s
            end if
        end if

        dist2 = max(s*(a*s + b*tt + 2.0d0*d) + tt*(b*s + c*tt + 2.0d0*e) &
            + dot_product(dv, dv), 0.0d0)
    end function point_tri_dist2

end module geometry_stl
