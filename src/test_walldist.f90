! Gate driver for the geometry-agnostic wall distance (IDDES phase T1b,
! gate b — docs/next_session_iddes.md): drives the PRODUCTION walldist
! machinery with a second indicator, a sphere, whose exact distance is the
! closed form |dist(x, centre) - R|. The sphere straddles the periodic x
! boundary so the minimum-image / periodic-image query path is exercised
! too. Run: mpirun -n 1 build_cpu/walldist_test

module sphere_body
    use, intrinsic :: iso_c_binding
    use :: init, only: dns_type
    use :: ibmm, only: ibm_type
    implicit none

    ! The test geometry: module state stands in for what isInBody reads
    ! from ibm/dns in production.
    real(C_DOUBLE) :: sph_centre(3), sph_radius, sph_leng(3)
    logical :: sph_periodic(3)

contains

    ! Minimum-image distance to the sphere centre (the indicator must be
    ! length-periodic along periodic directions, like isInBody).
    real(C_DOUBLE) function sphere_centre_dist(xIN) result(d)
        real(C_DOUBLE), intent(in) :: xIN(1:3)

        real(C_DOUBLE) :: dx(3)
        integer :: i

        dx = xIN - sph_centre
        do i = 1, 3
            if (sph_periodic(i)) dx(i) = dx(i) - sph_leng(i)*nint(dx(i)/sph_leng(i))
        end do
        d = sqrt(sum(dx*dx))
    end function sphere_centre_dist

    ! Same signature as isInBody (walldist body_indicator_i).
    logical function sphere_in(xIN, ibm, dns)
        real(C_DOUBLE), intent(in) :: xIN(1:3)
        type(ibm_type), intent(in) :: ibm
        type(dns_type), intent(in) :: dns

        sphere_in = sphere_centre_dist(xIN) < sph_radius
    end function sphere_in

end module sphere_body

program test_walldist
    use, intrinsic :: iso_c_binding
    use :: init, only: dns_type
    use :: ibmm, only: ibm_type
    use :: walldist, only: walldist_type, build_walldist, destroy_walldist, &
        walldist_distance
    use :: sphere_body
    implicit none

    integer, parameter :: n(3) = [64, 48, 56]
    real(C_DOUBLE), parameter :: leng(3) = [1.0d0, 0.9d0, 1.1d0]
    real(C_DOUBLE), parameter :: tols(3) = [1.0d-4, 1.0d-7, 1.0d-10]

    type(dns_type) :: dns
    type(ibm_type) :: ibm
    type(walldist_type) :: w
    real(C_DOUBLE), allocatable :: lineX(:), lineY(:), lineZ(:)
    real(C_DOUBLE) :: x(3), d, ref, err, me, maxerr(size(tols))
    integer :: i, j, k, t
    logical :: ok

    ! Sphere across the periodic x boundary; interior in the walled y/z.
    sph_centre = [0.12d0, 0.46d0, 0.55d0]
    sph_radius = 0.27d0
    sph_leng = leng
    sph_periodic = [.true., .false., .false.]

    allocate(lineX(0:n(1)), lineY(0:n(2)), lineZ(0:n(3)))
    do i = 0, n(1)
        lineX(i) = leng(1)*real(i, C_DOUBLE)/real(n(1), C_DOUBLE)
    end do
    do i = 0, n(2)
        lineY(i) = leng(2)*real(i, C_DOUBLE)/real(n(2), C_DOUBLE)
    end do
    do i = 0, n(3)
        lineZ(i) = leng(3)*real(i, C_DOUBLE)/real(n(3), C_DOUBLE)
    end do

    call build_walldist(w, sphere_in, ibm, dns, lineX, lineY, lineZ, n, leng, &
        sph_periodic, .true.)
    if (w%nCloud == 0) error stop "sphere produced no surface cloud"

    ! Probe an off-lattice grid of points (incommensurate with the sampling
    ! lattice) at each tolerance; the error must track the tolerance down.
    do t = 1, size(tols)
        me = 0.0d0
        !$omp parallel do collapse(2) private(i,j,k,x,d,ref,err) reduction(max:me)
        do k = 1, 23
            do j = 1, 23
                do i = 1, 23
                    x = [leng(1)*(i - 0.382d0)/23.0d0, &
                         leng(2)*(j - 0.271d0)/23.0d0, &
                         leng(3)*(k - 0.618d0)/23.0d0]
                    d = walldist_distance(w, sphere_in, ibm, dns, x, tols(t))
                    ref = abs(sphere_centre_dist(x) - sph_radius)
                    err = abs(d - ref)
                    me = max(me, err)
                end do
            end do
        end do
        !$omp end parallel do
        maxerr(t) = me
        print '(A,ES9.2,A,ES12.5)', " tol = ", tols(t), "  max|d - ref| = ", maxerr(t)
    end do

    call destroy_walldist(w)

    ! Gates: the distance converges with the polish tolerance and reaches
    ! ~1e-10 at the tightest setting (bisection floors slightly above it).
    ok = maxerr(size(tols)) <= 1.0d-9
    do t = 2, size(tols)
        ok = ok .and. (maxerr(t) <= maxerr(t-1) + 1.0d-12)
    end do
    if (ok) then
        print *, "PASS"
    else
        print *, "FAIL"
        error stop 1
    end if

end program test_walldist
