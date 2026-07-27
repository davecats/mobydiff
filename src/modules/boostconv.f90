!--------------------------!
!                          !
!    BoostConv module      !
!                          !
!--------------------------!
!
! Residual-recombination accelerator for steady fixed-point iterations
! (Citro & Palitta, "Residual Recombination Methods as Anderson-like
! Acceleration" -- the ROBUST Algorithm 3.1: skinny-QR residual basis
! with Givens downdating and a linear-dependency guard tau; see
! docs/next_session_boostconv.md and tutorials/naca/boostconv.pdf).
!
! The module is physics-blind: it operates on a PACKED device-resident
! state vector x(nDof) laid out as nVar equal blocks (the caller packs).
! One ACTIVE iteration = the caller's p solver steps: the caller calls
! boostconv_apply with the current packed state; the module forms the
! residual r_k = x - x_prev, replaces it by the boosted increment
! xi_k = r_k + W c_k with c_k = argmin ||r_k - V_k c|| over the last
! N stored residual differences (V kept as Q R), and overwrites x with
! x_prev + xi_k. All O(nDof) work runs as flat OpenMP-target loops on
! the device basis; the N x N algebra lives on the host. Inner products
! carry PER-VARIABLE weights 1/s_var^2 (dynamic RMS scaling, set by the
! caller each activation); if any scale drifts > 20 % from the value the
! current basis was built with, the basis is FLUSHED and restarted with
! the fresh scales (an incrementally-updated QR cannot survive a metric
! change). A divergence guard flushes the basis if the weighted residual
! norm grows over 3 consecutive activations.
module boostconv
    use, intrinsic :: iso_c_binding
    use :: comm, only: comm_type, comm_allreduce_sum
    implicit none
    private

    public :: boostconv_type, boostconv_init, boostconv_destroy
    public :: enter_boostconv_data, exit_boostconv_data
    public :: boostconv_apply

    type :: boostconv_type
        logical(C_BOOL) :: enabled = .false.
        integer(C_INT) :: capacity = 10      ! N: stored directions
        integer(C_INT) :: interval = 25      ! p: solver steps per activation
        real(C_DOUBLE) :: tau = 1.0d-3       ! rank guard
        integer(C_INT) :: nDof = 0
        integer(C_INT) :: nVar = 0
        integer(C_INT) :: varLen = 0         ! nDof/nVar
        integer(C_INT) :: nCols = 0          ! current basis size (0..capacity)
        integer(C_INT) :: nActive = 0        ! activation counter
        integer(C_INT) :: nGrow = 0          ! consecutive residual growths
        real(C_DOUBLE) :: rnorm_prev = -1.0d0
        ! device state (flat, nDof) and basis (nDof, capacity)
        real(C_DOUBLE), allocatable :: xwork(:)         ! caller pack target
        real(C_DOUBLE), allocatable :: xprev(:), rprev(:), xiprev(:)
        real(C_DOUBLE), allocatable :: rcur(:), xicur(:)
        real(C_DOUBLE), allocatable :: W(:,:), Q(:,:)
        ! host small algebra + weights
        real(C_DOUBLE), allocatable :: Rmat(:,:)         ! (capacity, capacity)
        real(C_DOUBLE), allocatable :: wvar(:)           ! (nVar) 1/s^2
        real(C_DOUBLE), allocatable :: svar0(:)          ! scales of the basis
    end type boostconv_type

contains

    subroutine boostconv_init(bc, nDof, nVar, capacity, interval, tau)
        type(boostconv_type), intent(inout) :: bc
        integer(C_INT), intent(in) :: nDof, nVar, capacity, interval
        real(C_DOUBLE), intent(in) :: tau

        bc%nDof = nDof
        bc%nVar = nVar
        bc%varLen = nDof/nVar
        bc%capacity = capacity
        bc%interval = interval
        bc%tau = tau
        allocate(bc%xwork(nDof))
        allocate(bc%xprev(nDof), bc%rprev(nDof), bc%xiprev(nDof))
        allocate(bc%rcur(nDof), bc%xicur(nDof))
        allocate(bc%W(nDof, capacity), bc%Q(nDof, capacity))
        allocate(bc%Rmat(capacity, capacity))
        allocate(bc%wvar(nVar), bc%svar0(nVar))
        bc%xwork = 0.0d0
        bc%xprev = 0.0d0; bc%rprev = 0.0d0; bc%xiprev = 0.0d0
        bc%rcur = 0.0d0; bc%xicur = 0.0d0
        bc%W = 0.0d0; bc%Q = 0.0d0; bc%Rmat = 0.0d0
        bc%wvar = 1.0d0; bc%svar0 = -1.0d0
        bc%enabled = .true.
    end subroutine boostconv_init

    subroutine boostconv_destroy(bc)
        type(boostconv_type), intent(inout) :: bc
        if (allocated(bc%W)) deallocate(bc%W, bc%Q, bc%xwork, bc%xprev, bc%rprev, &
            bc%xiprev, bc%rcur, bc%xicur, bc%Rmat, bc%wvar, bc%svar0)
        bc%enabled = .false.
    end subroutine boostconv_destroy

    subroutine enter_boostconv_data(bc)
        type(boostconv_type), intent(inout) :: bc
        if (.not. bc%enabled) return
        !$omp target enter data map(to: bc)
        !$omp target enter data map(to: bc%xwork, bc%xprev, bc%rprev, bc%xiprev, &
        !$omp& bc%rcur, bc%xicur, bc%W, bc%Q)
    end subroutine enter_boostconv_data

    subroutine exit_boostconv_data(bc)
        type(boostconv_type), intent(inout) :: bc
        if (.not. bc%enabled) return
        !$omp target exit data map(delete: bc%xwork, bc%xprev, bc%rprev, bc%xiprev, &
        !$omp& bc%rcur, bc%xicur, bc%W, bc%Q)
        !$omp target exit data map(delete: bc)
    end subroutine exit_boostconv_data

    ! Weighted global dot product of two device vectors: per-variable
    ! weights wvar (host array, small) applied blockwise; MPI-summed.
    function bdot(bc, a, b, c) result(s)
        type(boostconv_type), intent(in) :: bc
        real(C_DOUBLE), intent(in) :: a(:), b(:)
        type(comm_type), intent(in) :: c
        real(C_DOUBLE) :: s, sv(1), part
        integer :: v, i, lo, hi

        s = 0.0d0
        do v = 1, int(bc%nVar)
            lo = (v-1)*int(bc%varLen) + 1
            hi = v*int(bc%varLen)
            part = 0.0d0
            !$omp target teams distribute parallel do reduction(+:part) &
            !$omp& map(to: lo, hi) map(tofrom: part) private(i)
            do i = lo, hi
                part = part + a(i)*b(i)
            end do
            !$omp end target teams distribute parallel do
            s = s + bc%wvar(v)*part
        end do
        sv(1) = s
        call comm_allreduce_sum(c, sv)
        s = sv(1)
    end function bdot

    ! One residual-recombination activation. x is the packed CURRENT
    ! state (device); svar the caller's per-variable scales (host).
    ! On return x holds the boosted state (unchanged on the first two
    ! activations, which only prime the history).
    subroutine boostconv_apply(bc, x, svar, c, has_terminal)
        type(boostconv_type), intent(inout) :: bc
        real(C_DOUBLE), intent(inout) :: x(:)
        real(C_DOUBLE), intent(in) :: svar(:)
        type(comm_type), intent(in) :: c
        logical, intent(in) :: has_terminal

        integer :: i, j, m, n
        real(C_DOUBLE) :: rnorm, vnorm, qn, gc, gs, t1, t2
        real(C_DOUBLE), allocatable :: cvec(:), qtv(:)

        n = int(bc%nDof)

        ! dynamic-at-initialization scaling: the metric (per-variable
        ! weights) is captured when the basis is (re)started and FROZEN
        ! while it lives — an incrementally-updated QR cannot survive a
        ! metric change, and flushing on drift (the first attempt)
        ! thrashed the basis into uselessness during transients (V1).
        ! The growth guard below refreshes the metric when it flushes.
        if (bc%nCols == 0) bc%svar0 = svar
        do j = 1, int(bc%nVar)
            bc%wvar(j) = 1.0d0/max(bc%svar0(j)**2, 1.0d-30)
        end do

        if (bc%nActive == 0) then
            call dcopy_dev(x, bc%xprev, n)
            bc%nActive = 1
            return
        end if

        ! residual over the last interval: r = x - xprev
        !$omp target teams distribute parallel do private(i)
        do i = 1, n
            bc%rcur(i) = x(i) - bc%xprev(i)
        end do
        !$omp end target teams distribute parallel do

        rnorm = sqrt(max(bdot(bc, bc%rcur, bc%rcur, c), 0.0d0))
        if (bc%rnorm_prev > 0.0d0 .and. rnorm > bc%rnorm_prev) then
            bc%nGrow = bc%nGrow + 1
            if (bc%nGrow >= 3) then
                if (has_terminal) print '(a)', &
                    " [boostconv] residual grew 3 activations -- basis flushed"
                bc%nCols = 0
                bc%nGrow = 0
            end if
        else
            bc%nGrow = 0
        end if
        bc%rnorm_prev = rnorm

        if (bc%nActive == 1) then
            ! prime: xi_0 = r_0, no recombination possible yet
            call dcopy_dev(bc%rcur, bc%rprev, n)
            call dcopy_dev(bc%rcur, bc%xiprev, n)
            call dcopy_dev(x, bc%xprev, n)
            bc%nActive = 2
            return
        end if

        ! ---- new direction pair (paper eq. 2.4) ----
        ! v_new = r_{k-1} - r_k ; w_new = xi_{k-1} + r_k - r_{k-1}
        if (bc%nCols == int(bc%capacity)) then
            ! downdate: drop the oldest column. W shifts left; Q R gets
            ! the Givens treatment (R[:,2:] is upper Hessenberg).
            m = int(bc%capacity)
            do j = 2, m
                call ccopy_dev(bc%W, j, j-1, n)
            end do
            do j = 1, m-1
                bc%Rmat(:, j) = bc%Rmat(:, j+1)
            end do
            do j = 1, m-1
                t1 = bc%Rmat(j, j); t2 = bc%Rmat(j+1, j)
                qn = sqrt(t1*t1 + t2*t2)
                if (qn <= 0.0d0) cycle
                gc = t1/qn; gs = t2/qn
                do i = j, m-1
                    t1 = bc%Rmat(j, i); t2 = bc%Rmat(j+1, i)
                    bc%Rmat(j, i)   =  gc*t1 + gs*t2
                    bc%Rmat(j+1, i) = -gs*t1 + gc*t2
                end do
                call crot_dev(bc%Q, j, j+1, gc, gs, n)
            end do
            bc%nCols = m - 1
        end if

        m = int(bc%nCols)
        ! Gram-Schmidt the new v against Q(:,1:m) into rcur-scratch
        allocate(qtv(max(m,1)), cvec(max(m+1,1)))
        !$omp target teams distribute parallel do private(i)
        do i = 1, n
            bc%xicur(i) = bc%rprev(i) - bc%rcur(i)      ! v_new (scratch)
        end do
        !$omp end target teams distribute parallel do
        vnorm = sqrt(max(bdot(bc, bc%xicur, bc%xicur, c), 0.0d0))
        do j = 1, m
            qtv(j) = coldot(bc, bc%Q, j, bc%xicur, c)
        end do
        do j = 1, m
            call colaxpy_dev(bc%xicur, bc%Q, j, -qtv(j), n)
        end do
        qn = sqrt(max(bdot(bc, bc%xicur, bc%xicur, c), 0.0d0))

        if (qn >= bc%tau*max(vnorm, 1.0d-300) .and. m < int(bc%capacity)) then
            ! append: W(:,m+1) = w_new, Q(:,m+1) = q/||q||, R column
            call setw_dev(bc%W, m+1, bc%xiprev, bc%rcur, bc%rprev, n)
            call scalecol_dev(bc%xicur, bc%Q, m+1, 1.0d0/qn, n)
            bc%Rmat(1:m, m+1) = qtv(1:m)
            bc%Rmat(m+1, m+1) = qn
            bc%nCols = m + 1
        else if (has_terminal .and. qn < bc%tau*max(vnorm, 1.0d-300)) then
            print '(a)', " [boostconv] near-dependent direction discarded"
        end if

        ! ---- boosted increment: xi = r + W R^{-1} Q^T r ----
        m = int(bc%nCols)
        call dcopy_dev(bc%rcur, bc%xicur, n)
        if (m > 0) then
            do j = 1, m
                cvec(j) = coldot(bc, bc%Q, j, bc%rcur, c)
            end do
            ! back-substitute R c = Q^T r
            do j = m, 1, -1
                do i = j+1, m
                    cvec(j) = cvec(j) - bc%Rmat(j, i)*cvec(i)
                end do
                cvec(j) = cvec(j)/bc%Rmat(j, j)
            end do
            do j = 1, m
                call colaxpy_dev(bc%xicur, bc%W, j, cvec(j), n)
            end do
        end if

        ! SAFEGUARD (standard Anderson damping): a near-degenerate basis
        ! can produce wild extrapolations; cap ||xi|| at 10 ||r|| by
        ! rescaling the recombination part (xi -> r + s (xi - r)).
        t1 = sqrt(max(bdot(bc, bc%xicur, bc%xicur, c), 0.0d0))
        if (t1 > 10.0d0*rnorm .and. t1 > 0.0d0) then
            t2 = 10.0d0*rnorm/t1
            !$omp target teams distribute parallel do map(to: t2) private(i)
            do i = 1, n
                bc%xicur(i) = bc%rcur(i) + t2*(bc%xicur(i) - bc%rcur(i))
            end do
            !$omp end target teams distribute parallel do
            if (has_terminal) print '(a,es9.2)', &
                " [boostconv] correction capped, |xi|/|r| was ", t1/max(rnorm, 1.0d-300)
        end if

        ! overwrite the state: x = xprev + xi; roll the history
        !$omp target teams distribute parallel do private(i)
        do i = 1, n
            x(i) = bc%xprev(i) + bc%xicur(i)
        end do
        !$omp end target teams distribute parallel do
        call dcopy_dev(bc%rcur, bc%rprev, n)
        call dcopy_dev(bc%xicur, bc%xiprev, n)
        call dcopy_dev(x, bc%xprev, n)
        bc%nActive = bc%nActive + 1
        deallocate(qtv, cvec)
    end subroutine boostconv_apply

    ! ---- small device helpers (flat loops; b arrays are mapped) ----

    subroutine dcopy_dev(src, dst, n)
        real(C_DOUBLE), intent(in) :: src(:)
        real(C_DOUBLE), intent(inout) :: dst(:)
        integer, intent(in) :: n
        integer :: i
        !$omp target teams distribute parallel do private(i)
        do i = 1, n
            dst(i) = src(i)
        end do
        !$omp end target teams distribute parallel do
    end subroutine dcopy_dev

    subroutine ccopy_dev(A, jsrc, jdst, n)
        real(C_DOUBLE), intent(inout) :: A(:,:)
        integer, intent(in) :: jsrc, jdst, n
        integer :: i
        !$omp target teams distribute parallel do private(i)
        do i = 1, n
            A(i, jdst) = A(i, jsrc)
        end do
        !$omp end target teams distribute parallel do
    end subroutine ccopy_dev

    subroutine crot_dev(A, j1, j2, gc, gs, n)
        real(C_DOUBLE), intent(inout) :: A(:,:)
        integer, intent(in) :: j1, j2, n
        real(C_DOUBLE), intent(in) :: gc, gs
        integer :: i
        real(C_DOUBLE) :: t
        !$omp target teams distribute parallel do private(i, t)
        do i = 1, n
            t = A(i, j1)
            A(i, j1) =  gc*t + gs*A(i, j2)
            A(i, j2) = -gs*t + gc*A(i, j2)
        end do
        !$omp end target teams distribute parallel do
    end subroutine crot_dev

    function coldot(bc, A, j, v, c) result(s)
        type(boostconv_type), intent(in) :: bc
        real(C_DOUBLE), intent(in) :: A(:,:), v(:)
        integer, intent(in) :: j
        type(comm_type), intent(in) :: c
        real(C_DOUBLE) :: s, sv(1), part
        integer :: vv, i, lo, hi
        s = 0.0d0
        do vv = 1, int(bc%nVar)
            lo = (vv-1)*int(bc%varLen) + 1
            hi = vv*int(bc%varLen)
            part = 0.0d0
            !$omp target teams distribute parallel do reduction(+:part) &
            !$omp& map(to: lo, hi, j) map(tofrom: part) private(i)
            do i = lo, hi
                part = part + A(i, j)*v(i)
            end do
            !$omp end target teams distribute parallel do
            s = s + bc%wvar(vv)*part
        end do
        sv(1) = s
        call comm_allreduce_sum(c, sv)
        s = sv(1)
    end function coldot

    subroutine colaxpy_dev(v, A, j, alpha, n)
        real(C_DOUBLE), intent(inout) :: v(:)
        real(C_DOUBLE), intent(in) :: A(:,:)
        integer, intent(in) :: j, n
        real(C_DOUBLE), intent(in) :: alpha
        integer :: i
        !$omp target teams distribute parallel do map(to: alpha, j) private(i)
        do i = 1, n
            v(i) = v(i) + alpha*A(i, j)
        end do
        !$omp end target teams distribute parallel do
    end subroutine colaxpy_dev

    subroutine scalecol_dev(src, A, j, alpha, n)
        real(C_DOUBLE), intent(in) :: src(:)
        real(C_DOUBLE), intent(inout) :: A(:,:)
        integer, intent(in) :: j, n
        real(C_DOUBLE), intent(in) :: alpha
        integer :: i
        !$omp target teams distribute parallel do map(to: alpha, j) private(i)
        do i = 1, n
            A(i, j) = alpha*src(i)
        end do
        !$omp end target teams distribute parallel do
    end subroutine scalecol_dev

    ! W(:,j) = xiprev + rcur - rprev   (paper eq. 2.4, new W column)
    subroutine setw_dev(W, j, xiprev, rcur, rprev, n)
        real(C_DOUBLE), intent(inout) :: W(:,:)
        real(C_DOUBLE), intent(in) :: xiprev(:), rcur(:), rprev(:)
        integer, intent(in) :: j, n
        integer :: i
        !$omp target teams distribute parallel do map(to: j) private(i)
        do i = 1, n
            W(i, j) = xiprev(i) + rcur(i) - rprev(i)
        end do
        !$omp end target teams distribute parallel do
    end subroutine setw_dev

end module boostconv
