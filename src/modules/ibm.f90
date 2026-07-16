!--------------------------!
!                          !
!    Immersed Boundary     !
!         Module           !
!                          !
!--------------------------! 
! 
! authors: Dr.-Ing. Davide Gatti
!          B.Sc. Ahmet Cumhur
! 
! date:    28.04.26
! 


module ibmm
    use, intrinsic :: iso_c_binding
    use, intrinsic :: iso_fortran_env, only: int64
    use :: init, only: dns_type, grid_type, VAR_U, VAR_V, VAR_W, VAR_P, &
        is_face_staggered, face_at, cell_center_at
    use :: blocks, only: block_set_type, subdivide_node_line, FACE_PHYS, FACE_CLOSED
    use :: io, only: to_c_string, read_block_active, read_block_masks, read_mask_window
    use :: comm, only: comm_type, exchange_scalar_halos
    implicit none

    real(C_DOUBLE), parameter :: SOLID = 1.0d30
    real(C_DOUBLE), parameter :: DEFAULT_TOL = 1.0d-10
    integer(C_INT), parameter :: MAX_ITER = 200

    !========================
    ! IBM TYPE
    !========================
    type :: ibm_type
        integer :: n_wave_x, n_wave_z
        real(C_DOUBLE) :: amp_x, phase_x
        real(C_DOUBLE) :: amp_z, phase_z

        real(C_DOUBLE), allocatable :: coef(:,:,:,:,:) ! (0:nb+1,...,VAR_U:VAR_W,nBlocks)
        real(C_DOUBLE), allocatable :: mu(:,:,:,:,:)

        ! [ibm] band_filter: compressed list of near-body FLUID velocity
        ! DOFs (built by init_ibm_band, one entry per DOF, mapped to the
        ! device only when the filter is on — off allocates and maps
        ! nothing). bandDirs bit d (0/1/2 = x/y/z) is set when BOTH +-1
        ! same-component neighbours are fluid: the filter never reads a
        ! solid velocity.
        integer(C_INT) :: nBand = 0_C_INT
        integer(C_INT), allocatable :: bandI(:), bandJ(:), bandK(:)
        integer(C_INT), allocatable :: bandVar(:), bandBlk(:)
        integer(C_SIGNED_CHAR), allocatable :: bandDirs(:)

    end type ibm_type

    ! The one indicator signature every host-side geometry consumer takes
    ! (walldist, the classify_* routines, moby_prepare's coefficient
    ! tiles): exactly the analytic isInBody's. geometry_stl.f90 provides
    ! the same signature reading its own module state, so STL and analytic
    ! bodies flow through identical machinery. The DEVICE coefficient
    ! kernel (set_ibm_coeff) keeps calling the concrete isInBody directly:
    ! declare-target procedure arguments are not portable.
    abstract interface
        logical function body_indicator_i(xIN, ibm, dns)
            import :: C_DOUBLE, ibm_type, dns_type
            real(C_DOUBLE), intent(in) :: xIN(1:3)
            type(ibm_type), intent(in) :: ibm
            type(dns_type), intent(in) :: dns
        end function body_indicator_i
    end interface

    interface
        function fdm_h5_read_ibm_coeff_blocks(file_name, nbx, nby, nbz, n_blocks, id_start, &
                block_origin, block_level, lx, ly, lz, re, found, coef) &
                bind(C, name="fdm_h5_read_ibm_coeff_blocks") result(ierr)
            import :: C_CHAR, C_INT, C_DOUBLE
            character(kind=C_CHAR), intent(in) :: file_name(*)
            integer(C_INT), value :: nbx, nby, nbz, n_blocks, id_start
            integer(C_INT), intent(in) :: block_origin(*), block_level(*)
            real(C_DOUBLE), value :: lx, ly, lz, re
            integer(C_INT), intent(out) :: found
            real(C_DOUBLE), intent(inout) :: coef(*)
            integer(C_INT) :: ierr
        end function fdm_h5_read_ibm_coeff_blocks

        function fdm_h5_read_ibm_coeff(file_name, nbx, nby, nbz, n_blocks, block_origin, &
                global_nx, global_ny, global_nz, &
                lx, ly, lz, re, coef) bind(C, name="fdm_h5_read_ibm_coeff") result(ierr)
            import :: C_CHAR, C_INT, C_DOUBLE
            character(kind=C_CHAR), intent(in) :: file_name(*)
            integer(C_INT), value :: nbx, nby, nbz, n_blocks
            integer(C_INT), intent(in) :: block_origin(*)
            integer(C_INT), value :: global_nx, global_ny, global_nz
            real(C_DOUBLE), value :: lx, ly, lz, re
            real(C_DOUBLE), intent(inout) :: coef(*)
            integer(C_INT) :: ierr
        end function fdm_h5_read_ibm_coeff
    end interface

!$omp declare target(isInBody, distance3, bisection, add_neighbor_coeff)

contains

!========================
! INITIALIZE IBM
!========================
    ! The analytic wavy-wall geometry parameters, needed both by the
    ! coefficient kernels and by the block classification that runs before
    ! the block set exists.
    subroutine set_ibm_geometry_defaults(ibm)
        type(ibm_type), intent(inout) :: ibm

        ibm%n_wave_x = 1
        ibm%n_wave_z = 1
        ibm%amp_x = 2.5d-2
        ibm%amp_z = 2.5d-2
        ibm%phase_x = 0.0d0
        ibm%phase_z = 0.0d0
    end subroutine set_ibm_geometry_defaults

    subroutine init_ibm(ibm, blk)
        type(ibm_type), intent(inout) :: ibm
        type(block_set_type), intent(in) :: blk
        integer :: nx, ny, nz

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))

        call set_ibm_geometry_defaults(ibm)

        allocate(ibm%coef(0:nx+1,0:ny+1,0:nz+1,VAR_U:VAR_W,blk%nBlocks))
        allocate(ibm%mu(0:nx+1,0:ny+1,0:nz+1,VAR_U:VAR_W,blk%nBlocks))
        ibm%coef = 0.0d0
        ibm%mu = 1.0d0
    end subroutine init_ibm

    subroutine enter_ibm_data(ibm, dns)
        type(ibm_type), intent(inout) :: ibm
        type(dns_type), intent(in)    :: dns

#ifdef USE_OPENMP_OFFLOAD
        !$omp target enter data map(to: ibm)
        !$omp target enter data map(to: ibm%coef, ibm%mu)
#endif
    end subroutine enter_ibm_data

    subroutine exit_ibm_data(ibm, dns)
        type(ibm_type), intent(inout) :: ibm
        type(dns_type), intent(in)    :: dns

#ifdef USE_OPENMP_OFFLOAD
        if (allocated(ibm%bandI)) then
            !$omp target exit data map(delete: ibm%bandI, ibm%bandJ, ibm%bandK, &
            !$omp& ibm%bandVar, ibm%bandBlk, ibm%bandDirs)
        end if
        !$omp target exit data map(delete: ibm%coef, ibm%mu)
        !$omp target exit data map(delete: ibm)
#endif
    end subroutine exit_ibm_data

    ! [ibm] band_filter: build the compressed near-body band list. A DOF is
    ! in the band when a solid same-component DOF lies within band_width
    ! cells (box distance): the solid marker is dilated band_width times by
    ! one cell, with a halo exchange after every pass so the dilation
    ! reaches across block boundaries exactly (the halo layer is one cell
    ! deep). Runs on the DEVICE against the mapped coefficients (the
    ! analytic path fills coef only there); the list itself is built on the
    ! host and mapped once. Called only when the filter is on — off maps
    ! and allocates nothing.
    subroutine init_ibm_band(ibm, dns, blk, c, has_terminal)
        type(ibm_type), intent(inout) :: ibm
        type(dns_type), intent(in) :: dns
        type(block_set_type), intent(in) :: blk
        type(comm_type), intent(inout) :: c
        logical, intent(in) :: has_terminal

        real(C_DOUBLE), allocatable :: mark(:,:,:,:), tmp(:,:,:,:)
        real(C_DOUBLE), allocatable :: solidh(:,:,:,:), markh(:,:,:,:)
        integer(C_INT), allocatable :: bi(:), bj(:), bk(:), bv(:), bb(:)
        integer(C_SIGNED_CHAR), allocatable :: bd(:)
        real(C_DOUBLE) :: thresh
        integer :: nx, ny, nz, nBlocks, var, i, j, k, b, pass, n, cap
        integer :: dirs, istart, jstart, kstart
        logical :: fx, fy, fz

        nx = int(blk%nb(1)); ny = int(blk%nb(2)); nz = int(blk%nb(3))
        nBlocks = int(blk%nBlocks)
        thresh = 0.5d0*SOLID/dns%re

        allocate(mark(0:nx+1,0:ny+1,0:nz+1,nBlocks), tmp(0:nx+1,0:ny+1,0:nz+1,nBlocks))
        allocate(solidh(0:nx+1,0:ny+1,0:nz+1,nBlocks), markh(0:nx+1,0:ny+1,0:nz+1,nBlocks))
        !$omp target enter data map(alloc: mark, tmp)

        cap = 1024
        allocate(bi(cap), bj(cap), bk(cap), bv(cap), bb(cap), bd(cap))
        ibm%nBand = 0_C_INT

        do var = VAR_U, VAR_W
            ! solid marker (ghost-inclusive: coef carries ghost values on
            ! both the file and the analytic path)
            !$omp target teams distribute parallel do collapse(4) &
            !$omp& map(to: var, thresh, nx, ny, nz, nBlocks, ibm%coef) map(tofrom: mark) &
            !$omp& private(i,j,k,b)
            do b = 1, nBlocks
            do k = 0, nz+1
                do j = 0, ny+1
                    do i = 0, nx+1
                        mark(i,j,k,b) = merge(1.0d0, 0.0d0, ibm%coef(i,j,k,var,b) >= thresh)
                    end do
                end do
            end do
            end do
            !$omp end target teams distribute parallel do
            !$omp target update from(mark)
            solidh = mark

            ! band_width one-cell dilations, exchange-interleaved so the
            ! reach crosses block boundaries exactly
            do pass = 1, int(dns%ibm_band_width)
                !$omp target teams distribute parallel do collapse(4) &
                !$omp& map(to: nx, ny, nz, nBlocks, mark) map(tofrom: tmp) &
                !$omp& private(i,j,k,b)
                do b = 1, nBlocks
                do k = 1, nz
                    do j = 1, ny
                        do i = 1, nx
                            tmp(i,j,k,b) = max(mark(i,j,k,b), &
                                mark(i-1,j,k,b), mark(i+1,j,k,b), &
                                mark(i,j-1,k,b), mark(i,j+1,k,b), &
                                mark(i,j,k-1,b), mark(i,j,k+1,b))
                        end do
                    end do
                end do
                end do
                !$omp end target teams distribute parallel do
                !$omp target teams distribute parallel do collapse(4) &
                !$omp& map(to: nx, ny, nz, nBlocks, tmp) map(tofrom: mark) &
                !$omp& private(i,j,k,b)
                do b = 1, nBlocks
                do k = 1, nz
                    do j = 1, ny
                        do i = 1, nx
                            mark(i,j,k,b) = tmp(i,j,k,b)
                        end do
                    end do
                end do
                end do
                !$omp end target teams distribute parallel do
                call exchange_scalar_halos(c, mark, blk)
            end do
            !$omp target update from(mark)
            markh = mark

            ! host: compressed list of fluid band DOFs, skipping faces the
            ! predictor pins (physical/closed low faces, momentum_face_start)
            do b = 1, nBlocks
                istart = merge(2, 1, (blk%physLow(1,b) == FACE_PHYS .or. &
                    blk%physLow(1,b) == FACE_CLOSED) .and. var == VAR_U)
                jstart = merge(2, 1, (blk%physLow(2,b) == FACE_PHYS .or. &
                    blk%physLow(2,b) == FACE_CLOSED) .and. var == VAR_V)
                kstart = merge(2, 1, (blk%physLow(3,b) == FACE_PHYS .or. &
                    blk%physLow(3,b) == FACE_CLOSED) .and. var == VAR_W)
                do k = kstart, nz
                    do j = jstart, ny
                        do i = istart, nx
                            if (solidh(i,j,k,b) > 0.5d0) cycle          ! solid DOF
                            if (markh(i,j,k,b) < 0.5d0) cycle           ! not in band
                            fx = solidh(i-1,j,k,b) < 0.5d0 .and. solidh(i+1,j,k,b) < 0.5d0
                            fy = solidh(i,j-1,k,b) < 0.5d0 .and. solidh(i,j+1,k,b) < 0.5d0
                            fz = solidh(i,j,k-1,b) < 0.5d0 .and. solidh(i,j,k+1,b) < 0.5d0
                            if (.not. (fx .or. fy .or. fz)) cycle
                            dirs = 0
                            if (fx) dirs = dirs + 1
                            if (fy) dirs = dirs + 2
                            if (fz) dirs = dirs + 4
                            if (ibm%nBand >= cap) then
                                cap = 2*cap
                                call grow_int(bi, cap); call grow_int(bj, cap)
                                call grow_int(bk, cap); call grow_int(bv, cap)
                                call grow_int(bb, cap); call grow_byte(bd, cap)
                            end if
                            ibm%nBand = ibm%nBand + 1_C_INT
                            bi(ibm%nBand) = int(i, C_INT); bj(ibm%nBand) = int(j, C_INT)
                            bk(ibm%nBand) = int(k, C_INT); bv(ibm%nBand) = int(var, C_INT)
                            bb(ibm%nBand) = int(b, C_INT)
                            bd(ibm%nBand) = int(dirs, C_SIGNED_CHAR)
                        end do
                    end do
                end do
            end do
        end do

        !$omp target exit data map(delete: mark, tmp)
        deallocate(mark, tmp, solidh, markh)

        allocate(ibm%bandI(max(1, int(ibm%nBand))), ibm%bandJ(max(1, int(ibm%nBand))), &
            ibm%bandK(max(1, int(ibm%nBand))), ibm%bandVar(max(1, int(ibm%nBand))), &
            ibm%bandBlk(max(1, int(ibm%nBand))), ibm%bandDirs(max(1, int(ibm%nBand))))
        ibm%bandI(1:ibm%nBand) = bi(1:ibm%nBand); ibm%bandJ(1:ibm%nBand) = bj(1:ibm%nBand)
        ibm%bandK(1:ibm%nBand) = bk(1:ibm%nBand); ibm%bandVar(1:ibm%nBand) = bv(1:ibm%nBand)
        ibm%bandBlk(1:ibm%nBand) = bb(1:ibm%nBand); ibm%bandDirs(1:ibm%nBand) = bd(1:ibm%nBand)
        deallocate(bi, bj, bk, bv, bb, bd)
#ifdef USE_OPENMP_OFFLOAD
        !$omp target enter data map(to: ibm%bandI, ibm%bandJ, ibm%bandK, &
        !$omp& ibm%bandVar, ibm%bandBlk, ibm%bandDirs)
#endif
        if (has_terminal) print *, "IBM band filter: ", ibm%nBand, " band DOFs (width ", &
            dns%ibm_band_width, ", theta ", dns%ibm_band_theta, ")"
    end subroutine init_ibm_band

    subroutine grow_int(a, cap)
        integer(C_INT), allocatable, intent(inout) :: a(:)
        integer, intent(in) :: cap
        integer(C_INT), allocatable :: t(:)

        allocate(t(cap))
        t(1:size(a)) = a
        call move_alloc(t, a)
    end subroutine grow_int

    subroutine grow_byte(a, cap)
        integer(C_SIGNED_CHAR), allocatable, intent(inout) :: a(:)
        integer, intent(in) :: cap
        integer(C_SIGNED_CHAR), allocatable :: t(:)

        allocate(t(cap))
        t(1:size(a)) = a
        call move_alloc(t, a)
    end subroutine grow_byte

    subroutine read_ibm_coeff_file(ibm, dns, blk, has_terminal)
        type(ibm_type), intent(inout) :: ibm
        type(dns_type), intent(in) :: dns
        type(block_set_type), intent(in) :: blk
        logical, intent(in) :: has_terminal

        character(kind=C_CHAR,len=:), allocatable :: c_file_name
        integer(C_INT) :: ierr

        integer(C_INT) :: found

        if (len_trim(dns%ibm_coeff_file) == 0) return
        if (has_terminal) print *, "reading IBM coefficients: ", trim(dns%ibm_coeff_file)

        c_file_name = to_c_string(dns%ibm_coeff_file)
        ! Block-table layout first (refined runs); fall back to the legacy
        ! global ghost-layer layout, which holds level-0 data only.
        ierr = fdm_h5_read_ibm_coeff_blocks(c_file_name, blk%nb(1), blk%nb(2), blk%nb(3), &
            blk%nBlocks, blk%idStart, blk%origin, blk%level, &
            dns%leng(1), dns%leng(2), dns%leng(3), dns%re, found, ibm%coef)
        if (ierr == 2_C_INT) then
            if (has_terminal) print *, "error: coefficient file block table does not match", &
                " the solver's leaf table (stale file?): ", trim(dns%ibm_coeff_file)
            error stop
        end if
        if (ierr /= 0_C_INT) then
            if (has_terminal) print *, "error: could not read IBM coefficient file: ", trim(dns%ibm_coeff_file)
            error stop
        end if
        if (found /= 0_C_INT) return

        if (any(blk%level(1:blk%nBlocks) /= 0_C_INT)) then
            if (has_terminal) print *, "error: legacy coefficient file needs single-level blocks;", &
                " regenerate with mobygeom block-table"
            error stop
        end if
        ierr = fdm_h5_read_ibm_coeff(c_file_name, blk%nb(1), blk%nb(2), blk%nb(3), &
            blk%nBlocks, blk%origin, &
            dns%globalSize(1), dns%globalSize(2), dns%globalSize(3), &
            dns%leng(1), dns%leng(2), dns%leng(3), dns%re, ibm%coef)
        if (ierr /= 0_C_INT) then
            if (has_terminal) print *, "error: could not read IBM coefficient file: ", trim(dns%ibm_coeff_file)
            error stop
        end if
    end subroutine read_ibm_coeff_file


    ! Height of the analytic wavy bottom wall at streamwise position x.
    real(C_DOUBLE) function wavy_wall_height(x, ibm, dns) result(y_body)
!$omp declare target
        real(C_DOUBLE), intent(in) :: x
        type(ibm_type), intent(in) :: ibm
        type(dns_type), intent(in) :: dns

        real(C_DOUBLE), parameter :: pi = 3.141592653589793d0
        real(C_DOUBLE), parameter :: y_offset = 1.0d-2

        y_body = ibm%amp_x * 0.5d0 * &
                 (1.0d0 + sin(2.0d0*pi*real(ibm%n_wave_x,C_DOUBLE)*x/dns%leng(1) + ibm%phase_x)) + &
                 y_offset
    end function wavy_wall_height

    logical function isInBody(xIN, ibm, dns)
        implicit none

        real(C_DOUBLE), intent(in) :: xIN(1:3)
        type(ibm_type), intent(in) :: ibm
        type(dns_type), intent(in) :: dns

        isInBody = (xIN(2) < wavy_wall_height(xIN(1), ibm, dns))
    end function isInBody

    real(C_DOUBLE) function distance3(xA, xB) result(d)
        real(C_DOUBLE), intent(in) :: xA(1:3), xB(1:3)

        d = sqrt((xB(1)-xA(1))**2 + (xB(2)-xA(2))**2 + (xB(3)-xA(3))**2)
    end function distance3

    subroutine bisection(xAin,xB,ibm,dns)
        real(C_DOUBLE), intent(in) :: xAin(1:3)
        real(C_DOUBLE), intent(inout):: xB(1:3)
        type(ibm_type), intent(in) :: ibm
        type(dns_type), intent(in) :: dns

        real(C_DOUBLE) :: xA(1:3), xM(1:3)
        logical :: la, lm
        integer(C_INT) :: it

        xA = xAin
        xM = xA
        la = isInBody(xA, ibm, dns)

        do it = 1, MAX_ITER
            xM = 0.5d0*(xA + xB)
            if (distance3(xA, xB) < DEFAULT_TOL) exit

            lm = isInBody(xM, ibm, dns)
            if (lm .eqv. la) then
                xA = xM
            else
                xB = xM
            end if
        end do
        xB = xM
    end subroutine bisection

    subroutine add_neighbor_coeff(coeff, xA, xB, ibm, dns)
        real(C_DOUBLE), intent(inout) :: coeff
        real(C_DOUBLE), intent(in) :: xA(1:3)
        real(C_DOUBLE), intent(inout) :: xB(1:3)
        type(ibm_type), intent(in) :: ibm
        type(dns_type), intent(in) :: dns

        real(C_DOUBLE) :: d0, d

        if (isInBody(xB, ibm, dns)) then
            d0 = distance3(xA, xB)
            call bisection(xA, xB, ibm, dns)
            d = distance3(xA, xB)
            coeff = coeff + ((d0-d)/d)/d0**2
        end if
    end subroutine add_neighbor_coeff

    ! Phase 2 classification for the analytic IBM: a block is removable iff
    ! the block dilated by one halo cell is solid (isInBody) at cell centres
    ! and all three staggered locations. active is in x-fastest lattice
    ! raster order, 1 = keep.
    subroutine classify_active_blocks(active, dns, g, ibm, periodic, inside)
        integer(C_INT), intent(out) :: active(:)
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        type(ibm_type), intent(inout) :: ibm
        logical(C_BOOL), intent(in) :: periodic(1:3)
        procedure(body_indicator_i) :: inside

        integer :: nb, nTiles(3), gx, gy, gz, raster, d
        integer :: i, j, k, var, o(3)
        real(C_DOUBLE) :: xA(3)
        logical :: buried

        call set_ibm_geometry_defaults(ibm)

        nb = int(dns%block_nb)
        do d = 1, 3
            if (mod(int(dns%globalSize(d)), nb) /= 0) then
                error stop "[blocks] nb must divide the global grid in every direction"
            end if
            nTiles(d) = int(dns%globalSize(d))/nb
        end do
        if (size(active) /= product(nTiles)) error stop "block active mask size mismatch"

        raster = 0
        do gz = 0, nTiles(3) - 1
            do gy = 0, nTiles(2) - 1
                do gx = 0, nTiles(1) - 1
                    raster = raster + 1
                    o = [gx, gy, gz]*nb
                    buried = .true.
                    outer: do var = int(VAR_U), int(VAR_P)
                        ! Dilated window: 1-based cell indices o .. o+nb+1.
                        do k = o(3), o(3) + nb + 1
                            do j = o(2), o(2) + nb + 1
                                do i = o(1), o(1) + nb + 1
                                    xA(1) = location_coord(g%xNode, int(dns%globalSize(1)), &
                                        dns%leng(1), periodic(1), 1, var, i)
                                    xA(2) = location_coord(g%yNode, int(dns%globalSize(2)), &
                                        dns%leng(2), periodic(2), 2, var, j)
                                    xA(3) = location_coord(g%zNode, int(dns%globalSize(3)), &
                                        dns%leng(3), periodic(3), 3, var, k)
                                    if (.not. inside(xA, ibm, dns)) then
                                        buried = .false.
                                        exit outer
                                    end if
                                end do
                            end do
                        end do
                    end do outer
                    active(raster) = merge(0_C_INT, 1_C_INT, buried)
                end do
            end do
        end do
    end subroutine classify_active_blocks

    ! Coordinate of variable `var` at 1-based cell index `idx` along direction
    ! `dir`, matching slice_grid_direction's staggering. `n` is the cell count
    ! of the node line, so the same routine serves the global grid and the
    ! per-level refinement lines.
    real(C_DOUBLE) function location_coord(node, n, length, periodic, dir, var, idx) result(x)
        real(C_DOUBLE), intent(in) :: node(0:)
        integer, intent(in) :: n
        real(C_DOUBLE), intent(in) :: length
        logical(C_BOOL), intent(in) :: periodic
        integer, intent(in) :: dir, var, idx

        if (is_face_staggered(dir, var)) then
            x = face_at(node, n, length, idx - 1, periodic)
        else
            x = cell_center_at(node, n, length, idx, periodic)
        end if
    end function location_coord

    ! Per-level geometry masks for geometry-driven refinement (analytic
    ! IBM). For every level-l lattice cell, in x-fastest raster order:
    !   touch(c, l+1)  = the one-halo dilated block straddles the surface
    !                    (mixed solid/fluid over the 4 variable locations)
    !   buried(c, l+1) = the dilated block is solid everywhere (removable)
    ! Level lines are built here by midpoint subdivision, identical to the
    ! solver's per-level metric lines.
    subroutine classify_block_geometry(touch, buried, dns, g, ibm, periodic, nLevels, inside)
        integer(C_INT), intent(out) :: touch(:,:), buried(:,:)
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        type(ibm_type), intent(inout) :: ibm
        logical(C_BOOL), intent(in) :: periodic(1:3)
        integer, intent(in) :: nLevels
        procedure(body_indicator_i) :: inside

        integer :: nb, l, nTiles(3), gx, gy, gz, raster, d
        integer :: i, j, k, var, o(3), nf(3), nl(3)
        real(C_DOUBLE) :: xA(3)
        real(C_DOUBLE), allocatable :: lineX(:,:), lineY(:,:), lineZ(:,:)
        logical :: anySolid, anyFluid, isSolid

        call set_ibm_geometry_defaults(ibm)
        nb = int(dns%block_nb)

        ! Per-level lines follow [blocks] refine_dims: refined directions
        ! are midpoint-subdivided, fixed ones copy the global line
        ! (blocks.f90 build_level_lines is the solver-side counterpart).
        nf = int(dns%globalSize)*2**((nLevels - 1)*int(dns%block_refine_mask))
        allocate(lineX(0:nf(1), nLevels), lineY(0:nf(2), nLevels), lineZ(0:nf(3), nLevels))
        lineX(0:int(dns%globalSize(1)),1) = g%xNode
        lineY(0:int(dns%globalSize(2)),1) = g%yNode
        lineZ(0:int(dns%globalSize(3)),1) = g%zNode
        do l = 2, nLevels
            if (dns%block_refine_mask(1) == 1_C_INT) then
                call subdivide_node_line(lineX(0:int(dns%globalSize(1))*2**(l-2), l-1), &
                                         lineX(0:int(dns%globalSize(1))*2**(l-1), l))
            else
                lineX(0:int(dns%globalSize(1)), l) = lineX(0:int(dns%globalSize(1)), 1)
            end if
            if (dns%block_refine_mask(2) == 1_C_INT) then
                call subdivide_node_line(lineY(0:int(dns%globalSize(2))*2**(l-2), l-1), &
                                         lineY(0:int(dns%globalSize(2))*2**(l-1), l))
            else
                lineY(0:int(dns%globalSize(2)), l) = lineY(0:int(dns%globalSize(2)), 1)
            end if
            if (dns%block_refine_mask(3) == 1_C_INT) then
                call subdivide_node_line(lineZ(0:int(dns%globalSize(3))*2**(l-2), l-1), &
                                         lineZ(0:int(dns%globalSize(3))*2**(l-1), l))
            else
                lineZ(0:int(dns%globalSize(3)), l) = lineZ(0:int(dns%globalSize(3)), 1)
            end if
        end do

        do l = 1, nLevels
            nl = int(dns%globalSize)*2**((l-1)*int(dns%block_refine_mask))
            nTiles = nl/nb
            raster = 0
            do gz = 0, nTiles(3) - 1
                do gy = 0, nTiles(2) - 1
                    do gx = 0, nTiles(1) - 1
                        raster = raster + 1
                        o = [gx, gy, gz]*nb
                        anySolid = .false.
                        anyFluid = .false.
                        scan: do var = int(VAR_U), int(VAR_P)
                            do k = o(3), o(3) + nb + 1
                                do j = o(2), o(2) + nb + 1
                                    do i = o(1), o(1) + nb + 1
                                        xA(1) = location_coord(lineX(:,l), nl(1), dns%leng(1), &
                                            periodic(1), 1, var, i)
                                        xA(2) = location_coord(lineY(:,l), nl(2), dns%leng(2), &
                                            periodic(2), 2, var, j)
                                        xA(3) = location_coord(lineZ(:,l), nl(3), dns%leng(3), &
                                            periodic(3), 3, var, k)
                                        isSolid = inside(xA, ibm, dns)
                                        anySolid = anySolid .or. isSolid
                                        anyFluid = anyFluid .or. .not. isSolid
                                        if (anySolid .and. anyFluid) exit scan
                                    end do
                                end do
                            end do
                        end do scan
                        touch(raster, l) = merge(1_C_INT, 0_C_INT, anySolid .and. anyFluid)
                        buried(raster, l) = merge(1_C_INT, 0_C_INT, anySolid .and. .not. anyFluid)
                    end do
                end do
            end do
        end do

        deallocate(lineX, lineY, lineZ)
    end subroutine classify_block_geometry

    ! Produce the per-block keep mask for solid-block removal ([blocks]
    ! remove_solid): read it from the coefficient file (mobygeom
    ! block-active) or classify it from the analytic geometry. The mask is
    ! allocated here to the lattice size; the caller hands it to
    ! init_block_set. This is the geometry input to block-set construction,
    ! kept out of main so the construction dispatch stays thin.
    subroutine classify_active_mask(active, dns, g, ibm, periodic, has_terminal, inside)
        integer(C_INT), allocatable, intent(out) :: active(:)
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        type(ibm_type), intent(inout) :: ibm
        logical(C_BOOL), intent(in) :: periodic(1:3)
        logical, intent(in) :: has_terminal
        procedure(body_indicator_i) :: inside

        logical :: found

        if (any(mod(dns%globalSize, dns%block_nb) /= 0_C_INT)) then
            error stop "[blocks] nb must divide the global grid in every direction"
        end if
        if (dns%block_refine_nboxes > 0_C_INT) then
            error stop "solid-block removal with box refinement is unsupported; use refine_body"
        end if
        allocate(active(product(dns%globalSize/dns%block_nb)))
        if (len_trim(dns%ibm_coeff_file) > 0) then
            call read_block_active(active, found, dns, has_terminal)
            if (.not. found) then
                if (has_terminal) print *, &
                    "coefficient file has no block_active table; keeping all blocks"
                active = 1_C_INT
            end if
        else
            call classify_active_blocks(active, dns, g, ibm, periodic, inside)
        end if
    end subroutine classify_active_mask

    ! Produce the per-level touch/buried masks for geometry-driven
    ! refinement ([blocks] refine_body): read them from the coefficient
    ! file (mobygeom block-table) or classify them from the analytic
    ! geometry. Both arrays are allocated here for the finest lattice (one
    ! column per level); the caller hands them to init_block_set.
    subroutine classify_refinement_masks(touch, buried, maskLo, maskDims, dns, g, ibm, &
            periodic, has_terminal, inside)
        integer(C_INT), allocatable, intent(out) :: touch(:,:), buried(:,:)
        ! Per-level mask windows (block coords): deep-refinement rasters
        ! are stored/held windowed to the padded STL bbox; the analytic
        ! path and legacy files use full-lattice windows.
        integer(C_INT), allocatable, intent(out) :: maskLo(:,:), maskDims(:,:)
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        type(ibm_type), intent(inout) :: ibm
        logical(C_BOOL), intent(in) :: periodic(1:3)
        logical, intent(in) :: has_terminal
        procedure(body_indicator_i) :: inside

        integer :: level, maskCount, nLevels, maxCount
        integer(int64) :: finest
        logical :: found, has_win

        if (.not. dns%ibm_enabled) then
            error stop "[blocks] refine_body needs the IBM enabled"
        end if
        if (any(mod(dns%globalSize, dns%block_nb) /= 0_C_INT)) then
            error stop "[blocks] nb must divide the global grid in every direction"
        end if
        nLevels = int(dns%block_refine_levels) + 1
        allocate(maskLo(3, nLevels), maskDims(3, nLevels))
        if (len_trim(dns%ibm_coeff_file) > 0) then
            ! File-based geometry: masks computed by mobygeom block-table,
            ! WINDOWED on deep-refinement files (legacy = full rasters).
            do level = 0, nLevels - 1
                call read_mask_window(maskLo(:, level+1), maskDims(:, level+1), &
                    has_win, level, dns, has_terminal)
                if (.not. has_win) then
                    maskLo(:, level+1) = 0_C_INT
                    maskDims(:, level+1) = (dns%globalSize/dns%block_nb) &
                        *2**(int(level, C_INT)*dns%block_refine_mask)
                end if
            end do
            maxCount = 0
            do level = 0, nLevels - 1
                maxCount = max(maxCount, int(product(maskDims(:, level+1))))
            end do
            allocate(touch(max(1, maxCount), nLevels))
            allocate(buried(size(touch, 1), size(touch, 2)))
            do level = 0, nLevels - 1
                maskCount = int(product(maskDims(:, level+1)))
                call read_block_masks(touch(1:maskCount, level+1), &
                    buried(1:maskCount, level+1), level, maskCount, &
                    found, dns, has_terminal)
                if (.not. found) then
                    error stop "coefficient file has no refinement masks; run mobygeom block-table"
                end if
            end do
        else
            ! Analytic geometry: full-lattice rasters (classify_block_geometry
            ! is not windowed) -- guard the deep-refinement case where they
            ! no longer fit; the file path handles it.
            finest = int(product(int((dns%globalSize/dns%block_nb) &
                *2**(dns%block_refine_levels*dns%block_refine_mask), int64)), int64)
            if (finest > 200000000_int64) then
                error stop "[blocks] refine_body: analytic classification needs a " &
                    // "dense finest lattice too large for this depth; use the " &
                    // "mobygeom block-table file path"
            end if
            do level = 0, nLevels - 1
                maskLo(:, level+1) = 0_C_INT
                maskDims(:, level+1) = (dns%globalSize/dns%block_nb) &
                    *2**(int(level, C_INT)*dns%block_refine_mask)
            end do
            allocate(touch(int(finest), nLevels))
            allocate(buried(size(touch, 1), size(touch, 2)))
            call classify_block_geometry(touch, buried, dns, g, ibm, periodic, nLevels, inside)
        end if
    end subroutine classify_refinement_masks

    subroutine set_ibm_coeff(dns, blk, ibm, var)
        type(dns_type), intent(in) :: dns
        type(block_set_type), intent(in) :: blk
        type(ibm_type), intent(inout) :: ibm
        integer(C_INT), intent(in) :: var

        integer :: ix, iy, iz, b, nBlocks
        integer :: ilo, ihi, jlo, jhi, klo, khi
        real(C_DOUBLE) :: xA(1:3), xB(1:3), coeff
        real(C_DOUBLE) :: re_inv, solid_coef
        logical(C_BOOL) :: enabled

        if (var < VAR_U .or. var > VAR_W) error stop "invalid IBM coefficient variable"

        ilo = lbound(ibm%coef,1)
        ihi = ubound(ibm%coef,1)
        jlo = lbound(ibm%coef,2)
        jhi = ubound(ibm%coef,2)
        klo = lbound(ibm%coef,3)
        khi = ubound(ibm%coef,3)
        nBlocks = size(ibm%coef,5)

        enabled = dns%ibm_enabled
        re_inv = 1.0d0/dns%re
        solid_coef = SOLID*re_inv
#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do collapse(4) &
        !$omp& map(to: var, enabled, re_inv, solid_coef, ilo, ihi, jlo, jhi, klo, khi, nBlocks, dns, ibm, blk, &
        !$omp& blk%x, blk%y, blk%z) &
        !$omp& map(tofrom: ibm%coef) &
        !$omp& private(ix,iy,iz,b,xA,xB,coeff)
#endif
        do b = 1, nBlocks
        do iz = klo, khi
            do iy = jlo, jhi
                do ix = ilo, ihi
                    ibm%coef(ix,iy,iz,var,b) = 0.0d0
                    if (.not. enabled) cycle

                    xA(1) = blk%x(ix,var,b)
                    xA(2) = blk%y(iy,var,b)
                    xA(3) = blk%z(iz,var,b)
                    if (isInBody(xA, ibm, dns)) then
                        ibm%coef(ix,iy,iz,var,b) = solid_coef
#ifdef USE_IBM_SECONDORDER
                    else
                        coeff = 0.0d0

                        xB(1) = blk%x(ix-1,var,b); xB(2) = xA(2);             xB(3) = xA(3)
                        call add_neighbor_coeff(coeff, xA, xB, ibm, dns)
                        xB(1) = blk%x(ix+1,var,b); xB(2) = xA(2);             xB(3) = xA(3)
                        call add_neighbor_coeff(coeff, xA, xB, ibm, dns)
                        xB(1) = xA(1);             xB(2) = blk%y(iy-1,var,b); xB(3) = xA(3)
                        call add_neighbor_coeff(coeff, xA, xB, ibm, dns)
                        xB(1) = xA(1);             xB(2) = blk%y(iy+1,var,b); xB(3) = xA(3)
                        call add_neighbor_coeff(coeff, xA, xB, ibm, dns)
                        xB(1) = xA(1);             xB(2) = xA(2);             xB(3) = blk%z(iz-1,var,b)
                        call add_neighbor_coeff(coeff, xA, xB, ibm, dns)
                        xB(1) = xA(1);             xB(2) = xA(2);             xB(3) = blk%z(iz+1,var,b)
                        call add_neighbor_coeff(coeff, xA, xB, ibm, dns)

                        ibm%coef(ix,iy,iz,var,b) = coeff*re_inv
#endif
                    end if
                end do
            end do
        end do
        end do
#ifdef USE_OPENMP_OFFLOAD
        !$omp end target teams distribute parallel do
#endif
    end subroutine set_ibm_coeff

    ! Host twin of set_ibm_coeff for moby_prepare's STL path: the same
    ! graded sharp-interface arithmetic over ANY indicator (the device
    ! kernel above keeps calling the concrete isInBody -- declare-target
    ! procedure arguments are not portable). KEEP THE TWO IN LOCKSTEP.
    subroutine set_ibm_coeff_host(dns, blk, ibm, var, inside)
        type(dns_type), intent(in) :: dns
        type(block_set_type), intent(in) :: blk
        type(ibm_type), intent(inout) :: ibm
        integer(C_INT), intent(in) :: var
        procedure(body_indicator_i) :: inside

        integer :: ix, iy, iz, b, nBlocks
        integer :: ilo, ihi, jlo, jhi, klo, khi
        real(C_DOUBLE) :: xA(1:3), xB(1:3), coeff
        real(C_DOUBLE) :: re_inv, solid_coef

        if (var < VAR_U .or. var > VAR_W) error stop "invalid IBM coefficient variable"

        ilo = lbound(ibm%coef,1)
        ihi = ubound(ibm%coef,1)
        jlo = lbound(ibm%coef,2)
        jhi = ubound(ibm%coef,2)
        klo = lbound(ibm%coef,3)
        khi = ubound(ibm%coef,3)
        nBlocks = size(ibm%coef,5)

        re_inv = 1.0d0/dns%re
        solid_coef = SOLID*re_inv

        !$omp parallel do collapse(2) private(ix,iy,iz,b,xA,xB,coeff)
        do b = 1, nBlocks
        do iz = klo, khi
            do iy = jlo, jhi
                do ix = ilo, ihi
                    ibm%coef(ix,iy,iz,var,b) = 0.0d0
                    if (.not. dns%ibm_enabled) cycle

                    xA(1) = blk%x(ix,var,b)
                    xA(2) = blk%y(iy,var,b)
                    xA(3) = blk%z(iz,var,b)
                    if (inside(xA, ibm, dns)) then
                        ibm%coef(ix,iy,iz,var,b) = solid_coef
#ifdef USE_IBM_SECONDORDER
                    else
                        coeff = 0.0d0

                        xB(1) = blk%x(ix-1,var,b); xB(2) = xA(2);             xB(3) = xA(3)
                        call add_neighbor_coeff_host(coeff, xA, xB, ibm, dns, inside)
                        xB(1) = blk%x(ix+1,var,b); xB(2) = xA(2);             xB(3) = xA(3)
                        call add_neighbor_coeff_host(coeff, xA, xB, ibm, dns, inside)
                        xB(1) = xA(1);             xB(2) = blk%y(iy-1,var,b); xB(3) = xA(3)
                        call add_neighbor_coeff_host(coeff, xA, xB, ibm, dns, inside)
                        xB(1) = xA(1);             xB(2) = blk%y(iy+1,var,b); xB(3) = xA(3)
                        call add_neighbor_coeff_host(coeff, xA, xB, ibm, dns, inside)
                        xB(1) = xA(1);             xB(2) = xA(2);             xB(3) = blk%z(iz-1,var,b)
                        call add_neighbor_coeff_host(coeff, xA, xB, ibm, dns, inside)
                        xB(1) = xA(1);             xB(2) = xA(2);             xB(3) = blk%z(iz+1,var,b)
                        call add_neighbor_coeff_host(coeff, xA, xB, ibm, dns, inside)

                        ibm%coef(ix,iy,iz,var,b) = coeff*re_inv
#endif
                    end if
                end do
            end do
        end do
        end do
        !$omp end parallel do
    end subroutine set_ibm_coeff_host

    subroutine add_neighbor_coeff_host(coeff, xA, xB, ibm, dns, inside)
        real(C_DOUBLE), intent(inout) :: coeff
        real(C_DOUBLE), intent(in) :: xA(1:3)
        real(C_DOUBLE), intent(inout) :: xB(1:3)
        type(ibm_type), intent(in) :: ibm
        type(dns_type), intent(in) :: dns
        procedure(body_indicator_i) :: inside

        real(C_DOUBLE) :: d0, d

        if (inside(xB, ibm, dns)) then
            d0 = distance3(xA, xB)
            call bisection_host(xA, xB, ibm, dns, inside)
            d = distance3(xA, xB)
            coeff = coeff + ((d0-d)/d)/d0**2
        end if
    end subroutine add_neighbor_coeff_host

    subroutine bisection_host(xAin, xB, ibm, dns, inside)
        real(C_DOUBLE), intent(in) :: xAin(1:3)
        real(C_DOUBLE), intent(inout):: xB(1:3)
        type(ibm_type), intent(in) :: ibm
        type(dns_type), intent(in) :: dns
        procedure(body_indicator_i) :: inside

        real(C_DOUBLE) :: xA(1:3), xM(1:3)
        logical :: la, lm
        integer(C_INT) :: it

        xA = xAin
        xM = xA
        la = inside(xA, ibm, dns)

        do it = 1, MAX_ITER
            xM = 0.5d0*(xA + xB)
            if (distance3(xA, xB) < DEFAULT_TOL) exit

            lm = inside(xM, ibm, dns)
            if (lm .eqv. la) then
                xA = xM
            else
                xB = xM
            end if
        end do
        xB = xM
    end subroutine bisection_host

    subroutine update_ibm_mu(ibm, dt_gamma)
        type(ibm_type), intent(inout) :: ibm
        real(C_DOUBLE), intent(in) :: dt_gamma

        integer :: ix, iy, iz, var, b, nBlocks
        integer :: ilo, ihi, jlo, jhi, klo, khi

        nBlocks = size(ibm%coef,5)
        ilo = lbound(ibm%coef,1)
        ihi = ubound(ibm%coef,1)
        jlo = lbound(ibm%coef,2)
        jhi = ubound(ibm%coef,2)
        klo = lbound(ibm%coef,3)
        khi = ubound(ibm%coef,3)

#ifdef USE_OPENMP_OFFLOAD
        !$omp target teams distribute parallel do collapse(5) &
        !$omp& map(to: dt_gamma, ilo, ihi, jlo, jhi, klo, khi, nBlocks, ibm%coef) &
        !$omp& map(tofrom: ibm%mu) &
        !$omp& private(ix,iy,iz,var,b)
#endif
        do b = 1, nBlocks
        do var = VAR_U, VAR_W
            do iz = klo, khi
                do iy = jlo, jhi
                    do ix = ilo, ihi
                        ibm%mu(ix,iy,iz,var,b) = 1.0d0/(1.0d0 + dt_gamma*ibm%coef(ix,iy,iz,var,b))
                    end do
                end do
            end do
        end do
        end do
#ifdef USE_OPENMP_OFFLOAD
        !$omp end target teams distribute parallel do
#endif
    end subroutine update_ibm_mu

end module ibmm
