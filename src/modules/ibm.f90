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
    use :: init, only: dns_type, grid_type, VAR_U, VAR_W, VAR_P, &
        is_face_staggered, face_at, cell_center_at
    use :: blocks, only: block_set_type, subdivide_node_line
    use :: io, only: to_c_string
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

    end type ibm_type

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
        !$omp target exit data map(delete: ibm%coef, ibm%mu)
        !$omp target exit data map(delete: ibm)
#endif
    end subroutine exit_ibm_data

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


    logical function isInBody(xIN, ibm, dns)
        implicit none

        real(C_DOUBLE), intent(in) :: xIN(1:3)
        type(ibm_type), intent(in) :: ibm
        type(dns_type), intent(in) :: dns

        real(C_DOUBLE), parameter :: pi = 3.141592653589793d0
        real(C_DOUBLE), parameter :: y_offset = 1.0d-2
        real(C_DOUBLE) :: y_body

        y_body = ibm%amp_x * 0.5d0 * &
                 (1.0d0 + sin(2.0d0*pi*real(ibm%n_wave_x,C_DOUBLE)*xIN(1)/dns%leng(1) + ibm%phase_x)) + &
                 y_offset
        isInBody = (xIN(2) < y_body)
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
    subroutine classify_active_blocks(active, dns, g, ibm, periodic)
        integer(C_INT), intent(out) :: active(:)
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        type(ibm_type), intent(inout) :: ibm
        logical(C_BOOL), intent(in) :: periodic(1:3)

        integer :: nb, gnbt(3), gx, gy, gz, raster, d
        integer :: i, j, k, var, o(3)
        real(C_DOUBLE) :: xA(3)
        logical :: buried

        call set_ibm_geometry_defaults(ibm)

        nb = int(dns%block_nb)
        do d = 1, 3
            if (mod(int(dns%globalSize(d)), nb) /= 0) then
                error stop "[blocks] nb must divide the global grid in every direction"
            end if
            gnbt(d) = int(dns%globalSize(d))/nb
        end do
        if (size(active) /= product(gnbt)) error stop "block active mask size mismatch"

        raster = 0
        do gz = 0, gnbt(3) - 1
            do gy = 0, gnbt(2) - 1
                do gx = 0, gnbt(1) - 1
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
                                    if (.not. isInBody(xA, ibm, dns)) then
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
    ! `dir` on a node line of `n` cells (global grid or a refinement level),
    ! matching slice_grid_direction's staggering.
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
    subroutine classify_block_geometry(touch, buried, dns, g, ibm, periodic, nLevels)
        integer(C_INT), intent(out) :: touch(:,:), buried(:,:)
        type(dns_type), intent(in) :: dns
        type(grid_type), intent(in) :: g
        type(ibm_type), intent(inout) :: ibm
        logical(C_BOOL), intent(in) :: periodic(1:3)
        integer, intent(in) :: nLevels

        integer :: nb, l, gnbt(3), gx, gy, gz, raster, d
        integer :: i, j, k, var, o(3), nf(3), nl(3)
        real(C_DOUBLE) :: xA(3)
        real(C_DOUBLE), allocatable :: lineX(:,:), lineY(:,:), lineZ(:,:)
        logical :: anySolid, anyFluid, inside

        call set_ibm_geometry_defaults(ibm)
        nb = int(dns%block_nb)

        nf = int(dns%globalSize)*2**(nLevels - 1)
        allocate(lineX(0:nf(1), nLevels), lineY(0:nf(2), nLevels), lineZ(0:nf(3), nLevels))
        lineX(0:int(dns%globalSize(1)),1) = g%xNode
        lineY(0:int(dns%globalSize(2)),1) = g%yNode
        lineZ(0:int(dns%globalSize(3)),1) = g%zNode
        do l = 2, nLevels
            call subdivide_node_line(lineX(0:int(dns%globalSize(1))*2**(l-2), l-1), &
                                     lineX(0:int(dns%globalSize(1))*2**(l-1), l))
            call subdivide_node_line(lineY(0:int(dns%globalSize(2))*2**(l-2), l-1), &
                                     lineY(0:int(dns%globalSize(2))*2**(l-1), l))
            call subdivide_node_line(lineZ(0:int(dns%globalSize(3))*2**(l-2), l-1), &
                                     lineZ(0:int(dns%globalSize(3))*2**(l-1), l))
        end do

        do l = 1, nLevels
            nl = int(dns%globalSize)*2**(l-1)
            gnbt = nl/nb
            raster = 0
            do gz = 0, gnbt(3) - 1
                do gy = 0, gnbt(2) - 1
                    do gx = 0, gnbt(1) - 1
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
                                        inside = isInBody(xA, ibm, dns)
                                        anySolid = anySolid .or. inside
                                        anyFluid = anyFluid .or. .not. inside
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
