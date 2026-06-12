module io
    use, intrinsic :: iso_c_binding
    use :: init, only: dns_type, grid_type, VAR_U, VAR_V, VAR_W, VAR_P
    use :: blocks, only: block_set_type
    use :: boundary, only: boundary_type, NFACES
    use :: comm, only: comm_type
    implicit none

    interface
        function fdm_h5_write_field(file_name, nbx, nby, nbz, n_blocks, &
                n_blocks_global, id_start, block_origin, block_level, &
                rank, nranks, global_nx, global_ny, global_nz, &
                step, nsteps, lx, ly, lz, &
                re, dt, t_final, t_current, cfl, cflmax, pecletmax, dtmax, &
                forcing, pressure_niter, pressure_sor, &
                ibm_enabled, bc_count, periodic, bc_type, bc_value, grid_distribution, grid_stretch, &
                grid_natural_dyw_plus, x_node, y_node, z_node, q) &
                bind(C, name="fdm_h5_write_field") result(ierr)
            import :: C_CHAR, C_INT, C_DOUBLE
            character(kind=C_CHAR), intent(in) :: file_name(*)
            integer(C_INT), value :: nbx, nby, nbz, n_blocks, n_blocks_global, id_start
            integer(C_INT), value :: rank, nranks
            integer(C_INT), intent(in) :: block_origin(*), block_level(*)
            integer(C_INT), value :: global_nx, global_ny, global_nz
            integer(C_INT), value :: step, nsteps
            real(C_DOUBLE), value :: lx, ly, lz, re, dt, t_final, t_current
            real(C_DOUBLE), value :: cflmax, pecletmax, dtmax
            real(C_DOUBLE), intent(in) :: cfl(*)
            real(C_DOUBLE), intent(in) :: forcing(*)
            integer(C_INT), value :: pressure_niter
            real(C_DOUBLE), value :: pressure_sor
            integer(C_INT), value :: ibm_enabled
            integer(C_INT), value :: bc_count
            integer(C_INT), intent(in) :: periodic(*), bc_type(*), grid_distribution(*)
            real(C_DOUBLE), intent(in) :: bc_value(*), grid_stretch(*), grid_natural_dyw_plus(*)
            real(C_DOUBLE), intent(in) :: x_node(*), y_node(*), z_node(*)
            real(C_DOUBLE), intent(in) :: q(*)
            integer(C_INT) :: ierr
        end function fdm_h5_write_field

        function fdm_h5_write_grid(file_name, nx, ny, nz, lx, ly, lz, &
                periodic, grid_distribution, grid_stretch, grid_natural_dyw_plus, grid_natural_one_sided, &
                x_node, y_node, z_node, xu, yu, zu, xv, yv, zv, xw, yw, zw) &
                bind(C, name="fdm_h5_write_grid") result(ierr)
            import :: C_CHAR, C_INT, C_DOUBLE
            character(kind=C_CHAR), intent(in) :: file_name(*)
            integer(C_INT), value :: nx, ny, nz
            real(C_DOUBLE), value :: lx, ly, lz
            integer(C_INT), intent(in) :: periodic(*), grid_distribution(*), grid_natural_one_sided(*)
            real(C_DOUBLE), intent(in) :: grid_stretch(*), grid_natural_dyw_plus(*)
            real(C_DOUBLE), intent(in) :: x_node(*), y_node(*), z_node(*)
            real(C_DOUBLE), intent(in) :: xu(*), yu(*), zu(*)
            real(C_DOUBLE), intent(in) :: xv(*), yv(*), zv(*)
            real(C_DOUBLE), intent(in) :: xw(*), yw(*), zw(*)
            integer(C_INT) :: ierr
        end function fdm_h5_write_grid

        function fdm_h5_read_metadata(file_name, global_nx, global_ny, global_nz, &
                step, nsteps, lx, ly, lz, re, dt, t_final, t_current, cfl, &
                cflmax, pecletmax, dtmax, forcing, &
                pressure_niter, pressure_sor, ibm_enabled, bc_count, periodic, bc_type, bc_value, &
                grid_distribution, grid_stretch, grid_natural_dyw_plus) &
                bind(C, name="fdm_h5_read_metadata") result(ierr)
            import :: C_CHAR, C_INT, C_DOUBLE
            character(kind=C_CHAR), intent(in) :: file_name(*)
            integer(C_INT), intent(inout) :: global_nx, global_ny, global_nz
            integer(C_INT), intent(inout) :: step, nsteps
            real(C_DOUBLE), intent(inout) :: lx, ly, lz, re, dt, t_final, t_current
            real(C_DOUBLE), intent(inout) :: cfl(*)
            real(C_DOUBLE), intent(inout) :: cflmax, pecletmax, dtmax
            real(C_DOUBLE), intent(inout) :: forcing(*)
            integer(C_INT), intent(inout) :: pressure_niter
            real(C_DOUBLE), intent(inout) :: pressure_sor
            integer(C_INT), intent(inout) :: ibm_enabled
            integer(C_INT), value :: bc_count
            integer(C_INT), intent(inout) :: periodic(*), bc_type(*), grid_distribution(*)
            real(C_DOUBLE), intent(inout) :: bc_value(*), grid_stretch(*), grid_natural_dyw_plus(*)
            integer(C_INT) :: ierr
        end function fdm_h5_read_metadata

        function fdm_h5_read_block_active(file_name, n_lattice, block_nb, found, active) &
                bind(C, name="fdm_h5_read_block_active") result(ierr)
            import :: C_CHAR, C_INT
            character(kind=C_CHAR), intent(in) :: file_name(*)
            integer(C_INT), value :: n_lattice, block_nb
            integer(C_INT), intent(out) :: found
            integer(C_INT), intent(out) :: active(*)
            integer(C_INT) :: ierr
        end function fdm_h5_read_block_active

        function fdm_h5_read_field(file_name, nbx, nby, nbz, n_blocks, &
                n_blocks_global, id_start, block_origin, &
                global_nx, global_ny, global_nz, q) &
                bind(C, name="fdm_h5_read_field") result(ierr)
            import :: C_CHAR, C_INT, C_DOUBLE
            character(kind=C_CHAR), intent(in) :: file_name(*)
            integer(C_INT), value :: nbx, nby, nbz, n_blocks, n_blocks_global, id_start
            integer(C_INT), intent(in) :: block_origin(*)
            integer(C_INT), value :: global_nx, global_ny, global_nz
            real(C_DOUBLE), intent(out) :: q(*)
            integer(C_INT) :: ierr
        end function fdm_h5_read_field
    end interface

contains

logical function output_is_due(step, output_interval)
    integer, intent(in) :: step, output_interval

    if (output_interval <= 0) then
        output_is_due = .false.
    else
        output_is_due = modulo(step, output_interval) == 0
    end if
end function output_is_due

subroutine maybe_write_field(blk, dns, g, step, c, bc, pressure_niter, pressure_sor)
    type(block_set_type), intent(inout) :: blk
    type(dns_type), intent(in) :: dns
    type(grid_type), intent(in) :: g
    integer, intent(in) :: step
    type(comm_type), intent(in) :: c
    type(boundary_type), intent(in) :: bc
    integer(C_INT), intent(in) :: pressure_niter
    real(C_DOUBLE), intent(in) :: pressure_sor

    if (.not. output_is_due(step, dns%field_interval)) return
    call write_field(blk, dns, g, step, c, bc, pressure_niter, pressure_sor)
end subroutine maybe_write_field

subroutine write_field(blk, dns, g, step, c, bc, pressure_niter, pressure_sor)
    ! Parallel HDF5 call: all MPI ranks must enter this routine together.
    ! Global datasets, one hyperslab per block.
    type(block_set_type), intent(inout) :: blk
    type(dns_type), intent(in) :: dns
    type(grid_type), intent(in) :: g
    integer, intent(in) :: step
    type(comm_type), intent(in) :: c
    type(boundary_type), intent(in) :: bc
    integer(C_INT), intent(in) :: pressure_niter
    real(C_DOUBLE), intent(in) :: pressure_sor

    character(len=256) :: h5_file_name, xdmf_file_name
    character(kind=C_CHAR,len=:), allocatable :: c_file_name
    integer(C_INT) :: ierr
    integer(C_INT) :: periodic(1:3)
    integer(C_INT) :: ibm_enabled

    periodic = merge(1_C_INT, 0_C_INT, bc%isPeriodic)
    ibm_enabled = merge(1_C_INT, 0_C_INT, dns%ibm_enabled)

    call make_output_filename(trim(dns%field_prefix), step, ".h5", h5_file_name)
    call make_output_filename(trim(dns%field_prefix), step, ".xdmf", xdmf_file_name)
    if (c%has_terminal) then
        print *, "current time step: ", step, "   field filename: ", trim(h5_file_name)
    end if

#ifdef USE_OPENMP_OFFLOAD
    !$omp target update from(blk%q)
#endif

    c_file_name = to_c_string(h5_file_name)
    ierr = fdm_h5_write_field(c_file_name, blk%nb(1), blk%nb(2), blk%nb(3), &
        blk%nBlocks, blk%nBlocksGlobal, blk%idStart, blk%origin, blk%level, &
        int(c%world_rank, C_INT), int(c%world_size, C_INT), &
        dns%globalSize(1), dns%globalSize(2), dns%globalSize(3), &
        dns%step_current, dns%nsteps, &
        dns%leng(1), dns%leng(2), dns%leng(3), &
        dns%re, dns%dt, dns%t_final, dns%t_current, &
        dns%cfl, dns%cflmax, dns%pecletmax, dns%dtmax, &
        dns%forcing, &
        pressure_niter, pressure_sor, ibm_enabled, int(size(bc%faceBcType), C_INT), periodic, &
        bc%faceBcType(VAR_U:VAR_P,1:NFACES), bc%faceBcDefaultValue(VAR_U:VAR_P,1:NFACES), &
        g%distribution(1:3), g%stretch(1:3), g%natural_dyw_plus(1:3), &
        g%xNode(0:dns%globalSize(1)), g%yNode(0:dns%globalSize(2)), g%zNode(0:dns%globalSize(3)), &
        blk%q)
    if (ierr /= 0_C_INT) then
        if (c%has_terminal) print *, "error: could not write HDF5 field file: ", trim(h5_file_name)
        error stop
    end if

    ! No XDMF for the block-table layout; reassemble with
    ! tools/compare_fields.py --export-global for visualization.
end subroutine write_field

subroutine write_grid_export(dns, g, blk, bc, file_name, has_terminal)
    ! Serial preprocessing path (mobygrid): the single block covers the whole
    ! grid, so its slot-1 staggered coordinates are the global ones.
    type(dns_type), intent(in) :: dns
    type(grid_type), intent(in) :: g
    type(block_set_type), intent(in) :: blk
    type(boundary_type), intent(in) :: bc
    character(len=*), intent(in) :: file_name
    logical, intent(in), optional :: has_terminal

    character(kind=C_CHAR,len=:), allocatable :: c_file_name
    integer(C_INT) :: ierr
    integer(C_INT) :: nx, ny, nz
    integer(C_INT) :: periodic(1:3)
    integer(C_INT) :: natural_one_sided(1:3)
    logical :: terminal

    terminal = .true.
    if (present(has_terminal)) terminal = has_terminal

    nx = dns%globalSize(1)
    ny = dns%globalSize(2)
    nz = dns%globalSize(3)
    periodic = merge(1_C_INT, 0_C_INT, bc%isPeriodic)
    natural_one_sided = merge(1_C_INT, 0_C_INT, g%natural_one_sided)

    c_file_name = to_c_string(file_name)
    ierr = fdm_h5_write_grid(c_file_name, nx, ny, nz, dns%leng(1), dns%leng(2), dns%leng(3), &
        periodic, g%distribution(1:3), g%stretch(1:3), g%natural_dyw_plus(1:3), natural_one_sided, &
        g%xNode(0:int(nx)), g%yNode(0:int(ny)), g%zNode(0:int(nz)), &
        blk%x(0:int(nx)+1,VAR_U,1), blk%y(0:int(ny)+1,VAR_U,1), blk%z(0:int(nz)+1,VAR_U,1), &
        blk%x(0:int(nx)+1,VAR_V,1), blk%y(0:int(ny)+1,VAR_V,1), blk%z(0:int(nz)+1,VAR_V,1), &
        blk%x(0:int(nx)+1,VAR_W,1), blk%y(0:int(ny)+1,VAR_W,1), blk%z(0:int(nz)+1,VAR_W,1))
    if (ierr /= 0_C_INT) then
        if (terminal) print *, "error: could not write HDF5 grid file: ", trim(file_name)
        error stop
    end if
end subroutine write_grid_export

subroutine read_restart_metadata(dns, g, bc, pressure_niter, pressure_sor, file_name, c, &
        preserve_cflmax, preserve_pecletmax, preserve_dtmax, preserve_t_final)
    type(dns_type), intent(inout) :: dns
    type(grid_type), intent(inout) :: g
    type(boundary_type), intent(inout) :: bc
    integer(C_INT), intent(inout) :: pressure_niter
    real(C_DOUBLE), intent(inout) :: pressure_sor
    character(len=*), intent(in) :: file_name
    type(comm_type), intent(in) :: c
    logical, intent(in), optional :: preserve_cflmax, preserve_pecletmax, preserve_dtmax, preserve_t_final

    character(kind=C_CHAR,len=:), allocatable :: c_file_name
    integer(C_INT) :: ierr
    integer(C_INT) :: file_nsteps
    integer(C_INT) :: periodic(1:3)
    integer(C_INT) :: ibm_enabled
    integer :: dir
    real(C_DOUBLE) :: input_cflmax, input_pecletmax, input_dtmax, input_t_final

    file_nsteps = dns%nsteps
    periodic = merge(1_C_INT, 0_C_INT, bc%isPeriodic)
    ibm_enabled = merge(1_C_INT, 0_C_INT, dns%ibm_enabled)
    input_cflmax = dns%cflmax
    input_pecletmax = dns%pecletmax
    input_dtmax = dns%dtmax
    input_t_final = dns%t_final
    c_file_name = to_c_string(file_name)
    ierr = fdm_h5_read_metadata(c_file_name, dns%globalSize(1), dns%globalSize(2), dns%globalSize(3), &
        dns%step_current, file_nsteps, &
        dns%leng(1), dns%leng(2), dns%leng(3), dns%re, dns%dt, dns%t_final, dns%t_current, dns%cfl, &
        dns%cflmax, dns%pecletmax, dns%dtmax, dns%forcing, &
        pressure_niter, pressure_sor, ibm_enabled, int(size(bc%faceBcType), C_INT), periodic, &
        bc%faceBcType(VAR_U:VAR_P,1:NFACES), bc%faceBcDefaultValue(VAR_U:VAR_P,1:NFACES), &
        g%distribution(1:3), g%stretch(1:3), g%natural_dyw_plus(1:3))
    if (ierr /= 0_C_INT) then
        if (c%has_terminal) print *, "error: could not read restart metadata: ", trim(file_name)
        error stop
    end if

    if (present(preserve_cflmax)) then
        if (preserve_cflmax) dns%cflmax = input_cflmax
    end if
    if (present(preserve_pecletmax)) then
        if (preserve_pecletmax) dns%pecletmax = input_pecletmax
    end if
    if (present(preserve_dtmax)) then
        if (preserve_dtmax) dns%dtmax = input_dtmax
    end if
    if (present(preserve_t_final)) then
        if (preserve_t_final) dns%t_final = input_t_final
    end if

    if (dns%nsteps <= 0_C_INT) dns%nsteps = file_nsteps
    do dir = 1, 3
        bc%isPeriodic(dir) = periodic(dir) /= 0_C_INT
    end do
    dns%ibm_enabled = ibm_enabled /= 0_C_INT
end subroutine read_restart_metadata

! Per-block keep flags from the IBM coefficient file (mobygeom
! block-active). found is false when the file carries no table.
subroutine read_block_active(active, found, dns, has_terminal)
    integer(C_INT), intent(out) :: active(:)
    logical, intent(out) :: found
    type(dns_type), intent(in) :: dns
    logical, intent(in) :: has_terminal

    character(kind=C_CHAR,len=:), allocatable :: c_file_name
    integer(C_INT) :: ierr, c_found

    c_file_name = to_c_string(dns%ibm_coeff_file)
    ierr = fdm_h5_read_block_active(c_file_name, int(size(active), C_INT), &
        dns%block_nb, c_found, active)
    if (ierr /= 0_C_INT) then
        if (has_terminal) print *, "error: could not read block_active from: ", &
            trim(dns%ibm_coeff_file)
        error stop
    end if
    found = c_found /= 0_C_INT
end subroutine read_block_active

subroutine read_field(blk, dns, file_name, c)
    ! Parallel HDF5 call: all MPI ranks must enter this routine together.
    type(block_set_type), intent(inout) :: blk
    type(dns_type), intent(in) :: dns
    character(len=*), intent(in) :: file_name
    type(comm_type), intent(in) :: c

    character(kind=C_CHAR,len=:), allocatable :: c_file_name
    integer(C_INT) :: ierr

    c_file_name = to_c_string(file_name)
    ierr = fdm_h5_read_field(c_file_name, blk%nb(1), blk%nb(2), blk%nb(3), &
        blk%nBlocks, blk%nBlocksGlobal, blk%idStart, blk%origin, &
        dns%globalSize(1), dns%globalSize(2), dns%globalSize(3), &
        blk%q)
    if (ierr /= 0_C_INT) then
        if (c%has_terminal) print *, "error: could not read HDF5 field file: ", trim(file_name)
        error stop
    end if

#ifdef USE_OPENMP_OFFLOAD
    !$omp target update to(blk%q)
#endif
end subroutine read_field

subroutine write_xdmf(xdmf_file_name, h5_file_name, dns, g)
    character(len=*), intent(in) :: xdmf_file_name, h5_file_name
    type(dns_type), intent(in) :: dns
    type(grid_type), intent(in) :: g

    character(len=256) :: h5_ref
    integer :: io, slash
    integer(C_INT) :: nx, ny, nz

    nx = dns%globalSize(1)
    ny = dns%globalSize(2)
    nz = dns%globalSize(3)

    h5_ref = trim(h5_file_name)
    slash = scan(trim(h5_ref), "/", back=.true.)
    if (slash > 0) h5_ref = h5_ref(slash+1:)
    h5_ref = "./"//trim(h5_ref)

    open(newunit=io, file=trim(xdmf_file_name), status="replace", action="write")
    write(io,'(A)') '<?xml version="1.0" ?>'
    write(io,'(A)') '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
    write(io,'(A)') '<Xdmf Version="2.0">'
    write(io,'(A)') '  <Domain>'
    write(io,'(A)') '    <Grid Name="mobyDiff" GridType="Uniform">'
    write(io,'(A,ES20.12,A)') '      <Time Value="', dns%t_current, '"/>'
    write(io,'(A,I0,1X,I0,1X,I0,A)') &
        '      <Topology TopologyType="3DRectMesh" Dimensions="', nz+1_C_INT, ny+1_C_INT, nx+1_C_INT, '"/>'
    write(io,'(A)') '      <Geometry GeometryType="VXVYVZ">'
    call write_xdmf_coord(io, h5_ref, "/x", nx+1_C_INT)
    call write_xdmf_coord(io, h5_ref, "/y", ny+1_C_INT)
    call write_xdmf_coord(io, h5_ref, "/z", nz+1_C_INT)
    write(io,'(A)') '      </Geometry>'
    call write_xdmf_scalar(io, "un", h5_ref, "/un", nx, ny, nz)
    call write_xdmf_scalar(io, "vn", h5_ref, "/vn", nx, ny, nz)
    call write_xdmf_scalar(io, "wn", h5_ref, "/wn", nx, ny, nz)
    call write_xdmf_scalar(io, "pn", h5_ref, "/pn", nx, ny, nz)
    write(io,'(A)') '    </Grid>'
    write(io,'(A)') '  </Domain>'
    write(io,'(A)') '</Xdmf>'
    close(io)
end subroutine write_xdmf

subroutine write_xdmf_coord(io, h5_ref, dataset, n)
    integer, intent(in) :: io
    character(len=*), intent(in) :: h5_ref, dataset
    integer(C_INT), intent(in) :: n

    write(io,'(A,I0,A,A,A,A,A)') &
        '        <DataItem Dimensions="', n, &
        '" NumberType="Float" Precision="8" Format="HDF">', &
        trim(h5_ref), ':', trim(dataset), '</DataItem>'
end subroutine write_xdmf_coord

subroutine write_xdmf_scalar(io, name, h5_ref, dataset, nx, ny, nz)
    integer, intent(in) :: io
    character(len=*), intent(in) :: name, h5_ref, dataset
    integer(C_INT), intent(in) :: nx, ny, nz

    write(io,'(A,A,A)') '      <Attribute Name="', trim(name), '" AttributeType="Scalar" Center="Cell">'
    write(io,'(A,I0,1X,I0,1X,I0,A,A,A,A,A)') &
        '        <DataItem Dimensions="', nz, ny, nx, &
        '" NumberType="Float" Precision="8" Format="HDF">', &
        trim(h5_ref), ':', trim(dataset), '</DataItem>'
    write(io,'(A)') '      </Attribute>'
end subroutine write_xdmf_scalar


subroutine make_output_filename(prefix, step, extension, file_name)
    character(len=*), intent(in) :: prefix, extension
    integer, intent(in) :: step
    character(len=*), intent(out) :: file_name

    write(file_name,'(A,"_",I0,A)') trim(prefix), step, trim(extension)
end subroutine make_output_filename

function to_c_string(text) result(c_text)
    character(len=*), intent(in) :: text
    character(kind=C_CHAR,len=:), allocatable :: c_text
    integer :: i, n

    n = len_trim(text)
    allocate(character(kind=C_CHAR,len=n+1) :: c_text)
    do i = 1, n
        c_text(i:i) = text(i:i)
    end do
    c_text(n+1:n+1) = C_NULL_CHAR
end function to_c_string

end module io
