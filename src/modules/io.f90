module io
    use, intrinsic :: iso_c_binding
    use :: init, only: dns_type, grid_type, field_type, VAR_U, VAR_V, VAR_W, VAR_P
    use :: boundary, only: boundary_type, NFACES
    use :: comm, only: comm_type
    implicit none

    interface
        function fdm_h5_write_field(file_name, nx, ny, nz, rank, nranks, &
                global_nx, global_ny, global_nz, &
                local_i_first, local_i_last, local_j_first, local_j_last, &
                local_k_first, local_k_last, step, nsteps, lx, ly, lz, &
                re, dt, t_final, t_current, cfl, cflmax, pecletmax, dtmax, &
                forcing, pressure_niter, pressure_sor, &
                ibm_enabled, bc_count, periodic, bc_type, bc_value, grid_distribution, grid_stretch, &
                x_node, y_node, z_node, un, vn, wn, pn) &
                bind(C, name="fdm_h5_write_field") result(ierr)
            import :: C_CHAR, C_INT, C_DOUBLE
            character(kind=C_CHAR), intent(in) :: file_name(*)
            integer(C_INT), value :: nx, ny, nz, rank, nranks
            integer(C_INT), value :: global_nx, global_ny, global_nz
            integer(C_INT), value :: local_i_first, local_i_last
            integer(C_INT), value :: local_j_first, local_j_last
            integer(C_INT), value :: local_k_first, local_k_last
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
            real(C_DOUBLE), intent(in) :: bc_value(*), grid_stretch(*), x_node(*), y_node(*), z_node(*)
            real(C_DOUBLE), intent(in) :: un(*), vn(*), wn(*), pn(*)
            integer(C_INT) :: ierr
        end function fdm_h5_write_field

        function fdm_h5_write_grid(file_name, nx, ny, nz, lx, ly, lz, &
                periodic, grid_distribution, grid_stretch, grid_natural_one_sided, &
                x_node, y_node, z_node, xu, yu, zu, xv, yv, zv, xw, yw, zw) &
                bind(C, name="fdm_h5_write_grid") result(ierr)
            import :: C_CHAR, C_INT, C_DOUBLE
            character(kind=C_CHAR), intent(in) :: file_name(*)
            integer(C_INT), value :: nx, ny, nz
            real(C_DOUBLE), value :: lx, ly, lz
            integer(C_INT), intent(in) :: periodic(*), grid_distribution(*), grid_natural_one_sided(*)
            real(C_DOUBLE), intent(in) :: grid_stretch(*)
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
                grid_distribution, grid_stretch) &
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
            real(C_DOUBLE), intent(inout) :: bc_value(*), grid_stretch(*)
            integer(C_INT) :: ierr
        end function fdm_h5_read_metadata

        function fdm_h5_read_field(file_name, nx, ny, nz, &
                global_nx, global_ny, global_nz, &
                local_i_first, local_j_first, local_k_first, &
                un, vn, wn, pn) &
                bind(C, name="fdm_h5_read_field") result(ierr)
            import :: C_CHAR, C_INT, C_DOUBLE
            character(kind=C_CHAR), intent(in) :: file_name(*)
            integer(C_INT), value :: nx, ny, nz
            integer(C_INT), value :: global_nx, global_ny, global_nz
            integer(C_INT), value :: local_i_first, local_j_first, local_k_first
            real(C_DOUBLE), intent(out) :: un(*), vn(*), wn(*), pn(*)
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

subroutine maybe_write_field(f, dns, g, step, c, bc, pressure_niter, pressure_sor)
    type(field_type), intent(inout) :: f
    type(dns_type), intent(in) :: dns
    type(grid_type), intent(in) :: g
    integer, intent(in) :: step
    type(comm_type), intent(in) :: c
    type(boundary_type), intent(in) :: bc
    integer(C_INT), intent(in) :: pressure_niter
    real(C_DOUBLE), intent(in) :: pressure_sor

    if (.not. output_is_due(step, dns%field_interval)) return
    call write_field(f, dns, g, step, c, bc, pressure_niter, pressure_sor)
end subroutine maybe_write_field

subroutine write_field(f, dns, g, step, c, bc, pressure_niter, pressure_sor)
    ! Parallel HDF5 call: all MPI ranks must enter this routine together.
    type(field_type), intent(inout) :: f
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
    integer :: nx, ny, nz
    integer(C_INT) :: periodic(1:3)
    integer(C_INT) :: ibm_enabled

    nx = int(dns%localSize(1,2))
    ny = int(dns%localSize(2,2))
    nz = int(dns%localSize(3,2))
    periodic = merge(1_C_INT, 0_C_INT, bc%isPeriodic)
    ibm_enabled = merge(1_C_INT, 0_C_INT, dns%ibm_enabled)

    call make_output_filename(trim(dns%field_prefix), step, ".h5", h5_file_name)
    call make_output_filename(trim(dns%field_prefix), step, ".xdmf", xdmf_file_name)
    if (c%has_terminal) then
        print *, "current time step: ", step, "   field filename: ", trim(h5_file_name)
    end if

#ifdef USE_OPENMP_OFFLOAD
    !$omp target update from(f%q)
#endif

    c_file_name = to_c_string(h5_file_name)
    ierr = fdm_h5_write_field(c_file_name, int(nx, C_INT), int(ny, C_INT), int(nz, C_INT), &
        int(c%world_rank, C_INT), int(c%world_size, C_INT), &
        dns%globalSize(1), dns%globalSize(2), dns%globalSize(3), &
        dns%localSize(1,0), dns%localSize(1,1), dns%localSize(2,0), dns%localSize(2,1), &
        dns%localSize(3,0), dns%localSize(3,1), &
        dns%step_current, dns%nsteps, &
        dns%leng(1), dns%leng(2), dns%leng(3), &
        dns%re, dns%dt, dns%t_final, dns%t_current, &
        dns%cfl, dns%cflmax, dns%pecletmax, dns%dtmax, &
        dns%forcing, &
        pressure_niter, pressure_sor, ibm_enabled, int(size(bc%faceBcType), C_INT), periodic, &
        bc%faceBcType(VAR_U:VAR_P,1:NFACES), bc%faceBcDefaultValue(VAR_U:VAR_P,1:NFACES), &
        g%distribution(1:3), g%stretch(1:3), &
        g%xNode(0:dns%globalSize(1)), g%yNode(0:dns%globalSize(2)), g%zNode(0:dns%globalSize(3)), &
        f%q(0:nx+1,0:ny+1,0:nz+1,VAR_U), &
        f%q(0:nx+1,0:ny+1,0:nz+1,VAR_V), &
        f%q(0:nx+1,0:ny+1,0:nz+1,VAR_W), &
        f%q(0:nx+1,0:ny+1,0:nz+1,VAR_P))
    if (ierr /= 0_C_INT) then
        if (c%has_terminal) print *, "error: could not write HDF5 field file: ", trim(h5_file_name)
        error stop
    end if

    if (c%has_terminal) call write_xdmf(xdmf_file_name, h5_file_name, dns, g)
end subroutine write_field

subroutine write_grid_export(dns, g, bc, file_name, has_terminal)
    type(dns_type), intent(in) :: dns
    type(grid_type), intent(in) :: g
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
        periodic, g%distribution(1:3), g%stretch(1:3), natural_one_sided, &
        g%xNode(0:int(nx)), g%yNode(0:int(ny)), g%zNode(0:int(nz)), &
        g%x(0:int(nx)+1,VAR_U), g%y(0:int(ny)+1,VAR_U), g%z(0:int(nz)+1,VAR_U), &
        g%x(0:int(nx)+1,VAR_V), g%y(0:int(ny)+1,VAR_V), g%z(0:int(nz)+1,VAR_V), &
        g%x(0:int(nx)+1,VAR_W), g%y(0:int(ny)+1,VAR_W), g%z(0:int(nz)+1,VAR_W))
    if (ierr /= 0_C_INT) then
        if (terminal) print *, "error: could not write HDF5 grid file: ", trim(file_name)
        error stop
    end if
end subroutine write_grid_export

subroutine read_restart_metadata(dns, g, bc, pressure_niter, pressure_sor, file_name, c)
    type(dns_type), intent(inout) :: dns
    type(grid_type), intent(inout) :: g
    type(boundary_type), intent(inout) :: bc
    integer(C_INT), intent(inout) :: pressure_niter
    real(C_DOUBLE), intent(inout) :: pressure_sor
    character(len=*), intent(in) :: file_name
    type(comm_type), intent(in) :: c

    character(kind=C_CHAR,len=:), allocatable :: c_file_name
    integer(C_INT) :: ierr
    integer(C_INT) :: file_nsteps
    integer(C_INT) :: periodic(1:3)
    integer(C_INT) :: ibm_enabled
    integer :: dir

    file_nsteps = dns%nsteps
    periodic = merge(1_C_INT, 0_C_INT, bc%isPeriodic)
    ibm_enabled = merge(1_C_INT, 0_C_INT, dns%ibm_enabled)
    c_file_name = to_c_string(file_name)
    ierr = fdm_h5_read_metadata(c_file_name, dns%globalSize(1), dns%globalSize(2), dns%globalSize(3), &
        dns%step_current, file_nsteps, &
        dns%leng(1), dns%leng(2), dns%leng(3), dns%re, dns%dt, dns%t_final, dns%t_current, dns%cfl, &
        dns%cflmax, dns%pecletmax, dns%dtmax, dns%forcing, &
        pressure_niter, pressure_sor, ibm_enabled, int(size(bc%faceBcType), C_INT), periodic, &
        bc%faceBcType(VAR_U:VAR_P,1:NFACES), bc%faceBcDefaultValue(VAR_U:VAR_P,1:NFACES), &
        g%distribution(1:3), g%stretch(1:3))
    if (ierr /= 0_C_INT) then
        if (c%has_terminal) print *, "error: could not read restart metadata: ", trim(file_name)
        error stop
    end if

    if (dns%nsteps <= 0_C_INT) dns%nsteps = file_nsteps
    do dir = 1, 3
        bc%isPeriodic(dir) = periodic(dir) /= 0_C_INT
    end do
    dns%ibm_enabled = ibm_enabled /= 0_C_INT
end subroutine read_restart_metadata

subroutine read_field(f, dns, file_name, c)
    ! Parallel HDF5 call: all MPI ranks must enter this routine together.
    type(field_type), intent(inout) :: f
    type(dns_type), intent(in) :: dns
    character(len=*), intent(in) :: file_name
    type(comm_type), intent(in) :: c

    character(kind=C_CHAR,len=:), allocatable :: c_file_name
    integer(C_INT) :: ierr
    integer :: nx, ny, nz

    nx = int(dns%localSize(1,2))
    ny = int(dns%localSize(2,2))
    nz = int(dns%localSize(3,2))

    c_file_name = to_c_string(file_name)
    ierr = fdm_h5_read_field(c_file_name, int(nx, C_INT), int(ny, C_INT), int(nz, C_INT), &
        dns%globalSize(1), dns%globalSize(2), dns%globalSize(3), &
        dns%localSize(1,0), dns%localSize(2,0), dns%localSize(3,0), &
        f%q(0:nx+1,0:ny+1,0:nz+1,VAR_U), &
        f%q(0:nx+1,0:ny+1,0:nz+1,VAR_V), &
        f%q(0:nx+1,0:ny+1,0:nz+1,VAR_W), &
        f%q(0:nx+1,0:ny+1,0:nz+1,VAR_P))
    if (ierr /= 0_C_INT) then
        if (c%has_terminal) print *, "error: could not read HDF5 field file: ", trim(file_name)
        error stop
    end if

#ifdef USE_OPENMP_OFFLOAD
    !$omp target update to(f%q)
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
