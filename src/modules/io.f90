module io
    use, intrinsic :: iso_c_binding
    use :: init, only: grid_type, field_type
    implicit none

    interface
        function fdm_h5_write_field(file_name, nx, ny, nz, lx, ly, lz, dx, dy, dz, re, dt, t_current, &
                un, vn, wn, pn) bind(C, name="fdm_h5_write_field") result(ierr)
            import :: C_CHAR, C_INT, C_DOUBLE
            character(kind=C_CHAR), intent(in) :: file_name(*)
            integer(C_INT), value :: nx, ny, nz
            real(C_DOUBLE), value :: lx, ly, lz, dx, dy, dz, re, dt, t_current
            real(C_DOUBLE), intent(in) :: un(*), vn(*), wn(*), pn(*)
            integer(C_INT) :: ierr
        end function fdm_h5_write_field

        function fdm_h5_read_field(file_name, nx, ny, nz, un, vn, wn, pn) &
                bind(C, name="fdm_h5_read_field") result(ierr)
            import :: C_CHAR, C_INT, C_DOUBLE
            character(kind=C_CHAR), intent(in) :: file_name(*)
            integer(C_INT), value :: nx, ny, nz
            real(C_DOUBLE), intent(out) :: un(*), vn(*), wn(*), pn(*)
            integer(C_INT) :: ierr
        end function fdm_h5_read_field
    end interface

    type :: output_workspace_type
        real(C_DOUBLE), allocatable :: uc(:,:,:), vc(:,:,:), wc(:,:,:)
    end type output_workspace_type

contains

logical function output_is_due(step, output_interval)
    integer, intent(in) :: step, output_interval

    if (output_interval <= 0) then
        output_is_due = .false.
    else
        output_is_due = modulo(step, output_interval) == 0
    end if
end function output_is_due

subroutine init_output_workspace(out, g, output_interval)
    type(output_workspace_type), intent(inout) :: out
    type(grid_type), intent(in) :: g
    integer, intent(in) :: output_interval

    if (output_interval <= 0) return

    allocate(out%uc(1:g%nx,1:g%ny,1:g%nz))
    allocate(out%vc(1:g%nx,1:g%ny,1:g%nz))
    allocate(out%wc(1:g%nx,1:g%ny,1:g%nz))

    out%uc = 0.0d0
    out%vc = 0.0d0
    out%wc = 0.0d0

#ifdef USE_OPENMP_OFFLOAD
    !$omp target enter data map(to: out)
    !$omp target enter data map(alloc: &
    !$omp& out%uc(1:g%nx,1:g%ny,1:g%nz), &
    !$omp& out%vc(1:g%nx,1:g%ny,1:g%nz), &
    !$omp& out%wc(1:g%nx,1:g%ny,1:g%nz))
#endif
end subroutine init_output_workspace

subroutine destroy_output_workspace(out, g)
    type(output_workspace_type), intent(inout) :: out
    type(grid_type), intent(in) :: g

    if (.not. allocated(out%uc)) return

#ifdef USE_OPENMP_OFFLOAD
    !$omp target exit data map(delete: &
    !$omp& out%uc(1:g%nx,1:g%ny,1:g%nz), &
    !$omp& out%vc(1:g%nx,1:g%ny,1:g%nz), &
    !$omp& out%wc(1:g%nx,1:g%ny,1:g%nz))
    !$omp target exit data map(delete: out)
#endif

    deallocate(out%uc, out%vc, out%wc)
end subroutine destroy_output_workspace

subroutine maybe_write_vtk(out, f, g, step, output_interval, output_prefix)
    type(output_workspace_type), intent(inout) :: out
    type(field_type), intent(inout) :: f
    type(grid_type), intent(in) :: g
    integer, intent(in) :: step, output_interval
    character(len=*), intent(in) :: output_prefix

    if (.not. output_is_due(step, output_interval)) return
    call write_vtk(out, f, g, step, output_prefix)
end subroutine maybe_write_vtk

subroutine maybe_write_field(f, g, step, field_interval, field_prefix)
    type(field_type), intent(inout) :: f
    type(grid_type), intent(in) :: g
    integer, intent(in) :: step, field_interval
    character(len=*), intent(in) :: field_prefix

    if (.not. output_is_due(step, field_interval)) return
    call write_field(f, g, step, field_prefix)
end subroutine maybe_write_field

subroutine write_field(f, g, step, field_prefix)
    type(field_type), intent(inout) :: f
    type(grid_type), intent(in) :: g
    integer, intent(in) :: step
    character(len=*), intent(in) :: field_prefix

    character(len=256) :: file_name
    character(kind=C_CHAR,len=:), allocatable :: c_file_name
    integer(C_INT) :: ierr
    integer :: nx, ny, nz

    nx = g%nx
    ny = g%ny
    nz = g%nz

    write(file_name,'(A,"_",I0,".h5")') trim(field_prefix), step
    print *, "current time step: ", step, "   field filename: ", trim(file_name)

#ifdef USE_OPENMP_OFFLOAD
    !$omp target update from(f%un(0:nx+1,0:ny+1,0:nz+1), &
    !$omp& f%vn(0:nx+1,0:ny+1,0:nz+1), &
    !$omp& f%wn(0:nx+1,0:ny+1,0:nz+1), &
    !$omp& f%pn(0:nx+1,0:ny+1,0:nz+1))
#endif

    c_file_name = to_c_string(file_name)
    ierr = fdm_h5_write_field(c_file_name, int(nx, C_INT), int(ny, C_INT), int(nz, C_INT), &
        g%lx, g%ly, g%lz, g%dx, g%dy, g%dz, g%re, g%dt, g%t_current, f%un, f%vn, f%wn, f%pn)
    if (ierr /= 0_C_INT) then
        print *, "error: could not write HDF5 field file: ", trim(file_name)
        error stop
    end if
end subroutine write_field

subroutine read_field(f, g, file_name)
    type(field_type), intent(inout) :: f
    type(grid_type), intent(in) :: g
    character(len=*), intent(in) :: file_name

    character(kind=C_CHAR,len=:), allocatable :: c_file_name
    integer(C_INT) :: ierr
    integer :: nx, ny, nz

    nx = g%nx
    ny = g%ny
    nz = g%nz

    c_file_name = to_c_string(file_name)
    ierr = fdm_h5_read_field(c_file_name, int(nx, C_INT), int(ny, C_INT), int(nz, C_INT), &
        f%un, f%vn, f%wn, f%pn)
    if (ierr /= 0_C_INT) then
        print *, "error: could not read HDF5 field file: ", trim(file_name)
        error stop
    end if

#ifdef USE_OPENMP_OFFLOAD
    !$omp target update to(f%un(0:nx+1,0:ny+1,0:nz+1), &
    !$omp& f%vn(0:nx+1,0:ny+1,0:nz+1), &
    !$omp& f%wn(0:nx+1,0:ny+1,0:nz+1), &
    !$omp& f%pn(0:nx+1,0:ny+1,0:nz+1))
#endif
end subroutine read_field

subroutine write_vtk(out, f, g, step, output_prefix)
    type(output_workspace_type), intent(inout) :: out
    type(field_type), intent(inout) :: f
    type(grid_type), intent(in) :: g
    integer, intent(in) :: step
    character(len=*), intent(in) :: output_prefix

    character(len=256) :: file_name

    write(file_name,'(A,"_",I0,".vtk")') trim(output_prefix), step
    print *, "current time step: ", step, "   filename: ", trim(file_name), "   cfl:", g%cfl*g%dt

    call center_vel(out, f, g)
    call data_output(out, g, file_name)
end subroutine write_vtk

subroutine center_vel(out, f, g)
    type(output_workspace_type), intent(inout) :: out
    type(field_type), intent(inout) :: f
    type(grid_type), intent(in) :: g

    integer :: i, j, k, ip, jp, kp
    integer :: nx, ny, nz

    nx = g%nx
    ny = g%ny
    nz = g%nz

    !$omp target teams distribute parallel do collapse(3) &
    !$omp& map(to: f%un(0:nx+1,0:ny+1,0:nz+1), &
    !$omp& f%vn(0:nx+1,0:ny+1,0:nz+1), &
    !$omp& f%wn(0:nx+1,0:ny+1,0:nz+1)) &
    !$omp& map(tofrom: out%uc(1:nx,1:ny,1:nz), &
    !$omp& out%vc(1:nx,1:ny,1:nz), out%wc(1:nx,1:ny,1:nz)) &
    !$omp& private(i,j,k,ip,jp,kp)
    do i = 1, nx
        do j = 1, ny
            do k = 1, nz
                ip = i + 1
                jp = j + 1
                kp = k + 1

                out%uc(i,j,k) = 0.5d0 * (f%un(i,j,k) + f%un(ip,j,k))
                out%vc(i,j,k) = 0.5d0 * (f%vn(i,j,k) + f%vn(i,jp,k))
                out%wc(i,j,k) = 0.5d0 * (f%wn(i,j,k) + f%wn(i,j,kp))
            end do
        end do
    end do
    !$omp end target teams distribute parallel do

#ifdef USE_OPENMP_OFFLOAD
    !$omp target update from(out%uc(1:nx,1:ny,1:nz), &
    !$omp& out%vc(1:nx,1:ny,1:nz), out%wc(1:nx,1:ny,1:nz))
#endif
end subroutine center_vel

subroutine data_output(out, g, file_name)
    type(output_workspace_type), intent(in) :: out
    type(grid_type), intent(in) :: g
    character(len=*), intent(in) :: file_name

    integer :: io, npts, i, j, k

    npts = g%nx*g%ny*g%nz

    open(newunit=io, file=trim(file_name), status="replace", action="write")
    write(io,'(A)') "# vtk DataFile Version 3.0"
    write(io,'(A)') "3D velocity field"
    write(io,'(A)') "ASCII"
    write(io,'(A)') "DATASET STRUCTURED_POINTS"

    write(io,'(A,1X,I0,1X,I0,1X,I0)') "DIMENSIONS", g%nx+1, g%ny+1, g%nz+1
    write(io,'(A,1X,3ES20.12)') "ORIGIN", 0.0d0, 0.0d0, 0.0d0
    write(io,'(A,1X,3ES20.12)') "SPACING", g%dx, g%dy, g%dz

    write(io,'(A)') "FIELD FieldData 1"
    write(io,'(A)') "TIME 1 1 double"
    write(io,'(ES20.12)') g%t_current

    write(io,'(A,1X,I0)') "CELL_DATA", npts
    write(io,'(A)') "VECTORS U double"

    do k = 1, g%nz
        do j = 1, g%ny
            do i = 1, g%nx
                write(io,'(3ES20.12)') out%uc(i,j,k), out%vc(i,j,k), out%wc(i,j,k)
            end do
        end do
    end do

    close(io)
end subroutine data_output

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
