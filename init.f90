module init
    use, intrinsic :: iso_c_binding
    implicit none

    type :: grid_type
        integer :: nx, ny, nz
        integer :: nsteps
        real(C_DOUBLE) :: lx, ly, lz
        real(C_DOUBLE) :: dx, dy, dz
        real(C_DOUBLE) :: re, dt, t_final, t_current, cfl, cflmax, dtmax
    end type grid_type

    type :: field_type
        real(C_DOUBLE), allocatable :: un(:,:,:), us(:,:,:), oldrhsu(:,:,:)
        real(C_DOUBLE), allocatable :: vn(:,:,:), vs(:,:,:), oldrhsv(:,:,:)
        real(C_DOUBLE), allocatable :: wn(:,:,:), ws(:,:,:), oldrhsw(:,:,:)
        real(C_DOUBLE), allocatable :: pn(:,:,:), pc(:,:,:)
    end type field_type

contains

subroutine init_grid(g)
    type(grid_type), intent(inout) :: g

    g%nx = 200
    g%ny = 200
    g%nz = 200
    g%nsteps = 10000

    g%lx = 1.0d0
    g%ly = 1.0d0
    g%lz = 1.0d0
    g%re = 100.0d0

#ifdef USE_IBM
    g%dt = 1.0d-4
#else
    g%dt = 1.0d-3
#endif
    g%cflmax = 0.0d0
    g%dtmax = 1.0d-3
    g%t_current = 0.0d0
    g%cfl = 0.0d0

    call finalize_grid(g)
end subroutine init_grid

subroutine finalize_grid(g)
    type(grid_type), intent(inout) :: g

    g%dx = g%lx / real(g%nx, C_DOUBLE)
    g%dy = g%ly / real(g%ny, C_DOUBLE)
    g%dz = g%lz / real(g%nz, C_DOUBLE)
    g%t_final = g%dt * real(g%nsteps, C_DOUBLE)
end subroutine finalize_grid

subroutine init_field(f, g)
    type(field_type), intent(inout) :: f
    type(grid_type), intent(in)     :: g

    allocate(f%un(0:g%nx+1,0:g%ny+1,0:g%nz+1), f%us(0:g%nx+1,0:g%ny+1,0:g%nz+1), f%oldrhsu(1:g%nx,1:g%ny,1:g%nz))
    allocate(f%vn(0:g%nx+1,1:g%ny+1,0:g%nz+1), f%vs(0:g%nx+1,1:g%ny+1,0:g%nz+1), f%oldrhsv(1:g%nx,2:g%ny,1:g%nz))
    allocate(f%wn(0:g%nx+1,0:g%ny+1,0:g%nz+1), f%ws(0:g%nx+1,0:g%ny+1,0:g%nz+1), f%oldrhsw(1:g%nx,1:g%ny,1:g%nz))

    allocate(f%pn(0:g%nx+1,1:g%ny,0:g%nz+1), f%pc(0:g%nx+1,1:g%ny,0:g%nz+1))

    f%un = 0.0d0; f%us = 0.0d0; f%oldrhsu = 0.0d0
    f%vn = 0.0d0; f%vs = 0.0d0; f%oldrhsv = 0.0d0
    f%wn = 0.0d0; f%ws = 0.0d0; f%oldrhsw = 0.0d0
    f%pn = 0.0d0; f%pc = 0.0d0
end subroutine init_field

end module init
