module init
    use, intrinsic :: iso_c_binding
    implicit none

    integer(C_INT), parameter :: VAR_U = 1_C_INT
    integer(C_INT), parameter :: VAR_V = 2_C_INT
    integer(C_INT), parameter :: VAR_W = 3_C_INT
    integer(C_INT), parameter :: VAR_P = 4_C_INT
    integer(C_INT), parameter :: NVEL = 3_C_INT
    integer(C_INT), parameter :: NVAR = 4_C_INT
    integer(C_INT), parameter :: GRID_UNIFORM = 1_C_INT
    integer(C_INT), parameter :: GRID_COSINE  = 2_C_INT
    integer(C_INT), parameter :: GRID_TANH    = 3_C_INT

    type :: dns_type
        integer(C_INT) :: globalSize(1:3) = 0_C_INT
        integer(C_INT) :: localSize(1:3,0:2) = 0_C_INT
        integer(C_INT) :: mpiDims(1:3) = 0_C_INT
        integer(C_INT) :: gridDistribution(1:3) = GRID_UNIFORM
        integer(C_INT) :: step_current = 0_C_INT
        integer(C_INT) :: nsteps = 0_C_INT
        real(C_DOUBLE) :: leng(1:3) = 0.0d0
        real(C_DOUBLE) :: gridStretch(1:3) = 0.0d0
        real(C_DOUBLE) :: re = 0.0d0
        real(C_DOUBLE) :: dt = 0.0d0
        real(C_DOUBLE) :: t_final = 0.0d0
        real(C_DOUBLE) :: t_current = 0.0d0
        real(C_DOUBLE) :: cfl = 0.0d0
        real(C_DOUBLE) :: cflmax = 0.0d0
        real(C_DOUBLE) :: dtmax = 0.0d0
        real(C_DOUBLE) :: forcing(1:3) = 0.0d0
        logical(C_BOOL) :: ibm_enabled = .true.
        character(len=256) :: field_prefix = ""
        integer :: field_interval = 0
        character(len=256) :: restart_file = ""
    end type dns_type

    type :: grid_type
        real(C_DOUBLE) :: havg(1:3) = 0.0d0
        real(C_DOUBLE), allocatable :: xNode(:), yNode(:), zNode(:)
        real(C_DOUBLE), allocatable :: x(:,:), y(:,:), z(:,:)
        real(C_DOUBLE), allocatable :: d1x(:,:), d1y(:,:), d1z(:,:)
        real(C_DOUBLE), allocatable :: lapXm(:,:), lapX0(:,:), lapXp(:,:)
        real(C_DOUBLE), allocatable :: lapYm(:,:), lapY0(:,:), lapYp(:,:)
        real(C_DOUBLE), allocatable :: lapZm(:,:), lapZ0(:,:), lapZp(:,:)
    end type grid_type

    type :: field_type
        real(C_DOUBLE), allocatable :: q(:,:,:,:)      ! current u/v/w/p
        real(C_DOUBLE), allocatable :: qs(:,:,:,:)     ! tentative u/v/w
        real(C_DOUBLE), allocatable :: oldrhs(:,:,:,:) ! RK history for u/v/w
    end type field_type

contains

subroutine splash(has_terminal)
  logical, intent(in), optional :: has_terminal
  logical :: terminal

  terminal = .true.
  if (present(has_terminal)) terminal = has_terminal
  if (.not. terminal) return

  write(*,'(A)') "      __. - ~ ~ ~ - ."
  write(*,'(A)') "_   ,//           __  ' ,"
  write(*,'(A)') " \\  ||       __--  --    ,   ~~~~"
  write(*,'(A)') " , \\|\____---    o   \    ~~~    ~~~~"
  write(*,'(A)') ",   \ _            __/   ~~ ,  ~~~               mobyDiff"          
  write(*,'(A)') ",       \---/ / __--   ~~   ,~~                  commit: 7aa1c7b"
  write(*,'(A)') " ,          \/       ~~   ~~"
  write(*,'(A)') "  ,         ~~~ ~~~     ~~,"
  write(*,'(A)') "    ,    ~~~           , '"
  write(*,'(A)') "      ' - , _ _ _ ,  '"
  write(*,'(A)') ""
  write(*,'(A)') ""
end subroutine splash

subroutine init_grid(g, dns, periodic)
    type(grid_type), intent(inout) :: g
    type(dns_type), intent(inout)  :: dns
    logical(C_BOOL), intent(in)    :: periodic(1:3)

    integer :: nx, ny, nz

    nx = int(dns%localSize(1,2))
    ny = int(dns%localSize(2,2))
    nz = int(dns%localSize(3,2))

    call destroy_grid(g)

    g%havg(1) = dns%leng(1) / real(dns%globalSize(1), C_DOUBLE)
    g%havg(2) = dns%leng(2) / real(dns%globalSize(2), C_DOUBLE)
    g%havg(3) = dns%leng(3) / real(dns%globalSize(3), C_DOUBLE)

    allocate(g%xNode(0:int(dns%globalSize(1))))
    allocate(g%yNode(0:int(dns%globalSize(2))))
    allocate(g%zNode(0:int(dns%globalSize(3))))
    allocate(g%x(-1:nx+2,NVAR), g%d1x(0:nx+1,NVAR))
    allocate(g%y(-1:ny+2,NVAR), g%d1y(0:ny+1,NVAR))
    allocate(g%z(-1:nz+2,NVAR), g%d1z(0:nz+1,NVAR))
    allocate(g%lapXm(0:nx+1,NVAR), g%lapX0(0:nx+1,NVAR), g%lapXp(0:nx+1,NVAR))
    allocate(g%lapYm(0:ny+1,NVAR), g%lapY0(0:ny+1,NVAR), g%lapYp(0:ny+1,NVAR))
    allocate(g%lapZm(0:nz+1,NVAR), g%lapZ0(0:nz+1,NVAR), g%lapZp(0:nz+1,NVAR))

    call init_grid_direction(g%xNode, g%x, g%d1x, g%lapXm, g%lapX0, g%lapXp, &
        dns%globalSize(1), dns%localSize(1,0), nx, dns%leng(1), &
        dns%gridDistribution(1), dns%gridStretch(1), periodic(1), 1)
    call init_grid_direction(g%yNode, g%y, g%d1y, g%lapYm, g%lapY0, g%lapYp, &
        dns%globalSize(2), dns%localSize(2,0), ny, dns%leng(2), &
        dns%gridDistribution(2), dns%gridStretch(2), periodic(2), 2)
    call init_grid_direction(g%zNode, g%z, g%d1z, g%lapZm, g%lapZ0, g%lapZp, &
        dns%globalSize(3), dns%localSize(3,0), nz, dns%leng(3), &
        dns%gridDistribution(3), dns%gridStretch(3), periodic(3), 3)

    if (dns%nsteps <= 0_C_INT) then
        dns%nsteps = max(1_C_INT, int(ceiling((dns%t_final - dns%t_current)/dns%dt), C_INT))
    end if
    dns%t_final = dns%t_current + dns%dt * real(dns%nsteps, C_DOUBLE)
end subroutine init_grid

subroutine destroy_grid(g)
    type(grid_type), intent(inout) :: g

    if (allocated(g%xNode)) deallocate(g%xNode)
    if (allocated(g%yNode)) deallocate(g%yNode)
    if (allocated(g%zNode)) deallocate(g%zNode)
    if (allocated(g%x)) deallocate(g%x)
    if (allocated(g%y)) deallocate(g%y)
    if (allocated(g%z)) deallocate(g%z)
    if (allocated(g%d1x)) deallocate(g%d1x)
    if (allocated(g%d1y)) deallocate(g%d1y)
    if (allocated(g%d1z)) deallocate(g%d1z)
    if (allocated(g%lapXm)) deallocate(g%lapXm)
    if (allocated(g%lapX0)) deallocate(g%lapX0)
    if (allocated(g%lapXp)) deallocate(g%lapXp)
    if (allocated(g%lapYm)) deallocate(g%lapYm)
    if (allocated(g%lapY0)) deallocate(g%lapY0)
    if (allocated(g%lapYp)) deallocate(g%lapYp)
    if (allocated(g%lapZm)) deallocate(g%lapZm)
    if (allocated(g%lapZ0)) deallocate(g%lapZ0)
    if (allocated(g%lapZp)) deallocate(g%lapZp)
end subroutine destroy_grid

subroutine init_grid_direction(node, coord, d1, lapM, lap0, lapP, nGlobal, first, nLocal, &
        length, distribution, stretch, periodic, dir)
    real(C_DOUBLE), intent(inout) :: node(0:)
    real(C_DOUBLE), intent(inout) :: coord(-1:,:)
    real(C_DOUBLE), intent(inout) :: d1(0:,:)
    real(C_DOUBLE), intent(inout) :: lapM(0:,:), lap0(0:,:), lapP(0:,:)
    integer(C_INT), intent(in) :: nGlobal, first, distribution
    integer, intent(in) :: nLocal, dir
    real(C_DOUBLE), intent(in) :: length, stretch
    logical(C_BOOL), intent(in) :: periodic

    integer :: i, var, n, loCoord, hiCoord
    real(C_DOUBLE) :: s, hm, hp

    n = int(nGlobal)
    do i = 0, n
        s = real(i, C_DOUBLE) / real(n, C_DOUBLE)
        node(i) = distribution_coordinate(s, length, distribution, stretch)
    end do
    node(0) = 0.0d0
    node(n) = length


    loCoord = lbound(coord,1)
    hiCoord = ubound(coord,1)

    do var = VAR_U, VAR_P
        do i = loCoord, hiCoord
            if (is_face_staggered(dir, var)) then
                coord(i,var) = face_at(node, n, length, int(first) + i - 2, periodic)
            else
                coord(i,var) = cell_center_at(node, n, length, int(first) + i - 1, periodic)
            end if
        end do

        d1(:,var) = 0.0d0
        lapM(:,var) = 0.0d0
        lap0(:,var) = 0.0d0
        lapP(:,var) = 0.0d0

        do i = 0, nLocal+1
            if (is_face_staggered(dir, var)) then
                d1(i,var) = 1.0d0 / (cell_center_at(node, n, length, int(first) + i - 1, periodic) &
                                   - cell_center_at(node, n, length, int(first) + i - 2, periodic))
            else
                d1(i,var) = 1.0d0 / (face_at(node, n, length, int(first) + i - 1, periodic) &
                                   - face_at(node, n, length, int(first) + i - 2, periodic))
            end if
        end do

        do i = 1, nLocal
            hm = coord(i,var) - coord(i-1,var)
            hp = coord(i+1,var) - coord(i,var)
            lapM(i,var) = 2.0d0 / (hm * (hm + hp))
            lapP(i,var) = 2.0d0 / (hp * (hm + hp))
            lap0(i,var) = -(lapM(i,var) + lapP(i,var))
        end do
    end do
end subroutine init_grid_direction

logical function is_face_staggered(dir, var)
    integer, intent(in) :: dir, var

    is_face_staggered = (dir == 1 .and. var == VAR_U) .or. &
                        (dir == 2 .and. var == VAR_V) .or. &
                        (dir == 3 .and. var == VAR_W)
end function is_face_staggered

real(C_DOUBLE) function distribution_coordinate(s, length, distribution, stretch) result(x)
    real(C_DOUBLE), intent(in) :: s, length, stretch
    integer(C_INT), intent(in) :: distribution
    real(C_DOUBLE), parameter :: pi = 3.1415926535897932384626433832795d0
    real(C_DOUBLE) :: a

    select case (distribution)
    case (GRID_COSINE)
        x = 0.5d0 * length * (1.0d0 - cos(pi*s))
    case (GRID_TANH)
        a = max(stretch, 1.0d-12)
        x = 0.5d0 * length * (1.0d0 + tanh(a*(2.0d0*s - 1.0d0))/tanh(a))
    case default
        x = length*s
    end select
end function distribution_coordinate

real(C_DOUBLE) function face_at(node, n, length, idx, periodic) result(x)
    real(C_DOUBLE), intent(in) :: node(0:)
    integer, intent(in) :: n, idx
    real(C_DOUBLE), intent(in) :: length
    logical(C_BOOL), intent(in) :: periodic

    if (idx < 0) then
        if (periodic) then
            x = node(idx + n) - length
        else
            x = 2.0d0*node(0) - node(-idx)
        end if
    else if (idx > n) then
        if (periodic) then
            x = node(idx - n) + length
        else
            x = 2.0d0*node(n) - node(2*n - idx)
        end if
    else
        x = node(idx)
    end if
end function face_at

real(C_DOUBLE) function cell_center_at(node, n, length, idx, periodic) result(x)
    real(C_DOUBLE), intent(in) :: node(0:)
    integer, intent(in) :: n, idx
    real(C_DOUBLE), intent(in) :: length
    logical(C_BOOL), intent(in) :: periodic

    x = 0.5d0 * (face_at(node, n, length, idx - 1, periodic) + &
                 face_at(node, n, length, idx, periodic))
end function cell_center_at

subroutine init_field(f, dns)
    type(field_type), intent(inout) :: f
    type(dns_type), intent(in)      :: dns
    integer :: nx, ny, nz

    nx = int(dns%localSize(1,2))
    ny = int(dns%localSize(2,2))
    nz = int(dns%localSize(3,2))

    allocate(f%q(0:nx+1,0:ny+1,0:nz+1,NVAR))
    allocate(f%qs(0:nx+1,0:ny+1,0:nz+1,NVEL))
    allocate(f%oldrhs(1:nx,1:ny,1:nz,NVEL))

    f%q = 0.0d0
    f%qs = 0.0d0
    f%oldrhs = 0.0d0
end subroutine init_field

end module init
