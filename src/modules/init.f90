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
    integer(C_INT), parameter :: GRID_NATURAL = 4_C_INT
    integer(C_INT), parameter :: CFL_COURANT = 1_C_INT
    integer(C_INT), parameter :: CFL_PECLET  = 2_C_INT
    integer(C_INT), parameter :: NCFL = 2_C_INT

    ! Runtime/domain state shared by the solver modules.
    type :: dns_type
        integer(C_INT) :: globalSize(1:3) = 0_C_INT
        ! Per-direction local first index, last index, and count.
        integer(C_INT) :: localSize(1:3,0:2) = 0_C_INT
        integer(C_INT) :: step_current = 0_C_INT
        integer(C_INT) :: nsteps = 0_C_INT
        real(C_DOUBLE) :: leng(1:3) = 0.0d0
        real(C_DOUBLE) :: re = 0.0d0
        real(C_DOUBLE) :: dt = 0.0d0
        real(C_DOUBLE) :: t_final = 0.0d0
        real(C_DOUBLE) :: t_current = 0.0d0
        real(C_DOUBLE) :: cfl(1:NCFL) = 0.0d0
        real(C_DOUBLE) :: cflmax = 0.0d0
        real(C_DOUBLE) :: pecletmax = 0.0d0
        real(C_DOUBLE) :: peclet_rate = 0.0d0
        real(C_DOUBLE) :: dtmax = 0.0d0
        real(C_DOUBLE) :: forcing(1:3) = 0.0d0
        ! [flow] initial_u/v/w: uniform initial velocity (generic case).
        real(C_DOUBLE) :: initial_velocity(1:3) = 0.0d0
        real(C_DOUBLE) :: initial_noise = 0.0d0
        ! [blocks] nb: cubic block edge in cells; 0 = one block per rank box.
        integer(C_INT) :: block_nb = 0_C_INT
        ! [blocks] remove_solid: drop blocks buried inside the immersed body.
        logical(C_BOOL) :: block_remove_solid = .true.
        ! [blocks] refine = x0 x1 y0 y1 z0 z1: refine blocks intersecting this
        ! physical box (test option; lo > hi means unset).
        ! Up to 4 refinement boxes ([blocks] refine, repeatable key).
        real(C_DOUBLE) :: block_refine_box(6,4) = 0.0d0
        integer(C_INT) :: block_refine_nboxes = 0_C_INT
        ! [blocks] refine_levels: rounds of box refinement (max level).
        integer(C_INT) :: block_refine_levels = 1_C_INT
        ! [blocks] refine_body: refine blocks whose dilated region meets the
        ! immersed surface to the finest level (+1 block buffer), and remove
        ! buried blocks at every level (analytic IBM).
        logical(C_BOOL) :: block_refine_body = .false.
        logical(C_BOOL) :: ibm_enabled = .true.
        character(len=256) :: ibm_coeff_file = ""
        character(len=256) :: field_prefix = ""
        integer :: field_interval = 0
        character(len=256) :: restart_file = ""
    end type dns_type

    ! Grid generation parameters and the global node lines. The staggered
    ! coordinates and finite-difference metrics live per block in
    ! block_set_type, sliced from these lines (slice_grid_direction).
    type :: grid_type
        integer(C_INT) :: distribution(1:3) = GRID_UNIFORM
        ! Build the n-point line by midpoint subdivision of the (n/2)-point
        ! line generated with the same parameters - bitwise identical to
        ! one refinement level of the coarser line (blocks.f90 does the
        ! same subdivision), so a uniformly fine reference run can share
        ! its grid exactly with a refined run's fine level.
        logical(C_BOOL) :: subdivided(1:3) = .false.
        real(C_DOUBLE) :: stretch(1:3) = 0.0d0
        real(C_DOUBLE) :: natural_dyw_plus(1:3) = 0.05d0
        logical(C_BOOL) :: natural_one_sided(1:3) = .false.
        real(C_DOUBLE), allocatable :: xNode(:), yNode(:), zNode(:)
    end type grid_type

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

    call destroy_grid(g)

    allocate(g%xNode(0:int(dns%globalSize(1))))
    allocate(g%yNode(0:int(dns%globalSize(2))))
    allocate(g%zNode(0:int(dns%globalSize(3))))

    call build_node_line(g%xNode, dns%globalSize(1), dns%leng(1), &
        g%distribution(1), g%stretch(1), g%natural_one_sided(1), g%natural_dyw_plus(1), &
        g%subdivided(1))
    call build_node_line(g%yNode, dns%globalSize(2), dns%leng(2), &
        g%distribution(2), g%stretch(2), g%natural_one_sided(2), g%natural_dyw_plus(2), &
        g%subdivided(2))
    call build_node_line(g%zNode, dns%globalSize(3), dns%leng(3), &
        g%distribution(3), g%stretch(3), g%natural_one_sided(3), g%natural_dyw_plus(3), &
        g%subdivided(3))
end subroutine init_grid

subroutine destroy_grid(g)
    type(grid_type), intent(inout) :: g

    if (allocated(g%xNode)) deallocate(g%xNode)
    if (allocated(g%yNode)) deallocate(g%yNode)
    if (allocated(g%zNode)) deallocate(g%zNode)
end subroutine destroy_grid

recursive subroutine build_node_line(node, nGlobal, length, distribution, stretch, &
        natural_one_sided, natural_dyw_plus, subdivided)
    real(C_DOUBLE), intent(inout) :: node(0:)
    integer(C_INT), intent(in) :: nGlobal, distribution
    real(C_DOUBLE), intent(in) :: length, stretch, natural_dyw_plus
    logical(C_BOOL), intent(in) :: natural_one_sided
    logical(C_BOOL), intent(in), optional :: subdivided

    integer :: i, n
    real(C_DOUBLE) :: s
    real(C_DOUBLE), allocatable :: coarse(:)

    n = int(nGlobal)
    if (present(subdivided)) then
        if (subdivided) then
            ! Midpoint subdivision of the half-resolution line, exactly as
            ! blocks.f90 builds refinement-level lines.
            if (mod(n, 2) /= 0) error stop "subdivided grid needs an even point count"
            allocate(coarse(0:n/2))
            call build_node_line(coarse, int(n/2, C_INT), length, distribution, stretch, &
                natural_one_sided, natural_dyw_plus)
            do i = 0, n/2 - 1
                node(2*i) = coarse(i)
                node(2*i+1) = 0.5d0*(coarse(i) + coarse(i+1))
            end do
            node(n) = coarse(n/2)
            return
        end if
    end if
    do i = 0, n
        s = real(i, C_DOUBLE) / real(n, C_DOUBLE)
        node(i) = distribution_coordinate(s, length, distribution, stretch, &
            natural_one_sided, natural_dyw_plus, n)
    end do
    node(0) = 0.0d0
    node(n) = length
end subroutine build_node_line

! Sample the local, variable-staggered coordinates and the second-order
! finite-difference metrics for a window of nLocal cells starting at the
! 1-based global cell index `first` of a given global node line.
!
! This is shared by the rank-local grid setup (init_grid_direction above) and
! by the per-block metric setup in the blocks module, so the discrete
! operators are defined in exactly one place.
subroutine slice_grid_direction(node, coord, d1, lapM, lap0, lapP, nGlobal, first, nLocal, &
        length, periodic, dir)
    real(C_DOUBLE), intent(in) :: node(0:)
    real(C_DOUBLE), intent(inout) :: coord(-1:,:)
    real(C_DOUBLE), intent(inout) :: d1(0:,:)
    real(C_DOUBLE), intent(inout) :: lapM(0:,:), lap0(0:,:), lapP(0:,:)
    integer(C_INT), intent(in) :: nGlobal, first
    integer, intent(in) :: nLocal, dir
    real(C_DOUBLE), intent(in) :: length
    logical(C_BOOL), intent(in) :: periodic

    integer :: i, var, n, loCoord, hiCoord
    real(C_DOUBLE) :: hm, hp

    n = int(nGlobal)
    loCoord = lbound(coord,1)
    hiCoord = ubound(coord,1)

    do var = VAR_U, VAR_P
        ! Coordinates depend on both direction and variable because of staggering.
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

        ! First-derivative inverse spacings connect the opposite staggering.
        do i = 0, nLocal+1
            if (is_face_staggered(dir, var)) then
                d1(i,var) = 1.0d0 / (cell_center_at(node, n, length, int(first) + i - 1, periodic) &
                                   - cell_center_at(node, n, length, int(first) + i - 2, periodic))
            else
                d1(i,var) = 1.0d0 / (face_at(node, n, length, int(first) + i - 1, periodic) &
                                   - face_at(node, n, length, int(first) + i - 2, periodic))
            end if
        end do

        ! Three-point second-derivative stencil on nonuniform spacing.
        ! Extends to nLocal+1 (the high-side halo face) so the redundant
        ! top-face momentum computation has its diffusion coefficient
        ! lap*(nb+1); coord runs to nb+2, so the stencil closes. Index 0
        ! and the low side stay unused (the top-face stencil never reaches
        ! below v(0)).
        do i = 1, nLocal+1
            hm = coord(i,var) - coord(i-1,var)
            hp = coord(i+1,var) - coord(i,var)
            lapM(i,var) = 2.0d0 / (hm * (hm + hp))
            lapP(i,var) = 2.0d0 / (hp * (hm + hp))
            lap0(i,var) = -(lapM(i,var) + lapP(i,var))
        end do
    end do
end subroutine slice_grid_direction

logical function is_face_staggered(dir, var)
    integer, intent(in) :: dir, var

    is_face_staggered = (dir == 1 .and. var == VAR_U) .or. &
                        (dir == 2 .and. var == VAR_V) .or. &
                        (dir == 3 .and. var == VAR_W)
end function is_face_staggered

real(C_DOUBLE) function distribution_coordinate(s, length, distribution, stretch, natural_one_sided, &
        natural_dyw_plus, n) result(x)
    real(C_DOUBLE), intent(in) :: s, length, stretch, natural_dyw_plus
    integer(C_INT), intent(in) :: distribution
    logical(C_BOOL), intent(in) :: natural_one_sided
    integer, intent(in) :: n
    real(C_DOUBLE), parameter :: pi = 3.1415926535897932384626433832795d0
    real(C_DOUBLE) :: a

    select case (distribution)
    case (GRID_COSINE)
        x = 0.5d0 * length * (1.0d0 - cos(pi*s))
    case (GRID_TANH)
        a = max(stretch, 1.0d-12)
        x = 0.5d0 * length * (1.0d0 + tanh(a*(2.0d0*s - 1.0d0))/tanh(a))
    case (GRID_NATURAL)
        x = natural_channel_coordinate(s, length, stretch, natural_one_sided, natural_dyw_plus, n)
    case default
        x = length*s
    end select
end function distribution_coordinate

real(C_DOUBLE) function natural_channel_coordinate(s, length, blend_index, one_sided, dy_wall_plus, n) result(x)
    real(C_DOUBLE), intent(in) :: s, length, blend_index, dy_wall_plus
    logical(C_BOOL), intent(in) :: one_sided
    integer, intent(in) :: n

    real(C_DOUBLE) :: half_length, j, jmax, yplus, yplus_max

    if (.not. one_sided) then
        half_length = 0.5d0*length
        j = min(s, 1.0d0 - s) * real(n, C_DOUBLE)
        jmax = 0.5d0 * real(n, C_DOUBLE)
        yplus = natural_wall_coordinate(j, blend_index, dy_wall_plus)
        yplus_max = natural_wall_coordinate(jmax, blend_index, dy_wall_plus)
        if (yplus_max <= 0.0d0) then
            x = length*s
        else if (s <= 0.5d0) then
            x = half_length * yplus/yplus_max
        else
            x = length - half_length * yplus/yplus_max
        end if
    else
        j = s * real(n, C_DOUBLE)
        jmax = real(n, C_DOUBLE)
        yplus = natural_wall_coordinate(j, blend_index, dy_wall_plus)
        yplus_max = natural_wall_coordinate(jmax, blend_index, dy_wall_plus)
        if (yplus_max <= 0.0d0) then
            x = length*s
        else
            x = length * yplus/yplus_max
        end if
    end if
end function natural_channel_coordinate

real(C_DOUBLE) function natural_wall_coordinate(j, blend_index, dy_wall_plus) result(yplus)
    real(C_DOUBLE), intent(in) :: j, blend_index, dy_wall_plus

    real(C_DOUBLE), parameter :: alpha = 1.25d0
    real(C_DOUBLE), parameter :: c_eta = 0.8d0
    real(C_DOUBLE) :: jb, blend, outer, dy_wall

    if (j <= 0.0d0) then
        yplus = 0.0d0
        return
    end if

    jb = merge(blend_index, 40.0d0, blend_index > 0.0d0)
    dy_wall = merge(dy_wall_plus, 0.05d0, dy_wall_plus > 0.0d0)
    blend = (j/jb)**2
    outer = (0.75d0*alpha*c_eta*j)**(4.0d0/3.0d0)
    yplus = (dy_wall*j + outer*blend)/(1.0d0 + blend)
end function natural_wall_coordinate

real(C_DOUBLE) function face_at(node, n, length, idx, periodic) result(x)
    real(C_DOUBLE), intent(in) :: node(0:)
    integer, intent(in) :: n, idx
    real(C_DOUBLE), intent(in) :: length
    logical(C_BOOL), intent(in) :: periodic

    ! Extend the node line into halos by wrapping or mirroring.
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

end module init
