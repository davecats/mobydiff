module boundary
    use, intrinsic :: iso_c_binding
    use :: init, only: dns_type, VAR_U, VAR_V, VAR_W, VAR_P, VAR_S0, SOLID_FACE_THRESHOLD
    use :: blocks, only: block_set_type, FACE_PHYS
    implicit none

    integer(C_INT), parameter :: DIR_X = 1_C_INT
    integer(C_INT), parameter :: DIR_Y = 2_C_INT
    integer(C_INT), parameter :: DIR_Z = 3_C_INT
    integer(C_INT), parameter :: SIDE_MIN = 0_C_INT
    integer(C_INT), parameter :: SIDE_MAX = 1_C_INT
    integer, parameter :: NFACES = 6

    ! Domain-face patch types ([boundary] <dir>_<side>_patch = wall | patch |
    ! inlet | outlet). Meaningful on non-periodic faces only. UNSET (no
    ! declaration) falls back to the historical inference in
    ! domain_face_is_wall, so existing inis run unchanged. The patch type is
    ! the ONE user-facing face concept: resolve_face_bcs derives the
    ! per-variable BC rows of every declared face from it.
    integer(C_INT), parameter :: PATCH_UNSET = 0_C_INT
    integer(C_INT), parameter :: PATCH_GENERIC = 1_C_INT
    integer(C_INT), parameter :: PATCH_WALL = 2_C_INT
    integer(C_INT), parameter :: PATCH_INLET = 3_C_INT
    integer(C_INT), parameter :: PATCH_OUTLET = 4_C_INT

    ! Per-variable BC types (faceBcType values). DIRICHLET/NEUMANN come from
    ! the ini (_type keys); OUTFLOW is INTERNAL-only, derived by
    ! resolve_face_bcs for the normal velocity of an outlet face. Such a face
    ! is an UNKNOWN, not a boundary value: the momentum predictor advances
    ! it (step.f90 predict_outlet_faces), the pressure projection corrects
    ! it against the held outlet pressure, and it has NO boundary row.
    integer(C_INT), parameter :: BC_DIRICHLET = 0_C_INT
    integer(C_INT), parameter :: BC_NEUMANN = 1_C_INT
    integer(C_INT), parameter :: BC_OUTFLOW = 2_C_INT

    ! Boundary-value profiles ([boundary] <dir>_<side>_<var>_profile).
    ! CONSTANT (default) keeps the historical face-uniform value; PARABOLA
    ! scales it by prod over the face's NON-PERIODIC tangential directions of
    ! 4*s*(1-s), s = coord/leng (Poiseuille-type inlet; assumes the domain
    ! starts at 0). BLASIUS is the laminar flat-plate boundary layer, wall
    ! at y = 0, similarity variable eta = beta*y/theta(x) with theta =
    ! [boundary] blasius_theta the inlet momentum thickness and
    ! beta = 2 f''(0) the Blasius momentum-thickness constant. On an x FACE
    ! (inlet, station x = 0): u = value*f'(eta), v = value*beta*(eta f' -
    ! f)/(2 Re_theta) (the entrainment velocity). On a y FACE (top
    ! displacement BC): the v component only, the x-varying entrainment at
    ! the domain top y = leng(2), station x + x_v (x_v = Re_theta*theta/
    ! beta^2). The _value key of every blasius row is U_inf and
    ! Re_theta = U_inf*theta/nu from [flow] re. Evaluated once into
    ! pointBcValue at each variable's own staggered coordinate.
    integer(C_INT), parameter :: PROFILE_CONSTANT = 0_C_INT
    integer(C_INT), parameter :: PROFILE_PARABOLA = 1_C_INT
    integer(C_INT), parameter :: PROFILE_BLASIUS = 2_C_INT

    ! Per-face ghost modes for the generic cell-centred scalar BC applicator
    ! (apply_scalar_bc). MIRROR gives a zero face value (Dirichlet 0), COPY a
    ! zero normal gradient, VALUE a prescribed face value; NONE leaves the
    ! ghost untouched.
    integer(C_INT), parameter :: SCALAR_BC_NONE = 0_C_INT
    integer(C_INT), parameter :: SCALAR_BC_COPY = 1_C_INT
    integer(C_INT), parameter :: SCALAR_BC_MIRROR = 2_C_INT
    integer(C_INT), parameter :: SCALAR_BC_VALUE = 3_C_INT

    ! Resolved boundary-row kinds (bcKind). Every physical-face write of the
    ! solver is ONE affine operation, dst = w*src + C, and the BC TYPE is
    ! resolved at init into which cell pair (dst, src) it acts on, a
    ! per-(variable, face) weight w and a per-(variable, point) constant C:
    !   GHOST        tangential velocity, pressure, scalars: dst = the ghost
    !                cell, src = the adjacent interior cell.
    !                Dirichlet  w = -1, C = 2 v    (face midpoint carries v)
    !                Neumann    w = +1, C = dn v   (v = the normal derivative)
    !   FACE         normal velocity, Dirichlet: dst = the face dof, src =
    !                the neighbouring interior face, w = 0, C = v (the pin).
    !   NONE         no row. The normal velocity of an OUTLET is one: it is
    !                predicted and corrected like an interior face (see
    !                BC_OUTFLOW; the faces are flagged in faceOutlet).
    ! The kernels then carry no geometry and no type branches, and the
    ! velocity, pressure and scalar rows are columns of the same tables
    ! (scalar.f90 fills its columns through set_scalar_bc_rows).
    integer(C_INT), parameter :: BCK_NONE = 0_C_INT
    integer(C_INT), parameter :: BCK_GHOST = 1_C_INT
    integer(C_INT), parameter :: BCK_FACE = 2_C_INT

    type :: boundary_type
        logical(C_BOOL) :: isPeriodic(1:3)
        integer(C_INT) :: nTotal = 0_C_INT

        ! Active physical boundary points, block-local with the owning block
        ! slot. Periodic and block/MPI halos are handled by comm.f90.
        integer(C_INT), allocatable :: pointFace(:), slot(:), i(:), j(:), k(:)

        ! Face defaults seed the pointwise values (host-only, profiles
        ! evaluated per point); the kernels read the RESOLVED rows below.
        integer(C_INT) :: faceBcType(VAR_U:VAR_P,1:NFACES) = 0_C_INT
        real(C_DOUBLE) :: faceBcDefaultValue(VAR_U:VAR_P,1:NFACES) = 0.0d0
        real(C_DOUBLE), allocatable :: pointBcValue(:,:)

        ! Resolved affine rows (see BCK_*), one column per q variable
        ! 1..nVar (u, v, w, p, then the passive scalars at VAR_S0+is).
        integer(C_INT) :: nVar = 0_C_INT
        integer(C_INT), allocatable :: bcKind(:,:)   ! (nVar, NFACES)
        real(C_DOUBLE), allocatable :: bcW(:,:)      ! (nVar, NFACES)
        real(C_DOUBLE), allocatable :: bcC(:,:)      ! (nVar, nTotal)

        ! The declared OUTLET faces (1 = the face-normal velocity there is
        ! predicted, step.f90 predict_outlet_faces) and how many of this
        ! rank's boundary points lie on one. nOutlet = 0 on a rank, or in a
        ! case, without an outlet face -- nothing is launched then.
        integer(C_INT) :: faceOutlet(1:NFACES) = 0_C_INT
        integer(C_INT) :: nOutlet = 0_C_INT

        ! Which rows the ini set explicitly (_type/_value keys). Config is
        ! authority: read_restart_metadata keeps these rows over the restart
        ! file's, resolve_face_bcs seeds only unset rows of declared faces and
        ! hard-errors when an explicit type contradicts the declared patch.
        logical :: faceBcTypeSet(VAR_U:VAR_P,1:NFACES) = .false.
        logical :: faceBcValueSet(VAR_U:VAR_P,1:NFACES) = .false.

        ! Value profile per row (PROFILE_*), applied when seeding pointBcValue.
        integer(C_INT) :: faceBcProfile(VAR_U:VAR_P,1:NFACES) = PROFILE_CONSTANT

        ! Declared patch type per face (index = boundary_face_id).
        integer(C_INT) :: facePatchType(1:NFACES) = PATCH_UNSET

        ! Blasius profile state: inlet momentum thickness ([boundary]
        ! blasius_theta) and the f/f' similarity table (uniform eta step
        ! blasH, shooting-solved on first use; blasBeta = 2 f''(0) is the
        ! momentum-thickness constant, blasDisp = eta_max - f(eta_max) the
        ! displacement constant used beyond the table).
        real(C_DOUBLE) :: blasiusTheta = 1.0d0
        real(C_DOUBLE) :: blasH = 0.0d0, blasBeta = 0.0d0, blasDisp = 0.0d0
        real(C_DOUBLE), allocatable :: blasF(:), blasFp(:)
    end type boundary_type

contains

    subroutine init_bc(bc)
        type(boundary_type), intent(inout) :: bc

        call destroy_boundary_faces(bc)
        bc%isPeriodic(1:3) = .true.
        bc%faceBcType = BC_DIRICHLET
        bc%faceBcDefaultValue = 0.0d0
        bc%faceBcType(VAR_P,:) = BC_NEUMANN
        bc%faceBcTypeSet = .false.
        bc%faceBcValueSet = .false.
        bc%faceBcProfile = PROFILE_CONSTANT
        bc%facePatchType = PATCH_UNSET
    end subroutine init_bc

    ! A domain face is a wall iff its direction is non-periodic and the face
    ! is a wall patch: an explicit [boundary] <dir>_<side>_patch declaration
    ! wins; when absent, the historical inference applies (Dirichlet on both
    ! tangential velocity components = no-slip; Neumann tangential faces --
    ! free slip, symmetry -- carry no wall layer). NOTE the inference reads a
    ! Dirichlet velocity INLET as a wall; declare such a face `patch`.
    logical function domain_face_is_wall(bc, dir, side) result(is_wall)
        type(boundary_type), intent(in) :: bc
        integer, intent(in) :: dir, side

        integer :: face_id, var
        integer, parameter :: DIRICHLET = 0

        is_wall = .false.
        if (bc%isPeriodic(dir)) return

        face_id = boundary_face_id(dir, side)
        select case (bc%facePatchType(face_id))
        case (PATCH_WALL)
            is_wall = .true.
        case (PATCH_GENERIC, PATCH_INLET, PATCH_OUTLET)
            is_wall = .false.
        case default
            is_wall = .true.
            do var = int(VAR_U), int(VAR_W)
                if (var == dir) cycle   ! the normal component does not decide no-slip
                if (bc%faceBcType(var, face_id) /= DIRICHLET) is_wall = .false.
            end do
        end select
    end function domain_face_is_wall

    ! The declared patch type is the single face concept: derive the
    ! per-variable BC rows of every declared face, set-if-unset (explicit
    ! _type keys win when they agree; a direct contradiction with the
    ! declaration is a hard config error, so the consistency constraints are
    ! enforced by construction). PATCH_GENERIC and UNSET faces keep the raw
    ! per-variable keys untouched -- existing inis are bit-exact by
    ! construction. Runs with validate_patch_types: after parsing, with
    ! periodic_* final, BEFORE the first update_boundary_values.
    subroutine resolve_face_bcs(bc)
        type(boundary_type), intent(inout) :: bc
        integer :: dir, side, var, face_id

        do dir = 1, 3
            do side = 0, 1
                face_id = boundary_face_id(dir, side)
                select case (bc%facePatchType(face_id))
                case (PATCH_WALL, PATCH_INLET)
                    ! Velocity Dirichlet (wall: default 0 = no-slip; inlet:
                    ! the _value keys or the case supply the values),
                    ! pressure Neumann.
                    do var = int(VAR_U), int(VAR_W)
                        call resolve_bc_row(bc, var, face_id, BC_DIRICHLET)
                    end do
                    call resolve_bc_row(bc, int(VAR_P), face_id, BC_NEUMANN)
                case (PATCH_OUTLET)
                    ! Normal velocity: the internal outflow type (no boundary
                    ! row: the predictor and the projection own the face).
                    ! OUTFLOW cannot be written in an ini, so ANY explicit
                    ! normal type key contradicts the declaration. Tangential
                    ! Neumann 0, pressure Dirichlet (default 0 -- the held
                    ! outlet pressure).
                    do var = int(VAR_U), int(VAR_W)
                        call resolve_bc_row(bc, var, face_id, &
                            merge(BC_OUTFLOW, BC_NEUMANN, var == dir))
                    end do
                    call resolve_bc_row(bc, int(VAR_P), face_id, BC_DIRICHLET)
                end select
                ! A Neumann condition on the NORMAL velocity component is not a
                ! boundary condition the projection can honour: a row that
                ! rewrites the face from the corrected interior un-does the
                ! correction the SPD operator assumed there (numerics review
                ! F4, 2026-09-27). The outflow patch is the supported way to
                ! let flow leave -- the face is predicted, then corrected
                ! against the Dirichlet pressure -- so point there.
                if (.not. bc%isPeriodic(dir)) then
                    if (bc%faceBcType(dir, face_id) == BC_NEUMANN) then
                        print '(a,i0,a,i0,a)', " error: [boundary] Neumann type on the NORMAL velocity" // &
                            " component of face ", face_id, " (direction ", dir, &
                            "): declare the face `<dir>_<side>_patch = outlet` instead"
                        error stop "[boundary] Neumann normal velocity: use the outlet patch"
                    end if
                end if
            end do
        end do
    end subroutine resolve_face_bcs

    subroutine resolve_bc_row(bc, var, face_id, want)
        type(boundary_type), intent(inout) :: bc
        integer, intent(in) :: var, face_id
        integer(C_INT), intent(in) :: want

        if (bc%faceBcTypeSet(var, face_id)) then
            if (bc%faceBcType(var, face_id) /= want) then
                print '(a,i0,a,i0)', " error: explicit [boundary] _type key contradicts the " // &
                    "declared patch type on face ", face_id, ", variable ", var
                error stop "[boundary] BC type contradicts the declared patch type"
            end if
        else
            bc%faceBcType(var, face_id) = want
        end if
    end subroutine resolve_bc_row

    ! [boundary] <dir>_<side>_patch is meaningful on non-periodic faces only;
    ! declaring one on a periodic direction is a config error (checked here,
    ! after both the patch keys and the periodic_* flags are final).
    subroutine validate_patch_types(bc)
        type(boundary_type), intent(in) :: bc
        integer :: dir, side

        do dir = 1, 3
            if (.not. bc%isPeriodic(dir)) cycle
            do side = 0, 1
                if (bc%facePatchType(boundary_face_id(dir, side)) /= PATCH_UNSET) then
                    error stop "[boundary] patch type declared on a periodic direction"
                end if
            end do
        end do
    end subroutine validate_patch_types

    subroutine init_boundary_faces(bc, blk, dns)
        type(boundary_type), intent(inout) :: bc
        type(block_set_type), intent(in) :: blk
        type(dns_type), intent(in) :: dns
        integer :: nx, ny, nz
        integer :: b, dir, side, face_id, pos, total

        call validate_patch_types(bc)
        call resolve_face_bcs(bc)

        nx = int(blk%nb(1))
        ny = int(blk%nb(2))
        nz = int(blk%nb(3))

        call destroy_boundary_faces(bc)

        total = 0
        do b = 1, int(blk%nBlocks)
            do dir = 1, 3
                do side = 0, 1
                    if (block_face_is_physical(blk, b, dir, side)) then
                        total = total + boundary_face_n(dir, nx, ny, nz)
                    end if
                end do
            end do
        end do

        bc%nTotal = int(total, C_INT)
        if (total <= 0) return

        allocate(bc%pointFace(total), bc%slot(total), bc%i(total), bc%j(total), bc%k(total))
        allocate(bc%pointBcValue(VAR_U:VAR_P,total))
        bc%nVar = dns%nVar
        allocate(bc%bcKind(bc%nVar,NFACES), bc%bcW(bc%nVar,NFACES), bc%bcC(bc%nVar,total))
        bc%bcKind = BCK_NONE
        bc%bcW = 0.0d0
        bc%bcC = 0.0d0

        pos = 0
        do b = 1, int(blk%nBlocks)
            do dir = 1, 3
                do side = 0, 1
                    if (.not. block_face_is_physical(blk, b, dir, side)) cycle

                    face_id = boundary_face_id(dir, side)
                    call append_boundary_face_points(bc, face_id, b, pos, nx, ny, nz)
                end do
            end do
        end do

        call update_boundary_values(bc, blk, dns)

        ! The outlet faces: those carrying the OUTFLOW type on their normal
        ! component, and this rank's share of their points.
        do face_id = 1, NFACES
            bc%faceOutlet(face_id) = merge(1_C_INT, 0_C_INT, &
                bc%faceBcType(boundary_face_dir(face_id), face_id) == BC_OUTFLOW)
        end do
        bc%nOutlet = 0_C_INT
        do pos = 1, total
            bc%nOutlet = bc%nOutlet + bc%faceOutlet(int(bc%pointFace(pos)))
        end do
    end subroutine init_boundary_faces

    logical function block_face_is_physical(blk, b, dir, side)
        type(block_set_type), intent(in) :: blk
        integer, intent(in) :: b, dir, side

        ! FACE_PHYS only: FACE_CLOSED faces get zeroed halos, not boundary
        ! conditions.
        if (side == SIDE_MIN) then
            block_face_is_physical = blk%physLow(dir,b) == FACE_PHYS
        else
            block_face_is_physical = blk%physHigh(dir,b) == FACE_PHYS
        end if
    end function block_face_is_physical

    subroutine enter_boundary_data(bc)
        type(boundary_type), intent(inout) :: bc
        integer :: npts

#ifdef USE_OPENMP_OFFLOAD
        npts = int(bc%nTotal)
        !$omp target enter data map(to: bc)
        if (npts > 0) then
            !$omp target enter data map(to: bc%pointFace(1:npts), bc%slot(1:npts), &
            !$omp& bc%i(1:npts), bc%j(1:npts), bc%k(1:npts))
            !$omp target enter data map(to: bc%bcKind, bc%bcW, bc%bcC(1:bc%nVar,1:npts))
        end if
#endif
    end subroutine enter_boundary_data

    subroutine exit_boundary_data(bc)
        type(boundary_type), intent(inout) :: bc
        integer :: npts

#ifdef USE_OPENMP_OFFLOAD
        npts = int(bc%nTotal)
        if (npts > 0) then
            !$omp target exit data map(delete: bc%bcKind, bc%bcW, bc%bcC(1:bc%nVar,1:npts))
            !$omp target exit data map(delete: bc%pointFace(1:npts), bc%slot(1:npts), &
            !$omp& bc%i(1:npts), bc%j(1:npts), bc%k(1:npts))
        end if
        !$omp target exit data map(delete: bc)
#endif
    end subroutine exit_boundary_data

    subroutine destroy_boundary_faces(bc)
        type(boundary_type), intent(inout) :: bc

        if (allocated(bc%pointFace)) deallocate(bc%pointFace)
        if (allocated(bc%slot)) deallocate(bc%slot)
        if (allocated(bc%i)) deallocate(bc%i)
        if (allocated(bc%j)) deallocate(bc%j)
        if (allocated(bc%k)) deallocate(bc%k)
        if (allocated(bc%pointBcValue)) deallocate(bc%pointBcValue)
        if (allocated(bc%bcKind)) deallocate(bc%bcKind, bc%bcW, bc%bcC)
        bc%nTotal = 0_C_INT
        bc%nOutlet = 0_C_INT
        bc%faceOutlet = 0_C_INT
    end subroutine destroy_boundary_faces

    subroutine update_boundary_values(bc, blk, dns)
        type(boundary_type), intent(inout) :: bc
        type(block_set_type), intent(in) :: blk
        type(dns_type), intent(in) :: dns
        integer :: n, var, face_id, npts

        npts = int(bc%nTotal)
        if (npts <= 0) return

        if (any(bc%faceBcProfile == PROFILE_BLASIUS)) call prepare_blasius_profile(bc)

        do n = 1, npts
            face_id = int(bc%pointFace(n))
            do var = VAR_U, VAR_P
                if (bc%faceBcProfile(var,face_id) /= PROFILE_CONSTANT) then
                    bc%pointBcValue(var,n) = bc%faceBcDefaultValue(var,face_id) &
                        *profile_shape(bc, blk, dns, n, var, boundary_face_dir(face_id))
                else
                    bc%pointBcValue(var,n) = bc%faceBcDefaultValue(var,face_id)
                end if
            end do
        end do
        do var = VAR_U, VAR_P
            call resolve_affine_rows(bc, blk, var, bc%faceBcType(var,:), bc%pointBcValue(var,:))
        end do
    end subroutine update_boundary_values

    ! Boundary rows of a passive scalar stored in q at column `var`
    ! (scalar.f90, after its patch-derived types are final and before
    ! enter_boundary_data maps the tables): one type and one value per
    ! face, no profiles.
    subroutine set_scalar_bc_rows(bc, blk, var, bctype, value)
        type(boundary_type), intent(inout) :: bc
        type(block_set_type), intent(in) :: blk
        integer, intent(in) :: var
        integer(C_INT), intent(in) :: bctype(NFACES)
        real(C_DOUBLE), intent(in) :: value(NFACES)

        real(C_DOUBLE), allocatable :: pv(:)
        integer :: n, npts

        npts = int(bc%nTotal)
        if (npts <= 0) return
        allocate(pv(npts))
        do n = 1, npts
            pv(n) = value(int(bc%pointFace(n)))
        end do
        call resolve_affine_rows(bc, blk, var, bctype, pv)
    end subroutine set_scalar_bc_rows

    ! Resolve the BC type of column `var` on every face into the affine
    ! rows the kernels apply (see BCK_*): the kind and weight per face, the
    ! constant per point. `pointValue` is the row's boundary datum at each
    ! point (a Dirichlet value or a Neumann normal derivative). The Neumann
    ! constant is dn*v with dn the ghost-to-interior distance of the
    ! variable's OWN staggered coordinate along the face normal -- what the
    ! kernel used to compute per point; forming it here moves that product
    ! out of the kernel, which is why a nonzero Neumann datum is gated at
    ! nofma (the compiler used to contract q + dn*v into one FMA).
    subroutine resolve_affine_rows(bc, blk, var, bctype, pointValue)
        type(boundary_type), intent(inout) :: bc
        type(block_set_type), intent(in) :: blk
        integer, intent(in) :: var
        integer(C_INT), intent(in) :: bctype(NFACES)
        real(C_DOUBLE), intent(in) :: pointValue(:)

        integer :: f, n, npts, dir, side, b, ndir
        integer :: gi(3), ii(3)
        real(C_DOUBLE) :: dn

        npts = int(bc%nTotal)
        do f = 1, NFACES
            dir = boundary_face_dir(f)
            if (var == dir) then
                ! The normal velocity lives ON the face.
                select case (bctype(f))
                case (BC_DIRICHLET)
                    bc%bcKind(var,f) = BCK_FACE
                    bc%bcW(var,f) = 0.0d0
                case (BC_OUTFLOW)
                    ! An unknown, not a boundary value: no row.
                    bc%bcKind(var,f) = BCK_NONE
                    bc%bcW(var,f) = 0.0d0
                case default
                    ! Refused by resolve_face_bcs (numerics review F4).
                    error stop "boundary: Neumann normal velocity reached the row resolver"
                end select
            else
                select case (bctype(f))
                case (BC_DIRICHLET)
                    bc%bcKind(var,f) = BCK_GHOST
                    bc%bcW(var,f) = -1.0d0
                case (BC_NEUMANN, BC_OUTFLOW)
                    ! OUTFLOW never reaches a tangential/cell-centred row
                    ! (resolve_face_bcs gives those Neumann); kept as a
                    ! zero-gradient ghost should a caller pass it.
                    bc%bcKind(var,f) = BCK_GHOST
                    bc%bcW(var,f) = 1.0d0
                case default
                    error stop "boundary: unknown BC type in the row resolver"
                end select
            end if
        end do

        do n = 1, npts
            f = int(bc%pointFace(n))
            dir = boundary_face_dir(f)
            side = boundary_face_side(f)
            b = int(bc%slot(n))
            ndir = int(blk%nb(dir))
            gi = [int(bc%i(n)), int(bc%j(n)), int(bc%k(n))]
            ii = gi
            select case (bc%bcKind(var,f))
            case (BCK_GHOST)
                gi(dir) = merge(0, ndir + 1, side == SIDE_MIN)
                ii(dir) = merge(1, ndir, side == SIDE_MIN)
            case default
                gi(dir) = merge(1, ndir + 1, side == SIDE_MIN)
                ii(dir) = merge(2, ndir, side == SIDE_MIN)
            end select
            select case (dir)
            case (DIR_X)
                dn = blk%x(gi(1),var,b) - blk%x(ii(1),var,b)
            case (DIR_Y)
                dn = blk%y(gi(2),var,b) - blk%y(ii(2),var,b)
            case default
                dn = blk%z(gi(3),var,b) - blk%z(ii(3),var,b)
            end select
            if (bctype(f) == BC_DIRICHLET) then
                ! Ghost mirror 2v - q, or the pinned face value v itself.
                bc%bcC(var,n) = merge(2.0d0, 1.0d0, bc%bcKind(var,f) == BCK_GHOST)*pointValue(n)
            else
                ! Neumann data are stored as the normal derivative.
                bc%bcC(var,n) = dn*pointValue(n)
            end if
        end do
    end subroutine resolve_affine_rows

    ! Validate the blasius rows and build the similarity table once.
    ! Allowed: an x face (inlet) for u and v; a y face (top displacement)
    ! for v only. Host-only init path.
    subroutine prepare_blasius_profile(bc)
        type(boundary_type), intent(inout) :: bc
        integer :: var, face_id, fdir
        logical :: ok

        do face_id = 1, NFACES
            do var = VAR_U, VAR_P
                if (bc%faceBcProfile(var,face_id) /= PROFILE_BLASIUS) cycle
                fdir = boundary_face_dir(face_id)
                ok = (fdir == DIR_X .and. (var == VAR_U .or. var == VAR_V)) .or. &
                     (fdir == DIR_Y .and. var == VAR_V)
                if (.not. ok) then
                    error stop "[boundary] blasius profile: x faces (u,v) or y faces (v) only"
                end if
            end do
        end do
        if (.not. allocated(bc%blasFp)) call build_blasius_table(bc)
    end subroutine prepare_blasius_profile

    ! Solve the Blasius equation f''' = -f f''/2 (eta = y sqrt(U/(nu x)))
    ! by RK4 + secant shooting on f'(inf) = 1, tabulating f and f' on a
    ! uniform eta grid. beta = 2 f''(0) is exact by the momentum integral.
    subroutine build_blasius_table(bc)
        type(boundary_type), intent(inout) :: bc

        real(C_DOUBLE), parameter :: ETA_MAX = 12.0d0
        integer, parameter :: NTAB = 12000
        real(C_DOUBLE) :: fpp0, fpp0_prev, err, err_prev, fpp0_next
        integer :: it, ntop

        bc%blasH = ETA_MAX/real(NTAB, C_DOUBLE)
        allocate(bc%blasF(0:NTAB), bc%blasFp(0:NTAB))

        fpp0_prev = 0.33d0
        fpp0 = 0.34d0
        err_prev = blasius_integrate(bc, fpp0_prev, NTAB)
        do it = 1, 50
            err = blasius_integrate(bc, fpp0, NTAB)
            if (abs(err) < 1.0d-13 .or. err == err_prev) exit
            fpp0_next = fpp0 - err*(fpp0 - fpp0_prev)/(err - err_prev)
            fpp0_prev = fpp0
            err_prev = err
            fpp0 = fpp0_next
        end do

        ntop = NTAB
        bc%blasBeta = 2.0d0*fpp0
        bc%blasDisp = ETA_MAX - bc%blasF(ntop)
    end subroutine build_blasius_table

    ! One RK4 integration of the Blasius system from eta = 0 with f''(0) =
    ! fpp0, filling the f/f' tables; returns f'(eta_max) - 1 (the shooting
    ! residual).
    real(C_DOUBLE) function blasius_integrate(bc, fpp0, ntab) result(err)
        type(boundary_type), intent(inout) :: bc
        real(C_DOUBLE), intent(in) :: fpp0
        integer, intent(in) :: ntab

        real(C_DOUBLE) :: y(3), k1(3), k2(3), k3(3), k4(3), h
        integer :: i

        h = bc%blasH
        y = [0.0d0, 0.0d0, fpp0]
        bc%blasF(0) = 0.0d0
        bc%blasFp(0) = 0.0d0
        do i = 1, ntab
            k1 = blasius_rhs(y)
            k2 = blasius_rhs(y + 0.5d0*h*k1)
            k3 = blasius_rhs(y + 0.5d0*h*k2)
            k4 = blasius_rhs(y + h*k3)
            y = y + h*(k1 + 2.0d0*k2 + 2.0d0*k3 + k4)/6.0d0
            bc%blasF(i) = y(1)
            bc%blasFp(i) = y(2)
        end do
        err = y(2) - 1.0d0
    end function blasius_integrate

    pure function blasius_rhs(y) result(dy)
        real(C_DOUBLE), intent(in) :: y(3)
        real(C_DOUBLE) :: dy(3)

        dy = [y(2), y(3), -0.5d0*y(1)*y(3)]
    end function blasius_rhs

    ! Linearly interpolated f and f' at eta; beyond the table f' = 1 and
    ! f = eta - blasDisp (the outer asymptote).
    subroutine blasius_lookup(bc, eta, f, fp)
        type(boundary_type), intent(in) :: bc
        real(C_DOUBLE), intent(in) :: eta
        real(C_DOUBLE), intent(out) :: f, fp

        integer :: i, ntop
        real(C_DOUBLE) :: s

        ntop = ubound(bc%blasFp, 1)
        s = eta/bc%blasH
        if (s >= real(ntop, C_DOUBLE)) then
            f = eta - bc%blasDisp
            fp = 1.0d0
            return
        end if
        i = int(s)
        s = s - real(i, C_DOUBLE)
        f = (1.0d0 - s)*bc%blasF(i) + s*bc%blasF(i+1)
        fp = (1.0d0 - s)*bc%blasFp(i) + s*bc%blasFp(i+1)
    end subroutine blasius_lookup

    ! Profile factor at boundary point n for variable var.
    ! PARABOLA: product over the face's non-periodic tangential directions
    ! of 4*s*(1-s), with s the variable's own staggered coordinate
    ! normalized by the domain length (domain assumed to start at 0);
    ! periodic tangential directions contribute no shape (factor 1).
    ! BLASIUS: u = f'(eta), v = beta*(eta f' - f)/(2 Re_theta) at
    ! eta = beta*y/theta (see the PROFILE_* comment; the row value is U_inf,
    ! Re_theta = U_inf*theta*re with U_inf read from the face's u row).
    real(C_DOUBLE) function profile_shape(bc, blk, dns, n, var, face_dir) result(fac)
        type(boundary_type), intent(in) :: bc
        type(block_set_type), intent(in) :: blk
        type(dns_type), intent(in) :: dns
        integer, intent(in) :: n, var, face_dir

        integer :: d, b, face_id
        real(C_DOUBLE) :: coord, s, eta, f, fp, re_theta

        b = int(bc%slot(n))
        face_id = int(bc%pointFace(n))
        fac = 1.0d0

        select case (bc%faceBcProfile(var,face_id))
        case (PROFILE_PARABOLA)
            do d = 1, 3
                if (d == face_dir) cycle
                if (bc%isPeriodic(d)) cycle
                select case (d)
                case (1)
                    coord = blk%x(int(bc%i(n)), var, b)
                case (2)
                    coord = blk%y(int(bc%j(n)), var, b)
                case (3)
                    coord = blk%z(int(bc%k(n)), var, b)
                end select
                s = coord/dns%leng(d)
                fac = fac*4.0d0*s*(1.0d0 - s)
            end do
        case (PROFILE_BLASIUS)
            if (face_dir == DIR_X) then
                ! Inlet (x face): u = f'(eta), v = the entrainment component,
                ! eta = beta*y/theta_in (the station is x = 0). Re_theta from
                ! the face's u-row value = U_inf.
                eta = bc%blasBeta*blk%y(int(bc%j(n)), var, b)/bc%blasiusTheta
                call blasius_lookup(bc, eta, f, fp)
                if (var == VAR_U) then
                    fac = fp
                else
                    re_theta = bc%faceBcDefaultValue(VAR_U,face_id)*bc%blasiusTheta*dns%re
                    fac = bc%blasBeta*(eta*fp - f)/(2.0d0*re_theta)
                end if
            else
                ! Top (y face) displacement BC: the x-varying Blasius
                ! entrainment v at the domain top y = leng(2). Similarity
                ! station x + x_v with x_v = Re_theta,in*theta/beta^2;
                ! eta = y_top*sqrt(Re_theta,in/(theta*(x+x_v))); Re_theta,in
                ! from the v-row value = U_inf. v/U_inf =
                ! 0.5*sqrt(theta/(Re_theta,in*(x+x_v)))*(eta f' - f).
                re_theta = bc%faceBcDefaultValue(VAR_V,face_id)*bc%blasiusTheta*dns%re
                coord = blk%x(int(bc%i(n)), var, b) &
                      + re_theta*bc%blasiusTheta/(bc%blasBeta*bc%blasBeta)   ! x + x_v
                eta = dns%leng(2)*sqrt(re_theta/(bc%blasiusTheta*coord))
                call blasius_lookup(bc, eta, f, fp)
                fac = 0.5d0*sqrt(bc%blasiusTheta/(re_theta*coord))*(eta*fp - f)
            end if
        end select
    end function profile_shape

    subroutine append_boundary_face_points(bc, face_id, slot, pos, nx, ny, nz)
        type(boundary_type), intent(inout) :: bc
        integer, intent(in) :: face_id, slot, nx, ny, nz
        integer, intent(inout) :: pos
        integer :: i, j, k, dir, side

        dir = boundary_face_dir(face_id)
        side = boundary_face_side(face_id)

        select case (dir)
        case (DIR_X)
            i = merge(1, nx, side == SIDE_MIN)
            do k = 1, nz
                do j = 1, ny
                    pos = pos + 1
                    bc%pointFace(pos) = int(face_id, C_INT)
                    bc%slot(pos) = int(slot, C_INT)
                    bc%i(pos) = int(i, C_INT)
                    bc%j(pos) = int(j, C_INT)
                    bc%k(pos) = int(k, C_INT)
                end do
            end do
        case (DIR_Y)
            j = merge(1, ny, side == SIDE_MIN)
            do k = 1, nz
                do i = 1, nx
                    pos = pos + 1
                    bc%pointFace(pos) = int(face_id, C_INT)
                    bc%slot(pos) = int(slot, C_INT)
                    bc%i(pos) = int(i, C_INT)
                    bc%j(pos) = int(j, C_INT)
                    bc%k(pos) = int(k, C_INT)
                end do
            end do
        case (DIR_Z)
            k = merge(1, nz, side == SIDE_MIN)
            do j = 1, ny
                do i = 1, nx
                    pos = pos + 1
                    bc%pointFace(pos) = int(face_id, C_INT)
                    bc%slot(pos) = int(slot, C_INT)
                    bc%i(pos) = int(i, C_INT)
                    bc%j(pos) = int(j, C_INT)
                    bc%k(pos) = int(k, C_INT)
                end do
            end do
        case default
            error stop "invalid boundary direction"
        end select
    end subroutine append_boundary_face_points

    integer function boundary_face_n(dir, nx, ny, nz) result(n)
        integer, intent(in) :: dir, nx, ny, nz

        select case (dir)
        case (DIR_X)
            n = ny*nz
        case (DIR_Y)
            n = nx*nz
        case (DIR_Z)
            n = nx*ny
        case default
            error stop "invalid boundary direction"
        end select
    end function boundary_face_n

    integer function boundary_face_id(dir, side) result(face_id)
        integer, intent(in) :: dir, side

        face_id = 2*(dir - 1) + side + 1
    end function boundary_face_id

    integer function boundary_face_dir(face_id) result(dir)
        integer, intent(in) :: face_id

        dir = (face_id + 1)/2
    end function boundary_face_dir

    integer function boundary_face_side(face_id) result(side)
        integer, intent(in) :: face_id

        side = modulo(face_id - 1, 2)
    end function boundary_face_side

    ! Physical-face boundary writes of the q variables `vars` (default
    ! u, v, w, p; scalar_sync passes the scalar columns): every point's
    ! resolved row, dst = w*src + C (see BCK_*). Runs at init/restart and
    ! ONCE per substage, in the projection after its last correction (there
    ! is no call after the predictor: it commits only what it predicted and
    ! nothing reads a ghost before the projection ends), always
    ! BEFORE the halo exchange -- the exchange's tangential extension copies
    ! a neighbour block's ghosts into edge/corner halos, so the ghosts must
    ! be current when it reads them (every transfer extends that way, 2:1
    ! entries included: comm.f90 candidate_boxes).
    !
    ! The normal velocity of an OUTLET face has no row here (BCK_NONE): it
    ! is an unknown of the predictor and of the projection.
    subroutine apply_bc(blk, bc, vars)
        type(block_set_type), intent(inout) :: blk
        type(boundary_type), intent(in) :: bc
        integer(C_INT), intent(in), optional :: vars(:)

        integer :: n, npts, b, face_id, dir, side, v, nv, var, kind
        integer :: ghost_idx, interior_idx, face_idx, neighbor_idx
        integer :: gi(3), ii(3)
        integer(C_INT) :: local_n(1:3)
        integer(C_INT), allocatable :: vl(:)

        npts = int(bc%nTotal)
        if (npts <= 0) return

        if (present(vars)) then
            vl = vars
        else
            vl = [VAR_U, VAR_V, VAR_W, VAR_P]
        end if
        nv = size(vl)
        local_n = blk%nb(1:3)

        !$omp target teams distribute parallel do &
        !$omp& map(to: npts, nv, local_n(1:3), vl, &
        !$omp& bc%pointFace(1:npts), bc%slot(1:npts), bc%i(1:npts), bc%j(1:npts), bc%k(1:npts), &
        !$omp& bc%bcKind, bc%bcW, bc%bcC(1:bc%nVar,1:npts)) &
        !$omp& map(tofrom: blk%q) &
        !$omp& private(n,b,face_id,dir,side,v,var,kind, &
        !$omp& ghost_idx,interior_idx,face_idx,neighbor_idx,gi,ii)
        do n = 1, npts
            ! Each entry is one active physical boundary point of one block.
            face_id = int(bc%pointFace(n))
            b = int(bc%slot(n))
            dir = (face_id + 1)/2
            side = modulo(face_id - 1, 2)
            gi = [int(bc%i(n)), int(bc%j(n)), int(bc%k(n))]
            ii = gi

            ! The ghost/interior pair (cell-centred rows and tangential
            ! velocities) and the face/neighbour pair (the normal velocity,
            ! which lives on the boundary face in the staggered layout).
            if (side == SIDE_MIN) then
                ghost_idx = 0
                interior_idx = 1
                face_idx = 1
                neighbor_idx = 2
            else
                ghost_idx = int(local_n(dir)) + 1
                interior_idx = int(local_n(dir))
                face_idx = int(local_n(dir)) + 1
                neighbor_idx = int(local_n(dir))
            end if

            do v = 1, nv
                var = int(vl(v))
                kind = int(bc%bcKind(var,face_id))
                if (kind == BCK_NONE) cycle
                if (kind == BCK_GHOST) then
                    gi(dir) = ghost_idx
                    ii(dir) = interior_idx
                else
                    gi(dir) = face_idx
                    ii(dir) = neighbor_idx
                end if
                blk%q(gi(1),gi(2),gi(3),var,b) = &
                    bc%bcW(var,face_id)*blk%q(ii(1),ii(2),ii(3),var,b) + bc%bcC(var,n)
            end do
        end do
        !$omp end target teams distribute parallel do
    end subroutine apply_bc

    ! The INITIAL value of an outlet-normal face the run was not given: the
    ! zero-gradient value, face := interior neighbour. An outlet face is
    ! state (the predictor advances it from its own previous value), so this
    ! runs ONCE, before the first step -- on a cold start for every outlet
    ! face, on a restart only for the faces the file did not supply
    ! (`given(face)` true = leave the face alone: a low face is in the
    ! velocity dataset, a high one in io.f90 read_outlet_planes). The first
    ! projection then makes the boundary cell solenoidal.
    subroutine init_outlet_faces(blk, bc, given)
        type(block_set_type), intent(inout) :: blk
        type(boundary_type), intent(in) :: bc
        logical, intent(in) :: given(NFACES)

        integer :: n, npts, b, face_id, dir, side
        integer :: fi(3), ni(3)
        integer(C_INT) :: local_n(1:3), todo(NFACES)

        if (bc%nOutlet <= 0_C_INT) return
        npts = int(bc%nTotal)
        local_n = blk%nb(1:3)
        todo = merge(0_C_INT, bc%faceOutlet, given)

        ! The loop runs over ALL boundary points and skips the others, as
        ! apply_bc does: npts is then a loop bound as well as a map bound.
        !$omp target teams distribute parallel do &
        !$omp& map(to: npts, local_n(1:3), todo(1:NFACES), &
        !$omp& bc%pointFace(1:npts), bc%slot(1:npts), bc%i(1:npts), bc%j(1:npts), bc%k(1:npts)) &
        !$omp& map(tofrom: blk%q) &
        !$omp& private(n,b,face_id,dir,side,fi,ni)
        do n = 1, npts
            face_id = int(bc%pointFace(n))
            if (todo(face_id) == 0_C_INT) cycle
            b = int(bc%slot(n))
            dir = (face_id + 1)/2
            side = modulo(face_id - 1, 2)
            fi = [int(bc%i(n)), int(bc%j(n)), int(bc%k(n))]
            ni = fi
            if (side == SIDE_MIN) then
                fi(dir) = 1
                ni(dir) = 2
            else
                fi(dir) = int(local_n(dir)) + 1
                ni(dir) = int(local_n(dir))
            end if
            blk%q(fi(1),fi(2),fi(3),dir,b) = blk%q(ni(1),ni(2),ni(3),dir,b)
        end do
        !$omp end target teams distribute parallel do
    end subroutine init_outlet_faces

    ! Dirichlet VELOCITY data inside an immersed body
    ! (docs/next_session_body_at_outlet.md, B2). A Dirichlet face is pinned
    ! to its datum by apply_bc and no stage penalizes it (it is in no
    ! predicted range, and the projection does not correct it), so an inlet
    ! profile that is non-zero where the body crosses the inlet plane injects
    ! mass into solid cells, which pass it on through their penalized faces
    ! at O(u_in) under a pressure that grows every substage: a conduit
    ! through the body. The datum is therefore ZERO wherever the row's own
    ! staggered location is solid-centred -- the normal row at the face
    ! itself, a tangential row at its ghost position, where the mirror
    ! 2v - q with v = 0 makes the plane no-slip. Graded (fluid-centred)
    ! locations keep the datum: the same staircase every cut cell carries.
    ! Every Dirichlet velocity row, on every face, low and high alike (a wall
    ! datum is 0 already; a moving wall's would be wrong inside the body).
    ! Host code on the host coefficients, ONCE at init, before the rows are
    ! mapped: with no body it touches nothing. `nMasked` counts the NON-ZERO
    ! data this zeroed on this rank (a wall row inside the body is 0 already
    ! and would swamp the count), `nRows` the Dirichlet velocity rows examined.
    subroutine mask_dirichlet_velocity_in_solid(bc, blk, coef, nMasked, nRows)
        type(boundary_type), intent(inout) :: bc
        type(block_set_type), intent(in) :: blk
        real(C_DOUBLE), intent(in) :: coef(0:,0:,0:,1:,1:)
        integer, intent(out) :: nMasked, nRows

        integer :: n, b, f, dir, side, var, ndir
        integer :: gi(3)

        nMasked = 0
        nRows = 0
        do n = 1, int(bc%nTotal)
            f = int(bc%pointFace(n))
            dir = boundary_face_dir(f)
            side = boundary_face_side(f)
            b = int(bc%slot(n))
            ndir = int(blk%nb(dir))
            do var = VAR_U, VAR_W
                if (bc%bcKind(var,f) == BCK_NONE) cycle
                if (bc%bcKind(var,f) == BCK_GHOST .and. bc%bcW(var,f) /= -1.0d0) cycle   ! Neumann row
                gi = [int(bc%i(n)), int(bc%j(n)), int(bc%k(n))]
                if (bc%bcKind(var,f) == BCK_GHOST) then
                    gi(dir) = merge(0, ndir + 1, side == SIDE_MIN)
                else
                    gi(dir) = merge(1, ndir + 1, side == SIDE_MIN)
                end if
                nRows = nRows + 1
                if (abs(coef(gi(1),gi(2),gi(3),var,b)) > SOLID_FACE_THRESHOLD) then
                    if (bc%bcC(var,n) /= 0.0d0) nMasked = nMasked + 1
                    bc%bcC(var,n) = 0.0d0
                end if
            end do
        end do
    end subroutine mask_dirichlet_velocity_in_solid

    ! Physical-face ghosts of a standalone cell-centred scalar array (the
    ! RANS transport scalars, the projection's phi): the same affine write
    ! as apply_bc, with the row resolved per CALL from a per-face mode --
    ! COPY (zero normal gradient) w = 1, C = 0; MIRROR (face value 0)
    ! w = -1, C = 0; VALUE (prescribed face value) w = -1, C = 2 value;
    ! NONE leaves the ghost untouched. A separate routine only because a
    ! strided q slice cannot portably be passed under OpenMP target mapping.
    subroutine apply_scalar_bc(blk, bc, s, mode, value)
        type(block_set_type), intent(in) :: blk
        type(boundary_type), intent(in) :: bc
        real(C_DOUBLE), intent(inout) :: s(0:,0:,0:,1:)
        integer(C_INT), intent(in) :: mode(NFACES)
        real(C_DOUBLE), intent(in), optional :: value(NFACES)

        integer :: n, npts, b, face_id, dir, side, f
        integer :: gi(3), ii(3)
        integer(C_INT) :: local_n(1:3), on_l(NFACES)
        real(C_DOUBLE) :: w_l(NFACES), c_l(NFACES)

        npts = int(bc%nTotal)
        if (npts <= 0) return
        if (all(mode == SCALAR_BC_NONE)) return
        local_n = blk%nb(1:3)

        do f = 1, NFACES
            on_l(f) = merge(1_C_INT, 0_C_INT, mode(f) /= SCALAR_BC_NONE)
            w_l(f) = merge(1.0d0, -1.0d0, mode(f) == SCALAR_BC_COPY)
            c_l(f) = 0.0d0
            if (mode(f) == SCALAR_BC_VALUE .and. present(value)) c_l(f) = 2.0d0*value(f)
        end do

        !$omp target teams distribute parallel do &
        !$omp& map(to: npts, local_n(1:3), on_l(1:NFACES), w_l(1:NFACES), c_l(1:NFACES), &
        !$omp& bc%pointFace(1:npts), bc%slot(1:npts), bc%i(1:npts), bc%j(1:npts), bc%k(1:npts)) &
        !$omp& map(tofrom: s) &
        !$omp& private(n,b,face_id,dir,side,gi,ii)
        do n = 1, npts
            face_id = int(bc%pointFace(n))
            if (on_l(face_id) == 0_C_INT) cycle
            b = int(bc%slot(n))
            dir = (face_id + 1)/2
            side = modulo(face_id - 1, 2)
            gi = [int(bc%i(n)), int(bc%j(n)), int(bc%k(n))]
            ii = gi
            if (side == SIDE_MIN) then
                gi(dir) = 0
                ii(dir) = 1
            else
                gi(dir) = int(local_n(dir)) + 1
                ii(dir) = int(local_n(dir))
            end if
            s(gi(1),gi(2),gi(3),b) = w_l(face_id)*s(ii(1),ii(2),ii(3),b) + c_l(face_id)
        end do
        !$omp end target teams distribute parallel do
    end subroutine apply_scalar_bc

end module boundary
