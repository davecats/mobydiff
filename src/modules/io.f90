module io
    use, intrinsic :: iso_c_binding
    use :: init, only: dns_type, grid_type, VAR_U, VAR_V, VAR_W, VAR_P, NVAR, NVEL, &
        config_seen_type
    use :: blocks, only: block_set_type
    use :: boundary, only: boundary_type, NFACES
    use :: comm, only: comm_type
    implicit none

    ! Stride of one slot in the packed field-variable name table handed to
    ! field_hdf5.c (FDM_VAR_NAME_LEN there -- keep the two in lockstep).
    integer, parameter :: VAR_NAME_LEN = 32

    interface
        function fdm_h5_write_field(file_name, nbx, nby, nbz, n_blocks, &
                n_blocks_global, id_start, block_origin, block_level, &
                rank, nranks, global_nx, global_ny, global_nz, &
                step, nsteps, lx, ly, lz, &
                re, dt, t_final, t_current, cfl, cflmax, pecletmax, dtmax, &
                forcing, pressure_niter, pressure_sor, &
                ibm_enabled, bc_count, periodic, bc_type, bc_value, grid_distribution, grid_stretch, &
                grid_natural_dyw_plus, x_node, y_node, z_node, n_var, var_names, q) &
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
            ! n_var = NVAR + nScalar variables of q, named by a packed table
            ! of n_var NUL-terminated slots of VAR_NAME_LEN characters.
            integer(C_INT), value :: n_var
            character(kind=C_CHAR), intent(in) :: var_names(*)
            real(C_DOUBLE), intent(in) :: q(*)
            integer(C_INT) :: ierr
        end function fdm_h5_write_field

        ! Appends the LES eddy viscosity as a "nut" block dataset to an existing
        ! field file (collective). Called only when LES is active, so nut is a
        ! plain allocated array -- no c_loc (nvfortran ICEs on it).
        function fdm_h5_append_nut(file_name, nbx, nby, nbz, n_blocks, &
                n_blocks_global, id_start, nut) &
                bind(C, name="fdm_h5_append_nut") result(ierr)
            import :: C_CHAR, C_INT, C_DOUBLE
            character(kind=C_CHAR), intent(in) :: file_name(*)
            integer(C_INT), value :: nbx, nby, nbz, n_blocks, n_blocks_global, id_start
            real(C_DOUBLE), intent(in) :: nut(*)
            integer(C_INT) :: ierr
        end function fdm_h5_append_nut

        ! Named cell-centred block-layout scalar append/read (the RANS k and
        ! omega snapshot/restart datasets ride the nut machinery).
        function fdm_h5_append_scalar(file_name, name, nbx, nby, nbz, n_blocks, &
                n_blocks_global, id_start, s) &
                bind(C, name="fdm_h5_append_scalar") result(ierr)
            import :: C_CHAR, C_INT, C_DOUBLE
            character(kind=C_CHAR), intent(in) :: file_name(*), name(*)
            integer(C_INT), value :: nbx, nby, nbz, n_blocks, n_blocks_global, id_start
            real(C_DOUBLE), intent(in) :: s(*)
            integer(C_INT) :: ierr
        end function fdm_h5_append_scalar

        function fdm_h5_read_scalar(file_name, name, nbx, nby, nbz, n_blocks, &
                n_blocks_global, id_start, found, s) &
                bind(C, name="fdm_h5_read_scalar") result(ierr)
            import :: C_CHAR, C_INT, C_DOUBLE
            character(kind=C_CHAR), intent(in) :: file_name(*), name(*)
            integer(C_INT), value :: nbx, nby, nbz, n_blocks, n_blocks_global, id_start
            integer(C_INT), intent(out) :: found
            real(C_DOUBLE), intent(inout) :: s(*)
            integer(C_INT) :: ierr
        end function fdm_h5_read_scalar

        ! [blocks] refine_dims file-variant marker (xz quadtree mode): the
        ! per-direction refinement mask attribute. Written only when the
        ! mask is not all-ones so octree-mode files stay byte-identical;
        ! absent on read means the octree default (1,1,1). Collective.
        function fdm_h5_append_refine_dims(file_name, mask) &
                bind(C, name="fdm_h5_append_refine_dims") result(ierr)
            import :: C_CHAR, C_INT
            character(kind=C_CHAR), intent(in) :: file_name(*)
            integer(C_INT), intent(in) :: mask(*)
            integer(C_INT) :: ierr
        end function fdm_h5_append_refine_dims

        function fdm_h5_read_refine_dims(file_name, mask, has_blocks) &
                bind(C, name="fdm_h5_read_refine_dims") result(ierr)
            import :: C_CHAR, C_INT
            character(kind=C_CHAR), intent(in) :: file_name(*)
            integer(C_INT), intent(out) :: mask(*), has_blocks
            integer(C_INT) :: ierr
        end function fdm_h5_read_refine_dims

        ! Grid datasets in the case file (P3, mobygrid absorbed): node
        ! lines + the mobygrid-format attributes, so the case file doubles
        ! as the --grid-file of the retired mobygeom reference tooling.
        function fdm_h5_case_append_grid(file_name, nx, ny, nz, &
                periodic, grid_distribution, grid_stretch, grid_natural_dyw_plus, &
                grid_natural_one_sided, x_node, y_node, z_node, write_data) &
                bind(C, name="fdm_h5_case_append_grid") result(ierr)
            import :: C_CHAR, C_INT, C_DOUBLE
            character(kind=C_CHAR), intent(in) :: file_name(*)
            integer(C_INT), value :: nx, ny, nz, write_data
            integer(C_INT), intent(in) :: periodic(*), grid_distribution(*), grid_natural_one_sided(*)
            real(C_DOUBLE), intent(in) :: grid_stretch(*), grid_natural_dyw_plus(*)
            real(C_DOUBLE), intent(in) :: x_node(*), y_node(*), z_node(*)
            integer(C_INT) :: ierr
        end function fdm_h5_case_append_grid

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

        function fdm_h5_read_block_masks(file_name, level, n_raster, block_nb, found, touch, buried) &
                bind(C, name="fdm_h5_read_block_masks") result(ierr)
            import :: C_CHAR, C_INT
            character(kind=C_CHAR), intent(in) :: file_name(*)
            integer(C_INT), value :: level, n_raster
            integer(C_INT), intent(in) :: block_nb(*)
            integer(C_INT), intent(out) :: found
            integer(C_INT), intent(out) :: touch(*), buried(*)
            integer(C_INT) :: ierr
        end function fdm_h5_read_block_masks

        ! Per-level geometry-mask window attrs (deep-refinement files store
        ! the rasters windowed to the padded STL bbox); has_win = 0 on
        ! legacy full-raster files.
        function fdm_h5_read_mask_window(file_name, level, lo, dims, has_win) &
                bind(C, name="fdm_h5_read_mask_window") result(ierr)
            import :: C_CHAR, C_INT
            character(kind=C_CHAR), intent(in) :: file_name(*)
            integer(C_INT), value :: level
            integer(C_INT), intent(out) :: lo(*), dims(*), has_win
            integer(C_INT) :: ierr
        end function fdm_h5_read_mask_window

        function fdm_h5_read_block_active(file_name, n_lattice, block_nb, found, active) &
                bind(C, name="fdm_h5_read_block_active") result(ierr)
            import :: C_CHAR, C_INT
            character(kind=C_CHAR), intent(in) :: file_name(*)
            integer(C_INT), value :: n_lattice
            integer(C_INT), intent(in) :: block_nb(*)
            integer(C_INT), intent(out) :: found
            integer(C_INT), intent(out) :: active(*)
            integer(C_INT) :: ierr
        end function fdm_h5_read_block_active

        function fdm_h5_read_dwall_blocks(file_name, nbx, nby, nbz, n_blocks, id_start, &
                found, dwall) &
                bind(C, name="fdm_h5_read_dwall_blocks") result(ierr)
            import :: C_CHAR, C_INT, C_DOUBLE
            character(kind=C_CHAR), intent(in) :: file_name(*)
            integer(C_INT), value :: nbx, nby, nbz, n_blocks, id_start
            integer(C_INT), intent(out) :: found
            real(C_DOUBLE), intent(inout) :: dwall(*)
            integer(C_INT) :: ierr
        end function fdm_h5_read_dwall_blocks

        ! Case-file LAYOUT readers (step 7-2): the header attrs + the leaf
        ! table + the level-0 node lines, from which the solver builds its
        ! block set without the leaf builder. Serial opens: the table is
        ! global state every rank holds in full.
        function fdm_h5_case_read_layout(file_name, nx, ny, nz, lx, ly, lz, re, &
                block_nb, refine_mask, nb_auto, nb_ranks, n_blocks_global) &
                bind(C, name="fdm_h5_case_read_layout") result(ierr)
            import :: C_CHAR, C_INT, C_DOUBLE
            character(kind=C_CHAR), intent(in) :: file_name(*)
            integer(C_INT), intent(out) :: nx, ny, nz
            real(C_DOUBLE), intent(out) :: lx, ly, lz, re
            integer(C_INT), intent(out) :: block_nb(*), refine_mask(*)
            integer(C_INT), intent(out) :: nb_auto, nb_ranks, n_blocks_global
            integer(C_INT) :: ierr
        end function fdm_h5_case_read_layout

        function fdm_h5_case_read_inputs(file_name, buf, cap, found) &
                bind(C, name="fdm_h5_case_read_inputs") result(ierr)
            import :: C_CHAR, C_INT
            character(kind=C_CHAR), intent(in) :: file_name(*)
            character(kind=C_CHAR), intent(out) :: buf(*)
            integer(C_INT), value :: cap
            integer(C_INT), intent(out) :: found
            integer(C_INT) :: ierr
        end function fdm_h5_case_read_inputs

        function fdm_h5_case_read_blocks(file_name, n_blocks_global, rows) &
                bind(C, name="fdm_h5_case_read_blocks") result(ierr)
            import :: C_CHAR, C_INT
            character(kind=C_CHAR), intent(in) :: file_name(*)
            integer(C_INT), value :: n_blocks_global
            integer(C_INT), intent(out) :: rows(*)
            integer(C_INT) :: ierr
        end function fdm_h5_case_read_blocks

        function fdm_h5_case_read_grid(file_name, nx, ny, nz, x_node, y_node, z_node, found) &
                bind(C, name="fdm_h5_case_read_grid") result(ierr)
            import :: C_CHAR, C_INT, C_DOUBLE
            character(kind=C_CHAR), intent(in) :: file_name(*)
            integer(C_INT), value :: nx, ny, nz
            real(C_DOUBLE), intent(out) :: x_node(*), y_node(*), z_node(*)
            integer(C_INT), intent(out) :: found
            integer(C_INT) :: ierr
        end function fdm_h5_case_read_grid

        ! Restart guard: this rank's rows of a snapshot's blocks table
        ! against the solver's leaf table (0 match, 2 mismatch, 1 error; 0
        ! for a legacy file without a table).
        function fdm_h5_check_block_table(file_name, n_blocks, id_start, &
                block_origin, block_level) &
                bind(C, name="fdm_h5_check_block_table") result(ierr)
            import :: C_CHAR, C_INT
            character(kind=C_CHAR), intent(in) :: file_name(*)
            integer(C_INT), value :: n_blocks, id_start
            integer(C_INT), intent(in) :: block_origin(*), block_level(*)
            integer(C_INT) :: ierr
        end function fdm_h5_check_block_table

        function fdm_h5_write_rans_geometry(file_name, nbx, nby, nbz, n_blocks, &
                n_blocks_global, id_start, block_origin, block_level, &
                dwall, yeff, wallcell, xc, yc, zc) &
                bind(C, name="fdm_h5_write_rans_geometry") result(ierr)
            import :: C_CHAR, C_INT, C_DOUBLE
            character(kind=C_CHAR), intent(in) :: file_name(*)
            integer(C_INT), value :: nbx, nby, nbz, n_blocks, n_blocks_global, id_start
            integer(C_INT), intent(in) :: block_origin(*), block_level(*)
            real(C_DOUBLE), intent(in) :: dwall(*), yeff(*)
            integer(C_INT), intent(in) :: wallcell(*)
            real(C_DOUBLE), intent(in) :: xc(*), yc(*), zc(*)
            integer(C_INT) :: ierr
        end function fdm_h5_write_rans_geometry

        ! Case-file writers (moby_prepare): the exact inverses of the
        ! fdm_h5_read_* coefficient-file readers above. Collective over
        ! MPI_COMM_WORLD; write_data selects the one rank writing the
        ! lattice-global rasters.
        function fdm_h5_case_create(file_name, nx, ny, nz, lx, ly, lz, re, &
                block_nb, block_levels, refine_mask, nb_auto, nb_ranks, inputs, &
                n_blocks_global, id_start, n_blocks, block_origin, block_level) &
                bind(C, name="fdm_h5_case_create") result(ierr)
            import :: C_CHAR, C_INT, C_DOUBLE
            character(kind=C_CHAR), intent(in) :: file_name(*), inputs(*)
            integer(C_INT), value :: nx, ny, nz
            real(C_DOUBLE), value :: lx, ly, lz, re
            integer(C_INT), intent(in) :: block_nb(*), refine_mask(*)
            integer(C_INT), value :: block_levels, nb_auto, nb_ranks
            integer(C_INT), value :: n_blocks_global, id_start, n_blocks
            integer(C_INT), intent(in) :: block_origin(*), block_level(*)
            integer(C_INT) :: ierr
        end function fdm_h5_case_create

        function fdm_h5_case_append_masks(file_name, level, n_raster, touch, buried, &
                write_data) bind(C, name="fdm_h5_case_append_masks") result(ierr)
            import :: C_CHAR, C_INT
            character(kind=C_CHAR), intent(in) :: file_name(*)
            integer(C_INT), value :: level, n_raster, write_data
            integer(C_INT), intent(in) :: touch(*), buried(*)
            integer(C_INT) :: ierr
        end function fdm_h5_case_append_masks

        function fdm_h5_case_append_mask_window(file_name, level, lo, dims) &
                bind(C, name="fdm_h5_case_append_mask_window") result(ierr)
            import :: C_CHAR, C_INT
            character(kind=C_CHAR), intent(in) :: file_name(*)
            integer(C_INT), value :: level
            integer(C_INT), intent(in) :: lo(*), dims(*)
            integer(C_INT) :: ierr
        end function fdm_h5_case_append_mask_window

        function fdm_h5_case_append_active(file_name, n_lattice, active, write_data) &
                bind(C, name="fdm_h5_case_append_active") result(ierr)
            import :: C_CHAR, C_INT
            character(kind=C_CHAR), intent(in) :: file_name(*)
            integer(C_INT), value :: n_lattice, write_data
            integer(C_INT), intent(in) :: active(*)
            integer(C_INT) :: ierr
        end function fdm_h5_case_append_active

        ! n_comp = the coefficient array's component extent (3, or 4 when
        ! the case declares passive scalars): coef_blocks always holds the
        ! three staggered components, coef_p_blocks the cell-centred one.
        function fdm_h5_case_append_coef(file_name, nbx, nby, nbz, n_blocks, &
                n_blocks_global, id_start, n_comp, coef) &
                bind(C, name="fdm_h5_case_append_coef") result(ierr)
            import :: C_CHAR, C_INT, C_DOUBLE
            character(kind=C_CHAR), intent(in) :: file_name(*)
            integer(C_INT), value :: nbx, nby, nbz, n_blocks, n_blocks_global, id_start
            integer(C_INT), value :: n_comp
            real(C_DOUBLE), intent(in) :: coef(*)
            integer(C_INT) :: ierr
        end function fdm_h5_case_append_coef

        function fdm_h5_case_append_coef_p(file_name, nbx, nby, nbz, n_blocks, &
                n_blocks_global, id_start, n_comp, coef) &
                bind(C, name="fdm_h5_case_append_coef_p") result(ierr)
            import :: C_CHAR, C_INT, C_DOUBLE
            character(kind=C_CHAR), intent(in) :: file_name(*)
            integer(C_INT), value :: nbx, nby, nbz, n_blocks, n_blocks_global, id_start
            integer(C_INT), value :: n_comp
            real(C_DOUBLE), intent(in) :: coef(*)
            integer(C_INT) :: ierr
        end function fdm_h5_case_append_coef_p

        function fdm_h5_case_append_dwall(file_name, nbx, nby, nbz, n_blocks, &
                n_blocks_global, id_start, dwall) &
                bind(C, name="fdm_h5_case_append_dwall") result(ierr)
            import :: C_CHAR, C_INT, C_DOUBLE
            character(kind=C_CHAR), intent(in) :: file_name(*)
            integer(C_INT), value :: nbx, nby, nbz, n_blocks, n_blocks_global, id_start
            real(C_DOUBLE), intent(in) :: dwall(*)
            integer(C_INT) :: ierr
        end function fdm_h5_case_append_dwall

        function fdm_h5_read_field(file_name, nbx, nby, nbz, n_blocks, &
                n_blocks_global, id_start, block_origin, &
                global_nx, global_ny, global_nz, n_var, var_names, found, q) &
                bind(C, name="fdm_h5_read_field") result(ierr)
            import :: C_CHAR, C_INT, C_DOUBLE
            character(kind=C_CHAR), intent(in) :: file_name(*)
            integer(C_INT), value :: nbx, nby, nbz, n_blocks, n_blocks_global, id_start
            integer(C_INT), intent(in) :: block_origin(*)
            integer(C_INT), value :: global_nx, global_ny, global_nz
            integer(C_INT), value :: n_var
            character(kind=C_CHAR), intent(in) :: var_names(*)
            ! Per-variable presence flag; q is left untouched where 0.
            integer(C_INT), intent(out) :: found(*)
            real(C_DOUBLE), intent(inout) :: q(*)
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

subroutine maybe_write_field(blk, dns, g, step, c, bc, pressure_niter, pressure_sor, nut, &
        rans_k, rans_omg, rans_gam, rans_ret, iddes_fd, scalar_names, conj_vfrac)
    type(block_set_type), intent(inout) :: blk
    type(dns_type), intent(in) :: dns
    type(grid_type), intent(in) :: g
    integer, intent(in) :: step
    type(comm_type), intent(in) :: c
    type(boundary_type), intent(in) :: bc
    integer(C_INT), intent(in) :: pressure_niter
    real(C_DOUBLE), intent(in) :: pressure_sor
    real(C_DOUBLE), allocatable, intent(in) :: nut(:,:,:,:)   ! LES eddy viscosity (unallocated when LES off)
    real(C_DOUBLE), allocatable, intent(in), optional :: rans_k(:,:,:,:), rans_omg(:,:,:,:)
    real(C_DOUBLE), allocatable, intent(in), optional :: rans_gam(:,:,:,:), rans_ret(:,:,:,:)
    real(C_DOUBLE), allocatable, intent(in), optional :: iddes_fd(:,:,:,:)
    character(len=*), intent(in), optional :: scalar_names(:)
    real(C_DOUBLE), allocatable, intent(in), optional :: conj_vfrac(:,:,:,:)

    if (.not. output_is_due(step, dns%field_interval)) return
    call write_field(blk, dns, g, step, c, bc, pressure_niter, pressure_sor, nut, &
        rans_k, rans_omg, rans_gam, rans_ret, iddes_fd, scalar_names, conj_vfrac)
end subroutine maybe_write_field

! Packed variable-name table for the C field writer/reader: n_var slots of
! VAR_NAME_LEN characters, each NUL-terminated. Slots 1-4 are always
! un/vn/wn/pn, so a scalar-free file is byte-identical to the historical one.
function field_var_names(dns, scalar_names) result(table)
    type(dns_type), intent(in) :: dns
    character(len=*), intent(in), optional :: scalar_names(:)
    character(kind=C_CHAR,len=:), allocatable :: table

    character(len=*), parameter :: base(4) = [character(len=4) :: "un", "vn", "wn", "pn"]
    integer :: v, nv, pos, ln
    character(len=VAR_NAME_LEN) :: name

    nv = int(dns%nVar)
    allocate(character(kind=C_CHAR,len=VAR_NAME_LEN*nv) :: table)
    do v = 1, VAR_NAME_LEN*nv
        table(v:v) = C_NULL_CHAR
    end do
    do v = 1, nv
        if (v <= 4) then
            name = base(v)
        else if (present(scalar_names)) then
            name = scalar_names(v - 4)
        else
            write(name, '(a,i0)') "s", v - 4
        end if
        ln = len_trim(name)
        pos = (v - 1)*VAR_NAME_LEN
        table(pos+1:pos+ln) = name(1:ln)
    end do
end function field_var_names

subroutine write_field(blk, dns, g, step, c, bc, pressure_niter, pressure_sor, nut, &
        rans_k, rans_omg, rans_gam, rans_ret, iddes_fd, scalar_names, conj_vfrac)
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
    real(C_DOUBLE), allocatable, intent(in) :: nut(:,:,:,:)   ! LES eddy viscosity (unallocated when LES off)
    ! RANS transport scalars (T2); appended when present AND allocated, so
    ! LES/no-model output stays byte-identical. The T4 transition pair is
    ! additionally gated on the config flag: the arrays are 1-cell dummies
    ! (for uniform device maps) when transition is off.
    real(C_DOUBLE), allocatable, intent(in), optional :: rans_k(:,:,:,:), rans_omg(:,:,:,:)
    real(C_DOUBLE), allocatable, intent(in), optional :: rans_gam(:,:,:,:), rans_ret(:,:,:,:)
    ! IDDES DDES-shielding function: a 1-cell dummy unless model = iddes
    ! (size gate below), so non-iddes output is unchanged.
    real(C_DOUBLE), allocatable, intent(in), optional :: iddes_fd(:,:,:,:)
    ! Passive-scalar dataset names (scalar.f90); absent = the s1..sN default.
    character(len=*), intent(in), optional :: scalar_names(:)
    ! Conjugate fluid VOLUME FRACTION (C3): a 1-cell dummy unless a conjugate
    ! scalar is configured (size gate below), so every other file is
    ! unchanged. It is pure geometry, so it is written for the reader's
    ! benefit and never read back -- a restart rebuilds it from the case file
    ! exactly as a cold start does. What it buys: the cell capacity
    ! C = f + (1-f)C_s the scheme conserves is a solver-side quantity that no
    ! checker can reconstruct to the last bit from the geometry alone, so the
    ! conservation and interface-heat gates read the solver's own f.
    real(C_DOUBLE), allocatable, intent(in), optional :: conj_vfrac(:,:,:,:)

    character(len=256) :: h5_file_name
    character(kind=C_CHAR,len=:), allocatable :: c_file_name, var_names
    integer(C_INT) :: ierr
    integer(C_INT) :: periodic(1:3)
    integer(C_INT) :: ibm_enabled

    periodic = merge(1_C_INT, 0_C_INT, bc%isPeriodic)
    ibm_enabled = merge(1_C_INT, 0_C_INT, dns%ibm_enabled)

    call make_output_filename(trim(dns%field_prefix), step, ".h5", h5_file_name)
    if (c%has_terminal) then
        print *, "current time step: ", step, "   field filename: ", trim(h5_file_name)
    end if

#ifdef USE_OPENMP_OFFLOAD
    !$omp target update from(blk%q)
#endif

    c_file_name = to_c_string(h5_file_name)
    var_names = field_var_names(dns, scalar_names)
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
        dns%nVar, var_names, blk%q)
    if (ierr /= 0_C_INT) then
        if (c%has_terminal) print *, "error: could not write HDF5 field file: ", trim(h5_file_name)
        error stop
    end if

    ! xz-quadtree file variant: mark the per-direction refinement mask so
    ! readers (restart cross-check, compare_fields.py reassembly) scale
    ! block origins per direction. xyz files stay byte-identical.
    if (any(dns%block_refine_mask == 0_C_INT)) then
        ierr = fdm_h5_append_refine_dims(c_file_name, dns%block_refine_mask)
        if (ierr /= 0_C_INT) then
            if (c%has_terminal) print *, "error: could not append refine_dims to: ", &
                trim(h5_file_name)
            error stop
        end if
    end if

    ! Append the LES eddy viscosity as a "nut" dataset (only when LES is active,
    ! i.e. nut allocated, so the no-LES file is byte-identical). Collective; nut
    ! goes out on the host.
    if (allocated(nut)) then
#ifdef USE_OPENMP_OFFLOAD
        !$omp target update from(nut)
#endif
        ierr = fdm_h5_append_nut(c_file_name, blk%nb(1), blk%nb(2), blk%nb(3), &
            blk%nBlocks, blk%nBlocksGlobal, blk%idStart, nut)
        if (ierr /= 0_C_INT) then
            if (c%has_terminal) print *, "error: could not append nut to: ", trim(h5_file_name)
            error stop
        end if
    end if

    if (present(rans_k)) then
        if (allocated(rans_k)) call append_scalar_field(c_file_name, "k", rans_k, blk, &
            h5_file_name, c%has_terminal)
    end if
    if (present(rans_omg)) then
        if (allocated(rans_omg)) call append_scalar_field(c_file_name, "omega", rans_omg, blk, &
            h5_file_name, c%has_terminal)
    end if
    if (present(rans_gam) .and. dns%rans_transition) then
        if (allocated(rans_gam)) call append_scalar_field(c_file_name, "gamma", rans_gam, blk, &
            h5_file_name, c%has_terminal)
    end if
    if (present(rans_ret) .and. dns%rans_transition) then
        if (allocated(rans_ret)) call append_scalar_field(c_file_name, "rethetat", rans_ret, blk, &
            h5_file_name, c%has_terminal)
    end if
    if (present(iddes_fd)) then
        if (allocated(iddes_fd)) then
            if (size(iddes_fd) > 1) call append_scalar_field(c_file_name, "fd", iddes_fd, blk, &
                h5_file_name, c%has_terminal)
        end if
    end if
    if (present(conj_vfrac)) then
        if (allocated(conj_vfrac)) then
            if (size(conj_vfrac) > 1) call append_scalar_field(c_file_name, "vfrac", &
                conj_vfrac, blk, h5_file_name, c%has_terminal)
        end if
    end if

    ! No XDMF for the block-table layout; reassemble with
    ! tools/compare_fields.py --export-global for visualization.
end subroutine write_field

subroutine append_scalar_field(c_file_name, name, s, blk, h5_file_name, has_terminal)
    character(kind=C_CHAR,len=*), intent(in) :: c_file_name
    character(len=*), intent(in) :: name
    real(C_DOUBLE), allocatable, intent(in) :: s(:,:,:,:)
    type(block_set_type), intent(in) :: blk
    character(len=*), intent(in) :: h5_file_name
    logical, intent(in) :: has_terminal

    character(kind=C_CHAR,len=:), allocatable :: c_name
    integer(C_INT) :: ierr

#ifdef USE_OPENMP_OFFLOAD
    !$omp target update from(s)
#endif
    c_name = to_c_string(name)
    ierr = fdm_h5_append_scalar(c_file_name, c_name, blk%nb(1), blk%nb(2), blk%nb(3), &
        blk%nBlocks, blk%nBlocksGlobal, blk%idStart, s)
    if (ierr /= 0_C_INT) then
        if (has_terminal) print *, "error: could not append ", name, " to: ", trim(h5_file_name)
        error stop
    end if
end subroutine append_scalar_field

! Read one named RANS scalar from a restart file; found = .false. when the
! dataset is absent (old restart), which the caller handles by
! reinitializing and warning.
subroutine read_scalar_field(blk, name, file_name, s, found, has_terminal)
    type(block_set_type), intent(in) :: blk
    character(len=*), intent(in) :: name, file_name
    real(C_DOUBLE), intent(inout) :: s(:,:,:,:)
    logical, intent(out) :: found
    logical, intent(in) :: has_terminal

    character(kind=C_CHAR,len=:), allocatable :: c_file_name, c_name
    integer(C_INT) :: ierr, ifound

    c_file_name = to_c_string(file_name)
    c_name = to_c_string(name)
    ierr = fdm_h5_read_scalar(c_file_name, c_name, blk%nb(1), blk%nb(2), blk%nb(3), &
        blk%nBlocks, blk%nBlocksGlobal, blk%idStart, ifound, s)
    if (ierr /= 0_C_INT) then
        if (has_terminal) print *, "error: could not read ", name, " from: ", trim(file_name)
        error stop
    end if
    found = ifound /= 0_C_INT
end subroutine read_scalar_field

subroutine read_restart_metadata(dns, g, bc, pressure_niter, pressure_sor, file_name, c, seen)
    type(dns_type), intent(inout) :: dns
    type(grid_type), intent(inout) :: g
    type(boundary_type), intent(inout) :: bc
    integer(C_INT), intent(inout) :: pressure_niter
    real(C_DOUBLE), intent(inout) :: pressure_sor
    character(len=*), intent(in) :: file_name
    type(comm_type), intent(in) :: c
    ! Which values the ini set explicitly; config wins over the restart file for these.
    type(config_seen_type), intent(in) :: seen

    character(kind=C_CHAR,len=:), allocatable :: c_file_name
    integer(C_INT) :: ierr
    integer(C_INT) :: file_nsteps
    integer(C_INT) :: periodic(1:3)
    integer(C_INT) :: ibm_enabled
    integer :: dir
    real(C_DOUBLE) :: input_cflmax, input_pecletmax, input_dtmax, input_t_final
    integer(C_INT) :: input_pressure_niter
    real(C_DOUBLE) :: input_pressure_sor
    integer(C_INT) :: input_bc_type(VAR_U:VAR_P,1:NFACES)
    real(C_DOUBLE) :: input_bc_value(VAR_U:VAR_P,1:NFACES)

    input_bc_type = bc%faceBcType
    input_bc_value = bc%faceBcDefaultValue
    file_nsteps = dns%nsteps
    periodic = merge(1_C_INT, 0_C_INT, bc%isPeriodic)
    ibm_enabled = merge(1_C_INT, 0_C_INT, dns%ibm_enabled)
    input_cflmax = dns%cflmax
    input_pecletmax = dns%pecletmax
    input_dtmax = dns%dtmax
    input_t_final = dns%t_final
    input_pressure_niter = pressure_niter
    input_pressure_sor = pressure_sor
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

    ! Config is authority: for every value the ini set explicitly, keep the ini
    ! value rather than the one stored in the restart file. This matters most for
    ! the pressure solver -- a restart from an old red-black run stores sor=1.5,
    ! which DIVERGES under damped Jacobi (sor>0.8).
    if (seen%cflmax)         dns%cflmax     = input_cflmax
    if (seen%pecletmax)      dns%pecletmax  = input_pecletmax
    if (seen%dtmax)          dns%dtmax      = input_dtmax
    if (seen%t_final)        dns%t_final    = input_t_final
    if (seen%pressure_niter) pressure_niter = input_pressure_niter
    if (seen%pressure_sor)   pressure_sor   = input_pressure_sor
    ! BC rows the ini set explicitly also beat the restart file's (and the
    ! patch-type contradiction check must see the ini's values).
    where (bc%faceBcTypeSet)  bc%faceBcType         = input_bc_type
    where (bc%faceBcValueSet) bc%faceBcDefaultValue = input_bc_value

    if (dns%nsteps <= 0_C_INT) dns%nsteps = file_nsteps
    do dir = 1, 3
        bc%isPeriodic(dir) = periodic(dir) /= 0_C_INT
    end do
    dns%ibm_enabled = ibm_enabled /= 0_C_INT
end subroutine read_restart_metadata

! The case file's LAYOUT (step 7-2): header attributes checked against the
! ini -- the file is the single source of truth for grid, block size and
! leaf table, and an ini that disagrees is stale, which is a hard error
! naming the key rather than a silent override -- then the leaf table and
! the level-0 node lines. An explicit [blocks] nb must equal the file's; an
! unset one takes the file's (the nb rule). The node lines are compared
! BIT FOR BIT with the ones init_grid built from [grid.*]: one comparison
! covers every grid key at once, and the file's lines are what the solver
! keeps. Older block-table files (retired mobygeom) carry no node lines;
! the ini's are used and a note is printed.
subroutine read_case_layout(dns, g, file_name, inputs, nb, refMask, rows, has_layout, has_terminal)
    type(dns_type), intent(inout) :: dns
    type(grid_type), intent(inout) :: g
    character(len=*), intent(in) :: file_name
    ! The ini's input echo (config.f90 case_input_echo), compared line by
    ! line with the one stored in the file: the first differing line names
    ! the stale key. Files without the attribute (pre-step-7-3, or written
    ! by the retired mobygeom) fall back to the explicit checks below.
    character(len=*), intent(in) :: inputs
    integer(C_INT), intent(out) :: nb(3), refMask(3)
    integer(C_INT), allocatable, intent(out) :: rows(:,:)
    ! False for a LEGACY global-layout coefficient file, which carries no
    ! leaf table (the caller then still builds the layout from the ini;
    ! that path goes with the legacy reader in step 7-4).
    logical, intent(out) :: has_layout
    logical, intent(in) :: has_terminal

    character(kind=C_CHAR,len=:), allocatable :: c_file_name
    integer(C_INT) :: ierr, nx, ny, nz, nbAuto, nbRanks, nLeaves, found
    real(C_DOUBLE) :: lx, ly, lz, re
    real(C_DOUBLE), allocatable :: xn(:), yn(:), zn(:)
    integer, parameter :: ECHO_CAP = 16384
    character(kind=C_CHAR, len=ECHO_CAP) :: fileEcho

    c_file_name = to_c_string(file_name)
    ierr = fdm_h5_case_read_layout(c_file_name, nx, ny, nz, lx, ly, lz, re, &
        nb, refMask, nbAuto, nbRanks, nLeaves)
    if (ierr /= 0_C_INT) then
        if (has_terminal) print *, "error: could not read the case-file layout of: ", trim(file_name)
        error stop
    end if
    has_layout = nLeaves > 0_C_INT
    if (.not. has_layout) then
        if (has_terminal) print *, "note: legacy coefficient file without a leaf table;", &
            " block layout from the ini: ", trim(file_name)
        return
    end if
    ierr = fdm_h5_case_read_inputs(c_file_name, fileEcho, int(ECHO_CAP, C_INT), found)
    if (ierr /= 0_C_INT) then
        if (has_terminal) print *, "error: could not read the case_inputs echo of: ", trim(file_name)
        error stop
    end if
    if (found /= 0_C_INT) then
        call compare_case_inputs(fileEcho(1:index(fileEcho, C_NULL_CHAR) - 1), inputs, &
            file_name, has_terminal)
    else if (has_terminal) then
        print *, "note: case file carries no case_inputs echo (pre-step-7-3 file);", &
            " staleness checked on grid/re/nb/refine_dims only: ", trim(file_name)
    end if
    if (nx /= dns%globalSize(1) .or. ny /= dns%globalSize(2) .or. nz /= dns%globalSize(3)) then
        if (has_terminal) print *, "case file grid", nx, ny, nz, "differs from [grid] nx/ny/nz", &
            dns%globalSize, " -- stale case file: ", trim(file_name)
        error stop "case file / [grid] nx ny nz mismatch"
    end if
    if (lx /= dns%leng(1) .or. ly /= dns%leng(2) .or. lz /= dns%leng(3)) then
        if (has_terminal) print *, "case file lengths", lx, ly, lz, "differ from [grid] lx/ly/lz", &
            dns%leng, " -- stale case file: ", trim(file_name)
        error stop "case file / [grid] lx ly lz mismatch"
    end if
    if (re /= dns%re) then
        if (has_terminal) print *, "case file re", re, "differs from [flow] re", dns%re, &
            " -- stale case file (the coefficients carry 1/Re): ", trim(file_name)
        error stop "case file / [flow] re mismatch"
    end if
    if (any(refMask /= dns%block_refine_mask)) then
        if (has_terminal) print *, "case file refine_dims mask", refMask, &
            "differs from [blocks] refine_dims", dns%block_refine_mask, ": ", trim(file_name)
        error stop "case file / [blocks] refine_dims mismatch"
    end if
    if (all(dns%block_nb > 0_C_INT)) then
        if (any(nb /= dns%block_nb)) then
            if (has_terminal) print *, "case file block size", nb, "differs from [blocks] nb", &
                dns%block_nb, " -- stale case file: ", trim(file_name)
            error stop "case file / [blocks] nb mismatch"
        end if
    else
        dns%block_nb = nb
        dns%block_nb_auto = nbAuto /= 0_C_INT
        dns%block_nb_ranks = nbRanks
        if (has_terminal) then
            if (nbAuto /= 0_C_INT) then
                print '(A,3(1X,I0),A,I0,A)', " [blocks] nb from the case file:", nb, &
                    "  (one block per rank of a ", nbRanks, "-rank layout)"
            else
                print '(A,3(1X,I0))', " [blocks] nb from the case file:", nb
            end if
        end if
    end if

    allocate(rows(4, nLeaves))
    ierr = fdm_h5_case_read_blocks(c_file_name, nLeaves, rows)
    if (ierr /= 0_C_INT) then
        if (has_terminal) print *, "error: could not read the blocks table of: ", trim(file_name)
        error stop
    end if

    allocate(xn(0:int(nx)), yn(0:int(ny)), zn(0:int(nz)))
    ierr = fdm_h5_case_read_grid(c_file_name, nx, ny, nz, xn, yn, zn, found)
    if (ierr /= 0_C_INT) then
        if (has_terminal) print *, "error: could not read the node lines of: ", trim(file_name)
        error stop
    end if
    if (found /= 0_C_INT) then
        if (any(xn /= g%xNode) .or. any(yn /= g%yNode) .or. any(zn /= g%zNode)) then
            if (has_terminal) print *, "case file node lines differ from the [grid.x/y/z]", &
                " distribution -- stale case file: ", trim(file_name)
            error stop "case file / [grid.*] node-line mismatch"
        end if
        g%xNode = xn
        g%yNode = yn
        g%zNode = zn
    else if (has_terminal) then
        print *, "note: case file carries no node lines (pre-P3 block-table file);", &
            " using the [grid.*] lines: ", trim(file_name)
    end if
end subroutine read_case_layout

! Line-by-line comparison of two input echoes; the first difference is
! reported with its key and both values, then the run stops.
subroutine compare_case_inputs(fileEcho, iniEcho, file_name, has_terminal)
    character(len=*), intent(in) :: fileEcho, iniEcho, file_name
    logical, intent(in) :: has_terminal

    integer :: fa, fb, ia, ib, nf, ni
    character(len=:), allocatable :: fline, iline

    if (fileEcho == iniEcho) return
    fa = 1; ia = 1
    nf = len(fileEcho); ni = len(iniEcho)
    do
        if (fa > nf .and. ia > ni) exit
        fb = index(fileEcho(fa:), new_line("a"))
        ib = index(iniEcho(ia:), new_line("a"))
        if (fb == 0) then; fline = fileEcho(fa:nf); fa = nf + 1
        else; fline = fileEcho(fa:fa+fb-2); fa = fa + fb; end if
        if (ib == 0) then; iline = iniEcho(ia:ni); ia = ni + 1
        else; iline = iniEcho(ia:ia+ib-2); ia = ia + ib; end if
        if (fline /= iline) then
            if (has_terminal) then
                print *, "case file is STALE: ", trim(file_name)
                print *, "   file has: ", fline
                print *, "   ini  has: ", iline
                print *, "   re-run moby_prepare, or moby_solve --prepare"
            end if
            error stop "case file does not match the ini"
        end if
    end do
end subroutine compare_case_inputs

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

! Per-level refinement masks (mobygeom block-table); found is false when
! the file carries none.
subroutine read_block_masks(touch, buried, level, n_raster, found, dns, has_terminal)
    integer(C_INT), intent(out) :: touch(:), buried(:)
    integer, intent(in) :: level, n_raster
    logical, intent(out) :: found
    type(dns_type), intent(in) :: dns
    logical, intent(in) :: has_terminal

    character(kind=C_CHAR,len=:), allocatable :: c_file_name
    integer(C_INT) :: ierr, c_found

    c_file_name = to_c_string(dns%ibm_coeff_file)
    ierr = fdm_h5_read_block_masks(c_file_name, int(level, C_INT), int(n_raster, C_INT), &
        dns%block_nb, c_found, touch, buried)
    if (ierr /= 0_C_INT) then
        if (has_terminal) print *, "error: could not read block masks from: ", &
            trim(dns%ibm_coeff_file)
        error stop
    end if
    found = c_found /= 0_C_INT
end subroutine read_block_masks

! Per-level geometry-mask window (block coords) of a block-table file;
! has_win false on legacy full-raster files.
subroutine read_mask_window(lo, dims, has_win, level, dns, has_terminal)
    integer(C_INT), intent(out) :: lo(3), dims(3)
    logical, intent(out) :: has_win
    integer, intent(in) :: level
    type(dns_type), intent(in) :: dns
    logical, intent(in) :: has_terminal

    character(kind=C_CHAR,len=:), allocatable :: c_file_name
    integer(C_INT) :: ierr, c_has

    c_file_name = to_c_string(dns%ibm_coeff_file)
    ierr = fdm_h5_read_mask_window(c_file_name, int(level, C_INT), lo, dims, c_has)
    if (ierr /= 0_C_INT) then
        if (has_terminal) print *, "error: could not read mask windows from: ", &
            trim(dns%ibm_coeff_file)
        error stop
    end if
    has_win = c_has /= 0_C_INT
end subroutine read_mask_window

! Per-leaf wall-distance tiles from the case file (dataset dwall_blocks,
! row id = leaf id of the file's own blocks table). found is false when
! the file carries none.
subroutine read_dwall_blocks(dwall, found, dns, blk, has_terminal)
    real(C_DOUBLE), intent(inout) :: dwall(*)
    logical, intent(out) :: found
    type(dns_type), intent(in) :: dns
    type(block_set_type), intent(in) :: blk
    logical, intent(in) :: has_terminal

    character(kind=C_CHAR,len=:), allocatable :: c_file_name
    integer(C_INT) :: ierr, c_found

    c_file_name = to_c_string(dns%ibm_coeff_file)
    ierr = fdm_h5_read_dwall_blocks(c_file_name, blk%nb(1), blk%nb(2), blk%nb(3), &
        blk%nBlocks, blk%idStart, c_found, dwall)
    if (ierr /= 0_C_INT) then
        if (has_terminal) print *, "error: could not read dwall_blocks from: ", &
            trim(dns%ibm_coeff_file)
        error stop
    end if
    found = c_found /= 0_C_INT
end subroutine read_dwall_blocks

! RANS geometry diagnostic dump ([rans] dump_geometry): blocks table plus
! interior dwall/yeff/wallcell and the per-block cell-centre coordinates.
! Parallel HDF5 call: all MPI ranks must enter this routine together.
subroutine write_rans_geometry_file(file_name, blk, dwall, yeff, wallcell, &
        xc, yc, zc, has_terminal)
    character(len=*), intent(in) :: file_name
    type(block_set_type), intent(in) :: blk
    real(C_DOUBLE), intent(in) :: dwall(*), yeff(*)
    integer(C_INT), intent(in) :: wallcell(*)
    real(C_DOUBLE), intent(in) :: xc(*), yc(*), zc(*)
    logical, intent(in) :: has_terminal

    character(kind=C_CHAR,len=:), allocatable :: c_file_name
    integer(C_INT) :: ierr

    c_file_name = to_c_string(file_name)
    ierr = fdm_h5_write_rans_geometry(c_file_name, blk%nb(1), blk%nb(2), blk%nb(3), &
        blk%nBlocks, blk%nBlocksGlobal, blk%idStart, blk%origin, blk%level, &
        dwall, yeff, wallcell, xc, yc, zc)
    if (ierr /= 0_C_INT) then
        if (has_terminal) print *, "error: could not write RANS geometry file: ", &
            trim(file_name)
        error stop
    end if
end subroutine write_rans_geometry_file

! moby_prepare output (docs/prepare_solve_strategy.md P0): one case file in
! the block-table coefficient-file format the solver reads via [ibm]
! coeff_file -- header attributes + blocks leaf table + coef_blocks, plus
! the per-level refinement masks (refine_body), block_active (remove_solid)
! and dwall_blocks ([rans]). Parallel HDF5: all ranks enter together; each
! rank writes its own contiguous leaf-row range, rank 0 the lattice-global
! rasters (full rasters, no window attrs -- the analytic convention).
subroutine write_case_file(file_name, blk, dns, g, bc, c, nCoefComp, inputs, has_terminal, &
        coef, touch, buried, maskDims, active, dwall, maskLo)
    character(len=*), intent(in) :: file_name
    ! The input echo (config.f90 case_input_echo), stored for the read-side
    ! staleness check.
    character(len=*), intent(in) :: inputs
    type(block_set_type), intent(in) :: blk
    type(dns_type), intent(in) :: dns
    type(grid_type), intent(in) :: g
    type(boundary_type), intent(in) :: bc
    type(comm_type), intent(in) :: c
    ! Component extent of coef: 3, or 4 when the case declares passive
    ! scalars, in which case the cell-centred column is written as the
    ! separate OPTIONAL coef_p_blocks dataset (increment S3) and every
    ! other dataset in the file is untouched. A BODY-FREE case (step 7-1)
    ! passes no coef at all: the file then carries attrs + node lines +
    ! the leaf table and no coefficient/mask/dwall dataset.
    integer(C_INT), intent(in) :: nCoefComp
    logical, intent(in) :: has_terminal
    real(C_DOUBLE), intent(in), optional :: coef(*)
    integer(C_INT), intent(in), optional :: touch(:,:), buried(:,:), maskDims(:,:)
    integer(C_INT), intent(in), optional :: active(:)
    real(C_DOUBLE), intent(in), optional :: dwall(:,:,:,:)
    ! Per-level mask window origins (deep-refinement WINDOWED rasters):
    ! window attrs are written for any level whose raster does not cover
    ! the full lattice, so the readers reconstruct placement.
    integer(C_INT), intent(in), optional :: maskLo(:,:)

    character(kind=C_CHAR,len=:), allocatable :: c_file_name
    integer(C_INT) :: ierr, write_data, n_raster
    integer(C_INT) :: periodic(1:3), natural_one_sided(1:3)
    integer(C_INT) :: nb_auto
    integer :: level

    c_file_name = to_c_string(file_name)
    ! has_terminal in comm_type is world rank 0 -- the one raster writer.
    write_data = merge(1_C_INT, 0_C_INT, c%has_terminal)

    nb_auto = merge(1_C_INT, 0_C_INT, dns%block_nb_auto)
    ierr = fdm_h5_case_create(c_file_name, &
        dns%globalSize(1), dns%globalSize(2), dns%globalSize(3), &
        dns%leng(1), dns%leng(2), dns%leng(3), dns%re, &
        dns%block_nb, blk%nLevels - 1_C_INT, dns%block_refine_mask, &
        nb_auto, dns%block_nb_ranks, to_c_string(inputs), &
        blk%nBlocksGlobal, blk%idStart, blk%nBlocks, blk%origin, blk%level)
    call check_case_write(ierr, "create", file_name, has_terminal)

    ! Grid datasets (P3): the node lines + mobygrid-format attributes, so
    ! the case file replaces the retired mobygrid's grid handshake.
    periodic = merge(1_C_INT, 0_C_INT, bc%isPeriodic)
    natural_one_sided = merge(1_C_INT, 0_C_INT, g%natural_one_sided)
    ierr = fdm_h5_case_append_grid(c_file_name, &
        dns%globalSize(1), dns%globalSize(2), dns%globalSize(3), &
        periodic, g%distribution(1:3), g%stretch(1:3), g%natural_dyw_plus(1:3), &
        natural_one_sided, g%xNode(0:int(dns%globalSize(1))), &
        g%yNode(0:int(dns%globalSize(2))), g%zNode(0:int(dns%globalSize(3))), &
        write_data)
    call check_case_write(ierr, "grid datasets", file_name, has_terminal)

    if (present(touch)) then
        do level = 0, int(blk%nLevels) - 1
            n_raster = product(maskDims(:, level+1))
            ierr = fdm_h5_case_append_masks(c_file_name, int(level, C_INT), &
                n_raster, touch(1:n_raster, level+1), buried(1:n_raster, level+1), &
                write_data)
            call check_case_write(ierr, "refinement masks", file_name, has_terminal)
            if (present(maskLo)) then
                ! full-lattice rasters keep the legacy attr-free layout
                if (any(maskLo(:, level+1) /= 0_C_INT) .or. &
                    any(maskDims(:, level+1) /= (dns%globalSize/dns%block_nb) &
                        *2**(int(level, C_INT)*dns%block_refine_mask))) then
                    ierr = fdm_h5_case_append_mask_window(c_file_name, &
                        int(level, C_INT), maskLo(:, level+1), maskDims(:, level+1))
                    call check_case_write(ierr, "mask windows", file_name, has_terminal)
                end if
            end if
        end do
    end if
    if (present(active)) then
        ierr = fdm_h5_case_append_active(c_file_name, int(size(active), C_INT), &
            active, write_data)
        call check_case_write(ierr, "block_active", file_name, has_terminal)
    end if

    if (present(coef)) then
        ierr = fdm_h5_case_append_coef(c_file_name, blk%nb(1), blk%nb(2), blk%nb(3), &
            blk%nBlocks, blk%nBlocksGlobal, blk%idStart, nCoefComp, coef)
        call check_case_write(ierr, "coef_blocks", file_name, has_terminal)

        if (nCoefComp > 3_C_INT) then
            ierr = fdm_h5_case_append_coef_p(c_file_name, blk%nb(1), blk%nb(2), blk%nb(3), &
                blk%nBlocks, blk%nBlocksGlobal, blk%idStart, nCoefComp, coef)
            call check_case_write(ierr, "coef_p_blocks", file_name, has_terminal)
        end if
    end if

    if (present(dwall)) then
        ierr = fdm_h5_case_append_dwall(c_file_name, blk%nb(1), blk%nb(2), blk%nb(3), &
            blk%nBlocks, blk%nBlocksGlobal, blk%idStart, dwall)
        call check_case_write(ierr, "dwall_blocks", file_name, has_terminal)
    end if
end subroutine write_case_file

subroutine check_case_write(ierr, what, file_name, has_terminal)
    integer(C_INT), intent(in) :: ierr
    character(len=*), intent(in) :: what, file_name
    logical, intent(in) :: has_terminal

    if (ierr == 0_C_INT) return
    if (has_terminal) print *, "error: could not write ", what, &
        " to case file: ", trim(file_name)
    error stop
end subroutine check_case_write

subroutine read_field(blk, dns, file_name, c, scalar_names)
    ! Parallel HDF5 call: all MPI ranks must enter this routine together.
    type(block_set_type), intent(inout) :: blk
    type(dns_type), intent(in) :: dns
    character(len=*), intent(in) :: file_name
    type(comm_type), intent(in) :: c
    character(len=*), intent(in), optional :: scalar_names(:)

    character(kind=C_CHAR,len=:), allocatable :: c_file_name, var_names
    integer(C_INT) :: ierr, file_mask(1:3), has_blocks
    integer(C_INT), allocatable :: found(:)
    integer :: v

    c_file_name = to_c_string(file_name)
    ! The refine_dims variants store block origins in different index
    ! spaces (xz: y origins in GLOBAL cells) — a mixed BLOCK-layout
    ! restart must be a hard error, not a silent misplacement (the
    ! blocks-table row check alone could alias on y-symmetric layouts).
    ! Legacy global-3D files carry no block layout: any leaf table can
    ! slice them, so the check does not apply.
    ierr = fdm_h5_read_refine_dims(c_file_name, file_mask, has_blocks)
    if (ierr /= 0_C_INT) then
        if (c%has_terminal) print *, "error: could not read refine_dims from: ", trim(file_name)
        error stop
    end if
    if (has_blocks /= 0_C_INT .and. any(file_mask /= dns%block_refine_mask)) then
        if (c%has_terminal) print *, "restart file refine_dims mask", file_mask, &
            "does not match the configured [blocks] refine_dims", dns%block_refine_mask
        error stop "restart/config [blocks] refine_dims mismatch"
    end if
    ! A block-layout snapshot is sliced by leaf id, so its leaf table must
    ! BE the solver's (row for row over this rank's id range): a snapshot
    ! from a differently tiled or refined run would otherwise land silently
    ! on the wrong blocks. Until step 7-2 only the dataset extents were
    ! checked, which a same-count, different-layout table passes.
    ierr = fdm_h5_check_block_table(c_file_name, blk%nBlocks, blk%idStart, blk%origin, blk%level)
    if (ierr == 2_C_INT) then
        if (c%has_terminal) print *, "error: restart file block table does not match", &
            " the case layout (a snapshot of a different [blocks] tiling?): ", trim(file_name)
        error stop "restart/case block-table mismatch"
    end if
    if (ierr /= 0_C_INT) then
        if (c%has_terminal) print *, "error: could not read the blocks table of: ", trim(file_name)
        error stop
    end if
    var_names = field_var_names(dns, scalar_names)
    allocate(found(int(dns%nVar)))
    ierr = fdm_h5_read_field(c_file_name, blk%nb(1), blk%nb(2), blk%nb(3), &
        blk%nBlocks, blk%nBlocksGlobal, blk%idStart, blk%origin, &
        dns%globalSize(1), dns%globalSize(2), dns%globalSize(3), &
        dns%nVar, var_names, found, blk%q)
    if (ierr /= 0_C_INT) then
        if (c%has_terminal) print *, "error: could not read HDF5 field file: ", trim(file_name)
        error stop
    end if
    if (any(found(1:4) == 0_C_INT)) then
        if (c%has_terminal) print *, "error: restart file has no un/vn/wn/pn datasets: ", &
            trim(file_name)
        error stop
    end if
    ! A missing scalar dataset (a scalar added after the file was written, or
    ! an older restart) leaves the slot at its [scalar.N] initial condition --
    ! the RANS named-scalar precedent.
    do v = 5, int(dns%nVar)
        if (found(v) == 0_C_INT .and. c%has_terminal) then
            print *, "warning: restart file has no dataset for scalar", v - 4, &
                "; initialising it from [scalar.N] initial"
        end if
    end do
    deallocate(found)

#ifdef USE_OPENMP_OFFLOAD
    !$omp target update to(blk%q)
#endif
end subroutine read_field

subroutine read_force_file(f, blk, dns, file_name, has_terminal)
    ! Read a volumetric body force from an HDF5 field laid out like a
    ! velocity field: un/vn/wn = fx/fy/fz at the staggered points, pn ignored.
    ! Reuses the block-table field reader into a temporary q-shaped buffer,
    ! then copies the interior velocity components into the force array
    ! f(1:nb,1:nb,1:nb,NVEL,nBlocks). Host side; the caller maps f to device.
    ! Parallel HDF5 call: all MPI ranks must enter together.
    real(C_DOUBLE), intent(out) :: f(:,:,:,:,:)
    type(block_set_type), intent(in) :: blk
    type(dns_type), intent(in) :: dns
    character(len=*), intent(in) :: file_name
    logical, intent(in) :: has_terminal

    character(kind=C_CHAR,len=:), allocatable :: c_file_name
    real(C_DOUBLE), allocatable :: qtmp(:,:,:,:,:)
    integer(C_INT) :: ierr, found(NVAR)
    integer :: nx, ny, nz
    ! A velocity-layout file, NOT the scalar-extended one: this read stays at
    ! NVAR variables whatever dns%nVar is.
    character(kind=C_CHAR,len=VAR_NAME_LEN*NVAR) :: var_names

    if (len_trim(file_name) == 0) error stop "[force] type = file needs [force] file = path"

    nx = int(blk%nb(1)); ny = int(blk%nb(2)); nz = int(blk%nb(3))
    allocate(qtmp(0:nx+1, 0:ny+1, 0:nz+1, NVAR, blk%nBlocks))
    qtmp = 0.0d0

    var_names = repeat(C_NULL_CHAR, VAR_NAME_LEN*NVAR)
    var_names(1:2) = "un"
    var_names(VAR_NAME_LEN+1:VAR_NAME_LEN+2) = "vn"
    var_names(2*VAR_NAME_LEN+1:2*VAR_NAME_LEN+2) = "wn"
    var_names(3*VAR_NAME_LEN+1:3*VAR_NAME_LEN+2) = "pn"

    c_file_name = to_c_string(file_name)
    ierr = fdm_h5_read_field(c_file_name, blk%nb(1), blk%nb(2), blk%nb(3), &
        blk%nBlocks, blk%nBlocksGlobal, blk%idStart, blk%origin, &
        dns%globalSize(1), dns%globalSize(2), dns%globalSize(3), &
        NVAR, var_names, found, qtmp)
    if (ierr /= 0_C_INT .or. any(found(VAR_U:VAR_W) == 0_C_INT)) then
        if (has_terminal) print *, "error: could not read force field file: ", trim(file_name)
        error stop
    end if

    f(1:nx,1:ny,1:nz,VAR_U:VAR_W,:) = qtmp(1:nx,1:ny,1:nz,VAR_U:VAR_W,:)
    deallocate(qtmp)
end subroutine read_force_file

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
