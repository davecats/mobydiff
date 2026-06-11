#include <mpi.h>
#include <hdf5.h>
#include <stdlib.h>
#include <stdio.h>

#ifndef H5_HAVE_PARALLEL
#error "field_hdf5.c requires HDF5 built with parallel MPI-IO support"
#endif

static size_t linear_fortran(size_t i, size_t j, size_t k, size_t ni, size_t nj)
{
    return i + ni*(j + nj*k);
}

static size_t linear_hdf5(size_t i, size_t j, size_t k, size_t nj, size_t nk)
{
    return (i*nj + j)*nk + k;
}

static size_t linear_fortran4(size_t i, size_t j, size_t k, size_t v,
                              size_t ni, size_t nj, size_t nk)
{
    return i + ni*(j + nj*(k + nk*v));
}

static int double_mismatch(double a, double b)
{
    double diff = a - b;
    double scale = a;

    if (diff < 0.0) diff = -diff;
    if (scale < 0.0) scale = -scale;
    if (b < 0.0) {
        if (-b > scale) scale = -b;
    } else {
        if (b > scale) scale = b;
    }
    if (scale < 1.0) scale = 1.0;

    return diff > 1.0e-10*scale;
}

static hid_t create_parallel_file(const char *filename)
{
    hid_t plist = H5Pcreate(H5P_FILE_ACCESS);
    hid_t file;

    if (plist < 0) return -1;
    if (H5Pset_fapl_mpio(plist, MPI_COMM_WORLD, MPI_INFO_NULL) < 0) {
        H5Pclose(plist);
        return -1;
    }

    file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist);
    H5Pclose(plist);
    return file;
}

static hid_t open_parallel_file(const char *filename)
{
    hid_t plist = H5Pcreate(H5P_FILE_ACCESS);
    hid_t file;

    if (plist < 0) return -1;
    if (H5Pset_fapl_mpio(plist, MPI_COMM_WORLD, MPI_INFO_NULL) < 0) {
        H5Pclose(plist);
        return -1;
    }

    file = H5Fopen(filename, H5F_ACC_RDONLY, plist);
    H5Pclose(plist);
    return file;
}

static hid_t create_serial_file(const char *filename)
{
    return H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
}

static int write_attr_int(hid_t file, const char *name, int value)
{
    hid_t space = H5Screate(H5S_SCALAR);
    hid_t attr;
    herr_t status;

    if (space < 0) return 1;
    attr = H5Acreate2(file, name, H5T_NATIVE_INT, space, H5P_DEFAULT, H5P_DEFAULT);
    if (attr < 0) {
        H5Sclose(space);
        return 1;
    }

    status = H5Awrite(attr, H5T_NATIVE_INT, &value);
    H5Aclose(attr);
    H5Sclose(space);
    return status < 0;
}

static int write_attr_double(hid_t file, const char *name, double value)
{
    hid_t space = H5Screate(H5S_SCALAR);
    hid_t attr;
    herr_t status;

    if (space < 0) return 1;
    attr = H5Acreate2(file, name, H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT);
    if (attr < 0) {
        H5Sclose(space);
        return 1;
    }

    status = H5Awrite(attr, H5T_NATIVE_DOUBLE, &value);
    H5Aclose(attr);
    H5Sclose(space);
    return status < 0;
}

static int write_attr_int_array(hid_t file, const char *name, const int *values, hsize_t n)
{
    hid_t space = H5Screate_simple(1, &n, NULL);
    hid_t attr;
    herr_t status;

    if (space < 0) return 1;
    attr = H5Acreate2(file, name, H5T_NATIVE_INT, space, H5P_DEFAULT, H5P_DEFAULT);
    if (attr < 0) {
        H5Sclose(space);
        return 1;
    }

    status = H5Awrite(attr, H5T_NATIVE_INT, values);
    H5Aclose(attr);
    H5Sclose(space);
    return status < 0;
}

static int write_attr_double_array(hid_t file, const char *name, const double *values, hsize_t n)
{
    hid_t space = H5Screate_simple(1, &n, NULL);
    hid_t attr;
    herr_t status;

    if (space < 0) return 1;
    attr = H5Acreate2(file, name, H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT);
    if (attr < 0) {
        H5Sclose(space);
        return 1;
    }

    status = H5Awrite(attr, H5T_NATIVE_DOUBLE, values);
    H5Aclose(attr);
    H5Sclose(space);
    return status < 0;
}

static int read_attr_int(hid_t file, const char *name, int *value, int required)
{
    htri_t exists = H5Aexists(file, name);
    hid_t attr;
    herr_t status;

    if (exists <= 0) return required ? 1 : 0;

    attr = H5Aopen(file, name, H5P_DEFAULT);
    if (attr < 0) return 1;

    status = H5Aread(attr, H5T_NATIVE_INT, value);
    H5Aclose(attr);
    return status < 0;
}

static int read_attr_double(hid_t file, const char *name, double *value, int required)
{
    htri_t exists = H5Aexists(file, name);
    hid_t attr;
    herr_t status;

    if (exists <= 0) return required ? 1 : 0;

    attr = H5Aopen(file, name, H5P_DEFAULT);
    if (attr < 0) return 1;

    status = H5Aread(attr, H5T_NATIVE_DOUBLE, value);
    H5Aclose(attr);
    return status < 0;
}

static int read_attr_int_array(hid_t file, const char *name, int *values, hsize_t n, int required)
{
    htri_t exists = H5Aexists(file, name);
    hid_t attr;
    hid_t space;
    hsize_t dims[1] = {0};
    herr_t status;

    if (exists <= 0) return required ? 1 : 0;

    attr = H5Aopen(file, name, H5P_DEFAULT);
    if (attr < 0) return 1;

    space = H5Aget_space(attr);
    if (space < 0 || H5Sget_simple_extent_ndims(space) != 1) {
        if (space >= 0) H5Sclose(space);
        H5Aclose(attr);
        return 1;
    }

    H5Sget_simple_extent_dims(space, dims, NULL);
    if (dims[0] != n) {
        H5Sclose(space);
        H5Aclose(attr);
        return 1;
    }

    status = H5Aread(attr, H5T_NATIVE_INT, values);
    H5Sclose(space);
    H5Aclose(attr);
    return status < 0;
}

static int read_attr_double_array(hid_t file, const char *name, double *values, hsize_t n, int required)
{
    htri_t exists = H5Aexists(file, name);
    hid_t attr;
    hid_t space;
    hsize_t dims[1] = {0};
    herr_t status;

    if (exists <= 0) return required ? 1 : 0;

    attr = H5Aopen(file, name, H5P_DEFAULT);
    if (attr < 0) return 1;

    space = H5Aget_space(attr);
    if (space < 0 || H5Sget_simple_extent_ndims(space) != 1) {
        if (space >= 0) H5Sclose(space);
        H5Aclose(attr);
        return 1;
    }

    H5Sget_simple_extent_dims(space, dims, NULL);
    if (dims[0] != n) {
        H5Sclose(space);
        H5Aclose(attr);
        return 1;
    }

    status = H5Aread(attr, H5T_NATIVE_DOUBLE, values);
    H5Sclose(space);
    H5Aclose(attr);
    return status < 0;
}

/*
 * Global block table, rows indexed by global block id: zero-based cell
 * origin (3) and refinement level. Each rank writes its own contiguous id
 * range with an independent transfer under the collectively created
 * dataset.
 */
static int write_block_table(hid_t file, int n_blocks_global, int id_start,
                             int n_blocks, const int *block_origin,
                             const int *block_level)
{
    hsize_t dims[2] = {(hsize_t)n_blocks_global, 4};
    hsize_t start[2] = {(hsize_t)id_start, 0};
    hsize_t count[2] = {(hsize_t)n_blocks, 4};
    hid_t file_space = -1;
    hid_t mem_space = -1;
    hid_t dset = -1;
    hid_t xfer = -1;
    int *rows = NULL;
    herr_t status;

    rows = (int *)malloc((size_t)n_blocks*4*sizeof(int));
    if (rows == NULL) return 1;
    for (int b = 0; b < n_blocks; ++b) {
        rows[4*b + 0] = block_origin[3*b + 0];
        rows[4*b + 1] = block_origin[3*b + 1];
        rows[4*b + 2] = block_origin[3*b + 2];
        rows[4*b + 3] = block_level[b];
    }

    file_space = H5Screate_simple(2, dims, NULL);
    if (file_space < 0) {
        free(rows);
        return 1;
    }

    dset = H5Dcreate2(file, "blocks", H5T_NATIVE_INT, file_space,
                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dset < 0) {
        H5Sclose(file_space);
        free(rows);
        return 1;
    }

    mem_space = H5Screate_simple(2, count, NULL);
    xfer = H5Pcreate(H5P_DATASET_XFER);
    if (mem_space < 0 || xfer < 0 ||
        H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, count, NULL) < 0) {
        if (mem_space >= 0) H5Sclose(mem_space);
        if (xfer >= 0) H5Pclose(xfer);
        H5Dclose(dset);
        H5Sclose(file_space);
        free(rows);
        return 1;
    }
    H5Pset_dxpl_mpio(xfer, H5FD_MPIO_INDEPENDENT);

    status = H5Dwrite(dset, H5T_NATIVE_INT, mem_space, file_space, xfer, rows);

    H5Pclose(xfer);
    H5Sclose(mem_space);
    H5Dclose(dset);
    H5Sclose(file_space);
    free(rows);
    return status < 0;
}

/*
 * One hyperslab per block into the unchanged global dataset. Dataset
 * creation is collective; the per-block transfers are independent because
 * ranks own different block counts.
 */
static int write_global_dataset_blocks(hid_t file, const char *name,
                                       int nbx, int nby, int nbz,
                                       int n_blocks, const int *block_origin,
                                       int global_nx, int global_ny, int global_nz,
                                       const double *q, size_t block_stride)
{
    const size_t ni = (size_t)nbx + 2;
    const size_t nj = (size_t)nby + 2;
    const size_t n = (size_t)nbx*(size_t)nby*(size_t)nbz;
    hsize_t global_dims[3] = {(hsize_t)global_nz, (hsize_t)global_ny, (hsize_t)global_nx};
    hsize_t local_dims[3] = {(hsize_t)nbz, (hsize_t)nby, (hsize_t)nbx};
    hid_t file_space = -1;
    hid_t mem_space = -1;
    hid_t dset = -1;
    hid_t xfer = -1;
    double *buffer = NULL;
    herr_t status = 0;
    int ierr = 0;

    buffer = (double *)malloc(n*sizeof(double));
    if (buffer == NULL) return 1;

    file_space = H5Screate_simple(3, global_dims, NULL);
    if (file_space < 0) {
        free(buffer);
        return 1;
    }

    dset = H5Dcreate2(file, name, H5T_NATIVE_DOUBLE, file_space,
                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dset < 0) {
        H5Sclose(file_space);
        free(buffer);
        return 1;
    }

    mem_space = H5Screate_simple(3, local_dims, NULL);
    xfer = H5Pcreate(H5P_DATASET_XFER);
    if (mem_space < 0 || xfer < 0) {
        if (mem_space >= 0) H5Sclose(mem_space);
        if (xfer >= 0) H5Pclose(xfer);
        H5Dclose(dset);
        H5Sclose(file_space);
        free(buffer);
        return 1;
    }
    H5Pset_dxpl_mpio(xfer, H5FD_MPIO_INDEPENDENT);

    for (int b = 0; b < n_blocks; ++b) {
        const double *field = q + (size_t)b*block_stride;
        hsize_t start[3] = {
            (hsize_t)block_origin[3*b + 2],
            (hsize_t)block_origin[3*b + 1],
            (hsize_t)block_origin[3*b + 0]
        };

        for (size_t k = 1; k <= (size_t)nbz; ++k) {
            for (size_t j = 1; j <= (size_t)nby; ++j) {
                for (size_t i = 1; i <= (size_t)nbx; ++i) {
                    buffer[linear_hdf5(k-1, j-1, i-1, (size_t)nby, (size_t)nbx)] =
                        field[linear_fortran(i, j, k, ni, nj)];
                }
            }
        }

        if (H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, local_dims, NULL) < 0) {
            ierr = 1;
            break;
        }
        status = H5Dwrite(dset, H5T_NATIVE_DOUBLE, mem_space, file_space, xfer, buffer);
        if (status < 0) {
            ierr = 1;
            break;
        }
    }

    H5Pclose(xfer);
    H5Sclose(mem_space);
    H5Dclose(dset);
    H5Sclose(file_space);
    free(buffer);
    return ierr;
}

static int write_coord_dataset(hid_t file, const char *name, int n, int rank, const double *coord)
{
    hsize_t dims[1] = {(hsize_t)n};
    hid_t space = -1;
    hid_t dset = -1;
    herr_t status = 0;

    space = H5Screate_simple(1, dims, NULL);
    if (space < 0) return 1;

    dset = H5Dcreate2(file, name, H5T_NATIVE_DOUBLE, space,
                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dset < 0) {
        H5Sclose(space);
        return 1;
    }

    if (rank == 0) {
        status = H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, coord);
    }

    H5Dclose(dset);
    H5Sclose(space);
    return status < 0;
}

static int create_group(hid_t file, const char *name)
{
    hid_t group = H5Gcreate2(file, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (group < 0) return 1;
    return H5Gclose(group) < 0;
}

static int write_grid_staggered_coords(hid_t file, const char *var_name,
                                       int nx, int ny, int nz,
                                       const double *x, const double *y, const double *z)
{
    char group_name[64];
    char dataset_name[96];
    int ierr = 0;

    snprintf(group_name, sizeof(group_name), "staggered/%s", var_name);
    ierr |= create_group(file, group_name);

    snprintf(dataset_name, sizeof(dataset_name), "staggered/%s/x", var_name);
    ierr |= write_coord_dataset(file, dataset_name, nx + 2, 0, x);
    snprintf(dataset_name, sizeof(dataset_name), "staggered/%s/y", var_name);
    ierr |= write_coord_dataset(file, dataset_name, ny + 2, 0, y);
    snprintf(dataset_name, sizeof(dataset_name), "staggered/%s/z", var_name);
    ierr |= write_coord_dataset(file, dataset_name, nz + 2, 0, z);

    return ierr != 0;
}

int fdm_h5_write_grid(const char *filename, int nx, int ny, int nz,
                      double lx, double ly, double lz,
                      const int *periodic, const int *grid_distribution,
                      const double *grid_stretch, const double *grid_natural_dyw_plus,
                      const int *grid_natural_one_sided,
                      const double *x_node, const double *y_node, const double *z_node,
                      const double *xu, const double *yu, const double *zu,
                      const double *xv, const double *yv, const double *zv,
                      const double *xw, const double *yw, const double *zw)
{
    hid_t file;
    int ierr = 0;
    int staggered_counts[3] = {nx + 2, ny + 2, nz + 2};

    if (nx < 1 || ny < 1 || nz < 1) return 1;

    file = create_serial_file(filename);
    if (file < 0) return 1;

    ierr |= write_attr_int(file, "mobygrid_format", 1);
    ierr |= write_attr_int(file, "nx", nx);
    ierr |= write_attr_int(file, "ny", ny);
    ierr |= write_attr_int(file, "nz", nz);
    ierr |= write_attr_double(file, "lx", lx);
    ierr |= write_attr_double(file, "ly", ly);
    ierr |= write_attr_double(file, "lz", lz);
    ierr |= write_attr_int(file, "parallel_hdf5", 0);
    ierr |= write_attr_int_array(file, "periodic", periodic, 3);
    ierr |= write_attr_int_array(file, "grid_distribution", grid_distribution, 3);
    ierr |= write_attr_double_array(file, "grid_stretch", grid_stretch, 3);
    ierr |= write_attr_double_array(file, "grid_natural_dyw_plus", grid_natural_dyw_plus, 3);
    ierr |= write_attr_int_array(file, "grid_natural_one_sided", grid_natural_one_sided, 3);
    ierr |= write_attr_int_array(file, "staggered_counts", staggered_counts, 3);

    ierr |= write_coord_dataset(file, "x_nodes", nx + 1, 0, x_node);
    ierr |= write_coord_dataset(file, "y_nodes", ny + 1, 0, y_node);
    ierr |= write_coord_dataset(file, "z_nodes", nz + 1, 0, z_node);

    ierr |= create_group(file, "staggered");
    ierr |= write_grid_staggered_coords(file, "u", nx, ny, nz, xu, yu, zu);
    ierr |= write_grid_staggered_coords(file, "v", nx, ny, nz, xv, yv, zv);
    ierr |= write_grid_staggered_coords(file, "w", nx, ny, nz, xw, yw, zw);

    ierr |= H5Fclose(file) < 0;
    return ierr != 0;
}

static int read_global_dataset_blocks(hid_t file, const char *name,
                                      int nbx, int nby, int nbz,
                                      int n_blocks, const int *block_origin,
                                      int global_nx, int global_ny, int global_nz,
                                      double *q, size_t block_stride)
{
    const size_t ni = (size_t)nbx + 2;
    const size_t nj = (size_t)nby + 2;
    const size_t n = (size_t)nbx*(size_t)nby*(size_t)nbz;
    hsize_t expected_dims[3] = {(hsize_t)global_nz, (hsize_t)global_ny, (hsize_t)global_nx};
    hsize_t file_dims[3] = {0, 0, 0};
    hsize_t local_dims[3] = {(hsize_t)nbz, (hsize_t)nby, (hsize_t)nbx};
    hid_t dset = -1;
    hid_t file_space = -1;
    hid_t mem_space = -1;
    hid_t xfer = -1;
    double *buffer = NULL;
    herr_t status = 0;
    int ierr = 0;

    dset = H5Dopen2(file, name, H5P_DEFAULT);
    if (dset < 0) return 1;

    file_space = H5Dget_space(dset);
    if (file_space < 0 || H5Sget_simple_extent_ndims(file_space) != 3) {
        if (file_space >= 0) H5Sclose(file_space);
        H5Dclose(dset);
        return 1;
    }

    H5Sget_simple_extent_dims(file_space, file_dims, NULL);
    if (file_dims[0] != expected_dims[0] ||
        file_dims[1] != expected_dims[1] ||
        file_dims[2] != expected_dims[2]) {
        H5Sclose(file_space);
        H5Dclose(dset);
        return 1;
    }

    buffer = (double *)malloc(n*sizeof(double));
    mem_space = H5Screate_simple(3, local_dims, NULL);
    xfer = H5Pcreate(H5P_DATASET_XFER);
    if (buffer == NULL || mem_space < 0 || xfer < 0) {
        if (buffer != NULL) free(buffer);
        if (mem_space >= 0) H5Sclose(mem_space);
        if (xfer >= 0) H5Pclose(xfer);
        H5Sclose(file_space);
        H5Dclose(dset);
        return 1;
    }
    H5Pset_dxpl_mpio(xfer, H5FD_MPIO_INDEPENDENT);

    for (int b = 0; b < n_blocks; ++b) {
        double *field = q + (size_t)b*block_stride;
        hsize_t start[3] = {
            (hsize_t)block_origin[3*b + 2],
            (hsize_t)block_origin[3*b + 1],
            (hsize_t)block_origin[3*b + 0]
        };

        if (H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, local_dims, NULL) < 0) {
            ierr = 1;
            break;
        }
        status = H5Dread(dset, H5T_NATIVE_DOUBLE, mem_space, file_space, xfer, buffer);
        if (status < 0) {
            ierr = 1;
            break;
        }

        for (size_t k = 1; k <= (size_t)nbz; ++k) {
            for (size_t j = 1; j <= (size_t)nby; ++j) {
                for (size_t i = 1; i <= (size_t)nbx; ++i) {
                    field[linear_fortran(i, j, k, ni, nj)] =
                        buffer[linear_hdf5(k-1, j-1, i-1, (size_t)nby, (size_t)nbx)];
                }
            }
        }
    }

    H5Pclose(xfer);
    H5Sclose(mem_space);
    free(buffer);
    H5Sclose(file_space);
    H5Dclose(dset);
    return ierr;
}

int fdm_h5_write_field(const char *filename, int nbx, int nby, int nbz,
                       int n_blocks, int n_blocks_global, int id_start,
                       const int *block_origin, const int *block_level,
                       int rank, int nranks,
                       int global_nx, int global_ny, int global_nz,
                       int step, int nsteps,
                       double lx, double ly, double lz,
                       double re, double dt, double t_final, double t_current,
                       const double *cfl, double cflmax, double pecletmax, double dtmax,
                       const double *forcing,
                       int pressure_niter, double pressure_sor, int ibm_enabled, int bc_count,
                       const int *periodic, const int *bc_type, const double *bc_value,
                       const int *grid_distribution, const double *grid_stretch,
                       const double *grid_natural_dyw_plus,
                       const double *x_node, const double *y_node, const double *z_node,
                       const double *q)
{
    /* q is the solver's (0:nb+1,0:nb+1,0:nb+1, 4 vars, n_blocks) array. */
    const size_t var_stride = (size_t)(nbx + 2)*(size_t)(nby + 2)*(size_t)(nbz + 2);
    const size_t block_stride = var_stride*4;
    static const char *var_name[4] = {"un", "vn", "wn", "pn"};
    hid_t file;
    int ierr = 0;

    if (nbx < 1 || nby < 1 || nbz < 1 || n_blocks < 1) return 1;

    file = create_parallel_file(filename);
    if (file < 0) return 1;

    ierr |= write_attr_int(file, "nx", global_nx);
    ierr |= write_attr_int(file, "ny", global_ny);
    ierr |= write_attr_int(file, "nz", global_nz);
    ierr |= write_attr_int(file, "nranks", nranks);
    ierr |= write_attr_int(file, "parallel_hdf5", 1);
    ierr |= write_attr_int(file, "block_nb_x", nbx);
    ierr |= write_attr_int(file, "block_nb_y", nby);
    ierr |= write_attr_int(file, "block_nb_z", nbz);
    ierr |= write_attr_int(file, "n_blocks", n_blocks_global);
    ierr |= write_attr_int(file, "step", step);
    ierr |= write_attr_int(file, "nsteps", nsteps);
    ierr |= write_attr_double(file, "lx", lx);
    ierr |= write_attr_double(file, "ly", ly);
    ierr |= write_attr_double(file, "lz", lz);
    ierr |= write_attr_double(file, "re", re);
    ierr |= write_attr_double(file, "dt", dt);
    ierr |= write_attr_double(file, "t_final", t_final);
    ierr |= write_attr_double(file, "t_current", t_current);
    ierr |= write_attr_double_array(file, "cfl", cfl, 2);
    ierr |= write_attr_double(file, "cflmax", cflmax);
    ierr |= write_attr_double(file, "pecletmax", pecletmax);
    ierr |= write_attr_double(file, "dtmax", dtmax);
    ierr |= write_attr_double(file, "forcing_x", forcing[0]);
    ierr |= write_attr_double(file, "forcing_y", forcing[1]);
    ierr |= write_attr_double(file, "forcing_z", forcing[2]);
    ierr |= write_attr_int(file, "pressure_niter", pressure_niter);
    ierr |= write_attr_double(file, "pressure_sor", pressure_sor);
    ierr |= write_attr_int(file, "ibm_enabled", ibm_enabled);
    ierr |= write_attr_int_array(file, "periodic", periodic, 3);
    ierr |= write_attr_int_array(file, "bc_type", bc_type, (hsize_t)bc_count);
    ierr |= write_attr_double_array(file, "bc_value", bc_value, (hsize_t)bc_count);
    ierr |= write_attr_int_array(file, "grid_distribution", grid_distribution, 3);
    ierr |= write_attr_double_array(file, "grid_stretch", grid_stretch, 3);
    ierr |= write_attr_double_array(file, "grid_natural_dyw_plus", grid_natural_dyw_plus, 3);

    ierr |= write_block_table(file, n_blocks_global, id_start, n_blocks,
                              block_origin, block_level);

    for (int v = 0; v < 4; ++v) {
        ierr |= write_global_dataset_blocks(file, var_name[v], nbx, nby, nbz,
                                            n_blocks, block_origin,
                                            global_nx, global_ny, global_nz,
                                            q + (size_t)v*var_stride, block_stride);
    }
    ierr |= write_coord_dataset(file, "x", global_nx + 1, rank, x_node);
    ierr |= write_coord_dataset(file, "y", global_ny + 1, rank, y_node);
    ierr |= write_coord_dataset(file, "z", global_nz + 1, rank, z_node);

    ierr |= H5Fclose(file) < 0;
    return ierr != 0;
}

int fdm_h5_read_metadata(const char *filename,
                         int *global_nx, int *global_ny, int *global_nz,
                         int *step, int *nsteps,
                         double *lx, double *ly, double *lz,
                         double *re, double *dt, double *t_final, double *t_current,
                         double *cfl, double *cflmax, double *pecletmax, double *dtmax,
                         double *forcing,
                         int *pressure_niter, double *pressure_sor, int *ibm_enabled, int bc_count,
                         int *periodic, int *bc_type, double *bc_value,
                         int *grid_distribution, double *grid_stretch, double *grid_natural_dyw_plus)
{
    hid_t file;
    int file_nranks = 0;
    int file_parallel_hdf5 = 0;
    int ierr = 0;

    file = open_parallel_file(filename);
    if (file < 0) return 1;

    ierr |= read_attr_int(file, "nx", global_nx, 1);
    ierr |= read_attr_int(file, "ny", global_ny, 1);
    ierr |= read_attr_int(file, "nz", global_nz, 1);
    ierr |= read_attr_int(file, "nranks", &file_nranks, 1);
    ierr |= read_attr_int(file, "parallel_hdf5", &file_parallel_hdf5, 1);
    ierr |= file_nranks < 1;
    ierr |= file_parallel_hdf5 < 0 || file_parallel_hdf5 > 1;
    ierr |= read_attr_int(file, "step", step, 1);
    ierr |= read_attr_int(file, "nsteps", nsteps, 1);

    ierr |= read_attr_double(file, "lx", lx, 1);
    ierr |= read_attr_double(file, "ly", ly, 1);
    ierr |= read_attr_double(file, "lz", lz, 1);
    ierr |= read_attr_double(file, "re", re, 1);
    ierr |= read_attr_double(file, "dt", dt, 1);
    ierr |= read_attr_double(file, "t_final", t_final, 1);
    ierr |= read_attr_double(file, "t_current", t_current, 1);
    ierr |= read_attr_double_array(file, "cfl", cfl, 2, 1);
    ierr |= read_attr_double(file, "cflmax", cflmax, 1);
    ierr |= read_attr_double(file, "pecletmax", pecletmax, 1);
    ierr |= read_attr_double(file, "dtmax", dtmax, 1);
    ierr |= read_attr_double(file, "forcing_x", &forcing[0], 1);
    ierr |= read_attr_double(file, "forcing_y", &forcing[1], 1);
    ierr |= read_attr_double(file, "forcing_z", &forcing[2], 1);
    ierr |= read_attr_int(file, "pressure_niter", pressure_niter, 1);
    ierr |= read_attr_double(file, "pressure_sor", pressure_sor, 1);
    ierr |= read_attr_int(file, "ibm_enabled", ibm_enabled, 1);
    ierr |= read_attr_int_array(file, "periodic", periodic, 3, 1);
    ierr |= read_attr_int_array(file, "bc_type", bc_type, (hsize_t)bc_count, 1);
    ierr |= read_attr_double_array(file, "bc_value", bc_value, (hsize_t)bc_count, 1);
    ierr |= read_attr_int_array(file, "grid_distribution", grid_distribution, 3, 1);
    ierr |= read_attr_double_array(file, "grid_stretch", grid_stretch, 3, 1);
    for (int d = 0; d < 3; ++d) grid_natural_dyw_plus[d] = 0.05;
    ierr |= read_attr_double_array(file, "grid_natural_dyw_plus", grid_natural_dyw_plus, 3, 0);

    ierr |= H5Fclose(file) < 0;
    return ierr != 0;
}

int fdm_h5_read_field(const char *filename, int nbx, int nby, int nbz,
                      int n_blocks, const int *block_origin,
                      int global_nx, int global_ny, int global_nz,
                      double *q)
{
    const size_t var_stride = (size_t)(nbx + 2)*(size_t)(nby + 2)*(size_t)(nbz + 2);
    const size_t block_stride = var_stride*4;
    static const char *var_name[4] = {"un", "vn", "wn", "pn"};
    hid_t file;
    int ierr = 0;

    if (nbx < 1 || nby < 1 || nbz < 1 || n_blocks < 1) return 1;

    file = open_parallel_file(filename);
    if (file < 0) return 1;

    for (int v = 0; v < 4; ++v) {
        ierr |= read_global_dataset_blocks(file, var_name[v], nbx, nby, nbz,
                                           n_blocks, block_origin,
                                           global_nx, global_ny, global_nz,
                                           q + (size_t)v*var_stride, block_stride);
    }

    ierr |= H5Fclose(file) < 0;
    return ierr != 0;
}

/*
 * Optional per-block keep flags written by mobygeom block-active
 * (x-fastest lattice raster order). *found = 0 when the dataset is absent,
 * which is not an error: the solver then keeps every block.
 */
int fdm_h5_read_block_active(const char *filename, int n_lattice, int block_nb,
                             int *found, int *active)
{
    hsize_t dims[1] = {0};
    hid_t file = -1;
    hid_t dset = -1;
    hid_t space = -1;
    htri_t exists;
    int file_nb = 0;
    herr_t status;

    *found = 0;

    file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) return 1;

    exists = H5Lexists(file, "block_active", H5P_DEFAULT);
    if (exists <= 0) {
        H5Fclose(file);
        return 0;
    }

    if (read_attr_int(file, "block_nb", &file_nb, 1) != 0 || file_nb != block_nb) {
        H5Fclose(file);
        return 1;
    }

    dset = H5Dopen2(file, "block_active", H5P_DEFAULT);
    if (dset < 0) {
        H5Fclose(file);
        return 1;
    }
    space = H5Dget_space(dset);
    if (space < 0 || H5Sget_simple_extent_ndims(space) != 1) {
        if (space >= 0) H5Sclose(space);
        H5Dclose(dset);
        H5Fclose(file);
        return 1;
    }
    H5Sget_simple_extent_dims(space, dims, NULL);
    if (dims[0] != (hsize_t)n_lattice) {
        H5Sclose(space);
        H5Dclose(dset);
        H5Fclose(file);
        return 1;
    }

    status = H5Dread(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, active);
    H5Sclose(space);
    H5Dclose(dset);
    H5Fclose(file);
    if (status < 0) return 1;

    *found = 1;
    return 0;
}

int fdm_h5_read_ibm_coeff(const char *filename, int nbx, int nby, int nbz,
                          int n_blocks, const int *block_origin,
                          int global_nx, int global_ny, int global_nz,
                          double lx, double ly, double lz, double re,
                          double *coef)
{
    /* The coefficient file carries the global ghost layer, so a block's
     * window (halos included) starts at its zero-based cell origin. */
    const size_t ni = (size_t)nbx + 2;
    const size_t nj = (size_t)nby + 2;
    const size_t nk = (size_t)nbz + 2;
    const size_t n = ni*nj*nk*3;
    const size_t block_stride = n;
    hsize_t expected_dims[4] = {
        (hsize_t)global_nx + 2,
        (hsize_t)global_ny + 2,
        (hsize_t)global_nz + 2,
        3
    };
    hsize_t file_dims[4] = {0, 0, 0, 0};
    hsize_t local_dims[4] = {
        (hsize_t)nbx + 2,
        (hsize_t)nby + 2,
        (hsize_t)nbz + 2,
        3
    };
    hid_t file = -1;
    hid_t dset = -1;
    hid_t file_space = -1;
    hid_t mem_space = -1;
    double *buffer = NULL;
    int file_nx = 0, file_ny = 0, file_nz = 0;
    double file_lx = 0.0, file_ly = 0.0, file_lz = 0.0, file_re = 0.0;
    int ierr = 0;
    herr_t status = 0;

    if (nbx < 1 || nby < 1 || nbz < 1 || n_blocks < 1) return 1;

    file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) return 1;

    ierr |= read_attr_int(file, "nx", &file_nx, 1);
    ierr |= read_attr_int(file, "ny", &file_ny, 1);
    ierr |= read_attr_int(file, "nz", &file_nz, 1);
    ierr |= read_attr_double(file, "lx", &file_lx, 1);
    ierr |= read_attr_double(file, "ly", &file_ly, 1);
    ierr |= read_attr_double(file, "lz", &file_lz, 1);
    ierr |= read_attr_double(file, "re", &file_re, 1);
    ierr |= file_nx != global_nx || file_ny != global_ny || file_nz != global_nz;
    ierr |= double_mismatch(file_lx, lx) || double_mismatch(file_ly, ly) ||
            double_mismatch(file_lz, lz) || double_mismatch(file_re, re);

    dset = H5Dopen2(file, "coef", H5P_DEFAULT);
    if (dset < 0) {
        H5Fclose(file);
        return 1;
    }

    file_space = H5Dget_space(dset);
    if (file_space < 0 || H5Sget_simple_extent_ndims(file_space) != 4) {
        if (file_space >= 0) H5Sclose(file_space);
        H5Dclose(dset);
        H5Fclose(file);
        return 1;
    }

    H5Sget_simple_extent_dims(file_space, file_dims, NULL);
    for (int d = 0; d < 4; ++d) {
        ierr |= file_dims[d] != expected_dims[d];
    }

    buffer = (double *)malloc(n*sizeof(double));
    if (buffer == NULL) {
        H5Sclose(file_space);
        H5Dclose(dset);
        H5Fclose(file);
        return 1;
    }

    mem_space = H5Screate_simple(4, local_dims, NULL);
    if (mem_space < 0) {
        free(buffer);
        H5Sclose(file_space);
        H5Dclose(dset);
        H5Fclose(file);
        return 1;
    }

    for (int b = 0; b < n_blocks && ierr == 0; ++b) {
        double *block_coef = coef + (size_t)b*block_stride;
        hsize_t start[4] = {
            (hsize_t)block_origin[3*b + 0],
            (hsize_t)block_origin[3*b + 1],
            (hsize_t)block_origin[3*b + 2],
            0
        };

        if (H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, local_dims, NULL) < 0) {
            ierr = 1;
            break;
        }
        status = H5Dread(dset, H5T_NATIVE_DOUBLE, mem_space, file_space, H5P_DEFAULT, buffer);
        if (status < 0) {
            ierr = 1;
            break;
        }
        for (size_t v = 0; v < 3; ++v) {
            for (size_t k = 0; k < nk; ++k) {
                for (size_t j = 0; j < nj; ++j) {
                    for (size_t i = 0; i < ni; ++i) {
                        size_t h5_idx = (((i*nj) + j)*nk + k)*3 + v;
                        block_coef[linear_fortran4(i, j, k, v, ni, nj, nk)] = buffer[h5_idx];
                    }
                }
            }
        }
    }

    H5Sclose(mem_space);
    free(buffer);
    H5Sclose(file_space);
    H5Dclose(dset);
    ierr |= H5Fclose(file) < 0;
    return ierr != 0 || status < 0;
}

static int write_dataset1(hid_t file, const char *name, hsize_t n, const double *values)
{
    hid_t space = -1;
    hid_t dset = -1;
    herr_t status;

    space = H5Screate_simple(1, &n, NULL);
    if (space < 0) return 1;

    dset = H5Dcreate2(file, name, H5T_NATIVE_DOUBLE, space,
                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dset < 0) {
        H5Sclose(space);
        return 1;
    }

    status = H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, values);
    H5Dclose(dset);
    H5Sclose(space);
    return status < 0;
}

static int write_dataset2_fortran(hid_t file, const char *name,
                                  int nwall, int nstat, const double *values)
{
    hsize_t dims[2] = {(hsize_t)nwall, (hsize_t)nstat};
    size_t n = (size_t)nwall*(size_t)nstat;
    hid_t space = -1;
    hid_t dset = -1;
    double *buffer = NULL;
    herr_t status;

    buffer = (double *)malloc(n*sizeof(double));
    if (buffer == NULL) return 1;

    for (int i = 0; i < nwall; ++i) {
        for (int s = 0; s < nstat; ++s) {
            buffer[(size_t)i*(size_t)nstat + (size_t)s] =
                values[(size_t)s + (size_t)nstat*(size_t)i];
        }
    }

    space = H5Screate_simple(2, dims, NULL);
    if (space < 0) {
        free(buffer);
        return 1;
    }

    dset = H5Dcreate2(file, name, H5T_NATIVE_DOUBLE, space,
                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dset < 0) {
        H5Sclose(space);
        free(buffer);
        return 1;
    }

    status = H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);

    H5Dclose(dset);
    H5Sclose(space);
    free(buffer);
    return status < 0;
}

static int read_dataset1(hid_t file, const char *name, hsize_t n, double *values)
{
    hsize_t dims[1] = {0};
    hid_t dset = -1;
    hid_t space = -1;
    herr_t status;

    dset = H5Dopen2(file, name, H5P_DEFAULT);
    if (dset < 0) return 1;

    space = H5Dget_space(dset);
    if (space < 0 || H5Sget_simple_extent_ndims(space) != 1) {
        if (space >= 0) H5Sclose(space);
        H5Dclose(dset);
        return 1;
    }

    H5Sget_simple_extent_dims(space, dims, NULL);
    if (dims[0] != n) {
        H5Sclose(space);
        H5Dclose(dset);
        return 1;
    }

    status = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, values);
    H5Sclose(space);
    H5Dclose(dset);
    return status < 0;
}

static int read_dataset2_fortran(hid_t file, const char *name,
                                 int nwall, int nstat, double *values)
{
    hsize_t expected[2] = {(hsize_t)nwall, (hsize_t)nstat};
    hsize_t dims[2] = {0, 0};
    size_t n = (size_t)nwall*(size_t)nstat;
    hid_t dset = -1;
    hid_t space = -1;
    double *buffer = NULL;
    herr_t status;

    dset = H5Dopen2(file, name, H5P_DEFAULT);
    if (dset < 0) return 1;

    space = H5Dget_space(dset);
    if (space < 0 || H5Sget_simple_extent_ndims(space) != 2) {
        if (space >= 0) H5Sclose(space);
        H5Dclose(dset);
        return 1;
    }

    H5Sget_simple_extent_dims(space, dims, NULL);
    if (dims[0] != expected[0] || dims[1] != expected[1]) {
        H5Sclose(space);
        H5Dclose(dset);
        return 1;
    }

    buffer = (double *)malloc(n*sizeof(double));
    if (buffer == NULL) {
        H5Sclose(space);
        H5Dclose(dset);
        return 1;
    }

    status = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);
    if (status >= 0) {
        for (int i = 0; i < nwall; ++i) {
            for (int s = 0; s < nstat; ++s) {
                values[(size_t)s + (size_t)nstat*(size_t)i] =
                    buffer[(size_t)i*(size_t)nstat + (size_t)s];
            }
        }
    }

    free(buffer);
    H5Sclose(space);
    H5Dclose(dset);
    return status < 0;
}

int fdm_h5_write_channel_stats(const char *filename, int nwall, int nstat,
                               int step, double t_current, int wall_dir, double re,
                               const double *forcing, const double *coord,
                               const double *profile, const double *raw_sum,
                               const double *count)
{
    hid_t file;
    int ierr = 0;

    if (nwall < 1 || nstat < 1) return 1;

    file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file < 0) return 1;

    ierr |= write_attr_int(file, "nwall", nwall);
    ierr |= write_attr_int(file, "nstat", nstat);
    ierr |= write_attr_int(file, "sample_weighting", 1);
    ierr |= write_attr_int(file, "step", step);
    ierr |= write_attr_int(file, "wall_dir", wall_dir);
    ierr |= write_attr_double(file, "t_current", t_current);
    ierr |= write_attr_double(file, "re", re);
    ierr |= write_attr_double(file, "forcing_x", forcing[0]);
    ierr |= write_attr_double(file, "forcing_y", forcing[1]);
    ierr |= write_attr_double(file, "forcing_z", forcing[2]);
    ierr |= write_dataset1(file, "coord", (hsize_t)nwall, coord);
    ierr |= write_dataset1(file, "count", (hsize_t)nwall, count);
    ierr |= write_dataset2_fortran(file, "profile", nwall, nstat, profile);
    ierr |= write_dataset2_fortran(file, "raw_sum", nwall, nstat, raw_sum);

    ierr |= H5Fclose(file) < 0;
    return ierr != 0;
}

int fdm_h5_read_channel_stats(const char *filename, int nwall, int nstat,
                              int *step, double *t_current,
                              double *raw_sum, double *count)
{
    hid_t file;
    int file_nwall = nwall;
    int file_nstat = nstat;
    int sample_weighting = 0;
    int ierr = 0;

    file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) return 1;

    ierr |= read_attr_int(file, "nwall", &file_nwall, 1);
    ierr |= read_attr_int(file, "nstat", &file_nstat, 1);
    ierr |= read_attr_int(file, "sample_weighting", &sample_weighting, 1);
    ierr |= read_attr_int(file, "step", step, 1);
    ierr |= read_attr_double(file, "t_current", t_current, 1);

    if (file_nwall != nwall || file_nstat != nstat || sample_weighting != 1) {
        H5Fclose(file);
        return 1;
    }

    ierr |= read_dataset1(file, "count", (hsize_t)nwall, count);
    ierr |= read_dataset2_fortran(file, "raw_sum", nwall, nstat, raw_sum);

    ierr |= H5Fclose(file) < 0;
    return ierr != 0;
}
