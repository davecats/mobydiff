#include <mpi.h>
#include <hdf5.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#ifndef H5_HAVE_PARALLEL
#error "field_hdf5.c requires HDF5 built with parallel MPI-IO support"
#endif

/* Field-variable name table (fdm_h5_write_field / fdm_h5_read_field): the
 * Fortran side packs n_var NUL-terminated names of this fixed stride, so the
 * variable count is a runtime value (u,v,w,p plus the passive scalars,
 * docs/next_session_scalar.md) instead of the hardcoded 4. */
#define FDM_VAR_NAME_LEN 32

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

/* Variable-length-free fixed string attribute (the case-input echo). */
static int write_attr_string(hid_t file, const char *name, const char *text)
{
    hid_t space = H5Screate(H5S_SCALAR);
    hid_t type = H5Tcopy(H5T_C_S1);
    hid_t attr = -1;
    size_t n = strlen(text);
    herr_t status = -1;

    if (space < 0 || type < 0 || H5Tset_size(type, n + 1) < 0) {
        if (space >= 0) H5Sclose(space);
        if (type >= 0) H5Tclose(type);
        return 1;
    }
    attr = H5Acreate2(file, name, type, space, H5P_DEFAULT, H5P_DEFAULT);
    if (attr >= 0) status = H5Awrite(attr, type, text);
    if (attr >= 0) H5Aclose(attr);
    H5Tclose(type);
    H5Sclose(space);
    return status < 0;
}

/* *found = 0 when the attribute is absent (older case files). The text is
 * NUL-terminated into buf; longer texts are a read error. */
static int read_attr_string(hid_t file, const char *name, char *buf, size_t cap, int *found)
{
    htri_t exists = H5Aexists(file, name);
    hid_t attr = -1, type = -1;
    size_t n;
    herr_t status = -1;

    *found = 0;
    if (exists <= 0) return 0;
    attr = H5Aopen(file, name, H5P_DEFAULT);
    type = attr >= 0 ? H5Aget_type(attr) : -1;
    if (attr < 0 || type < 0) {
        if (attr >= 0) H5Aclose(attr);
        return 1;
    }
    n = H5Tget_size(type);
    if (n + 1 <= cap) {
        status = H5Aread(attr, type, buf);
        buf[n] = '\0';
    }
    H5Tclose(type);
    H5Aclose(attr);
    if (status < 0) return 1;
    *found = 1;
    return 0;
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
 * Block-table dataset layout (strategy doc Section 8): one row of
 * (nbz, nby, nbx) interior cells per global block id. Each rank writes its
 * contiguous id range with one independent transfer under the collectively
 * created dataset.
 */
static int write_block_dataset(hid_t file, const char *name,
                               int nbx, int nby, int nbz,
                               int n_blocks, int n_blocks_global, int id_start,
                               const double *q, size_t block_stride)
{
    const size_t ni = (size_t)nbx + 2;
    const size_t nj = (size_t)nby + 2;
    const size_t n = (size_t)nbx*(size_t)nby*(size_t)nbz;
    hsize_t global_dims[4] = {(hsize_t)n_blocks_global, (hsize_t)nbz, (hsize_t)nby, (hsize_t)nbx};
    hsize_t local_dims[4] = {(hsize_t)n_blocks, (hsize_t)nbz, (hsize_t)nby, (hsize_t)nbx};
    hsize_t start[4] = {(hsize_t)id_start, 0, 0, 0};
    hid_t file_space = -1;
    hid_t mem_space = -1;
    hid_t dset = -1;
    hid_t xfer = -1;
    double *buffer = NULL;
    herr_t status = 0;
    int ierr = 0;

    buffer = (double *)malloc((size_t)n_blocks*n*sizeof(double));
    if (buffer == NULL) return 1;

    for (int b = 0; b < n_blocks; ++b) {
        const double *field = q + (size_t)b*block_stride;
        double *row = buffer + (size_t)b*n;
        for (size_t k = 1; k <= (size_t)nbz; ++k) {
            for (size_t j = 1; j <= (size_t)nby; ++j) {
                for (size_t i = 1; i <= (size_t)nbx; ++i) {
                    row[linear_hdf5(k-1, j-1, i-1, (size_t)nby, (size_t)nbx)] =
                        field[linear_fortran(i, j, k, ni, nj)];
                }
            }
        }
    }

    file_space = H5Screate_simple(4, global_dims, NULL);
    if (file_space < 0) {
        free(buffer);
        return 1;
    }
    dset = H5Dcreate2(file, name, H5T_NATIVE_DOUBLE, file_space,
                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    mem_space = H5Screate_simple(4, local_dims, NULL);
    xfer = H5Pcreate(H5P_DATASET_XFER);
    if (dset < 0 || mem_space < 0 || xfer < 0 ||
        H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, local_dims, NULL) < 0) {
        if (dset >= 0) H5Dclose(dset);
        if (mem_space >= 0) H5Sclose(mem_space);
        if (xfer >= 0) H5Pclose(xfer);
        H5Sclose(file_space);
        free(buffer);
        return 1;
    }
    H5Pset_dxpl_mpio(xfer, H5FD_MPIO_INDEPENDENT);

    status = H5Dwrite(dset, H5T_NATIVE_DOUBLE, mem_space, file_space, xfer, buffer);
    if (status < 0) ierr = 1;

    H5Pclose(xfer);
    H5Sclose(mem_space);
    H5Dclose(dset);
    H5Sclose(file_space);
    free(buffer);
    return ierr;
}

/*
 * Integer variant of write_block_dataset for per-cell markers held
 * interior-only ((nbx,nby,nbz) per block, Fortran order, no ghost layer).
 */
static int write_block_int_dataset(hid_t file, const char *name,
                                   int nbx, int nby, int nbz,
                                   int n_blocks, int n_blocks_global, int id_start,
                                   const int *vals)
{
    const size_t n = (size_t)nbx*(size_t)nby*(size_t)nbz;
    hsize_t global_dims[4] = {(hsize_t)n_blocks_global, (hsize_t)nbz, (hsize_t)nby, (hsize_t)nbx};
    hsize_t local_dims[4] = {(hsize_t)n_blocks, (hsize_t)nbz, (hsize_t)nby, (hsize_t)nbx};
    hsize_t start[4] = {(hsize_t)id_start, 0, 0, 0};
    hid_t file_space = -1;
    hid_t mem_space = -1;
    hid_t dset = -1;
    hid_t xfer = -1;
    int *buffer = NULL;
    herr_t status = 0;
    int ierr = 0;

    buffer = (int *)malloc((size_t)n_blocks*n*sizeof(int));
    if (buffer == NULL) return 1;

    for (int b = 0; b < n_blocks; ++b) {
        const int *field = vals + (size_t)b*n;
        int *row = buffer + (size_t)b*n;
        for (size_t k = 0; k < (size_t)nbz; ++k) {
            for (size_t j = 0; j < (size_t)nby; ++j) {
                for (size_t i = 0; i < (size_t)nbx; ++i) {
                    row[linear_hdf5(k, j, i, (size_t)nby, (size_t)nbx)] =
                        field[linear_fortran(i, j, k, (size_t)nbx, (size_t)nby)];
                }
            }
        }
    }

    file_space = H5Screate_simple(4, global_dims, NULL);
    dset = file_space >= 0 ? H5Dcreate2(file, name, H5T_NATIVE_INT, file_space,
                                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) : -1;
    mem_space = H5Screate_simple(4, local_dims, NULL);
    xfer = H5Pcreate(H5P_DATASET_XFER);
    if (file_space < 0 || dset < 0 || mem_space < 0 || xfer < 0 ||
        H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, local_dims, NULL) < 0) {
        ierr = 1;
    } else {
        H5Pset_dxpl_mpio(xfer, H5FD_MPIO_INDEPENDENT);
        status = H5Dwrite(dset, H5T_NATIVE_INT, mem_space, file_space, xfer, buffer);
        if (status < 0) ierr = 1;
    }

    if (xfer >= 0) H5Pclose(xfer);
    if (mem_space >= 0) H5Sclose(mem_space);
    if (dset >= 0) H5Dclose(dset);
    if (file_space >= 0) H5Sclose(file_space);
    free(buffer);
    return ierr;
}

/*
 * One row of n doubles per global block id (e.g. per-block 1D coordinate
 * lines); each rank writes its contiguous id range independently.
 */
static int write_block_rows(hid_t file, const char *name, int n,
                            int n_blocks, int n_blocks_global, int id_start,
                            const double *rows)
{
    hsize_t global_dims[2] = {(hsize_t)n_blocks_global, (hsize_t)n};
    hsize_t local_dims[2] = {(hsize_t)n_blocks, (hsize_t)n};
    hsize_t start[2] = {(hsize_t)id_start, 0};
    hid_t file_space = -1;
    hid_t mem_space = -1;
    hid_t dset = -1;
    hid_t xfer = -1;
    herr_t status = 0;
    int ierr = 0;

    file_space = H5Screate_simple(2, global_dims, NULL);
    dset = file_space >= 0 ? H5Dcreate2(file, name, H5T_NATIVE_DOUBLE, file_space,
                                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) : -1;
    mem_space = H5Screate_simple(2, local_dims, NULL);
    xfer = H5Pcreate(H5P_DATASET_XFER);
    if (file_space < 0 || dset < 0 || mem_space < 0 || xfer < 0 ||
        H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, local_dims, NULL) < 0) {
        ierr = 1;
    } else {
        H5Pset_dxpl_mpio(xfer, H5FD_MPIO_INDEPENDENT);
        status = H5Dwrite(dset, H5T_NATIVE_DOUBLE, mem_space, file_space, xfer, rows);
        if (status < 0) ierr = 1;
    }

    if (xfer >= 0) H5Pclose(xfer);
    if (mem_space >= 0) H5Sclose(mem_space);
    if (dset >= 0) H5Dclose(dset);
    if (file_space >= 0) H5Sclose(file_space);
    return ierr;
}

/*
 * The file row that holds each of this rank's blocks, matched on
 * (origin, level) against the file's own blocks table.
 *
 * A block-layout dataset is sliced by ROW, and the row order is the writer's
 * leaf order -- which need not be the reader's: a snapshot written under one
 * key-bit order (blocks.f90 leaf_key) restarts under another, and so does an
 * initial condition generated by a tool. Matching on (origin, level) makes
 * the row order of a field file irrelevant, and it is also the layout guard:
 * a block the file does not hold means a different tiling (return 2).
 *
 * *row_of is malloc'ed (n_blocks entries, the caller frees); *reordered says
 * whether the map differs from the identity id_start + b. A file with no
 * blocks table (a legacy global-3D restart) has no rows to map: identity.
 * 0 = ok, 1 = read error, 2 = layout mismatch.
 */
typedef struct { int key[4]; int row; } block_row_t;   /* level, z, y, x */

static int block_row_cmp(const void *pa, const void *pb)
{
    const block_row_t *a = (const block_row_t *)pa;
    const block_row_t *b = (const block_row_t *)pb;
    for (int i = 0; i < 4; ++i) {
        if (a->key[i] != b->key[i]) return a->key[i] < b->key[i] ? -1 : 1;
    }
    return 0;
}

static int block_row_map(hid_t file, int n_blocks, int n_blocks_global, int id_start,
                         const int *block_origin, const int *block_level,
                         int **row_of, int *reordered)
{
    hid_t bset = -1, bspace = -1;
    hsize_t dims[2] = {0, 0};
    int *rows = NULL;
    block_row_t *sorted = NULL;
    int *map = (int *)malloc((size_t)n_blocks*sizeof(int));
    int same = 1;
    int ierr = 0;

    *row_of = NULL;
    *reordered = 0;
    if (map == NULL) return 1;
    for (int b = 0; b < n_blocks; ++b) map[b] = id_start + b;

    if (H5Lexists(file, "blocks", H5P_DEFAULT) <= 0) {
        *row_of = map;
        return 0;
    }

    bset = H5Dopen2(file, "blocks", H5P_DEFAULT);
    bspace = bset >= 0 ? H5Dget_space(bset) : -1;
    if (bset < 0 || bspace < 0 || H5Sget_simple_extent_ndims(bspace) != 2) {
        ierr = 1;
    } else {
        H5Sget_simple_extent_dims(bspace, dims, NULL);
        if (dims[1] != 4) ierr = 1;
        else if (dims[0] != (hsize_t)n_blocks_global) ierr = 2;
    }
    if (!ierr) {
        rows = (int *)malloc((size_t)n_blocks_global*4*sizeof(int));
        if (rows == NULL ||
            H5Dread(bset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, rows) < 0) ierr = 1;
    }
    if (bspace >= 0) H5Sclose(bspace);
    if (bset >= 0) H5Dclose(bset);

    /* The common case: the file is in this run's order. */
    for (int b = 0; b < n_blocks && !ierr && same; ++b) {
        const int *r = rows + 4*(size_t)(id_start + b);
        same = r[0] == block_origin[3*b+0] && r[1] == block_origin[3*b+1] &&
               r[2] == block_origin[3*b+2] && r[3] == block_level[b];
    }

    if (!ierr && !same) {
        sorted = (block_row_t *)malloc((size_t)n_blocks_global*sizeof(block_row_t));
        if (sorted == NULL) {
            ierr = 1;
        } else {
            for (int i = 0; i < n_blocks_global; ++i) {
                sorted[i].key[0] = rows[4*(size_t)i+3];
                sorted[i].key[1] = rows[4*(size_t)i+2];
                sorted[i].key[2] = rows[4*(size_t)i+1];
                sorted[i].key[3] = rows[4*(size_t)i+0];
                sorted[i].row = i;
            }
            qsort(sorted, (size_t)n_blocks_global, sizeof(block_row_t), block_row_cmp);
            for (int b = 0; b < n_blocks; ++b) {
                block_row_t want;
                const block_row_t *hit;
                want.key[0] = block_level[b];
                want.key[1] = block_origin[3*b+2];
                want.key[2] = block_origin[3*b+1];
                want.key[3] = block_origin[3*b+0];
                want.row = 0;
                hit = (const block_row_t *)bsearch(&want, sorted, (size_t)n_blocks_global,
                                                   sizeof(block_row_t), block_row_cmp);
                if (hit == NULL) {
                    ierr = 2;
                    break;
                }
                map[b] = hit->row;
            }
            *reordered = 1;
        }
    }

    if (sorted != NULL) free(sorted);
    if (rows != NULL) free(rows);
    if (ierr) {
        free(map);
        return ierr;
    }
    *row_of = map;
    return 0;
}

/* One (n_blocks_global, nbz, nby, nbx) dataset into this rank's blocks;
 * row_of[b] is the file row of local block b (block_row_map). Consecutive
 * rows are read as one hyperslab, so a file in this run's own order is a
 * single read, exactly as before the map existed. */
static int read_block_dataset(hid_t file, hid_t dset, hid_t file_space,
                              int nbx, int nby, int nbz,
                              int n_blocks, int n_blocks_global, const int *row_of,
                              double *q, size_t block_stride)
{
    const size_t ni = (size_t)nbx + 2;
    const size_t nj = (size_t)nby + 2;
    const size_t n = (size_t)nbx*(size_t)nby*(size_t)nbz;
    hsize_t expected[4] = {(hsize_t)n_blocks_global, (hsize_t)nbz, (hsize_t)nby, (hsize_t)nbx};
    hsize_t file_dims[4] = {0, 0, 0, 0};
    hid_t xfer = -1;
    double *buffer = NULL;
    int ierr = 0;

    (void)file;
    H5Sget_simple_extent_dims(file_space, file_dims, NULL);
    for (int d = 0; d < 4; ++d) {
        if (file_dims[d] != expected[d]) return 1;
    }

    buffer = (double *)malloc((size_t)n_blocks*n*sizeof(double));
    xfer = H5Pcreate(H5P_DATASET_XFER);
    if (buffer == NULL || xfer < 0) {
        if (buffer != NULL) free(buffer);
        if (xfer >= 0) H5Pclose(xfer);
        return 1;
    }
    H5Pset_dxpl_mpio(xfer, H5FD_MPIO_INDEPENDENT);

    for (int b0 = 0; b0 < n_blocks && !ierr; ) {
        int b1 = b0;
        hsize_t start[4] = {(hsize_t)row_of[b0], 0, 0, 0};
        hsize_t count[4] = {0, (hsize_t)nbz, (hsize_t)nby, (hsize_t)nbx};
        hid_t mem_space;

        while (b1 + 1 < n_blocks && row_of[b1+1] == row_of[b1] + 1) ++b1;
        count[0] = (hsize_t)(b1 - b0 + 1);
        mem_space = H5Screate_simple(4, count, NULL);
        if (mem_space < 0 ||
            H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, count, NULL) < 0 ||
            H5Dread(dset, H5T_NATIVE_DOUBLE, mem_space, file_space, xfer,
                    buffer + (size_t)b0*n) < 0) ierr = 1;
        if (mem_space >= 0) H5Sclose(mem_space);
        b0 = b1 + 1;
    }

    if (!ierr) {
        for (int b = 0; b < n_blocks; ++b) {
            double *field = q + (size_t)b*block_stride;
            const double *row = buffer + (size_t)b*n;
            for (size_t k = 1; k <= (size_t)nbz; ++k) {
                for (size_t j = 1; j <= (size_t)nby; ++j) {
                    for (size_t i = 1; i <= (size_t)nbx; ++i) {
                        field[linear_fortran(i, j, k, ni, nj)] =
                            row[linear_hdf5(k-1, j-1, i-1, (size_t)nby, (size_t)nbx)];
                    }
                }
            }
        }
    }

    H5Pclose(xfer);
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

static int read_global_dataset_blocks(hid_t file, hid_t dset, hid_t file_space,
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
    hid_t mem_space = -1;
    hid_t xfer = -1;
    double *buffer = NULL;
    herr_t status = 0;
    int ierr = 0;

    if (H5Sget_simple_extent_ndims(file_space) != 3) return 1;

    H5Sget_simple_extent_dims(file_space, file_dims, NULL);
    if (file_dims[0] != expected_dims[0] ||
        file_dims[1] != expected_dims[1] ||
        file_dims[2] != expected_dims[2]) {
        return 1;
    }

    buffer = (double *)malloc(n*sizeof(double));
    mem_space = H5Screate_simple(3, local_dims, NULL);
    xfer = H5Pcreate(H5P_DATASET_XFER);
    if (buffer == NULL || mem_space < 0 || xfer < 0) {
        if (buffer != NULL) free(buffer);
        if (mem_space >= 0) H5Sclose(mem_space);
        if (xfer >= 0) H5Pclose(xfer);
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
                       int n_var, const char *var_names,
                       const double *q)
{
    /* q is the solver's (0:nb+1,0:nb+1,0:nb+1, n_var, n_blocks) array;
     * var_names holds n_var NUL-terminated names of stride FDM_VAR_NAME_LEN
     * ("un","vn","wn","pn" plus one per passive scalar). */
    const size_t var_stride = (size_t)(nbx + 2)*(size_t)(nby + 2)*(size_t)(nbz + 2);
    const size_t block_stride = var_stride*(size_t)n_var;
    hid_t file;
    int ierr = 0;

    if (nbx < 1 || nby < 1 || nbz < 1 || n_blocks < 1 || n_var < 1) return 1;

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

    for (int v = 0; v < n_var; ++v) {
        ierr |= write_block_dataset(file, var_names + (size_t)v*FDM_VAR_NAME_LEN,
                                    nbx, nby, nbz,
                                    n_blocks, n_blocks_global, id_start,
                                    q + (size_t)v*var_stride, block_stride);
    }
    ierr |= write_coord_dataset(file, "x", global_nx + 1, rank, x_node);
    ierr |= write_coord_dataset(file, "y", global_ny + 1, rank, y_node);
    ierr |= write_coord_dataset(file, "z", global_nz + 1, rank, z_node);

    ierr |= H5Fclose(file) < 0;
    return ierr != 0;
}

/* Append the LES eddy viscosity as a "nut" block dataset to a field file that
 * fdm_h5_write_field already produced. Kept separate so the no-LES output is
 * byte-identical (write_field untouched) and so the Fortran side can pass nut as
 * a plain array (no c_loc, which nvfortran ICEs on for allocatable components).
 * Collective: all ranks must call it together (LES on/off is a global setting).
 * nut has the same haloed per-block layout as one velocity component. */
int fdm_h5_append_scalar(const char *filename, const char *name,
                         int nbx, int nby, int nbz,
                         int n_blocks, int n_blocks_global, int id_start,
                         const double *s)
{
    const size_t var_stride = (size_t)(nbx + 2)*(size_t)(nby + 2)*(size_t)(nbz + 2);
    hid_t plist, file;
    int ierr = 0;

    if (nbx < 1 || nby < 1 || nbz < 1 || n_blocks < 1) return 1;

    plist = H5Pcreate(H5P_FILE_ACCESS);
    if (plist < 0) return 1;
    if (H5Pset_fapl_mpio(plist, MPI_COMM_WORLD, MPI_INFO_NULL) < 0) {
        H5Pclose(plist);
        return 1;
    }
    file = H5Fopen(filename, H5F_ACC_RDWR, plist);
    H5Pclose(plist);
    if (file < 0) return 1;

    ierr |= write_block_dataset(file, name, nbx, nby, nbz,
                                n_blocks, n_blocks_global, id_start,
                                s, var_stride);
    ierr |= H5Fclose(file) < 0;
    return ierr != 0;
}

int fdm_h5_append_nut(const char *filename, int nbx, int nby, int nbz,
                      int n_blocks, int n_blocks_global, int id_start,
                      const double *nut)
{
    return fdm_h5_append_scalar(filename, "nut", nbx, nby, nbz,
                                n_blocks, n_blocks_global, id_start, nut);
}

/* [blocks] refine_dims marker of the xz-quadtree file variant: the
 * per-direction refinement mask as a 3-int attribute. Written only when
 * the mask is not all-ones, so octree-mode files stay byte-identical;
 * absent on read means the octree default {1,1,1}. Collective. */
int fdm_h5_append_refine_dims(const char *filename, const int *mask)
{
    hid_t plist, file;
    int ierr = 0;

    plist = H5Pcreate(H5P_FILE_ACCESS);
    if (plist < 0) return 1;
    if (H5Pset_fapl_mpio(plist, MPI_COMM_WORLD, MPI_INFO_NULL) < 0) {
        H5Pclose(plist);
        return 1;
    }
    file = H5Fopen(filename, H5F_ACC_RDWR, plist);
    H5Pclose(plist);
    if (file < 0) return 1;

    ierr |= write_attr_int_array(file, "refine_dims", mask, 3);
    ierr |= H5Fclose(file) < 0;
    return ierr != 0;
}

/* has_blocks reports whether the file carries a block table at all:
 * legacy global-3D restarts have no block layout, so the reader skips the
 * refine_dims cross-check for them (any leaf table can slice a global
 * field). */
int fdm_h5_read_refine_dims(const char *filename, int *mask, int *has_blocks)
{
    hid_t file;
    int ierr = 0;

    mask[0] = mask[1] = mask[2] = 1;
    *has_blocks = 0;
    file = open_parallel_file(filename);
    if (file < 0) return 1;
    *has_blocks = H5Lexists(file, "blocks", H5P_DEFAULT) > 0;
    ierr |= read_attr_int_array(file, "refine_dims", mask, 3, 0);
    ierr |= H5Fclose(file) < 0;
    return ierr != 0;
}

/* Read a named cell-centred block-layout scalar (the RANS k/omega restart
 * datasets). *found = 0 when the dataset is absent, which is not an error:
 * the solver then reinitializes the scalar and warns (old restart files).
 * Collective. Only the block-table (4D) layout is supported -- legacy
 * global 3D restarts never carried these datasets. */
int fdm_h5_read_scalar(const char *filename, const char *name,
                       int nbx, int nby, int nbz,
                       int n_blocks, int n_blocks_global, int id_start,
                       const int *block_origin, const int *block_level,
                       int *found, double *s)
{
    const size_t var_stride = (size_t)(nbx + 2)*(size_t)(nby + 2)*(size_t)(nbz + 2);
    hid_t file, dset, file_space;
    htri_t exists;
    int *row_of = NULL;
    int reordered = 0;
    int ierr = 0;

    *found = 0;
    if (nbx < 1 || nby < 1 || nbz < 1 || n_blocks < 1) return 1;

    file = open_parallel_file(filename);
    if (file < 0) return 1;

    exists = H5Lexists(file, name, H5P_DEFAULT);
    if (exists <= 0) {
        H5Fclose(file);
        return 0;
    }

    dset = H5Dopen2(file, name, H5P_DEFAULT);
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
    ierr |= block_row_map(file, n_blocks, n_blocks_global, id_start,
                          block_origin, block_level, &row_of, &reordered) != 0;
    if (!ierr) {
        ierr |= read_block_dataset(file, dset, file_space, nbx, nby, nbz,
                                   n_blocks, n_blocks_global, row_of,
                                   s, var_stride);
    }
    if (row_of != NULL) free(row_of);
    H5Sclose(file_space);
    H5Dclose(dset);
    ierr |= H5Fclose(file) < 0;
    if (!ierr) *found = 1;
    return ierr != 0;
}

/*
 * One row of n doubles per block, appended to a field file that
 * fdm_h5_write_field already produced: the (n_blocks_global, n) layout of
 * write_block_rows. Used for the outlet-normal velocity on a HIGH domain
 * face (io.f90 write_outlet_planes): that face is state, and index nb+1 is
 * not part of the per-variable block datasets. Written only by cases that
 * have such a face, so every other file is unchanged. Collective.
 */
int fdm_h5_append_block_rows(const char *filename, const char *name, int n,
                             int n_blocks, int n_blocks_global, int id_start,
                             const double *rows)
{
    hid_t plist, file;
    int ierr = 0;

    if (n < 1 || n_blocks < 1) return 1;

    plist = H5Pcreate(H5P_FILE_ACCESS);
    if (plist < 0) return 1;
    if (H5Pset_fapl_mpio(plist, MPI_COMM_WORLD, MPI_INFO_NULL) < 0) {
        H5Pclose(plist);
        return 1;
    }
    file = H5Fopen(filename, H5F_ACC_RDWR, plist);
    H5Pclose(plist);
    if (file < 0) return 1;

    ierr |= write_block_rows(file, name, n, n_blocks, n_blocks_global, id_start, rows);
    ierr |= H5Fclose(file) < 0;
    return ierr != 0;
}

/*
 * Read back a per-block row dataset (fdm_h5_append_block_rows). *found = 0
 * when the file does not carry it, which is not an error (an older
 * snapshot, a generated initial condition): the caller falls back. Rows are
 * matched to this rank's blocks on (origin, level), like every block-layout
 * dataset. Collective.
 */
int fdm_h5_read_block_rows(const char *filename, const char *name, int n,
                           int n_blocks, int n_blocks_global, int id_start,
                           const int *block_origin, const int *block_level,
                           int *found, double *rows)
{
    hid_t file, dset, file_space, xfer;
    hsize_t dims[2] = {0, 0};
    int *row_of = NULL;
    int reordered = 0;
    int ierr = 0;

    *found = 0;
    if (n < 1 || n_blocks < 1) return 1;

    file = open_parallel_file(filename);
    if (file < 0) return 1;
    if (H5Lexists(file, name, H5P_DEFAULT) <= 0) {
        H5Fclose(file);
        return 0;
    }
    dset = H5Dopen2(file, name, H5P_DEFAULT);
    if (dset < 0) {
        H5Fclose(file);
        return 1;
    }
    file_space = H5Dget_space(dset);
    xfer = H5Pcreate(H5P_DATASET_XFER);
    if (file_space < 0 || xfer < 0 || H5Sget_simple_extent_ndims(file_space) != 2) {
        ierr = 1;
    } else {
        H5Sget_simple_extent_dims(file_space, dims, NULL);
        if (dims[0] != (hsize_t)n_blocks_global || dims[1] != (hsize_t)n) ierr = 1;
    }
    if (!ierr) {
        ierr |= block_row_map(file, n_blocks, n_blocks_global, id_start,
                              block_origin, block_level, &row_of, &reordered) != 0;
    }
    if (!ierr) H5Pset_dxpl_mpio(xfer, H5FD_MPIO_INDEPENDENT);

    /* Consecutive file rows as one hyperslab (read_block_dataset's rule). */
    for (int b0 = 0; b0 < n_blocks && !ierr; ) {
        int b1 = b0;
        hsize_t start[2] = {(hsize_t)row_of[b0], 0};
        hsize_t count[2] = {0, (hsize_t)n};
        hid_t mem_space;

        while (b1 + 1 < n_blocks && row_of[b1+1] == row_of[b1] + 1) ++b1;
        count[0] = (hsize_t)(b1 - b0 + 1);
        mem_space = H5Screate_simple(2, count, NULL);
        if (mem_space < 0 ||
            H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, count, NULL) < 0 ||
            H5Dread(dset, H5T_NATIVE_DOUBLE, mem_space, file_space, xfer,
                    rows + (size_t)b0*(size_t)n) < 0) ierr = 1;
        if (mem_space >= 0) H5Sclose(mem_space);
        b0 = b1 + 1;
    }

    if (row_of != NULL) free(row_of);
    if (xfer >= 0) H5Pclose(xfer);
    if (file_space >= 0) H5Sclose(file_space);
    H5Dclose(dset);
    ierr |= H5Fclose(file) < 0;
    if (!ierr) *found = 1;
    return ierr != 0;
}

/*
 * RANS geometry diagnostic file (rans.f90, [rans] dump_geometry): the leaf
 * block table, interior dwall/yeff/wallcell rows in the block-table layout
 * and the per-block cell-centre coordinate lines. Self-contained file --
 * the field-output layout is not touched. Collective: all ranks enter.
 * dwall/yeff carry the solver's ghost layout ((nb+2)^3 per block); the
 * interior is extracted here. wallcell and xc/yc/zc are interior-only.
 */
int fdm_h5_write_rans_geometry(const char *filename, int nbx, int nby, int nbz,
                               int n_blocks, int n_blocks_global, int id_start,
                               const int *block_origin, const int *block_level,
                               const double *dwall, const double *yeff,
                               const int *wallcell,
                               const double *xc, const double *yc, const double *zc)
{
    const size_t block_stride = (size_t)(nbx + 2)*(size_t)(nby + 2)*(size_t)(nbz + 2);
    hid_t file;
    int ierr = 0;

    if (nbx < 1 || nby < 1 || nbz < 1 || n_blocks < 1) return 1;

    file = create_parallel_file(filename);
    if (file < 0) return 1;

    ierr |= write_attr_int(file, "block_nb_x", nbx);
    ierr |= write_attr_int(file, "block_nb_y", nby);
    ierr |= write_attr_int(file, "block_nb_z", nbz);
    ierr |= write_attr_int(file, "n_blocks", n_blocks_global);

    ierr |= write_block_table(file, n_blocks_global, id_start, n_blocks,
                              block_origin, block_level);
    ierr |= write_block_dataset(file, "dwall", nbx, nby, nbz,
                                n_blocks, n_blocks_global, id_start,
                                dwall, block_stride);
    ierr |= write_block_dataset(file, "yeff", nbx, nby, nbz,
                                n_blocks, n_blocks_global, id_start,
                                yeff, block_stride);
    ierr |= write_block_int_dataset(file, "wallcell", nbx, nby, nbz,
                                    n_blocks, n_blocks_global, id_start, wallcell);
    ierr |= write_block_rows(file, "xc", nbx, n_blocks, n_blocks_global, id_start, xc);
    ierr |= write_block_rows(file, "yc", nby, n_blocks, n_blocks_global, id_start, yc);
    ierr |= write_block_rows(file, "zc", nbz, n_blocks, n_blocks_global, id_start, zc);

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

/* found[v] = 0 marks a dataset the file does not carry; the caller decides
 * whether that is fatal (u,v,w,p) or a reinitialise-and-warn (a scalar added
 * after the restart file was written -- the RANS named-scalar precedent). */
int fdm_h5_read_field(const char *filename, int nbx, int nby, int nbz,
                      int n_blocks, int n_blocks_global, int id_start,
                      const int *block_origin, const int *block_level,
                      int global_nx, int global_ny, int global_nz,
                      int n_var, const char *var_names, int *found,
                      double *q)
{
    const size_t var_stride = (size_t)(nbx + 2)*(size_t)(nby + 2)*(size_t)(nbz + 2);
    const size_t block_stride = var_stride*(size_t)n_var;
    hid_t file;
    int *row_of = NULL;
    int reordered = 0;
    int ierr = 0;

    if (nbx < 1 || nby < 1 || nbz < 1 || n_blocks < 1 || n_var < 1) return 1;

    file = open_parallel_file(filename);
    if (file < 0) return 1;
    /* Block-layout datasets are read through the row map; a legacy global-3D
     * file has no blocks table and gets the (unused) identity. */
    if (block_row_map(file, n_blocks, n_blocks_global, id_start,
                      block_origin, block_level, &row_of, &reordered) != 0) {
        H5Fclose(file);
        return 1;
    }

    for (int v = 0; v < n_var; ++v) {
        const char *name = var_names + (size_t)v*FDM_VAR_NAME_LEN;
        hid_t dset;
        hid_t file_space = -1;
        int ndims = 0;

        found[v] = 0;
        if (H5Lexists(file, name, H5P_DEFAULT) <= 0) continue;
        dset = H5Dopen2(file, name, H5P_DEFAULT);
        if (dset < 0) {
            ierr = 1;
            break;
        }
        found[v] = 1;
        file_space = H5Dget_space(dset);
        if (file_space < 0) {
            H5Dclose(dset);
            ierr = 1;
            break;
        }
        ndims = H5Sget_simple_extent_ndims(file_space);
        if (ndims == 4) {
            ierr |= read_block_dataset(file, dset, file_space, nbx, nby, nbz,
                                       n_blocks, n_blocks_global, row_of,
                                       q + (size_t)v*var_stride, block_stride);
        } else {
            /* Legacy global 3D layout (pre block-table restarts): one
             * hyperslab per block, valid only for level-0 blocks whose
             * origins are level-0 cells. */
            ierr |= read_global_dataset_blocks(file, dset, file_space, nbx, nby, nbz,
                                               n_blocks, block_origin,
                                               global_nx, global_ny, global_nz,
                                               q + (size_t)v*var_stride, block_stride);
        }
        H5Sclose(file_space);
        H5Dclose(dset);
    }

    free(row_of);
    ierr |= H5Fclose(file) < 0;
    return ierr != 0;
}

/*
 * Restart guard: every block of this rank must have a row in the file's
 * blocks table (block_row_map). 0 = it has, 2 = the file is a snapshot of a
 * different tiling, 1 = read error; *reordered = 1 when the file's row order
 * is not this run's (the readers then map rows, see block_row_map). A file
 * without a blocks table (a legacy global-3D restart) has no layout to
 * disagree with: 0.
 */
int fdm_h5_check_block_table(const char *filename, int n_blocks, int n_blocks_global,
                             int id_start, const int *block_origin, const int *block_level,
                             int *reordered)
{
    hid_t file;
    int *row_of = NULL;
    int ierr;

    *reordered = 0;
    file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) return 1;
    ierr = block_row_map(file, n_blocks, n_blocks_global, id_start,
                         block_origin, block_level, &row_of, reordered);
    if (row_of != NULL) free(row_of);
    if (H5Fclose(file) < 0 && ierr == 0) ierr = 1;
    return ierr;
}

/*
 * Case-file layout header (step 7-2, docs/next_session_prepare_everything.md):
 * the attributes fdm_h5_case_create wrote plus the row count of the blocks
 * table. This is what lets the solver take its block set FROM the file
 * instead of rebuilding it from the ini. block_nb_xyz is the per-direction
 * block size (the pre-per-direction scalar block_nb is accepted); the level
 * count is NOT taken from the block_levels attribute, which the retired
 * mobygeom wrote as "number of levels" and moby_prepare as "finest level"
 * -- the reader derives it from the table itself. nb_auto/nb_ranks record a
 * DERIVED nb (the nb rule) and the rank count it was derived for; absent on
 * older files (0). key_order (capacity key_cap) receives the block_key_order
 * record, n_key_bits its length (0 = not recorded). */
int fdm_h5_case_read_layout(const char *filename, int *nx, int *ny, int *nz,
                            double *lx, double *ly, double *lz, double *re,
                            int *block_nb, int *refine_mask,
                            int *nb_auto, int *nb_ranks,
                            int key_cap, int *n_key_bits, int *key_order,
                            int *n_blocks_global)
{
    hid_t file, dset = -1, space = -1;
    hsize_t dims[2] = {0, 0};
    int nb_scalar = 0;
    int ierr = 0;

    file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) return 1;

    /* No blocks table: a legacy global-layout coefficient file (the retired
     * mobygeom stl-ibm-coeff), which carries no layout at all and none of
     * the block attributes -- reported as n_blocks_global = 0, not as an
     * error, until step 7-4 retires it. */
    *n_blocks_global = 0;
    *n_key_bits = 0;
    if (H5Lexists(file, "blocks", H5P_DEFAULT) <= 0) {
        ierr |= H5Fclose(file) < 0;
        return ierr != 0;
    }

    ierr |= read_attr_int(file, "nx", nx, 1);
    ierr |= read_attr_int(file, "ny", ny, 1);
    ierr |= read_attr_int(file, "nz", nz, 1);
    ierr |= read_attr_double(file, "lx", lx, 1);
    ierr |= read_attr_double(file, "ly", ly, 1);
    ierr |= read_attr_double(file, "lz", lz, 1);
    ierr |= read_attr_double(file, "re", re, 1);
    block_nb[0] = block_nb[1] = block_nb[2] = 0;
    if (read_attr_int_array(file, "block_nb_xyz", block_nb, 3, 0) != 0 ||
        block_nb[0] <= 0) {
        ierr |= read_attr_int(file, "block_nb", &nb_scalar, 1);
        block_nb[0] = block_nb[1] = block_nb[2] = nb_scalar;
    }
    refine_mask[0] = refine_mask[1] = refine_mask[2] = 1;
    ierr |= read_attr_int_array(file, "refine_dims", refine_mask, 3, 0);
    *nb_auto = 0;
    *nb_ranks = 0;
    ierr |= read_attr_int(file, "block_nb_auto", nb_auto, 0);
    ierr |= read_attr_int(file, "block_nb_ranks", nb_ranks, 0);
    /* The leaf key's bit order (blocks.f90 keyOrder): directions 1..3 from
     * the most to the least significant key bit. Absent on a file written
     * before 2026-10-01 (and always in xz mode): *n_key_bits = 0, which the
     * reader takes as the legacy interleave. */
    if (H5Aexists(file, "block_key_order") > 0) {
        hid_t attr = H5Aopen(file, "block_key_order", H5P_DEFAULT);
        hid_t aspace = attr >= 0 ? H5Aget_space(attr) : -1;
        hsize_t adims[1] = {0};
        if (attr < 0 || aspace < 0 || H5Sget_simple_extent_ndims(aspace) != 1) {
            ierr = 1;
        } else {
            H5Sget_simple_extent_dims(aspace, adims, NULL);
            if (adims[0] < 1 || adims[0] > (hsize_t)key_cap ||
                H5Aread(attr, H5T_NATIVE_INT, key_order) < 0) ierr = 1;
            else *n_key_bits = (int)adims[0];
        }
        if (aspace >= 0) H5Sclose(aspace);
        if (attr >= 0) H5Aclose(attr);
    }

    {
        dset = H5Dopen2(file, "blocks", H5P_DEFAULT);
        space = dset >= 0 ? H5Dget_space(dset) : -1;
        if (dset < 0 || space < 0 || H5Sget_simple_extent_ndims(space) != 2) {
            ierr = 1;
        } else {
            H5Sget_simple_extent_dims(space, dims, NULL);
            if (dims[1] != 4 || dims[0] < 1) ierr = 1;
            *n_blocks_global = (int)dims[0];
        }
        if (space >= 0) H5Sclose(space);
        if (dset >= 0) H5Dclose(dset);
    }
    ierr |= H5Fclose(file) < 0;
    return ierr != 0;
}

/* The case-input echo attribute (see fdm_h5_case_create); *found = 0 on a
 * file written before step 7-3 or by the retired mobygeom. */
int fdm_h5_case_read_inputs(const char *filename, char *buf, int cap, int *found)
{
    hid_t file;
    int ierr;

    file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) return 1;
    ierr = read_attr_string(file, "case_inputs", buf, (size_t)cap, found);
    ierr |= H5Fclose(file) < 0;
    return ierr != 0;
}

/* The whole blocks table, (n_blocks_global, 4) row-major: origin x,y,z in
 * level-l cells and the level, row i = leaf id i (Morton order). Every rank
 * reads all of it -- the table is global state, exactly as the builder
 * held it. */
int fdm_h5_case_read_blocks(const char *filename, int n_blocks_global, int *rows)
{
    hid_t file, dset = -1, space = -1;
    hsize_t dims[2] = {0, 0};
    int ierr = 0;

    file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) return 1;
    dset = H5Dopen2(file, "blocks", H5P_DEFAULT);
    space = dset >= 0 ? H5Dget_space(dset) : -1;
    if (dset < 0 || space < 0 || H5Sget_simple_extent_ndims(space) != 2) {
        ierr = 1;
    } else {
        H5Sget_simple_extent_dims(space, dims, NULL);
        if (dims[0] != (hsize_t)n_blocks_global || dims[1] != 4) {
            ierr = 1;
        } else {
            ierr |= H5Dread(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, rows) < 0;
        }
    }
    if (space >= 0) H5Sclose(space);
    if (dset >= 0) H5Dclose(dset);
    ierr |= H5Fclose(file) < 0;
    return ierr != 0;
}

static int read_dline(hid_t file, const char *name, int n, double *values)
{
    hid_t dset = -1, space = -1;
    hsize_t dims[1] = {0};
    int ierr = 0;

    dset = H5Dopen2(file, name, H5P_DEFAULT);
    space = dset >= 0 ? H5Dget_space(dset) : -1;
    if (dset < 0 || space < 0 || H5Sget_simple_extent_ndims(space) != 1) {
        ierr = 1;
    } else {
        H5Sget_simple_extent_dims(space, dims, NULL);
        if (dims[0] != (hsize_t)n) {
            ierr = 1;
        } else {
            ierr |= H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, values) < 0;
        }
    }
    if (space >= 0) H5Sclose(space);
    if (dset >= 0) H5Dclose(dset);
    return ierr;
}

/* The level-0 node lines a case file carries (fdm_h5_case_append_grid).
 * *found = 0 when the file has none (a retired-mobygeom block-table file
 * built before P3), which is not an error. */
int fdm_h5_case_read_grid(const char *filename, int nx, int ny, int nz,
                          double *x_node, double *y_node, double *z_node, int *found)
{
    hid_t file;
    int ierr = 0;

    *found = 0;
    file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) return 1;
    if (H5Lexists(file, "x_nodes", H5P_DEFAULT) > 0) {
        ierr |= read_dline(file, "x_nodes", nx + 1, x_node);
        ierr |= read_dline(file, "y_nodes", ny + 1, y_node);
        ierr |= read_dline(file, "z_nodes", nz + 1, z_node);
        if (!ierr) *found = 1;
    }
    ierr |= H5Fclose(file) < 0;
    return ierr != 0;
}

/*
 * Block-table coefficients: coef_blocks rows are (nb+2)^3 ghost windows x 3
 * staggered vars at each leaf's level, row id = leaf id of the file's own
 * blocks table -- which the solver's leaf table IS since step 7-2
 * (fdm_h5_case_read_blocks), so there is nothing to cross-check here.
 * *found = 0 when the file has no coef_blocks (legacy global layout).
 */
int fdm_h5_read_ibm_coeff_blocks(const char *filename, int nbx, int nby, int nbz,
                                 int n_blocks, int id_start,
                                 double lx, double ly, double lz, double re,
                                 int n_comp, int *found, double *coef)
{
    const size_t ni = (size_t)nbx + 2;
    const size_t nj = (size_t)nby + 2;
    const size_t nk = (size_t)nbz + 2;
    const size_t n = ni*nj*nk*3;
    /* The DESTINATION array carries n_comp components (4 with the
     * cell-centred scalar column, increment S3); the dataset always
     * carries the three staggered ones. */
    const size_t dst_stride = ni*nj*nk*(size_t)n_comp;
    hsize_t local_dims[5] = {(hsize_t)n_blocks, ni, nj, nk, 3};
    hsize_t start[5] = {(hsize_t)id_start, 0, 0, 0, 0};
    hsize_t file_dims[5] = {0, 0, 0, 0, 0};
    hid_t file = -1, dset = -1, file_space = -1, mem_space = -1;
    double *buffer = NULL;
    double file_lx = 0.0, file_ly = 0.0, file_lz = 0.0, file_re = 0.0;
    int ierr = 0;
    herr_t status = -1;

    *found = 0;
    file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) return 1;

    if (H5Lexists(file, "coef_blocks", H5P_DEFAULT) <= 0) {
        H5Fclose(file);
        return 0;
    }

    ierr |= read_attr_double(file, "lx", &file_lx, 1);
    ierr |= read_attr_double(file, "ly", &file_ly, 1);
    ierr |= read_attr_double(file, "lz", &file_lz, 1);
    ierr |= read_attr_double(file, "re", &file_re, 1);
    ierr |= double_mismatch(file_lx, lx) || double_mismatch(file_ly, ly) ||
            double_mismatch(file_lz, lz) || double_mismatch(file_re, re);
    if (ierr) {
        H5Fclose(file);
        return ierr;
    }

    dset = H5Dopen2(file, "coef_blocks", H5P_DEFAULT);
    file_space = dset >= 0 ? H5Dget_space(dset) : -1;
    if (dset < 0 || file_space < 0 || H5Sget_simple_extent_ndims(file_space) != 5) {
        if (file_space >= 0) H5Sclose(file_space);
        if (dset >= 0) H5Dclose(dset);
        H5Fclose(file);
        return 1;
    }
    H5Sget_simple_extent_dims(file_space, file_dims, NULL);
    if (file_dims[1] != ni || file_dims[2] != nj || file_dims[3] != nk || file_dims[4] != 3) {
        H5Sclose(file_space);
        H5Dclose(dset);
        H5Fclose(file);
        return 1;
    }

    buffer = (double *)malloc((size_t)n_blocks*n*sizeof(double));
    mem_space = H5Screate_simple(5, local_dims, NULL);
    if (buffer == NULL || mem_space < 0 ||
        H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, local_dims, NULL) < 0) {
        if (buffer != NULL) free(buffer);
        if (mem_space >= 0) H5Sclose(mem_space);
        H5Sclose(file_space);
        H5Dclose(dset);
        H5Fclose(file);
        return 1;
    }

    status = H5Dread(dset, H5T_NATIVE_DOUBLE, mem_space, file_space, H5P_DEFAULT, buffer);
    if (status >= 0) {
        for (int b = 0; b < n_blocks; ++b) {
            double *block_coef = coef + (size_t)b*dst_stride;
            const double *row = buffer + (size_t)b*n;
            for (size_t v = 0; v < 3; ++v) {
                for (size_t k = 0; k < nk; ++k) {
                    for (size_t j = 0; j < nj; ++j) {
                        for (size_t i = 0; i < ni; ++i) {
                            size_t h5_idx = (((i*nj) + j)*nk + k)*3 + v;
                            block_coef[linear_fortran4(i, j, k, v, ni, nj, nk)] = row[h5_idx];
                        }
                    }
                }
            }
        }
    } else {
        ierr = 1;
    }

    H5Sclose(mem_space);
    free(buffer);
    H5Sclose(file_space);
    H5Dclose(dset);
    H5Fclose(file);
    if (!ierr) *found = 1;
    return ierr;
}

/*
 * ONE ghost-inclusive per-leaf tile dataset (n_blocks_global, nb+2, nb+2,
 * nb+2), transposed into a Fortran-ordered destination whose per-block
 * stride and component-plane offset the caller supplies -- shared by
 * dwall_blocks and coef_p_blocks (which lands in the VAR_P plane of the
 * 4-component coefficient array). Row id = leaf id of the file's own blocks
 * table (see fdm_h5_read_ibm_coeff_blocks).
 * *found = 0 when the file carries no such dataset.
 */
static int read_leaf_tiles(const char *filename, const char *name,
                           int nbx, int nby, int nbz, int n_blocks, int id_start,
                           size_t dst_stride, size_t dst_offset,
                           int *found, double *dst)
{
    const size_t ni = (size_t)nbx + 2;
    const size_t nj = (size_t)nby + 2;
    const size_t nk = (size_t)nbz + 2;
    const size_t n = ni*nj*nk;
    hsize_t local_dims[4] = {(hsize_t)n_blocks, ni, nj, nk};
    hsize_t start[4] = {(hsize_t)id_start, 0, 0, 0};
    hsize_t file_dims[4] = {0, 0, 0, 0};
    hid_t file = -1, dset = -1, file_space = -1, mem_space = -1;
    double *buffer = NULL;
    int ierr = 0;
    herr_t status = -1;

    *found = 0;
    file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) return 1;

    if (H5Lexists(file, name, H5P_DEFAULT) <= 0) {
        H5Fclose(file);
        return 0;
    }

    dset = H5Dopen2(file, name, H5P_DEFAULT);
    file_space = dset >= 0 ? H5Dget_space(dset) : -1;
    if (dset < 0 || file_space < 0 || H5Sget_simple_extent_ndims(file_space) != 4) {
        if (file_space >= 0) H5Sclose(file_space);
        if (dset >= 0) H5Dclose(dset);
        H5Fclose(file);
        return 1;
    }
    H5Sget_simple_extent_dims(file_space, file_dims, NULL);
    if (file_dims[1] != ni || file_dims[2] != nj || file_dims[3] != nk) {
        H5Sclose(file_space);
        H5Dclose(dset);
        H5Fclose(file);
        return 1;
    }

    buffer = (double *)malloc((size_t)n_blocks*n*sizeof(double));
    mem_space = H5Screate_simple(4, local_dims, NULL);
    if (buffer == NULL || mem_space < 0 ||
        H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, local_dims, NULL) < 0) {
        if (buffer != NULL) free(buffer);
        if (mem_space >= 0) H5Sclose(mem_space);
        H5Sclose(file_space);
        H5Dclose(dset);
        H5Fclose(file);
        return 1;
    }

    status = H5Dread(dset, H5T_NATIVE_DOUBLE, mem_space, file_space, H5P_DEFAULT, buffer);
    if (status >= 0) {
        for (int b = 0; b < n_blocks; ++b) {
            double *block_dst = dst + (size_t)b*dst_stride + dst_offset;
            const double *row = buffer + (size_t)b*n;
            for (size_t k = 0; k < nk; ++k) {
                for (size_t j = 0; j < nj; ++j) {
                    for (size_t i = 0; i < ni; ++i) {
                        block_dst[linear_fortran(i, j, k, ni, nj)] =
                            row[(i*nj + j)*nk + k];
                    }
                }
            }
        }
    } else {
        ierr = 1;
    }

    H5Sclose(mem_space);
    free(buffer);
    H5Sclose(file_space);
    H5Dclose(dset);
    H5Fclose(file);
    if (!ierr) *found = 1;
    return ierr;
}

/* Per-leaf wall-distance tiles (mobygeom block-table): dwall_blocks rows
 * are (nb+2)^3 ghost windows of the cell-centred distance to the immersed
 * surface, evaluated at each leaf's level. */
int fdm_h5_read_dwall_blocks(const char *filename, int nbx, int nby, int nbz,
                             int n_blocks, int id_start,
                             int *found, double *dwall)
{
    const size_t n = ((size_t)nbx + 2)*((size_t)nby + 2)*((size_t)nbz + 2);

    return read_leaf_tiles(filename, "dwall_blocks", nbx, nby, nbz, n_blocks,
                           id_start, n, 0, found, dwall);
}

/* Per-leaf CELL-CENTRED (pressure-position) IBM coefficient tiles, the
 * passive scalars' penalization coefficient (increment S3). Optional: only
 * a case file prepared with [scalar] carries coef_p_blocks, and it lands in
 * component 3 (0-based) of the solver's 4-component coefficient array. */
int fdm_h5_read_ibm_coeff_p_blocks(const char *filename, int nbx, int nby, int nbz,
                                   int n_blocks, int id_start,
                                   int n_comp, int *found, double *coef)
{
    const size_t n = ((size_t)nbx + 2)*((size_t)nby + 2)*((size_t)nbz + 2);

    if (n_comp < 4) return 1;
    return read_leaf_tiles(filename, "coef_p_blocks", nbx, nby, nbz, n_blocks,
                           id_start, n*(size_t)n_comp, n*3, found, coef);
}

/*
 * Case-file writers (moby_prepare, docs/prepare_solve_strategy.md P0): the
 * preprocessing output IS the block-table coefficient file the solver
 * already reads through the fdm_h5_read_* functions above, so each writer
 * below is the exact inverse of its reader (dataset names, shapes, index
 * order, attributes). All functions are collective over MPI_COMM_WORLD;
 * per-leaf datasets are written as independent contiguous row ranges,
 * lattice-global rasters once by the rank called with write_data = 1.
 */

static hid_t open_parallel_rdwr(const char *filename)
{
    hid_t plist = H5Pcreate(H5P_FILE_ACCESS);
    hid_t file;

    if (plist < 0) return -1;
    if (H5Pset_fapl_mpio(plist, MPI_COMM_WORLD, MPI_INFO_NULL) < 0) {
        H5Pclose(plist);
        return -1;
    }
    file = H5Fopen(filename, H5F_ACC_RDWR, plist);
    H5Pclose(plist);
    return file;
}

/* Create the case file: the header attributes every reader checks plus the
 * blocks leaf table. refine_mask follows the fdm_h5_append_refine_dims
 * convention (written only when not all-ones). */
int fdm_h5_case_create(const char *filename,
                       int nx, int ny, int nz,
                       double lx, double ly, double lz, double re,
                       const int *block_nb, int block_levels, const int *refine_mask,
                       int nb_auto, int nb_ranks, const char *inputs,
                       int n_key_bits, const int *key_order,
                       int n_blocks_global, int id_start, int n_blocks,
                       const int *block_origin, const int *block_level)
{
    hid_t file;
    int ierr = 0;

    if (n_blocks < 1 || n_blocks_global < 1) return 1;

    file = create_parallel_file(filename);
    if (file < 0) return 1;

    ierr |= write_attr_int(file, "nx", nx);
    ierr |= write_attr_int(file, "ny", ny);
    ierr |= write_attr_int(file, "nz", nz);
    ierr |= write_attr_double(file, "lx", lx);
    ierr |= write_attr_double(file, "ly", ly);
    ierr |= write_attr_double(file, "lz", lz);
    ierr |= write_attr_double(file, "re", re);
    /* Per-direction block size. The legacy scalar is kept for cubic layouts so
     * case files stay readable by older builds and by the retired-mobygeom
     * reference tooling; a non-cubic layout has no meaningful scalar. */
    ierr |= write_attr_int_array(file, "block_nb_xyz", block_nb, 3);
    if (block_nb[0] == block_nb[1] && block_nb[1] == block_nb[2]) {
        ierr |= write_attr_int(file, "block_nb", block_nb[0]);
    }
    ierr |= write_attr_int(file, "block_levels", block_levels);
    /* The nb rule (step 7-1): a DERIVED block size is a property of the file,
     * tagged with the rank count it was derived for, so a solve on another
     * rank count can tell "one block per rank" from an explicit choice. */
    ierr |= write_attr_int(file, "block_nb_auto", nb_auto);
    ierr |= write_attr_int(file, "block_nb_ranks", nb_ranks);
    /* The INPUT ECHO (step 7-3): every ini value this file is a function
     * of, one key=value per line (config.f90 case_input_echo). The reader
     * rebuilds the same text from its ini and names the first line that
     * differs -- staleness attribute by attribute, readable with h5dump. */
    ierr |= write_attr_string(file, "case_inputs", inputs);
    if (refine_mask[0] != 1 || refine_mask[1] != 1 || refine_mask[2] != 1) {
        ierr |= write_attr_int_array(file, "refine_dims", refine_mask, 3);
    }
    /* The ROW ORDER of the blocks table: the leaf key's bit order (blocks.f90
     * keyOrder), directions 1..3 from the most to the least significant bit.
     * Not written in xz mode (n_key_bits = 0), which has one fixed key. */
    if (n_key_bits > 0) {
        ierr |= write_attr_int_array(file, "block_key_order", key_order, (hsize_t)n_key_bits);
    }

    ierr |= write_block_table(file, n_blocks_global, id_start, n_blocks,
                              block_origin, block_level);

    ierr |= H5Fclose(file) < 0;
    return ierr != 0;
}

/* One lattice-global int raster (block_active, block_touch_l*, ...):
 * collectively created, written by the write_data rank only. */
static int case_append_raster(const char *filename, const char *name,
                              int n_raster, const int *values, int write_data)
{
    hsize_t dims[1] = {(hsize_t)n_raster};
    hid_t file = -1, space = -1, dset = -1, xfer = -1;
    int ierr = 0;

    file = open_parallel_rdwr(filename);
    if (file < 0) return 1;

    space = H5Screate_simple(1, dims, NULL);
    dset = space >= 0 ? H5Dcreate2(file, name, H5T_NATIVE_INT, space,
                                   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) : -1;
    xfer = H5Pcreate(H5P_DATASET_XFER);
    if (space < 0 || dset < 0 || xfer < 0) {
        ierr = 1;
    } else if (write_data) {
        H5Pset_dxpl_mpio(xfer, H5FD_MPIO_INDEPENDENT);
        ierr |= H5Dwrite(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, xfer, values) < 0;
    }

    if (xfer >= 0) H5Pclose(xfer);
    if (dset >= 0) H5Dclose(dset);
    if (space >= 0) H5Sclose(space);
    ierr |= H5Fclose(file) < 0;
    return ierr != 0;
}

int fdm_h5_case_append_masks(const char *filename, int level, int n_raster,
                             const int *touch, const int *buried, int write_data)
{
    char name[64];
    int ierr = 0;

    snprintf(name, sizeof(name), "block_touch_l%d", level);
    ierr |= case_append_raster(filename, name, n_raster, touch, write_data);
    snprintf(name, sizeof(name), "block_buried_l%d", level);
    ierr |= case_append_raster(filename, name, n_raster, buried, write_data);
    return ierr != 0;
}

/* WINDOWED per-level mask attributes (deep-refinement case files: the
 * rasters cover only the window; the reader treats absent attrs as full
 * rasters). Collective like the other case appends. */
int fdm_h5_case_append_mask_window(const char *filename, int level,
                                   const int *lo, const int *dims)
{
    hid_t plist, file;
    char name[64];
    int ierr = 0;

    plist = H5Pcreate(H5P_FILE_ACCESS);
    if (plist < 0) return 1;
    if (H5Pset_fapl_mpio(plist, MPI_COMM_WORLD, MPI_INFO_NULL) < 0) {
        H5Pclose(plist);
        return 1;
    }
    file = H5Fopen(filename, H5F_ACC_RDWR, plist);
    H5Pclose(plist);
    if (file < 0) return 1;
    snprintf(name, sizeof(name), "mask_win_lo_l%d", level);
    ierr |= write_attr_int_array(file, name, lo, 3);
    snprintf(name, sizeof(name), "mask_win_dims_l%d", level);
    ierr |= write_attr_int_array(file, name, dims, 3);
    ierr |= H5Fclose(file) < 0;
    return ierr != 0;
}

int fdm_h5_case_append_active(const char *filename, int n_lattice,
                              const int *active, int write_data)
{
    return case_append_raster(filename, "block_active", n_lattice,
                              active, write_data);
}

/* coef_blocks: (n_blocks_global, nb+2, nb+2, nb+2, 3) ghost-inclusive
 * per-leaf staggered-coefficient tiles; the transpose is the exact inverse
 * of fdm_h5_read_ibm_coeff_blocks. */
int fdm_h5_case_append_coef(const char *filename, int nbx, int nby, int nbz,
                            int n_blocks, int n_blocks_global, int id_start,
                            int n_comp, const double *coef)
{
    const size_t ni = (size_t)nbx + 2;
    const size_t nj = (size_t)nby + 2;
    const size_t nk = (size_t)nbz + 2;
    const size_t n = ni*nj*nk*3;
    /* Source stride: the coefficient array may carry a fourth,
     * cell-centred component (coef_p_blocks writes it separately). */
    const size_t src_stride = ni*nj*nk*(size_t)n_comp;
    hsize_t global_dims[5] = {(hsize_t)n_blocks_global, ni, nj, nk, 3};
    hsize_t local_dims[5] = {(hsize_t)n_blocks, ni, nj, nk, 3};
    hsize_t start[5] = {(hsize_t)id_start, 0, 0, 0, 0};
    hid_t file = -1, file_space = -1, mem_space = -1, dset = -1, xfer = -1;
    double *buffer = NULL;
    int ierr = 0;

    if (nbx < 1 || nby < 1 || nbz < 1 || n_blocks < 1) return 1;

    file = open_parallel_rdwr(filename);
    if (file < 0) return 1;

    buffer = (double *)malloc((size_t)n_blocks*n*sizeof(double));
    if (buffer == NULL) {
        H5Fclose(file);
        return 1;
    }
    for (int b = 0; b < n_blocks; ++b) {
        const double *block_coef = coef + (size_t)b*src_stride;
        double *row = buffer + (size_t)b*n;
        for (size_t v = 0; v < 3; ++v) {
            for (size_t k = 0; k < nk; ++k) {
                for (size_t j = 0; j < nj; ++j) {
                    for (size_t i = 0; i < ni; ++i) {
                        size_t h5_idx = (((i*nj) + j)*nk + k)*3 + v;
                        row[h5_idx] = block_coef[linear_fortran4(i, j, k, v, ni, nj, nk)];
                    }
                }
            }
        }
    }

    file_space = H5Screate_simple(5, global_dims, NULL);
    dset = file_space >= 0 ? H5Dcreate2(file, "coef_blocks", H5T_NATIVE_DOUBLE,
                                        file_space, H5P_DEFAULT, H5P_DEFAULT,
                                        H5P_DEFAULT) : -1;
    mem_space = H5Screate_simple(5, local_dims, NULL);
    xfer = H5Pcreate(H5P_DATASET_XFER);
    if (file_space < 0 || dset < 0 || mem_space < 0 || xfer < 0 ||
        H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, local_dims, NULL) < 0) {
        ierr = 1;
    } else {
        H5Pset_dxpl_mpio(xfer, H5FD_MPIO_INDEPENDENT);
        ierr |= H5Dwrite(dset, H5T_NATIVE_DOUBLE, mem_space, file_space, xfer, buffer) < 0;
    }

    if (xfer >= 0) H5Pclose(xfer);
    if (mem_space >= 0) H5Sclose(mem_space);
    if (dset >= 0) H5Dclose(dset);
    if (file_space >= 0) H5Sclose(file_space);
    free(buffer);
    ierr |= H5Fclose(file) < 0;
    return ierr != 0;
}

/* Grid datasets in the case file (P3, mobygrid absorbed): the node lines
 * plus the mobygrid-format attributes, so the case file doubles as the
 * --grid-file of the retired mobygeom reference tooling (which needs
 * mobygrid_format = 1, {x,y,z}_nodes and validates the grid attributes
 * when present; the /staggered inspection group was never read).
 * Collective open/creates; data written by the write_data rank. */
static int append_dline(hid_t file, const char *name, int n,
                        const double *values, int write_data)
{
    hsize_t dims[1] = {(hsize_t)n};
    hid_t space = -1, dset = -1, xfer = -1;
    int ierr = 0;

    space = H5Screate_simple(1, dims, NULL);
    dset = space >= 0 ? H5Dcreate2(file, name, H5T_NATIVE_DOUBLE, space,
                                   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT) : -1;
    xfer = H5Pcreate(H5P_DATASET_XFER);
    if (space < 0 || dset < 0 || xfer < 0) {
        ierr = 1;
    } else if (write_data) {
        H5Pset_dxpl_mpio(xfer, H5FD_MPIO_INDEPENDENT);
        ierr |= H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, xfer, values) < 0;
    }
    if (xfer >= 0) H5Pclose(xfer);
    if (dset >= 0) H5Dclose(dset);
    if (space >= 0) H5Sclose(space);
    return ierr != 0;
}

int fdm_h5_case_append_grid(const char *filename, int nx, int ny, int nz,
                            const int *periodic, const int *grid_distribution,
                            const double *grid_stretch,
                            const double *grid_natural_dyw_plus,
                            const int *grid_natural_one_sided,
                            const double *x_node, const double *y_node,
                            const double *z_node, int write_data)
{
    hid_t file;
    int ierr = 0;

    file = open_parallel_rdwr(filename);
    if (file < 0) return 1;

    ierr |= write_attr_int(file, "mobygrid_format", 1);
    ierr |= write_attr_int_array(file, "periodic", periodic, 3);
    ierr |= write_attr_int_array(file, "grid_distribution", grid_distribution, 3);
    ierr |= write_attr_double_array(file, "grid_stretch", grid_stretch, 3);
    ierr |= write_attr_double_array(file, "grid_natural_dyw_plus", grid_natural_dyw_plus, 3);
    ierr |= write_attr_int_array(file, "grid_natural_one_sided", grid_natural_one_sided, 3);

    ierr |= append_dline(file, "x_nodes", nx + 1, x_node, write_data);
    ierr |= append_dline(file, "y_nodes", ny + 1, y_node, write_data);
    ierr |= append_dline(file, "z_nodes", nz + 1, z_node, write_data);

    ierr |= H5Fclose(file) < 0;
    return ierr != 0;
}

/* One ghost-inclusive per-leaf tile dataset (n_blocks_global, nb+2, nb+2,
 * nb+2) from a Fortran-ordered source with the caller's per-block stride
 * and component-plane offset -- the exact inverse of read_leaf_tiles, and
 * shared by dwall_blocks and coef_p_blocks. */
static int case_append_leaf_tiles(const char *filename, const char *name,
                                  int nbx, int nby, int nbz, int n_blocks,
                                  int n_blocks_global, int id_start,
                                  size_t src_stride, size_t src_offset,
                                  const double *src)
{
    const size_t ni = (size_t)nbx + 2;
    const size_t nj = (size_t)nby + 2;
    const size_t nk = (size_t)nbz + 2;
    const size_t n = ni*nj*nk;
    hsize_t global_dims[4] = {(hsize_t)n_blocks_global, ni, nj, nk};
    hsize_t local_dims[4] = {(hsize_t)n_blocks, ni, nj, nk};
    hsize_t start[4] = {(hsize_t)id_start, 0, 0, 0};
    hid_t file = -1, file_space = -1, mem_space = -1, dset = -1, xfer = -1;
    double *buffer = NULL;
    int ierr = 0;

    if (nbx < 1 || nby < 1 || nbz < 1 || n_blocks < 1) return 1;

    file = open_parallel_rdwr(filename);
    if (file < 0) return 1;

    buffer = (double *)malloc((size_t)n_blocks*n*sizeof(double));
    if (buffer == NULL) {
        H5Fclose(file);
        return 1;
    }
    for (int b = 0; b < n_blocks; ++b) {
        const double *block_src = src + (size_t)b*src_stride + src_offset;
        double *row = buffer + (size_t)b*n;
        for (size_t k = 0; k < nk; ++k) {
            for (size_t j = 0; j < nj; ++j) {
                for (size_t i = 0; i < ni; ++i) {
                    row[(i*nj + j)*nk + k] = block_src[linear_fortran(i, j, k, ni, nj)];
                }
            }
        }
    }

    file_space = H5Screate_simple(4, global_dims, NULL);
    dset = file_space >= 0 ? H5Dcreate2(file, name, H5T_NATIVE_DOUBLE,
                                        file_space, H5P_DEFAULT, H5P_DEFAULT,
                                        H5P_DEFAULT) : -1;
    mem_space = H5Screate_simple(4, local_dims, NULL);
    xfer = H5Pcreate(H5P_DATASET_XFER);
    if (file_space < 0 || dset < 0 || mem_space < 0 || xfer < 0 ||
        H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, local_dims, NULL) < 0) {
        ierr = 1;
    } else {
        H5Pset_dxpl_mpio(xfer, H5FD_MPIO_INDEPENDENT);
        ierr |= H5Dwrite(dset, H5T_NATIVE_DOUBLE, mem_space, file_space, xfer, buffer) < 0;
    }

    if (xfer >= 0) H5Pclose(xfer);
    if (mem_space >= 0) H5Sclose(mem_space);
    if (dset >= 0) H5Dclose(dset);
    if (file_space >= 0) H5Sclose(file_space);
    free(buffer);
    ierr |= H5Fclose(file) < 0;
    return ierr != 0;
}

/* dwall_blocks: ghost-inclusive cell-centred raw body-distance tiles;
 * inverse of fdm_h5_read_dwall_blocks. */
int fdm_h5_case_append_dwall(const char *filename, int nbx, int nby, int nbz,
                             int n_blocks, int n_blocks_global, int id_start,
                             const double *dwall)
{
    const size_t n = ((size_t)nbx + 2)*((size_t)nby + 2)*((size_t)nbz + 2);

    return case_append_leaf_tiles(filename, "dwall_blocks", nbx, nby, nbz, n_blocks,
                                  n_blocks_global, id_start, n, 0, dwall);
}

/* coef_p_blocks: the OPTIONAL cell-centred (pressure-position) coefficient
 * tiles the passive scalars penalise with, written only when the prepare
 * ini declares [scalar] -- so a case file for a scalar-free run is
 * byte-identical to what P0-P3 produced. Source component 3 (0-based) of
 * the 4-component coefficient array; inverse of
 * fdm_h5_read_ibm_coeff_p_blocks. */
int fdm_h5_case_append_coef_p(const char *filename, int nbx, int nby, int nbz,
                              int n_blocks, int n_blocks_global, int id_start,
                              int n_comp, const double *coef)
{
    const size_t n = ((size_t)nbx + 2)*((size_t)nby + 2)*((size_t)nbz + 2);

    if (n_comp < 4) return 1;
    return case_append_leaf_tiles(filename, "coef_p_blocks", nbx, nby, nbz, n_blocks,
                                  n_blocks_global, id_start,
                                  n*(size_t)n_comp, n*3, coef);
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

/* Boundary-layer statistics: profiles on the global (x,y) plane, averaged
 * in the spanwise (z) direction and time. The stat tables are stored as
 * (nx*ny, nstat) with the (x,y) cell flattened y-fastest (row = (gx-1)*ny +
 * gy), reusing the 1D/2D helpers. Serial write from rank 0 (the caller has
 * already allreduce-summed the accumulators). */
int fdm_h5_write_bl_stats(const char *filename, int nx, int ny, int nstat,
                          int step, double t_current, double re,
                          const double *xcoord, const double *ycoord,
                          const double *profile, const double *raw_sum,
                          const double *count)
{
    hid_t file;
    int npt = nx*ny;
    int ierr = 0;

    if (nx < 1 || ny < 1 || nstat < 1) return 1;

    file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file < 0) return 1;

    ierr |= write_attr_int(file, "nx", nx);
    ierr |= write_attr_int(file, "ny", ny);
    ierr |= write_attr_int(file, "nstat", nstat);
    ierr |= write_attr_int(file, "sample_weighting", 1);
    ierr |= write_attr_int(file, "step", step);
    ierr |= write_attr_double(file, "t_current", t_current);
    ierr |= write_attr_double(file, "re", re);
    ierr |= write_dataset1(file, "xcoord", (hsize_t)nx, xcoord);
    ierr |= write_dataset1(file, "ycoord", (hsize_t)ny, ycoord);
    ierr |= write_dataset1(file, "count", (hsize_t)npt, count);
    ierr |= write_dataset2_fortran(file, "profile", npt, nstat, profile);
    ierr |= write_dataset2_fortran(file, "raw_sum", npt, nstat, raw_sum);

    ierr |= H5Fclose(file) < 0;
    return ierr != 0;
}

int fdm_h5_read_bl_stats(const char *filename, int nx, int ny, int nstat,
                         int *step, double *t_current,
                         double *raw_sum, double *count)
{
    hid_t file;
    int file_nx = nx, file_ny = ny, file_nstat = nstat;
    int sample_weighting = 0;
    int npt = nx*ny;
    int ierr = 0;

    file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) return 1;

    ierr |= read_attr_int(file, "nx", &file_nx, 1);
    ierr |= read_attr_int(file, "ny", &file_ny, 1);
    ierr |= read_attr_int(file, "nstat", &file_nstat, 1);
    ierr |= read_attr_int(file, "sample_weighting", &sample_weighting, 1);
    ierr |= read_attr_int(file, "step", step, 1);
    ierr |= read_attr_double(file, "t_current", t_current, 1);

    if (file_nx != nx || file_ny != ny || file_nstat != nstat || sample_weighting != 1) {
        H5Fclose(file);
        return 1;
    }

    ierr |= read_dataset1(file, "count", (hsize_t)npt, count);
    ierr |= read_dataset2_fortran(file, "raw_sum", npt, nstat, raw_sum);

    ierr |= H5Fclose(file) < 0;
    return ierr != 0;
}
