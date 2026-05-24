#include <mpi.h>
#include <hdf5.h>
#include <stdlib.h>

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

static int write_decomposition(hid_t file, int rank, int nranks,
                               int local_i_first, int local_i_last,
                               int local_j_first, int local_j_last,
                               int local_k_first, int local_k_last)
{
    hsize_t dims[2] = {(hsize_t)nranks, 6};
    hsize_t start[2] = {(hsize_t)rank, 0};
    hsize_t count[2] = {1, 6};
    int range[6] = {
        local_i_first, local_i_last,
        local_j_first, local_j_last,
        local_k_first, local_k_last
    };
    hid_t file_space = -1;
    hid_t mem_space = -1;
    hid_t dset = -1;
    hid_t xfer = -1;
    herr_t status;

    file_space = H5Screate_simple(2, dims, NULL);
    if (file_space < 0) return 1;

    dset = H5Dcreate2(file, "rank_local_range", H5T_NATIVE_INT, file_space,
                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dset < 0) {
        H5Sclose(file_space);
        return 1;
    }

    if (H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, count, NULL) < 0) {
        H5Dclose(dset);
        H5Sclose(file_space);
        return 1;
    }

    mem_space = H5Screate_simple(2, count, NULL);
    if (mem_space < 0) {
        H5Dclose(dset);
        H5Sclose(file_space);
        return 1;
    }

    xfer = H5Pcreate(H5P_DATASET_XFER);
    if (xfer < 0) {
        H5Sclose(mem_space);
        H5Dclose(dset);
        H5Sclose(file_space);
        return 1;
    }
    H5Pset_dxpl_mpio(xfer, H5FD_MPIO_COLLECTIVE);

    status = H5Dwrite(dset, H5T_NATIVE_INT, mem_space, file_space, xfer, range);

    H5Pclose(xfer);
    H5Sclose(mem_space);
    H5Dclose(dset);
    H5Sclose(file_space);
    return status < 0;
}

static int write_global_dataset3(hid_t file, const char *name,
                                 int nx, int ny, int nz,
                                 int global_nx, int global_ny, int global_nz,
                                 int local_i_first, int local_j_first, int local_k_first,
                                 const double *field)
{
    const size_t ni = (size_t)nx + 2;
    const size_t nj = (size_t)ny + 2;
    const size_t n = (size_t)nx*(size_t)ny*(size_t)nz;
    hsize_t global_dims[3] = {(hsize_t)global_nz, (hsize_t)global_ny, (hsize_t)global_nx};
    hsize_t local_dims[3] = {(hsize_t)nz, (hsize_t)ny, (hsize_t)nx};
    hsize_t start[3] = {
        (hsize_t)local_k_first - 1,
        (hsize_t)local_j_first - 1,
        (hsize_t)local_i_first - 1
    };
    hid_t file_space = -1;
    hid_t mem_space = -1;
    hid_t dset = -1;
    hid_t xfer = -1;
    double *buffer = NULL;
    herr_t status;

    buffer = (double *)malloc(n*sizeof(double));
    if (buffer == NULL) return 1;

    for (size_t k = 1; k <= (size_t)nz; ++k) {
        for (size_t j = 1; j <= (size_t)ny; ++j) {
            for (size_t i = 1; i <= (size_t)nx; ++i) {
                buffer[linear_hdf5(k-1, j-1, i-1, (size_t)ny, (size_t)nx)] =
                    field[linear_fortran(i, j, k, ni, nj)];
            }
        }
    }

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

    if (H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, local_dims, NULL) < 0) {
        H5Dclose(dset);
        H5Sclose(file_space);
        free(buffer);
        return 1;
    }

    mem_space = H5Screate_simple(3, local_dims, NULL);
    if (mem_space < 0) {
        H5Dclose(dset);
        H5Sclose(file_space);
        free(buffer);
        return 1;
    }

    xfer = H5Pcreate(H5P_DATASET_XFER);
    if (xfer < 0) {
        H5Sclose(mem_space);
        H5Dclose(dset);
        H5Sclose(file_space);
        free(buffer);
        return 1;
    }
    H5Pset_dxpl_mpio(xfer, H5FD_MPIO_COLLECTIVE);

    status = H5Dwrite(dset, H5T_NATIVE_DOUBLE, mem_space, file_space, xfer, buffer);

    H5Pclose(xfer);
    H5Sclose(mem_space);
    H5Dclose(dset);
    H5Sclose(file_space);
    free(buffer);
    return status < 0;
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

static int read_global_dataset3(hid_t file, const char *name,
                                int nx, int ny, int nz,
                                int global_nx, int global_ny, int global_nz,
                                int local_i_first, int local_j_first, int local_k_first,
                                double *field)
{
    const size_t ni = (size_t)nx + 2;
    const size_t nj = (size_t)ny + 2;
    const size_t n = (size_t)nx*(size_t)ny*(size_t)nz;
    hsize_t expected_dims[3] = {(hsize_t)global_nz, (hsize_t)global_ny, (hsize_t)global_nx};
    hsize_t file_dims[3] = {0, 0, 0};
    hsize_t local_dims[3] = {(hsize_t)nz, (hsize_t)ny, (hsize_t)nx};
    hsize_t start[3] = {
        (hsize_t)local_k_first - 1,
        (hsize_t)local_j_first - 1,
        (hsize_t)local_i_first - 1
    };
    hid_t dset = -1;
    hid_t file_space = -1;
    hid_t mem_space = -1;
    hid_t xfer = -1;
    double *buffer = NULL;
    herr_t status;

    dset = H5Dopen2(file, name, H5P_DEFAULT);
    if (dset < 0) return 1;

    file_space = H5Dget_space(dset);
    if (file_space < 0) {
        H5Dclose(dset);
        return 1;
    }

    if (H5Sget_simple_extent_ndims(file_space) != 3) {
        H5Sclose(file_space);
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
    if (buffer == NULL) {
        H5Sclose(file_space);
        H5Dclose(dset);
        return 1;
    }

    if (H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, local_dims, NULL) < 0) {
        free(buffer);
        H5Sclose(file_space);
        H5Dclose(dset);
        return 1;
    }

    mem_space = H5Screate_simple(3, local_dims, NULL);
    if (mem_space < 0) {
        free(buffer);
        H5Sclose(file_space);
        H5Dclose(dset);
        return 1;
    }

    xfer = H5Pcreate(H5P_DATASET_XFER);
    if (xfer < 0) {
        H5Sclose(mem_space);
        free(buffer);
        H5Sclose(file_space);
        H5Dclose(dset);
        return 1;
    }
    H5Pset_dxpl_mpio(xfer, H5FD_MPIO_COLLECTIVE);

    status = H5Dread(dset, H5T_NATIVE_DOUBLE, mem_space, file_space, xfer, buffer);
    if (status >= 0) {
        for (size_t k = 1; k <= (size_t)nz; ++k) {
            for (size_t j = 1; j <= (size_t)ny; ++j) {
                for (size_t i = 1; i <= (size_t)nx; ++i) {
                    field[linear_fortran(i, j, k, ni, nj)] =
                        buffer[linear_hdf5(k-1, j-1, i-1, (size_t)ny, (size_t)nx)];
                }
            }
        }
    }

    H5Pclose(xfer);
    H5Sclose(mem_space);
    free(buffer);
    H5Sclose(file_space);
    H5Dclose(dset);
    return status < 0;
}

int fdm_h5_write_field(const char *filename, int nx, int ny, int nz,
                       int rank, int nranks,
                       int global_nx, int global_ny, int global_nz,
                       int local_i_first, int local_i_last,
                       int local_j_first, int local_j_last,
                       int local_k_first, int local_k_last,
                       int step, int nsteps,
                       double lx, double ly, double lz,
                       double dx, double dy, double dz,
                       double re, double dt, double t_final, double t_current,
                       double cfl, double cflmax, double dtmax,
                       const double *forcing,
                       int pressure_niter, double pressure_sor, int ibm_enabled, int bc_count,
                       const int *periodic, const int *bc_type, const double *bc_value,
                       const int *grid_distribution, const double *grid_stretch,
                       const double *x_node, const double *y_node, const double *z_node,
                       const double *un, const double *vn,
                       const double *wn, const double *pn)
{
    hid_t file;
    int ierr = 0;

    if (nx < 1 || ny < 1 || nz < 1) return 1;

    file = create_parallel_file(filename);
    if (file < 0) return 1;

    ierr |= write_attr_int(file, "nx", global_nx);
    ierr |= write_attr_int(file, "ny", global_ny);
    ierr |= write_attr_int(file, "nz", global_nz);
    ierr |= write_attr_int(file, "nranks", nranks);
    ierr |= write_attr_int(file, "parallel_hdf5", 1);
    ierr |= write_attr_int(file, "step", step);
    ierr |= write_attr_int(file, "nsteps", nsteps);
    ierr |= write_attr_double(file, "lx", lx);
    ierr |= write_attr_double(file, "ly", ly);
    ierr |= write_attr_double(file, "lz", lz);
    ierr |= write_attr_double(file, "dx", dx);
    ierr |= write_attr_double(file, "dy", dy);
    ierr |= write_attr_double(file, "dz", dz);
    ierr |= write_attr_double(file, "re", re);
    ierr |= write_attr_double(file, "dt", dt);
    ierr |= write_attr_double(file, "t_final", t_final);
    ierr |= write_attr_double(file, "t_current", t_current);
    ierr |= write_attr_double(file, "cfl", cfl);
    ierr |= write_attr_double(file, "cflmax", cflmax);
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

    ierr |= write_decomposition(file, rank, nranks,
                                local_i_first, local_i_last,
                                local_j_first, local_j_last,
                                local_k_first, local_k_last);

    ierr |= write_global_dataset3(file, "un", nx, ny, nz,
                                  global_nx, global_ny, global_nz,
                                  local_i_first, local_j_first, local_k_first, un);
    ierr |= write_global_dataset3(file, "vn", nx, ny, nz,
                                  global_nx, global_ny, global_nz,
                                  local_i_first, local_j_first, local_k_first, vn);
    ierr |= write_global_dataset3(file, "wn", nx, ny, nz,
                                  global_nx, global_ny, global_nz,
                                  local_i_first, local_j_first, local_k_first, wn);
    ierr |= write_global_dataset3(file, "pn", nx, ny, nz,
                                  global_nx, global_ny, global_nz,
                                  local_i_first, local_j_first, local_k_first, pn);
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
                         double *cfl, double *cflmax, double *dtmax,
                         double *forcing,
                         int *pressure_niter, double *pressure_sor, int *ibm_enabled, int bc_count,
                         int *periodic, int *bc_type, double *bc_value,
                         int *grid_distribution, double *grid_stretch)
{
    hid_t file;
    int ierr = 0;

    file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) return 1;

    ierr |= read_attr_int(file, "nx", global_nx, 1);
    ierr |= read_attr_int(file, "ny", global_ny, 1);
    ierr |= read_attr_int(file, "nz", global_nz, 1);
    ierr |= read_attr_int(file, "step", step, 0);
    ierr |= read_attr_int(file, "nsteps", nsteps, 0);

    ierr |= read_attr_double(file, "lx", lx, 1);
    ierr |= read_attr_double(file, "ly", ly, 1);
    ierr |= read_attr_double(file, "lz", lz, 1);
    ierr |= read_attr_double(file, "re", re, 1);
    ierr |= read_attr_double(file, "dt", dt, 1);
    ierr |= read_attr_double(file, "t_final", t_final, 0);
    ierr |= read_attr_double(file, "t_current", t_current, 1);
    ierr |= read_attr_double(file, "cfl", cfl, 0);
    ierr |= read_attr_double(file, "cflmax", cflmax, 0);
    ierr |= read_attr_double(file, "dtmax", dtmax, 0);
    ierr |= read_attr_double(file, "forcing_x", &forcing[0], 0);
    ierr |= read_attr_double(file, "forcing_y", &forcing[1], 0);
    ierr |= read_attr_double(file, "forcing_z", &forcing[2], 0);
    ierr |= read_attr_int(file, "pressure_niter", pressure_niter, 0);
    ierr |= read_attr_double(file, "pressure_sor", pressure_sor, 0);
    ierr |= read_attr_int(file, "ibm_enabled", ibm_enabled, 1);
    ierr |= read_attr_int_array(file, "periodic", periodic, 3, 0);
    ierr |= read_attr_int_array(file, "bc_type", bc_type, (hsize_t)bc_count, 0);
    ierr |= read_attr_double_array(file, "bc_value", bc_value, (hsize_t)bc_count, 0);
    ierr |= read_attr_int_array(file, "grid_distribution", grid_distribution, 3, 1);
    ierr |= read_attr_double_array(file, "grid_stretch", grid_stretch, 3, 1);

    ierr |= H5Fclose(file) < 0;
    return ierr != 0;
}

int fdm_h5_read_field(const char *filename, int nx, int ny, int nz,
                      int global_nx, int global_ny, int global_nz,
                      int local_i_first, int local_j_first, int local_k_first,
                      double *un, double *vn, double *wn, double *pn)
{
    hid_t file;
    int ierr = 0;

    if (nx < 1 || ny < 1 || nz < 1) return 1;

    file = open_parallel_file(filename);
    if (file < 0) return 1;

    ierr |= read_global_dataset3(file, "un", nx, ny, nz,
                                 global_nx, global_ny, global_nz,
                                 local_i_first, local_j_first, local_k_first, un);
    ierr |= read_global_dataset3(file, "vn", nx, ny, nz,
                                 global_nx, global_ny, global_nz,
                                 local_i_first, local_j_first, local_k_first, vn);
    ierr |= read_global_dataset3(file, "wn", nx, ny, nz,
                                 global_nx, global_ny, global_nz,
                                 local_i_first, local_j_first, local_k_first, wn);
    ierr |= read_global_dataset3(file, "pn", nx, ny, nz,
                                 global_nx, global_ny, global_nz,
                                 local_i_first, local_j_first, local_k_first, pn);

    ierr |= H5Fclose(file) < 0;
    return ierr != 0;
}
