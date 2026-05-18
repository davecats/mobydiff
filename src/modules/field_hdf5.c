#include <hdf5.h>
#include <stdlib.h>

static size_t linear_fortran(size_t i, size_t j, size_t k, size_t ni, size_t nj)
{
    return i + ni*(j + nj*k);
}

static size_t linear_hdf5(size_t i, size_t j, size_t k, size_t nj, size_t nk)
{
    return (i*nj + j)*nk + k;
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

static int write_dataset3(hid_t file, const char *name, int nx, int ny, int nz, const double *field)
{
    const size_t ni = (size_t)nx + 2;
    const size_t nj = (size_t)ny + 2;
    const size_t nk = (size_t)nz + 2;
    const size_t n = ni*nj*nk;
    hsize_t dims[3] = {(hsize_t)ni, (hsize_t)nj, (hsize_t)nk};
    hid_t space = -1;
    hid_t dset = -1;
    double *buffer = NULL;
    herr_t status;

    buffer = (double *)malloc(n*sizeof(double));
    if (buffer == NULL) return 1;

    for (size_t k = 0; k < nk; ++k) {
        for (size_t j = 0; j < nj; ++j) {
            for (size_t i = 0; i < ni; ++i) {
                buffer[linear_hdf5(i, j, k, nj, nk)] = field[linear_fortran(i, j, k, ni, nj)];
            }
        }
    }

    space = H5Screate_simple(3, dims, NULL);
    if (space < 0) {
        free(buffer);
        return 1;
    }
    dset = H5Dcreate2(file, name, H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
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

static int read_dataset3(hid_t file, const char *name, int nx, int ny, int nz, double *field)
{
    const size_t ni = (size_t)nx + 2;
    const size_t nj = (size_t)ny + 2;
    const size_t nk = (size_t)nz + 2;
    const size_t n = ni*nj*nk;
    hsize_t dims[3] = {0, 0, 0};
    hid_t dset = -1;
    hid_t space = -1;
    double *buffer = NULL;
    herr_t status;

    dset = H5Dopen2(file, name, H5P_DEFAULT);
    if (dset < 0) return 1;

    space = H5Dget_space(dset);
    if (space < 0) {
        H5Dclose(dset);
        return 1;
    }

    if (H5Sget_simple_extent_ndims(space) != 3) {
        H5Sclose(space);
        H5Dclose(dset);
        return 1;
    }
    H5Sget_simple_extent_dims(space, dims, NULL);
    if (dims[0] != (hsize_t)ni || dims[1] != (hsize_t)nj || dims[2] != (hsize_t)nk) {
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
        for (size_t k = 0; k < nk; ++k) {
            for (size_t j = 0; j < nj; ++j) {
                for (size_t i = 0; i < ni; ++i) {
                    field[linear_fortran(i, j, k, ni, nj)] = buffer[linear_hdf5(i, j, k, nj, nk)];
                }
            }
        }
    }

    free(buffer);
    H5Sclose(space);
    H5Dclose(dset);
    return status < 0;
}

int fdm_h5_write_field(const char *filename, int nx, int ny, int nz,
                       double lx, double ly, double lz,
                       double dx, double dy, double dz,
                       double re, double dt, double t_current,
                       const double *un, const double *vn,
                       const double *wn, const double *pn)
{
    hid_t file;
    int ierr = 0;

    if (nx < 1 || ny < 1 || nz < 1) return 1;

    file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file < 0) return 1;

    ierr |= write_attr_int(file, "nx", nx);
    ierr |= write_attr_int(file, "ny", ny);
    ierr |= write_attr_int(file, "nz", nz);
    ierr |= write_attr_double(file, "lx", lx);
    ierr |= write_attr_double(file, "ly", ly);
    ierr |= write_attr_double(file, "lz", lz);
    ierr |= write_attr_double(file, "dx", dx);
    ierr |= write_attr_double(file, "dy", dy);
    ierr |= write_attr_double(file, "dz", dz);
    ierr |= write_attr_double(file, "re", re);
    ierr |= write_attr_double(file, "dt", dt);
    ierr |= write_attr_double(file, "t_current", t_current);

    ierr |= write_dataset3(file, "un", nx, ny, nz, un);
    ierr |= write_dataset3(file, "vn", nx, ny, nz, vn);
    ierr |= write_dataset3(file, "wn", nx, ny, nz, wn);
    ierr |= write_dataset3(file, "pn", nx, ny, nz, pn);

    ierr |= H5Fclose(file) < 0;
    return ierr != 0;
}

int fdm_h5_read_field(const char *filename, int nx, int ny, int nz,
                      double *un, double *vn, double *wn, double *pn)
{
    hid_t file;
    int ierr = 0;

    if (nx < 1 || ny < 1 || nz < 1) return 1;

    file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) return 1;

    ierr |= read_dataset3(file, "un", nx, ny, nz, un);
    ierr |= read_dataset3(file, "vn", nx, ny, nz, vn);
    ierr |= read_dataset3(file, "wn", nx, ny, nz, wn);
    ierr |= read_dataset3(file, "pn", nx, ny, nz, pn);

    ierr |= H5Fclose(file) < 0;
    return ierr != 0;
}
