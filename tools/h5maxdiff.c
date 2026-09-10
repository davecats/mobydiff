/* Max |a-b| per dataset between two HDF5 field files.
 *
 *     h5maxdiff a.h5 b.h5 [dataset ...]
 *
 * Exists because HoreKa has neither h5diff nor h5py, and the project's
 * bit-exactness argument is "max_abs 0 on the fields" -- which needs a number,
 * not a byte comparison (HDF5 object headers carry modification times, so two
 * runs writing identical data still differ byte-wise).
 *
 * Datasets default to the field-file set; name them explicitly to compare
 * others. Absent datasets are skipped with a note rather than failing, so the
 * same call works on no-LES and RANS outputs. Version-independent on purpose:
 * no H5Lvisit / H5Oget_info, whose signatures moved between 1.10 and 1.12.
 *
 * Build:  gcc -O2 -o h5maxdiff h5maxdiff.c -I$HDF5_ROOT/include \
 *              -L$HDF5_ROOT/lib -lhdf5 -Wl,-rpath,$HDF5_ROOT/lib
 * Exit 0 if every dataset present in both matches EXACTLY, 1 otherwise.
 */
#include <hdf5.h>
#include <stdio.h>
#include <stdlib.h>

static const char *DEFAULT_SETS[] = {
    "un", "vn", "wn", "pn", "nut", "k", "omega", "gamma", "rethetat", "fd", NULL
};

static double *slurp(hid_t f, const char *name, hsize_t *n)
{
    hid_t d = H5Dopen2(f, name, H5P_DEFAULT);
    if (d < 0) return NULL;
    hid_t s = H5Dget_space(d);
    hssize_t np = H5Sget_simple_extent_npoints(s);
    double *buf = (np > 0) ? malloc((size_t)np * sizeof(double)) : NULL;
    if (buf && H5Dread(d, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf) < 0) {
        free(buf); buf = NULL;
    }
    *n = (hsize_t)np;
    H5Sclose(s); H5Dclose(d);
    return buf;
}

int main(int argc, char **argv)
{
    if (argc < 3) { fprintf(stderr, "usage: h5maxdiff a.h5 b.h5 [dataset ...]\n"); return 2; }
    H5Eset_auto2(H5E_DEFAULT, NULL, NULL);
    hid_t fa = H5Fopen(argv[1], H5F_ACC_RDONLY, H5P_DEFAULT);
    hid_t fb = H5Fopen(argv[2], H5F_ACC_RDONLY, H5P_DEFAULT);
    if (fa < 0 || fb < 0) { fprintf(stderr, "cannot open input files\n"); return 2; }

    const char **sets = DEFAULT_SETS;
    const char *argsets[64];
    if (argc > 3) {
        int n = 0;
        for (int i = 3; i < argc && n < 63; i++) argsets[n++] = argv[i];
        argsets[n] = NULL;
        sets = argsets;
    }

    printf("comparing %s vs %s\n", argv[1], argv[2]);
    int bad = 0, compared = 0;
    double worst = 0.0;
    for (int i = 0; sets[i]; i++) {
        hsize_t na = 0, nb = 0;
        double *a = slurp(fa, sets[i], &na);
        double *b = slurp(fb, sets[i], &nb);
        if (!a && !b) { free(a); free(b); continue; }          /* absent in both: fine */
        if (!a || !b || na != nb) {
            printf("  %-12s MISSING or SHAPE MISMATCH (%llu vs %llu)\n", sets[i],
                   (unsigned long long)na, (unsigned long long)nb);
            bad = 1; free(a); free(b); continue;
        }
        double m = 0.0;
        hsize_t ndiff = 0;
        for (hsize_t j = 0; j < na; j++) {
            double d = a[j] - b[j];
            if (d != 0.0) { ndiff++; if (d < 0) d = -d; if (d > m) m = d; }
        }
        printf("  %-12s n=%-12llu max_abs=%.17g%s\n", sets[i],
               (unsigned long long)na, m, ndiff ? "   ** DIFFERS **" : "");
        if (m > worst) worst = m;
        if (ndiff) bad = 1;
        compared++;
        free(a); free(b);
    }
    if (!compared && !bad) { printf("FAIL: no common datasets compared\n"); bad = 1; }
    printf("%s worst max_abs = %.17g over %d dataset(s)\n",
           bad ? "FAIL:" : "OK:", worst, compared);
    H5Fclose(fa); H5Fclose(fb);
    return bad;
}
