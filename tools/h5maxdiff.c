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
 * Block-layout field files hold one ROW per leaf block, in the writer's leaf
 * order, with the leaves listed in the file's own `blocks` table (origin
 * x,y,z + level). Two runs of the same case from case files with different
 * block orders (blocks.f90 leaf_key, the block_key_order record) write the
 * same fields in different row orders, so when both files carry a `blocks`
 * table and the tables differ, the rows of b are MATCHED to those of a on
 * (origin, level) before comparing -- and a table that cannot be matched is
 * a different tiling, reported as such.
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

/* The blocks table (n rows of origin x,y,z + level), or NULL. */
static int *slurp_blocks(hid_t f, hsize_t *n)
{
    hsize_t dims[2] = {0, 0};
    int *rows = NULL;
    hid_t d, s;

    *n = 0;
    if (H5Lexists(f, "blocks", H5P_DEFAULT) <= 0) return NULL;
    d = H5Dopen2(f, "blocks", H5P_DEFAULT);
    if (d < 0) return NULL;
    s = H5Dget_space(d);
    if (H5Sget_simple_extent_ndims(s) == 2) {
        H5Sget_simple_extent_dims(s, dims, NULL);
        if (dims[1] == 4 && dims[0] > 0) {
            rows = malloc((size_t)dims[0] * 4 * sizeof(int));
            if (rows && H5Dread(d, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, rows) < 0) {
                free(rows); rows = NULL;
            }
        }
    }
    if (rows) *n = dims[0];
    H5Sclose(s); H5Dclose(d);
    return rows;
}

typedef struct { int key[4]; hsize_t row; } leaf_t;

static int leaf_cmp(const void *pa, const void *pb)
{
    const leaf_t *a = pa, *b = pb;
    for (int i = 0; i < 4; i++)
        if (a->key[i] != b->key[i]) return a->key[i] < b->key[i] ? -1 : 1;
    return 0;
}

/* perm[i] = the row of b holding the leaf of a's row i; NULL when the two
 * tables are identical row for row (or absent). *bad = 1 when both files
 * have a table and the leaves do not match as SETS. */
static hsize_t *row_match(hid_t fa, hid_t fb, hsize_t *nrows, int *bad)
{
    hsize_t na = 0, nb = 0, i;
    int *ra = slurp_blocks(fa, &na);
    int *rb = slurp_blocks(fb, &nb);
    hsize_t *perm = NULL;
    leaf_t *sorted = NULL;
    int same = 1;

    *nrows = 0;
    if (!ra || !rb) { free(ra); free(rb); return NULL; }
    if (na != nb) {
        printf("  blocks       DIFFERENT TILINGS (%llu vs %llu leaves)\n",
               (unsigned long long)na, (unsigned long long)nb);
        *bad = 1; free(ra); free(rb); return NULL;
    }
    for (i = 0; i < 4 * na && same; i++) same = ra[i] == rb[i];
    if (same) { free(ra); free(rb); return NULL; }

    sorted = malloc((size_t)nb * sizeof(leaf_t));
    perm = malloc((size_t)na * sizeof(hsize_t));
    for (i = 0; i < nb; i++) {
        for (int c = 0; c < 4; c++) sorted[i].key[c] = rb[4 * i + c];
        sorted[i].row = i;
    }
    qsort(sorted, (size_t)nb, sizeof(leaf_t), leaf_cmp);
    for (i = 0; i < na; i++) {
        leaf_t want;
        const leaf_t *hit;
        for (int c = 0; c < 4; c++) want.key[c] = ra[4 * i + c];
        hit = bsearch(&want, sorted, (size_t)nb, sizeof(leaf_t), leaf_cmp);
        if (!hit) {
            printf("  blocks       DIFFERENT TILINGS (a leaf of the first file is not in the second)\n");
            *bad = 1; free(perm); perm = NULL; break;
        }
        perm[i] = hit->row;
    }
    if (perm) {
        printf("  (block row orders differ: rows matched on origin + level)\n");
        *nrows = na;
    }
    free(sorted); free(ra); free(rb);
    return perm;
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
    hsize_t nrows = 0;
    hsize_t *perm = row_match(fa, fb, &nrows, &bad);
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
        /* One row per leaf when the dataset is block-layout (its point
         * count is a multiple of the leaf count -- every field dataset is). */
        hsize_t per = (perm && nrows && na % nrows == 0) ? na / nrows : 0;
        for (hsize_t j = 0; j < na; j++) {
            double d = a[j] - (per ? b[perm[j / per] * per + j % per] : b[j]);
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
    free(perm);
    H5Fclose(fa); H5Fclose(fb);
    return bad;
}
