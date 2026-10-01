/* The commit a binary was built from, for the solver banner.
 *
 * CMakeLists.txt asks git at CONFIGURE time (`git describe --always --dirty --exclude *`,
 * so an uncommitted tree says so) and defines MOBY_COMMIT on THIS file
 * alone: a new commit recompiles these few lines and relinks, instead of
 * rebuilding every module the way a define on init.f90 would. compile.sh
 * re-runs cmake on every build, which keeps it current; a bare
 * `cmake --build` keeps the value of the last configure.
 */
#include <string.h>

#ifndef MOBY_COMMIT
#define MOBY_COMMIT "unknown"
#endif

/* Fill the blank-padded Fortran character buffer buf(1:len). */
void moby_build_commit(char *buf, int len)
{
    size_t n = strlen(MOBY_COMMIT);

    if (len <= 0) return;
    memset(buf, ' ', (size_t)len);
    memcpy(buf, MOBY_COMMIT, n < (size_t)len ? n : (size_t)len);
}
