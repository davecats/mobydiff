#!/usr/bin/env python3
"""Global volume-weighted divergence sum Sum(vol * div) from a MOBY_DIVDUMP file.

This is the continuity null-space metric the projection can NEVER fix: the
Jacobi projection preserves the MEAN of div(qs) (the composite operator is
conservative, Sum vol*DG phi = 0), so the global Sum(vol*div(qs)) after the
predictor+sync equals the net 2:1-interface flux mismatch. On a single grid it
telescopes to ~0; a nonzero value is the un-conserved interface flux. It must
stay at ROUND-OFF after any momentum-interface change (Axis 2 of the gate).

vol of a cell at level l (uniform cube) = (lx/(nx*2^l))^3. Pass the divpre file
(D of qs before projection). Reports the signed global sum and the sum of |.|."""
import sys
import h5py
import numpy as np


def main():
    fn = sys.argv[1]
    h5 = h5py.File(fn, "r")
    lx = float(h5.attrs["lx"]); nx = int(h5.attrs["nx"])
    blocks = h5["blocks"][...]
    D = h5["pn"][...]                     # (nblk, nb, nb, nb), the divergence
    h5.close()
    total = 0.0
    absum = 0.0
    for bid, (ox, oy, oz, lev) in enumerate(blocks):
        vol = (lx / (nx * 2 ** lev)) ** 3
        s = np.nansum(D[bid])
        total += vol * s
        absum += vol * np.nansum(np.abs(D[bid]))
    print(f"file={fn}")
    print(f"  Sum(vol*div)     = {total:+.6e}")
    print(f"  Sum(vol*|div|)   = {absum:.6e}")


if __name__ == "__main__":
    main()
