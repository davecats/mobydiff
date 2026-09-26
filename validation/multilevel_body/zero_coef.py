#!/usr/bin/env python3
"""Zero the penalization coefficients of a prepared block-table case file.

The uniform-flow gate needs a file that carries the REAL touch-driven 3-level
refinement -- the interfaces under test -- while exerting no force, so that a
constant velocity field is preserved exactly and any deviation is the
interface transfers' doing and nothing else.

moby_prepare already handles the other half: `[blocks] keep_buried = true`
zeroes the buried masks, so no leaf is removed and no closed face appears in
the flow, and the `blocks` table prepare writes is consistent with those zeroed
masks (the solver cross-checks it row by row at read). Only the coefficients
are left, and nothing in the config can ask for zero ones.

Replaces the mobygeom-importing make_uniform_twin.py, which had to rebuild the
leaf table in Python because it zeroed the masks after the fact.

Usage: zero_coef.py <case.h5>
"""
import sys

import h5py
import numpy as np

# Every coefficient dataset a prepared case file can carry: the staggered
# velocity tiles, and the cell-centred column that appears only when passive
# scalars are configured.
COEF = ("coef_blocks", "coef_p_blocks")


def main() -> int:
    if len(sys.argv) != 2:
        print(__doc__)
        return 2
    path = sys.argv[1]
    touched = []
    with h5py.File(path, "r+") as f:
        for name in COEF:
            if name not in f:
                continue
            d = f[name]
            # Written in place: the shape and the leaf ordering are prepare's
            # and must not change -- the solver matches them against the
            # blocks table.
            d[...] = np.zeros(d.shape, dtype=d.dtype)
            touched.append(f"{name}{d.shape}")
        if not touched:
            print(f"{path}: no coefficient dataset found -- not a prepared "
                  f"block-table case file?", file=sys.stderr)
            return 1
    print(f"   zeroed {', '.join(touched)} in {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
