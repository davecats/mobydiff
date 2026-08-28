#!/usr/bin/env python3
"""Seed a snapshot with the EXACT two-material steady profile (C1 gate 1).

The gate's claim is that the piecewise-linear two-material solution is an
exact FIXED POINT of the discrete operator -- which it is only if the cut
face's series resistance, and hence the level-set weight w, is exact. Testing
a fixed point by starting AT it is both sharper and vastly cheaper than
waiting out a transient: any error in w changes k_face and the profile leaves
immediately, at a rate set by the error rather than by the slowest eigenmode.
(The cold-start convergence run in the `converge` group is the companion
statement that the fixed point is also the one the solver reaches.)

    ./seed_slab_ic.py base.h5 out.h5 --y-wall Y --kappa K
"""

from __future__ import annotations

import argparse
import os
import shutil
import sys

import numpy as np
import h5py

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scalar"))
from scalar_tools import BlockGeometry                       # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("base")
    ap.add_argument("out")
    ap.add_argument("--y-wall", type=float, required=True)
    ap.add_argument("--kappa", type=float, required=True)
    ap.add_argument("--contact", type=float, default=0.0)
    ap.add_argument("--name", default="theta")
    a = ap.parse_args()

    shutil.copyfile(a.base, a.out)
    with h5py.File(a.out, "r+") as h5:
        geo = BlockGeometry(h5)
        ly = geo.leng[1]
        q = 1.0 / (a.y_wall / a.kappa + a.contact + (ly - a.y_wall))
        th = h5[a.name]
        for bid in range(geo.n_blocks):
            _, y, _, _ = geo.mesh(bid)
            yy = np.broadcast_to(y, th[bid].shape)
            th[bid] = np.where(yy <= a.y_wall, q * yy / a.kappa,
                               q * (a.y_wall / a.kappa + a.contact + (yy - a.y_wall)))
        # Nothing else moves: the velocity is identically zero in this case.
    print(f"{a.out}: exact profile seeded, kappa_s = {a.kappa:g}, q = {q:.16e}")


if __name__ == "__main__":
    main()
