#!/usr/bin/env python3
"""Seed a snapshot with a C2 manufactured solution (the C1 seed_slab_ic idiom).

    ./seed_manufactured.py plane    base.h5 out.h5 --theta 30 --kappa 10 ...
    ./seed_manufactured.py cylinder base.h5 out.h5 --radius 0.2 --kappa 10 ...

Starting AT the exact solution is both sharper and vastly cheaper than
waiting out a transient: the field leaves at a rate set by the scheme's own
error rather than by the slowest eigenmode, and the indicator printed at
step 0 is then a measurement of the SCHEME on a prescribed field.
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
from check_oblique import Plane                              # noqa: E402
from check_cylinder import Shell                             # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("plane")
    p.add_argument("base")
    p.add_argument("out")
    p.add_argument("--theta", type=float, required=True)
    p.add_argument("--kappa", type=float, required=True)
    p.add_argument("--x0", type=float, default=0.5)
    p.add_argument("--y0", type=float, default=0.5)
    p.add_argument("--q-n", type=float, default=1.0, dest="q_n")
    p.add_argument("--amp", type=float, default=0.0)
    p.add_argument("--name", default="theta")

    q = sub.add_parser("cylinder")
    q.add_argument("base")
    q.add_argument("out")
    q.add_argument("--radius", type=float, required=True)
    q.add_argument("--kappa", type=float, required=True)
    q.add_argument("--centre", type=float, nargs=2, default=(0.5, 0.5))
    q.add_argument("--source", type=float, default=1.0)
    q.add_argument("--name", default="theta")

    a = ap.parse_args()
    if a.cmd == "plane":
        field = Plane(a.theta, a.x0, a.y0, a.kappa, a.q_n, a.amp)
    else:
        field = Shell(a.centre[0], a.centre[1], a.radius, a.kappa, a.source)

    shutil.copyfile(a.base, a.out)
    with h5py.File(a.out, "r+") as h5:
        geo = BlockGeometry(h5)
        th = h5[a.name]
        for bid in range(geo.n_blocks):
            x, y, z, _ = geo.mesh(bid)
            th[bid] = field.temperature(*np.broadcast_arrays(x, y, z))
    print(f"{a.out}: {a.cmd} manufactured solution seeded (kappa_s = {a.kappa:g})")


if __name__ == "__main__":
    main()
