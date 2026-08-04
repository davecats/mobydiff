#!/usr/bin/env python3
"""Write a manufactured initial condition into a copy of a solver snapshot.

The solver has no analytic non-uniform scalar IC (by design: [scalar.N]
initial / init_profile cover the production cases), so the transport gates
seed their fields through the RESTART path -- which also exercises the S0
scalar restart round-trip.

  make_scalar_ic.py in.h5 out.h5 --scalar s1 --wave 1 1 0 --amp 1.0
                                 [--velocity U V W] [--zero-pressure]

--wave gives the integer wavenumbers (k_d = 2 pi n_d / L_d) of
s = amp * sin(k.x); --velocity overwrites u,v,w with a uniform field (the
staggered face values are the same constant, so the discrete divergence is
exactly zero).
"""

from __future__ import annotations

import argparse
import shutil
import sys

import h5py
import numpy as np

from scalar_tools import BlockGeometry


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("infile")
    ap.add_argument("outfile")
    ap.add_argument("--scalar", default="s1")
    ap.add_argument("--wave", type=int, nargs=3, default=[1, 0, 0])
    ap.add_argument("--amp", type=float, default=1.0)
    ap.add_argument("--offset", type=float, default=0.0,
                    help="constant added to the wave (a non-zero mean makes "
                         "int s dV a real quantity for the conservation gates)")
    ap.add_argument("--velocity", type=float, nargs=3, default=None)
    ap.add_argument("--zero-pressure", action="store_true")
    a = ap.parse_args()

    if a.infile != a.outfile:
        shutil.copyfile(a.infile, a.outfile)

    with h5py.File(a.outfile, "r+") as f:
        geo = BlockGeometry(f)
        k = [2.0 * np.pi * n / L for n, L in zip(a.wave, geo.leng)]
        if a.scalar not in f:
            print(f"error: {a.outfile} has no dataset '{a.scalar}'", file=sys.stderr)
            return 1
        dset = f[a.scalar]
        for bid in range(geo.n_blocks):
            x, y, z, _ = geo.mesh(bid)
            dset[bid] = a.offset + a.amp * np.sin(k[0] * x + k[1] * y + k[2] * z)
        if a.velocity is not None:
            for name, val in zip(("un", "vn", "wn"), a.velocity):
                f[name][...] = val
        if a.zero_pressure:
            f["pn"][...] = 0.0

    print(f"wrote {a.outfile}: {a.scalar} = {a.offset} + {a.amp} sin(k.x), n = {a.wave}"
          + (f", u = {a.velocity}" if a.velocity else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
