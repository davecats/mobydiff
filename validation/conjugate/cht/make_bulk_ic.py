#!/usr/bin/env python3
"""Seed the BULK-HEATING conjugate channel from an already-developed field.

The velocity does not have to be re-developed: this reads a snapshot of the
antisymmetric campaign on the SAME grid, keeps u, v, w, p untouched, and
rewrites only the six scalars with the steady profile of the bulk-heating
problem. That saves the whole develop phase for the momentum field.

THE SEED IS EXACT WHERE IT MATTERS. With a uniform source S in the fluid and
both outer solid faces at zero, the steady mean balance dJ/dy = S gives

    J_wall = S h                    (all the heat generated in a half-channel
                                     leaves through that wall)
    theta_i = J_wall R_s ,          R_s = d/(kappa_s D_f)

so the SOLID profile -- linear from 0 at the outer face to theta_i at the
interface -- is known in closed form, with no measured fluid resistance in it.
That is the difference from the antisymmetric case, where the solid mean had
to be inferred from the run and the high-capacity walls could never relax to
it: here there is nothing to relax to, because the seed is already right.

The FLUID part uses the Reynolds analogy, which is unusually well founded in
this configuration: the total stress and the total heat flux BOTH fall
linearly to zero at the centreline, so theta+ and U+ satisfy the same
transport problem up to Pr, and with S = 1 and forcing_x = 1 one gets
theta_tau = u_tau = 1 and simply theta - theta_i ~ U(y). It is still only a
seed -- the Pr-dependent sublayer offset is deliberately wrong and the run
finds it in a few time units.

    ./make_bulk_ic.py stat_80000.h5 IC_bulk.h5 [--source 1.0]
"""

from __future__ import annotations

import argparse
import os
import shutil
import sys

import h5py
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "scalar"))
from scalar_tools import BlockGeometry                       # noqa: E402

Y_LO, Y_HI = 0.6, 2.6
RE, PR = 180.0, 0.71
SWEEP = [("k1", 1.0), ("k2", 1.0), ("k3", 1.0),
         ("a1", 0.1), ("a2", 10.0), ("a3", 100.0)]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("infile", help="a developed snapshot on THIS grid")
    ap.add_argument("out")
    ap.add_argument("--source", type=float, default=1.0)
    ap.add_argument("--interfaces", type=float, nargs=2, default=(Y_LO, Y_HI))
    a = ap.parse_args()

    ylo, yhi = a.interfaces
    d = ylo                      # slab thickness (the domain is symmetric)
    h = 0.5 * (yhi - ylo)        # half the fluid gap
    df = 1.0 / (RE * PR)
    j_wall = a.source * h        # exact: dJ/dy = S over a half-channel

    shutil.copyfile(a.infile, a.out)
    with h5py.File(a.out, "r+") as f:
        geo = BlockGeometry(f)
        ynode = f["y"][...]
        yc = 0.5 * (ynode[:-1] + ynode[1:])
        blocks = f["blocks"][...]
        nb = geo.nb[1]
        ny = yc.size

        # the developed field's own mean streamwise velocity, per row
        un = f["un"][...]
        umean = np.zeros(ny)
        cnt = np.zeros(ny)
        for bid, (ox, oy, oz, _) in enumerate(blocks):
            umean[oy:oy + nb] += un[bid].mean(axis=(0, 2))
            cnt[oy:oy + nb] += 1.0
        umean /= np.maximum(cnt, 1.0)

        fl = (yc > ylo) & (yc < yhi)
        lo, hi = yc < ylo, yc > yhi
        print(f"   source S = {a.source:g}   h = {h:g}  ->  J_wall = {j_wall:.6f} "
              f"(exact, the same for every scalar)   theta_tau = {j_wall:.6f}")
        for nm, ks in SWEEP:
            if nm not in f:
                continue
            rs = d / (ks * df)
            ti = j_wall * rs                       # exact interface temperature
            th = np.empty(ny)
            th[lo] = ti * (yc[lo] - (ylo - d)) / d          # linear, 0 -> ti
            th[hi] = ti * ((yhi + d) - yc[hi]) / d          # mirror
            # Reynolds analogy in the fluid, on the run's own mean velocity
            th[fl] = ti + umean[fl] * j_wall
            rows = np.zeros(f[nm].shape)
            for bid, (ox, oy, oz, _) in enumerate(blocks):
                rows[bid] = np.broadcast_to(th[oy:oy + nb][None, :, None],
                                            (nb, nb, nb))
            f[nm][...] = rows
            print(f"   {nm}: kappa_s = {ks:<6g} R_s = {rs:8.2f}  theta_interface"
                  f" = {ti:10.4f}   theta_centre = {th[fl].max():10.4f}")
        f.attrs["step"] = np.int32(0)
        f.attrs["t_current"] = 0.0
    print(f"{a.out}: scalars re-seeded for bulk heating; u, v, w, p untouched")


if __name__ == "__main__":
    main()
