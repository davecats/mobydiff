#!/usr/bin/env python3
"""Re-seed the SOLID mean profile from the developed fluid field.

WHY THIS EXISTS. With a wall one half-height thick the solid carries most of
the total thermal resistance (R_s = d/(kappa_s D_f) = 128 against R_f ~ 25 for
kappa_s = 1), and its mean profile relaxes on the SOLID conduction time
d^2/alpha_s. For alpha_s = alpha_f that is t ~ 80; for the high-capacity
scalars of this sweep (C_s = 100, 10^4 at kappa_s = 1) it is t ~ 8e3 and
8e5 -- unreachable, and not worth reaching, because the capacity cannot change
the steady state at all (that is C1's gate 1c, measured to 1e-13).

So the solid mean is not something to wait for, it is something to SET. The
fluid side equilibrates on fluid time scales (t ~ 5), so after the develop
phase the run knows its own fluid resistance

    R_f = (theta_interface - theta_centre) / J        (measured, per side)

and the exact steady state follows in closed form:

    J = 2/(2 R_s + 2 R_f_half),   theta_i = 1 - J R_s,

with a LINEAR profile through each slab. This replaces the solid MEAN by that
profile and leaves the solid FLUCTUATIONS untouched -- they equilibrate on the
penetration-depth time, which is short, and they are the quantity being
measured. The fluid is not touched at all.

This is the standard move for conjugate DNS with a high-capacity wall: Tiselj
et al. note that a source term in the solid changes only the mean and has no
influence on the temperature fluctuations. Here the mean is corrected rather
than driven, which is the same statement.

    ./reseed_solid.py in.h5 out.h5
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

Y_LO, Y_HI = 1.0, 3.0
RE, PR = 180.0, 0.71
SWEEP = [("k1", 1.0), ("k2", 1.0), ("k3", 1.0),
         ("a1", 0.1), ("a2", 10.0), ("a3", 100.0)]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("infile")
    ap.add_argument("out")
    ap.add_argument("--walls", type=float, nargs=2, default=(1.0, -1.0))
    a = ap.parse_args()

    shutil.copyfile(a.infile, a.out)
    d = Y_LO                                   # slab thickness (= ly - Y_HI too)
    df = 1.0 / (RE * PR)
    with h5py.File(a.out, "r+") as f:
        geo = BlockGeometry(f)
        ynode = f["y"][...]
        yc = 0.5 * (ynode[:-1] + ynode[1:])
        blocks = f["blocks"][...]
        nb = geo.nb[1]
        ny = yc.size

        for nm, ks in SWEEP:
            if nm not in f:
                continue
            data = f[nm][...]
            # global x-z averaged mean profile
            mean = np.zeros(ny)
            cnt = np.zeros(ny)
            for bid, (ox, oy, oz, _) in enumerate(blocks):
                mean[oy:oy + nb] += data[bid].mean(axis=(0, 2))
                cnt[oy:oy + nb] += 1.0
            mean /= np.maximum(cnt, 1.0)

            fl = (yc > Y_LO) & (yc < Y_HI)
            yf = yc[fl]
            mf = mean[fl]
            # the fluid's own resistance, from the two near-interface values
            # and the centreline, measured on the developed field
            i0, i1 = 0, mf.size - 1
            mid = mf.size // 2
            dth_half = 0.5 * (abs(mf[i0] - mf[mid]) + abs(mf[i1] - mf[mid]))
            # J from the near-wall gradient (molecular, at the interface)
            grad0 = (mf[i0] - a.walls[0] * 0.0) * 0.0     # placeholder, see below
            # J is measured directly: the wall-normal conduction flux at the
            # first fluid cell face is not stored here, so use the fluid
            # resistance and the series relation instead.
            rf_half = dth_half / max(abs(mf[i0] - mf[mid]), 1e-30)  # = 1, kept explicit
            # R_f_half in flux units: dtheta_half / J. J is unknown, so take
            # the analogy-free route: the fluid resistance is a property of the
            # flow, so measure it as dtheta_half divided by the flux implied by
            # the CURRENT solid profile (which is linear and known).
            sol = yc < Y_LO
            if sol.sum() >= 2:
                ys, ms = yc[sol], mean[sol]
                slope = np.polyfit(ys, ms, 1)[0]           # d<theta>/dy
                jcur = -ks * df * slope                    # conduction flux, +y
            else:
                jcur = np.nan
            rf_half = dth_half / max(abs(jcur), 1e-30)
            rs = d / (ks * df)
            j = 2.0 / (2.0 * rs + 2.0 * rf_half)
            ti = a.walls[0] - j * rs

            # rebuild: linear in each slab, fluctuations preserved
            newmean = mean.copy()
            lo = yc < Y_LO
            hi = yc > Y_HI
            newmean[lo] = a.walls[0] - (a.walls[0] - ti) * (yc[lo] - (Y_LO - d)) / d
            newmean[hi] = a.walls[1] + (a.walls[0] - ti) * ((Y_HI + d) - yc[hi]) / d
            shift = newmean - mean
            for bid, (ox, oy, oz, _) in enumerate(blocks):
                data[bid] += shift[oy:oy + nb][None, :, None]
            f[nm][...] = data
            print(f"   {nm}: kappa_s = {ks:<6g} R_s = {rs:8.2f}  R_f/2 = {rf_half:7.2f}"
                  f"  J_measured = {jcur:.6e} -> J_steady = {j:.6e}"
                  f"   theta_i {mean[lo][-1]:.4f} -> {ti:.4f}")
    print(f"{a.out}: solid mean re-seeded, fluctuations untouched")


if __name__ == "__main__":
    main()
