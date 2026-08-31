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

    R_f/2 = (theta_interface - theta_centre) / J      (measured, ONCE)

R_f is a property of the FLOW, identical for every scalar whatever its solid
is, so it is measured from all of them and shared -- not taken per scalar.
That distinction is not cosmetic: the first campaign used each scalar's own
estimate, which for a high-capacity wall is contaminated by its own frozen
transient, and it left C_s = 100 and 10^4 running 12 % and 20 % off their
steady flux for the entire window.

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

Y_LO, Y_HI = 0.6, 2.6
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

        # PASS 1: the fluid resistance. R_f/2 = dtheta_half / J is a property
        # of the FLOW -- the same number for every scalar, whatever its solid
        # is -- so it is measured once and shared. The first campaign got this
        # wrong: it used each scalar's OWN estimate, and for the high-capacity
        # walls that estimate is contaminated by their frozen transient, which
        # left k2 and k3 12 % and 20 % off their steady flux for the whole run.
        # J is taken at the interface FACE with the scheme's own coefficient
        # (w = 1/2 here, so the plain harmonic mean), which is a near-interface
        # quantity and therefore equilibrated.
        means, rfs = {}, []
        for nm, ks in SWEEP:
            if nm not in f:
                continue
            data = f[nm][...]
            mean = np.zeros(ny)
            cnt = np.zeros(ny)
            for bid, (ox, oy, oz, _) in enumerate(blocks):
                mean[oy:oy + nb] += data[bid].mean(axis=(0, 2))
                cnt[oy:oy + nb] += 1.0
            mean /= np.maximum(cnt, 1.0)
            means[nm] = mean

            fl = (yc > Y_LO) & (yc < Y_HI)
            idx = np.flatnonzero(fl)
            i0, i1, mid = idx[0], idx[-1], idx[idx.size // 2]
            dy = yc[i0] - yc[i0 - 1]
            kf = 2.0 * ks / (ks + 1.0) * df          # harmonic mean at w = 1/2
            j_lo = kf * (mean[i0] - mean[i0 - 1]) / dy
            j_hi = kf * (mean[i1 + 1] - mean[i1]) / dy
            jint = 0.5 * abs(j_lo + j_hi)
            dth_half = 0.5 * (abs(mean[i0] - mean[mid]) + abs(mean[i1] - mean[mid]))
            rfs.append(dth_half / max(jint, 1e-30))
        rf_half = float(np.median(rfs))
        print(f"   fluid resistance R_f/2: per-scalar {np.array2string(np.array(rfs), precision=2)}"
              f"  -> median {rf_half:.2f} (spread {100*(max(rfs)/min(rfs)-1):.1f} %)")

        # PASS 2: every solid gets the exact steady linear profile implied by
        # that ONE fluid resistance. Fluctuations are untouched.
        for nm, ks in SWEEP:
            if nm not in f:
                continue
            data = f[nm][...]
            mean = means[nm]
            rs = d / (ks * df)
            j = 2.0 / (2.0 * rs + 2.0 * rf_half)
            ti = a.walls[0] - j * rs

            newmean = mean.copy()
            lo = yc < Y_LO
            hi = yc > Y_HI
            newmean[lo] = a.walls[0] - (a.walls[0] - ti) * (yc[lo] - (Y_LO - d)) / d
            newmean[hi] = a.walls[1] + (a.walls[0] - ti) * ((Y_HI + d) - yc[hi]) / d
            shift = newmean - mean
            for bid, (ox, oy, oz, _) in enumerate(blocks):
                data[bid] += shift[oy:oy + nb][None, :, None]
            f[nm][...] = data
            print(f"   {nm}: kappa_s = {ks:<6g} R_s = {rs:8.2f}  J_steady = {j:.6e}"
                  f"   theta_i {mean[lo][-1]:+.4f} -> {ti:+.4f}"
                  f"   (moved {abs(ti - mean[lo][-1]):.2e})")
    print(f"{a.out}: solid mean re-seeded, fluctuations untouched")


if __name__ == "__main__":
    main()
