#!/usr/bin/env python3
"""Leg 1 of the three-leg plan: put every solid on its EXACT steady mean.

    ./reseed_solid_exact.py Fstat_640000.h5 IC_clean.h5

WHY NOT A CAPACITY-ACCELERATED RUN. The intent of leg 1 is to reach the solid's
steady state without waiting d^2/alpha_s for it, and lowering C_s does that
because the steady state is capacity-INDEPENDENT (C1 gate 1c proves it
exactly). But here the same end is reached in closed form and in one step: at
steady state the solid profile is LINEAR with slope J/(kappa_s D_f), and the
only free quantity is its level -- which, with an imposed flux at BOTH outer
faces, is a pure-Neumann null-space mode and therefore carries no physics.
So the transient is removed by replacing each solid's MEAN profile with the
exact one anchored at its own fluid-side interface value. Fluctuations are
untouched: subtracting a y-dependent MEAN leaves theta' exactly as it was, and
the operator is linear.

Running leg 1 instead would cost ~76 h and would still leave the kappa_s = 1,
C_s = 10^4 scalar unequilibrated (alpha_s = 10^-4 -> d^2/alpha_s ~ 10^6 time
units). This is why the earlier campaigns' note says the thick-wall solid must
be SEEDED, not waited for.

ALSO CORRECTED HERE: the source constant. `source` was set from the PREVIOUS
run's bulk velocity, so generation exceeded the imposed outflux by 0.30 %,
which drove a slow level drift (predicted -7.5e-4 per time unit, measured
-9.3e-4). The corrected value is written to stdout for the ini.
"""

from __future__ import annotations

import argparse
import os
import shutil
import sys

import h5py
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "scalar"))
from scalar_tools import BlockGeometry                       # noqa: E402

RE, PR = 149.0, 0.71
YLO, YHI = 1.0, 3.0
SWEEP = [("k1", 1.0), ("k2", 1.0), ("k3", 1.0),
         ("a1", 0.1), ("a2", 10.0), ("a3", 100.0)]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("src")
    ap.add_argument("out")
    a = ap.parse_args()

    shutil.copyfile(a.src, a.out)
    df = 1.0 / (RE * PR)
    with h5py.File(a.out, "r+") as f:
        geo = BlockGeometry(f)
        nb = geo.nb[1]
        yn = f["y"][...]
        yc = 0.5 * (yn[:-1] + yn[1:])
        bl = f["blocks"][...]
        ny = yc.size
        fl = (yc > YLO) & (yc < YHI)
        d = YLO                                   # slab thickness

        for nm, ks in SWEEP:
            if nm not in f:
                continue
            # current plane mean, per row
            mean = np.zeros(ny)
            cnt = np.zeros(ny)
            for b, (ox, oy, oz, _) in enumerate(bl):
                mean[oy:oy + nb] += f[nm][b].mean(axis=(0, 2))
                cnt[oy:oy + nb] += 1.0
            mean /= np.maximum(cnt, 1.0)

            # the EXACT steady solid mean: linear, slope J/(kappa_s D_f) with
            # J = 1 by the imposed-flux BC, anchored on the FLUID-side
            # interface value so the interface stays continuous.
            target = mean.copy()
            j_lo = int(np.argmax(fl))              # first fluid row
            j_hi = int(len(fl) - 1 - np.argmax(fl[::-1]))
            slope = 1.0 / (ks * df)
            lo = yc < YLO
            hi = yc > YHI
            target[lo] = mean[j_lo] - (YLO - yc[lo]) * slope
            target[hi] = mean[j_hi] - (yc[hi] - YHI) * slope
            shift = target - mean                  # zero in the fluid
            rows = np.zeros(f[nm].shape)
            src = f[nm][...]
            for b, (ox, oy, oz, _) in enumerate(bl):
                rows[b] = src[b] + shift[oy:oy + nb][None, :, None]
            f[nm][...] = rows
            print(f"   {nm}: kappa_s {ks:<6g}  solid-mean shift "
                  f"{shift[lo].min():+9.4f} .. {shift[lo].max():+9.4f}"
                  f"   (fluctuations untouched)")
        f.attrs["step"] = np.int32(0)
        f.attrs["t_current"] = 0.0

        # the corrected source constant, from THIS field's own mean velocity
        u = np.zeros(ny); c = np.zeros(ny)
        for b, (ox, oy, oz, _) in enumerate(bl):
            u[oy:oy + nb] += f["un"][b].mean(axis=(0, 2)); c[oy:oy + nb] += 1.0
        u /= np.maximum(c, 1.0)
        half = fl & (yc < 0.5 * (YLO + YHI))
        I = float((u[half] * np.diff(yn)[half]).sum())
    print(f"\n{a.out}: written")
    print(f"   int(<u> dy) = {I:.4f}  ->  set [scalar.N] source = {1.0/I:.6f}"
          f"   so generation exactly balances the imposed outflux")


if __name__ == "__main__":
    main()
