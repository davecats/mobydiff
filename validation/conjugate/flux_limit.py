#!/usr/bin/env python3
"""The kappa_s -> infinity / kappa_s -> 0 limits, as a RATE (C1 gate 2).

The field comparison against the S3 `dirichlet` / `adiabatic` twins can only
be an O(h) statement -- those modes place their effective boundary on the
staircase while the conjugate interface sits at its true position. The sharp
half of the limit is the interface flux itself:

    kappa_s -> infinity :  q -> q_Dirichlet = 1/(L - y_wall)   like 1/kappa_s
    kappa_s -> 0        :  q -> 0                              like kappa_s

Both are read off the SOLVED field: q = (T(L) - T_wall_top)/(L - y_wall) is
the fluid-side gradient of the converged linear profile, taken between the
top Dirichlet face and the first fluid cell centre, so it needs nothing from
the solver but theta.

    ./flux_limit.py --y-wall Y --high a.h5 b.h5 --high-kappa 1e3 1e5 ...
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import h5py

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scalar"))
from scalar_tools import BlockGeometry                       # noqa: E402


def fluid_flux(path, y_wall, name="theta"):
    """Fluid-side heat flux of the converged linear profile (D_fluid = 1)."""
    with h5py.File(path, "r") as h5:
        geo = BlockGeometry(h5)
        th = h5[name][...]
        ly = geo.leng[1]
        ys, ts = [], []
        for bid in range(geo.n_blocks):
            (_, _), (yc, _), (_, _) = geo.block_axes(bid)
            for j, y in enumerate(yc):
                if y > y_wall:
                    ys.append(y)
                    ts.append(float(th[bid][:, j, :].mean()))
    ys = np.asarray(ys)
    ts = np.asarray(ts)
    o = np.argsort(ys)
    ys, ts = ys[o], ts[o]
    # The top Dirichlet face is at y = ly with theta = 1; use the outermost
    # two fluid rows, which is exact for a linear profile.
    return (ts[-1] - ts[-2]) / (ys[-1] - ys[-2])


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--y-wall", type=float, required=True)
    ap.add_argument("--high", nargs="*", default=[])
    ap.add_argument("--high-kappa", nargs="*", type=float, default=[])
    ap.add_argument("--low", nargs="*", default=[])
    ap.add_argument("--low-kappa", nargs="*", type=float, default=[])
    ap.add_argument("--ly", type=float, default=1.0)
    a = ap.parse_args()

    ok = True
    if a.high:
        q_inf = 1.0 / (a.ly - a.y_wall)
        print(f"      kappa -> inf: q_Dirichlet = {q_inf:.12f}")
        errs = []
        for path, ka in zip(a.high, a.high_kappa):
            q = fluid_flux(path, a.y_wall)
            errs.append(abs(q - q_inf))
            print(f"      kappa_s = {ka:9.1e}   q = {q:.12f}   |q - q_inf| = {errs[-1]:.3e}")
        # Each decade of kappa_s must cut the residual by ~that decade.
        for i in range(1, len(errs)):
            r = errs[i - 1] / max(errs[i], 1e-300)
            k = a.high_kappa[i] / a.high_kappa[i - 1]
            print(f"      residual ratio {r:.3e} for a kappa ratio {k:.0e}"
                  f"  -> order {np.log(r)/np.log(k):.3f}")
            if not 0.7 <= np.log(r) / np.log(k) <= 1.3:
                ok = False
    if a.low:
        print("      kappa -> 0: q must vanish like kappa_s")
        qs = []
        for path, ka in zip(a.low, a.low_kappa):
            q = abs(fluid_flux(path, a.y_wall))
            qs.append(q)
            print(f"      kappa_s = {ka:9.1e}   |q| = {q:.6e}   |q|/kappa_s = {q/ka:.6e}")
        for i in range(1, len(qs)):
            r = qs[i - 1] / max(qs[i], 1e-300)
            k = a.low_kappa[i - 1] / a.low_kappa[i]
            print(f"      flux ratio {r:.3e} for a kappa ratio {k:.0e}"
                  f"  -> order {np.log(r)/np.log(k):.3f}")
            if not 0.7 <= np.log(r) / np.log(k) <= 1.3:
                ok = False
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
