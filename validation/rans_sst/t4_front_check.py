#!/usr/bin/env python3
"""T4 STEP-0 evidence (docs/next_session_iddes.md; the first-order-upwind
deviation comment in rans.f90): on the channel gates — the only transition
cases the solver can host, since the bc machinery has no inflow/outflow —
the gamma front is WALL-NORMAL and the mean wall-normal velocity is ~0, so
first-order upwind convection cannot smear it. This script quantifies that
on a gate field: the upwind numerical diffusivity across the front,
D_num = |v| dy/2, against the front's physical diffusivity
D_phys = nu + <nut> (sigma_f = 1). PASS if D_num/D_phys < 1e-2 at every
front row (in practice it is orders of magnitude smaller).

Usage: t4_front_check.py <field.h5>
"""
import sys

import h5py
import numpy as np

from rans_channel_check import load_raw, profiles_single_level


def main() -> int:
    path = sys.argv[1]
    blocks, data, y_nodes, re, _ly, _ny = load_raw(path)
    if data["gam"] is None:
        raise SystemExit("no gamma dataset — not a transition run")
    with h5py.File(path, "r") as f:
        data["absv"] = np.abs(f["vn"][...])
    yc, prof = profiles_single_level(blocks, data, y_nodes)
    with h5py.File(path, "r") as f:
        vmax_global = np.abs(f["vn"][...]).max()
    dy = np.diff(y_nodes)

    g = prof["gam"]
    lo, hi = g.min(), g.max()
    if hi - lo < 1e-3:
        print(f"gamma is uniform ({lo:.3f}..{hi:.3f}) — no front to smear; PASS")
        return 0
    front = (g > lo + 0.1*(hi - lo)) & (g < lo + 0.9*(hi - lo))
    if not front.any():
        print("gamma front narrower than one row — treated as sharp; PASS")
        return 0

    nu = 1.0/re
    nut_row = prof["nut"] if prof["nut"] is not None else np.zeros_like(g)
    # Conservative bound: the GLOBAL max |v| (not the row mean) drives the
    # cross-front upwind diffusivity.
    dnum = vmax_global*dy/2.0
    dphys = nu + nut_row                  # gamma diffusivity (sigma_f = 1)
    ratio = dnum[front]/dphys[front]

    print(f"gamma profile: min {lo:.3f} max {hi:.3f}; front spans "
          f"{front.sum()} rows, y in [{yc[front].min():.3f}, {yc[front].max():.3f}]")
    print(f"max |v| over field = {vmax_global:.3e}; at the front rows: "
          f"mean|v| max = {prof['absv'][front].max():.3e}")
    print(f"upwind numerical diffusivity across the front D_num = max|v| dy/2: "
          f"max {dnum[front].max():.3e}")
    print(f"physical gamma diffusivity there  D_phys = nu + nut: "
          f"min {dphys[front].min():.3e}")
    print(f"D_num/D_phys at the front: max {ratio.max():.3e}  (tol 1e-2)")
    ok = ratio.max() < 1e-2
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
