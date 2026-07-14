#!/usr/bin/env python3
"""A3 INCREMENT 3 gates: SD7003 at Re 6e4, aoa 4, gamma-Re_thetat ON.

  check_sd7003.py <snapshot.h5> <forces.txt>

- LSB presence: a reversed-flow (u < -0.005) patch adjacent to the
  suction surface;
- transition location x_t/c: first chord station where the near-wall
  gamma band crosses 0.5 (gate 0.5 +- 0.1 at tu = 0.1 %);
- gamma-front chordwise smearing: stations from gamma 0.1 to 0.9,
  REPORTED in level-4 cells (the first-order-upwind measurement that
  decides the separate TVD increment — do not gate on it);
- C_L/C_D tail means vs the published gamma-Re_thetat scatter
  (C_L ~ 0.58-0.62, C_D ~ 0.021-0.023; gate at +-15 %).

Geometry assumptions: LE at (4.5, 6.0), chord 1, suction side above
y = 6.0; window resolution = the finest level (slice_field replication).
"""
import argparse
import sys

import numpy as np

sys.path.insert(0, "../naca0012")
from slice_field import load_window  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("h5")
    ap.add_argument("forces")
    ap.add_argument("--tail", type=float, default=0.2)
    a = ap.parse_args()

    names = ["un", "gamma"]
    out, xc, yc = load_window(a.h5, names, 4.4, 5.75, 5.95, 6.35)
    u, gam = out["un"], out["gamma"]
    solid = np.isnan(u) | (np.abs(u) < 1e-25)

    nx = xc.size
    xs, cf_sign, gam_bl, rev_any = [], [], [], []
    for i in range(nx):
        col_solid = np.where(solid[:, i])[0]
        if col_solid.size == 0:
            continue  # ahead of the LE / behind the TE
        j_surf = col_solid.max()  # topmost solid cell = suction surface
        if j_surf + 12 >= yc.size:
            continue
        xs.append(xc[i])
        cf_sign.append(np.sign(u[j_surf + 1, i]))
        rev_any.append(bool((u[j_surf + 1:j_surf + 9, i] < -0.005).any()))
        gam_bl.append(np.nanmean(gam[j_surf + 2:j_surf + 8, i]))
    xs = np.array(xs); gam_bl = np.array(gam_bl); rev_any = np.array(rev_any)
    xoc = xs - 4.5

    # LSB: contiguous reversed-flow region
    if rev_any.any():
        i0, i1 = np.where(rev_any)[0][[0, -1]]
        print(f"LSB: reversed flow from x/c = {xoc[i0]:.3f} to {xoc[i1]:.3f} "
              f"({i1 - i0 + 1} stations)")
        lsb = True
    else:
        print("LSB: NO reversed-flow patch found")
        lsb = False

    # transition location + front width from the near-wall gamma band
    def cross(level):
        idx = np.where(gam_bl >= level)[0]
        return xoc[idx[0]] if idx.size else np.nan
    x10, x50, x90 = cross(0.1), cross(0.5), cross(0.9)
    h4 = 12.0 / 512 / 16
    width_cells = (x90 - x10) / h4 if np.isfinite(x90 - x10) else np.nan
    print(f"gamma front: x/c(0.1) = {x10:.3f}, x_t/c(0.5) = {x50:.3f}, "
          f"x/c(0.9) = {x90:.3f}")
    print(f"gamma-front chordwise smearing: {width_cells:.1f} level-4 cells "
          f"({x90 - x10:.4f} c) [REPORT ONLY — triggers the TVD increment]")

    d = np.loadtxt(a.forces, skiprows=1)
    tl = d[int((1.0 - a.tail) * d.shape[0]):, :]
    cl, cd = float(tl[:, 2].mean()), float(tl[:, 3].mean())
    print(f"C_L = {cl:.4f} (published ~0.58-0.62), C_D = {cd:.4f} (~0.021-0.023)")

    ok = lsb
    if not (0.4 <= x50 <= 0.6):
        print(f"FAIL: x_t/c = {x50:.3f} outside 0.5 +- 0.1"); ok = False
    if not (0.49 <= cl <= 0.72):
        print(f"FAIL: C_L outside the +-15 % band around 0.60"); ok = False
    if not (0.017 <= cd <= 0.027):
        print(f"FAIL: C_D outside the +-15 % band around 0.022"); ok = False
    print("sd7003 gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
