#!/usr/bin/env python3
"""A3 INCREMENT 3 gates: SD7003 at Re 6e4, aoa 4, gamma-Re_thetat ON.

  check_sd7003.py <snapshot.h5> <forces.txt>

- LSB presence: a reversed-flow (u < -0.005) patch adjacent to the
  suction surface;
- transition location x_t/c: TURBULENCE onset — the first chord station
  where the near-wall k band exceeds 100 k_inf (gate 0.5 +- 0.1 at
  tu = 0.1 %). NOTE gamma itself is 1 in the FREESTREAM (Langtry-Menter
  convention), so a gamma threshold on a thin-BL band reads the
  freestream and reports the LE — the first version of this checker did
  exactly that;
- gamma-front chordwise smearing: stations across the k-onset rise
  (k from 10 to 1000 k_inf along the near-wall band), REPORTED in
  level-4 cells (the first-order-upwind measurement that decides the
  separate TVD increment — do not gate on it);
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
    ap.add_argument("--span", choices=("z", "y"), default="z",
                    help="span axis (y = the refine_dims xz orientation: the "
                         "suction side is then above LIFT z = 6.0)")
    ap.add_argument("--lmax", type=int, default=4,
                    help="finest refinement level (smearing reported in its cells)")
    a = ap.parse_args()

    names = ["un", "k", "gamma"]
    out, xc, yc = load_window(a.h5, names, 4.4, 5.75, 5.95, 6.35, span=a.span)
    u, kf, gam = out["un"], out["k"], out["gamma"]
    solid = np.isnan(u) | (np.abs(u) < 1e-25)
    k_inf = 1.5e-6   # tu = 0.1 %, U_inf = 1

    nx = xc.size
    xs, u1s, k_bl, rev_any = [], [], [], []
    for i in range(nx):
        col_solid = np.where(solid[:, i])[0]
        if col_solid.size == 0:
            continue  # ahead of the LE / behind the TE
        j_surf = col_solid.max()  # topmost solid cell = suction surface
        if j_surf + 12 >= yc.size:
            continue
        xs.append(xc[i])
        u1s.append(u[j_surf + 1, i])
        rev_any.append(bool((u[j_surf + 1:j_surf + 9, i] < -0.005).any()))
        k_bl.append(np.nanmax(kf[j_surf + 1:j_surf + 10, i]))
    xs = np.array(xs); k_bl = np.array(k_bl); rev_any = np.array(rev_any)
    u1s = np.array(u1s)
    xoc = xs - 4.5

    # LSB: separation from the first-fluid-cell u sign (Cf proxy)
    sep = u1s < -1e-3
    if sep.any():
        i0 = np.where(sep)[0][0]
        # reattachment: first station past i0 where u1 stays positive
        pos = np.where((xoc > xoc[i0]) & (u1s > 1e-3))[0]
        x_r = xoc[pos[0]] if pos.size else np.nan
        print(f"LSB: separation x_s/c = {xoc[i0]:.3f}, reattachment x_r/c = {x_r:.3f} "
              f"(published x_s ~ 0.22-0.30, x_r ~ 0.65-0.70)")
        lsb = True
    else:
        print("LSB: NO separated region found")
        lsb = False

    # transition = turbulence onset along the near-wall k band
    def cross(mult):
        idx = np.where(k_bl >= mult * k_inf)[0]
        return xoc[idx[0]] if idx.size else np.nan
    x10, x50, x90 = cross(10.0), cross(100.0), cross(1000.0)
    hf = 12.0 / 512 / 2**a.lmax
    width_cells = (x90 - x10) / hf if np.isfinite(x90 - x10) else np.nan
    print(f"k onset: x/c(10 k_inf) = {x10:.3f}, x_t/c(100 k_inf) = {x50:.3f}, "
          f"x/c(1000 k_inf) = {x90:.3f}")
    print(f"transition-front chordwise smearing: {width_cells:.1f} level-{a.lmax} cells "
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
