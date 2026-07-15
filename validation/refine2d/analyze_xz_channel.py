#!/usr/bin/env python3
"""R2D-2 turbulent gate: developed-stats comparison of the xz wall-band
channel (runs/xz, coarse 128x64x128 core + x,z-refined wall bands, ONE
shared stretched y line) against the committed uniform128 control (the
core's exact resolution everywhere) and the uniform-256 reference.

Gates:
  1. NO interface band: u'/v'/w' smooth across the band edge y rows
     (jump ratio of each profile at the interface row vs the neighbour
     rows ~ 1; the reflux-band signature of the y-refined campaign was
     1.56 excess / 0.31 kink).
  2. Core (0.7 < y < 1.3, identical resolution): mean fluctuation
     ratios xz/uniform128 ~ 1, mean-U and -<u'v'> close.
  3. Bulk/mean profile sane vs the reference.

Usage: analyze_xz_channel.py [--dir ../channel_interface/developed/runs]
"""
import argparse
import glob
import os

import h5py
import numpy as np

U, V, W, UU, VV, WW, UV = range(7)


def read_stats(base, name):
    p = os.path.join(base, name, "stats", "channel_stats.h5")
    ys, pr = [], []
    for lf in [p] + sorted(glob.glob(p.replace(".h5", "") + "_l*.h5")):
        if not os.path.exists(lf):
            continue
        with h5py.File(lf) as f:
            c = f["coord"][...]; P = f["profile"][...]; n = f["count"][...]
        m = n > 0
        ys.append(c[m]); pr.append(P[m])
    y = np.concatenate(ys); P = np.concatenate(pr, 0)
    o = np.argsort(y)
    y, P = y[o], P[o]
    # The xz case samples the SAME y rows in the core (l0) and band (l1)
    # files near the interface; average duplicate rows (weighted equally:
    # both are full x-z plane averages of their own region -- keep them
    # separate instead: duplicates only occur if a y row spans both
    # levels, which cannot happen for whole-row bands).
    return y, P


def comp(P):
    return (np.sqrt(np.clip(P[:, UU] - P[:, U]**2, 0, None)),
            np.sqrt(np.clip(P[:, VV] - P[:, V]**2, 0, None)),
            np.sqrt(np.clip(P[:, WW] - P[:, W]**2, 0, None)),
            -(P[:, UV] - P[:, U]*P[:, V]))


def jump_ratio(y, f, y_if):
    """Smoothness of profile f across the interface: value at the row
    nearest y_if vs the linear interpolation through its 2nd neighbours
    on each side. 1 = perfectly smooth."""
    i = int(np.argmin(np.abs(y - y_if)))
    i0, i1 = i - 2, i + 2
    fit = f[i0] + (f[i1] - f[i0])*(y[i] - y[i0])/(y[i1] - y[i0])
    return f[i]/fit


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "..", "channel_interface", "developed", "runs"))
    ap.add_argument("--case", default="xz")
    a = ap.parse_args()

    runs = {}
    for name in (a.case, "uniform128", "reference"):
        try:
            runs[name] = read_stats(a.dir, name)
        except (OSError, ValueError):
            print(f"({name} stats not found -- skipped)")
    y_if = (0.6205, 1.3795)

    yx, Px = runs[a.case]
    cx = comp(Px)
    print(f"{a.case}: {yx.size} stat rows, U_bulk = {np.trapezoid(Px[:, U], yx)/2.0:.4f}, "
          f"U_max = {Px[:, U].max():.4f}")

    # Gate 1: interface smoothness of the fluctuation profiles.
    print("\ninterface jump ratios (1 = smooth; y-refined campaign band was "
          "1.56 u' excess / 0.31 v' kink):")
    ok = True
    for lab, f in zip(("u'", "v'", "w'", "-u'v'"), cx):
        r = [jump_ratio(yx, f, yi) for yi in y_if]
        flag = all(abs(v - 1.0) < 0.10 for v in r)
        ok &= flag
        print(f"  {lab:6s} lower {r[0]:.3f}  upper {r[1]:.3f}  "
              f"{'ok' if flag else 'BAND?'}")

    # Gate 2: core ratios vs the uniform128 control (identical resolution).
    if "uniform128" in runs:
        yc, Pc = runs["uniform128"]
        cc = comp(Pc)
        core = (yc > 0.7) & (yc < 1.3)
        print("\ncore (0.7 < y < 1.3) ratios vs uniform128 (same resolution):")
        for lab, fx, fc in zip(("u'", "v'", "w'", "-u'v'"), cx, cc):
            r = float(np.mean(np.interp(yc[core], yx, fx)/fc[core]))
            print(f"  {lab:6s} {r:.4f}")
        du = float(np.mean(np.abs(np.interp(yc[core], yx, Px[:, U]) - Pc[core, U])
                           /Pc[core, U]))
        print(f"  mean-U rel diff {du:.4%}")

    if "reference" in runs:
        yr, Pr = runs["reference"]
        print(f"\nreference (uniform 256): U_bulk = {np.trapezoid(Pr[:, U], yr)/2.0:.4f}, "
              f"U_max = {Pr[:, U].max():.4f}")

    print("\ninterface-smoothness gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
