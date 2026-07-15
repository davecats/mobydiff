#!/usr/bin/env python3
"""LE-fan strip metric of the R1 tables (validation/naca0012 README):
rms of the y-high-passed u along vertical strips in the LE band
y in [5.95, 6.05], at FIXED physical distances upstream of the nose
(x = 4.5 - d, d = 0.006/0.012/0.023/0.047/0.094 c), on a slice_field.py
npz (finest-lattice painted window). High-pass = subtract the 5-pixel
moving average along y.

  fan_metric.py slice.npz [slice2.npz ...] [--nose 4.5]
"""
import argparse

import numpy as np

DISTS = (0.006, 0.012, 0.023, 0.047, 0.094)


def strips(path, nose, band=(5.95, 6.05), hp=5):
    d = np.load(path)
    xc, yc, u = d["xc"], d["yc"], d["un"]
    jb = (yc > band[0]) & (yc < band[1])
    out = []
    for dist in DISTS:
        i = int(np.argmin(np.abs(xc - (nose - dist))))
        col = u[jb, i].astype(np.float64)
        ker = np.ones(hp)/hp
        low = np.convolve(col, ker, mode="same")
        e = (col - low)[hp//2:-(hp//2)]
        out.append(float(np.sqrt(np.nanmean(e*e))))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("npz", nargs="+")
    ap.add_argument("--nose", type=float, default=4.5)
    a = ap.parse_args()
    print("case                     " + "".join(f"{d:>9.3f}c" for d in DISTS))
    for p in a.npz:
        s = strips(p, a.nose)
        print(f"{p:24s} " + "".join(f"{v:9.4f} " for v in s))


if __name__ == "__main__":
    main()
