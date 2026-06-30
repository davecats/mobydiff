#!/usr/bin/env python3
"""2:1 interface-NORMAL velocity asymmetry metric (no LES needed).

The embedded core patch is centred on the channel centreline, so its two y-faces
are mirror-symmetric: physics demands EQUAL v' on the lower and upper face. The
Phase-3c "low-block-owns-face" convention makes the COARSE-OWNS face (the box's
lower y-face: coarse below, fine above) carry an excess v' (the fine cells get
their interface-normal velocity by prolong-injection of the under-resolved coarse
value), while the FINE-OWNS face (upper) is clean. This script measures the
excess directly from a developed-flow snapshot set, reassembled onto the finest
lattice: per y-face it reports the time-avg v'_rms mean and peakedness (max/mean)
over the fine cells adjacent to the face, plus a mid-box interior reference.

A clean (fixed) scheme gives lower == upper (both ~ the interior peakedness).

Usage:
  python3 vface_asym.py 'runs/gate5/patch/stats/channel_field_*.h5'
"""
import os, sys, glob
import numpy as np
import h5py

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "tools"))
from compare_fields import block_geometry  # noqa: E402


def main():
    pattern = sys.argv[1] if len(sys.argv) > 1 else \
        os.path.join(HERE, "runs/gate5/patch/stats/channel_field_*.h5")
    snaps = sorted(glob.glob(pattern))
    if not snaps:
        sys.exit(f"no snapshots matched {pattern}")
    with h5py.File(snaps[0]) as f:
        bl, nb, lmax, shape = block_geometry(f)
    bl = np.array(bl)
    nz, ny, nx = shape

    fine = bl[bl[:, 3] > 0]
    if len(fine) == 0:
        sys.exit("no refined leaves -- this metric is for the refined patch run")
    x0, x1 = fine[:, 0].min(), fine[:, 0].max() + nb[0]
    y0, y1 = fine[:, 1].min(), fine[:, 1].max() + nb[1]
    z0, z1 = fine[:, 2].min(), fine[:, 2].max() + nb[2]
    print(f"box finest cells x[{x0}:{x1}] y[{y0}:{y1}] z[{z0}:{z1}]  ({len(snaps)} snaps)")

    s1 = np.zeros(shape); s2 = np.zeros(shape); n = 0
    for s in snaps:
        with h5py.File(s) as f:
            D = f["vn"][...]
        G = np.full(shape, np.nan)
        for bid, (ox, oy, oz, l) in enumerate(bl):
            scl = 2 ** (lmax - int(l)); ax0, ay0, az0 = ox * scl, oy * scl, oz * scl
            ex = np.repeat(np.repeat(np.repeat(D[bid], scl, 0), scl, 1), scl, 2)
            G[az0:az0 + nb[2] * scl, ay0:ay0 + nb[1] * scl, ax0:ax0 + nb[0] * scl] = ex
        s1 += G; s2 += G * G; n += 1
    mean = s1 / n
    rms = np.sqrt(np.maximum(s2 / n - mean ** 2, 0.0))

    def stat(yrow, label):
        plane = rms[z0:z1, yrow, x0:x1]
        mx, mn = np.nanmax(plane), np.nanmean(plane)
        print(f"  {label:28s} y={yrow:3d}  v'_rms mean={mn:.4f}  max/mean={mx/mn:.3f}")
        return mn, mx / mn

    print("v' (normal) rms over the box's y-face planes (fine cells):")
    lo_m, lo_r = stat(y0, "LOWER face (coarse-owns)")
    up_m, up_r = stat(y1 - 1, "UPPER face (fine-owns)")
    stat((y0 + y1) // 2, "interior (mid-box)")
    print(f"  asymmetry: lower/upper v'_rms mean ratio = {lo_m / up_m:.3f} "
          f"(1.000 = symmetric/fixed); peakedness {lo_r:.3f} vs {up_r:.3f}")


if __name__ == "__main__":
    main()
