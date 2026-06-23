#!/usr/bin/env python3
"""Visualise a Taylor-Green vortex field AND its error vs the exact solution.

For a TGV output (MOBY_TGV runs: uniform / refined / slab_x ...), draws an x-y
cross-section of each variable (u, v, w, p) next to its error (scheme - exact),
evaluated at the staggered location where each variable is stored. Block
outlines are drawn; 2:1 coarse-fine interfaces (where a fine block borders a
coarse one) are highlighted in green so the interface artifact is obvious.

Exact decaying TGV (k = 2*pi/lx, nu = 1/Re, F = exp(-2 nu k^2 t)):
  u = -cos(kx) sin(ky) F   (x-face, y-centre)
  v =  sin(kx) cos(ky) F   (x-centre, y-face)
  w = 0
  p = -1/4 (cos 2kx + cos 2ky) F^2   (cell centre)

Usage:
  python3 tools/plot_tgv_error.py FIELD.h5 [--z Z] [--out out.png]
"""
import argparse
import os
import sys

import h5py
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_field_section import Field, subdivide  # noqa: E402

VARS = [("u", "un"), ("v", "vn"), ("w", "wn"), ("p", "pn")]


def exact(label, xe, ye, k, F):
    """Exact TGV value at the staggered location for `label` over a block whose
    x/y node lines are xe/ye (length nb+1). Returns plane[y, x] (nb, nb)."""
    xc = 0.5 * (xe[:-1] + xe[1:])     # cell centres
    yc = 0.5 * (ye[:-1] + ye[1:])
    xf = xe[:-1]                      # low faces (storage convention: q[i] at low face)
    yf = ye[:-1]
    if label == "u":                 # x-face, y-centre
        X, Y = np.meshgrid(xf, yc)
        return -np.cos(k * X) * np.sin(k * Y) * F
    if label == "v":                 # x-centre, y-face
        X, Y = np.meshgrid(xc, yf)
        return np.sin(k * X) * np.cos(k * Y) * F
    if label == "w":
        X, Y = np.meshgrid(xc, yc)
        return np.zeros_like(X)
    X, Y = np.meshgrid(xc, yc)        # p: cell centre
    return -0.25 * (np.cos(2 * k * X) + np.cos(2 * k * Y)) * F * F


def interface_segments(f, z):
    """Line segments along every 2:1 coarse-fine block boundary in the z-plane."""
    nb = f.nb
    boxes = []
    for bid, (ox, oy, oz, lev) in enumerate(f.blocks):
        zl = f.linesZ[lev]
        if not (zl[oz] <= z < zl[oz + nb]):
            continue
        boxes.append((f.linesX[lev][ox], f.linesX[lev][ox + nb],
                      f.linesY[lev][oy], f.linesY[lev][oy + nb], int(lev)))
    eps = 1e-9
    segs = []
    for ax0, ax1, ay0, ay1, al in boxes:
        for bx0, bx1, by0, by1, bl in boxes:
            if bl == al:
                continue
            if abs(ax1 - bx0) < eps and min(ay1, by1) - max(ay0, by0) > eps:
                segs.append(((ax1, max(ay0, by0)), (ax1, min(ay1, by1))))
            if abs(ay1 - by0) < eps and min(ax1, bx1) - max(ax0, bx0) > eps:
                segs.append(((max(ax0, bx0), ay1), (min(ax1, bx1), ay1)))
    return segs


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("field")
    ap.add_argument("--z", type=float, default=None)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    f = Field(args.field)
    h = f.h5
    lx = float(h.attrs["lx"])
    k = 2 * np.pi / lx
    nu = 1.0 / float(h.attrs["re"])
    t = float(h.attrs["t_current"])
    F = np.exp(-2 * nu * k * k * t)
    z = args.z if args.z is not None else 0.5 * f.lz()
    out = args.out or (os.path.splitext(args.field)[0] + "_tgv.png")
    iface = interface_segments(f, z)

    # which block boundaries are 2:1 interfaces (a level-0 block touching level-1)
    fig, axes = plt.subplots(4, 2, figsize=(13, 18))
    for row, (label, dset) in enumerate(VARS):
        patches = f.patches(dset, z)
        fld_vals, err_vals, exs = [], [], []
        for xe, ye, plane in patches:
            ex = exact(label, xe, ye, k, F)
            exs.append((xe, ye, plane, ex))
            fld_vals.append(plane.ravel())
            err_vals.append((plane - ex).ravel())
        fld_vals = np.concatenate(fld_vals)
        err_vals = np.concatenate(err_vals)
        fa = np.abs(fld_vals - (fld_vals.mean() if label == "p" else 0)).max() or 1.0
        ea = np.abs(err_vals).max() or 1.0
        foff = fld_vals.mean() if label == "p" else 0.0
        for col, (a, off, which) in enumerate(
                [(fa, foff, "field"), (ea, 0.0, "error")]):
            ax = axes[row, col]
            for xe, ye, plane, ex in exs:
                data = plane - off if which == "field" else plane - ex
                ax.pcolormesh(xe, ye, data, cmap="RdBu_r", vmin=-a, vmax=a,
                              shading="flat")
                ax.add_patch(Rectangle((xe[0], ye[0]), xe[-1]-xe[0], ye[-1]-ye[0],
                                       fill=False, ec="0.5", lw=0.3))
            for (p0, p1) in iface:        # highlight 2:1 coarse-fine interfaces
                ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color="lime", lw=1.6)
            ax.set_title(f"{label}  {which}   (|max|={a:.2e})", fontsize=10)
            ax.set_aspect("equal")
            ax.set_xlim(0, lx); ax.set_ylim(0, float(h.attrs["ly"]))
            cb = fig.colorbar(ax.collections[0], ax=ax, fraction=0.046, pad=0.04)
            cb.ax.tick_params(labelsize=7)
    fig.suptitle(f"{os.path.basename(args.field)}   t={t:.3f}  F={F:.3f}  "
                 f"Re={float(h.attrs['re']):.0f}   (left: field, right: error vs exact)",
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.99])
    fig.savefig(out, dpi=110)
    print("wrote", out)
    # quick text summary: max error per variable
    for label, dset in VARS:
        patches = f.patches(dset, z)
        me = max(np.abs(p - exact(label, xe, ye, k, F)).max()
                 for xe, ye, p in patches)
        print(f"  max|{label}-exact| (z-plane) = {me:.3e}")


if __name__ == "__main__":
    main()
