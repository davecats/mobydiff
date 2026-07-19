#!/usr/bin/env python3
"""Span-averaged lift-velocity (w) field around the C10 airfoil with the
refinement-block outlines overlaid, plus an LE zoom.

  plot_c10_field.py <field.h5> [--out c10_field.png]
"""
import argparse

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from cv_forces import load_plane


def block_outlines(path, x0, x1, z0, z1, lmin=1):
    with h5py.File(path, "r") as f:
        blocks = f["blocks"][...]
        nb = int(f.attrs["block_nb_x"])
        nx = int(f.attrs["nx"]); nz = int(f.attrs["nz"])
        lx = float(f.attrs["lx"]); lz = float(f.attrs["lz"])
    rects = []
    for ox, oy, oz, lev in blocks:
        if lev < lmin:
            continue
        hx = lx/(nx*2**int(lev)); hz = lz/(nz*2**int(lev))
        bx, bz = ox*hx, oz*hz
        if bx > x1 or bx + nb*hx < x0 or bz > z1 or bz + nb*hz < z0:
            continue
        rects.append((bx, bz, nb*hx, nb*hz, int(lev)))
    return rects


def panel(ax, path, x0, x1, z0, z1, vmax):
    xc, zc, u, w, p, nut, hx, hz = load_plane(path, x0, x1, z0, z1)
    im = ax.pcolormesh(xc, zc, w, cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                       shading="nearest")
    for bx, bz, dx, dz, lev in block_outlines(path, x0, x1, z0, z1):
        ax.add_patch(plt.Rectangle((bx, bz), dx, dz, fill=False,
                                   ec="k", lw=0.25, alpha=0.35))
    ax.set_xlim(x0, x1); ax.set_ylim(z0, z1)
    ax.set_aspect("equal"); ax.set_xlabel("x/c"); ax.set_ylabel("z/c")
    return im


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("h5")
    ap.add_argument("--out", default="c10_field.png")
    a = ap.parse_args()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    im1 = panel(ax1, a.h5, 47.0, 56.0, 44.5, 51.5, 0.15)
    ax1.set_title("w (lift velocity), block outlines")
    im2 = panel(ax2, a.h5, 49.85, 50.45, 47.85, 48.25, 0.30)
    ax2.set_title("LE zoom")
    for im, ax in ((im1, ax1), (im2, ax2)):
        fig.colorbar(im, ax=ax, shrink=0.75)
    fig.tight_layout()
    fig.savefig(a.out, dpi=150)
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
