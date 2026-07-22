#!/usr/bin/env python3
"""Nose-closeup frames of the v1 trip-run instability: u and w around
the LE for every snapshot of a prefix, plus a w contact sheet.

  nose_closeup_series.py [--prefix c11_aoa5tv] [--window 49.92 50.28 47.86 48.14]

Colour ranges are FIXED across frames so growth is visible; each panel
is annotated with max|w| in the window (the late frames saturate the
scale deliberately).
"""
import argparse
import glob
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from cv_forces import load_plane


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prefix", default="c11_aoa5tv")
    ap.add_argument("--window", type=float, nargs=4,
                    default=[49.92, 50.28, 47.86, 48.14])
    a = ap.parse_args()
    x0, x1, z0, z1 = a.window

    files = sorted(glob.glob(f"{a.prefix}_*.h5"),
                   key=lambda p: int(re.search(r"_(\d+)\.h5$", p).group(1)))
    frames = []
    for p in files:
        step = int(re.search(r"_(\d+)\.h5$", p).group(1))
        xc, zc, u, w, pr, nut, hx, hz = load_plane(p, x0, x1, z0, z1)
        frames.append((step, xc, zc, u, w))
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.2))
        im1 = ax1.pcolormesh(xc, zc, u, cmap="RdBu_r", vmin=-0.5, vmax=1.7,
                             shading="nearest")
        im2 = ax2.pcolormesh(xc, zc, w, cmap="RdBu_r", vmin=-0.5, vmax=0.5,
                             shading="nearest")
        mw = np.nanmax(np.abs(w))
        ax1.set_title(f"u   step {step}")
        ax2.set_title(f"w   step {step}   max|w| = {mw:.3g}")
        for im, ax in ((im1, ax1), (im2, ax2)):
            ax.set_aspect("equal"); ax.set_xlabel("x/c"); ax.set_ylabel("z/c")
            fig.colorbar(im, ax=ax, shrink=0.8)
        fig.tight_layout()
        out = f"nose_{a.prefix}_{step}.png"
        fig.savefig(out, dpi=130)
        plt.close(fig)
        print(f"{out}  max|w| = {mw:.4g}")

    n = len(frames)
    if n > 1:
        cols = 3
        rows = (n + cols - 1)//cols
        fig, axes = plt.subplots(rows, cols, figsize=(5.2*cols, 4.0*rows))
        for ax, (step, xc, zc, u, w) in zip(np.ravel(axes), frames):
            ax.pcolormesh(xc, zc, w, cmap="RdBu_r", vmin=-0.5, vmax=0.5,
                          shading="nearest")
            ax.set_aspect("equal")
            ax.set_title(f"step {step}  max|w| {np.nanmax(np.abs(w)):.3g}",
                         fontsize=9)
            ax.set_xticks([]); ax.set_yticks([])
        for ax in np.ravel(axes)[n:]:
            ax.axis("off")
        fig.suptitle(f"{a.prefix}: w nose closeup every 250 steps")
        fig.tight_layout()
        fig.savefig(f"nose_{a.prefix}_sheet.png", dpi=120)
        print(f"nose_{a.prefix}_sheet.png ({n} frames)")


if __name__ == "__main__":
    main()
