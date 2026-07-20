#!/usr/bin/env python3
"""Plot the u, v, p fields of a boundary-layer snapshot over the x-y plane.

  plot_fields.py <field.h5> [--out blasius2d.png] [--ymax 40]

Fields are z-averaged (quasi-2D) and drawn as pcolormesh on the staggered
node coordinates. y is clipped to --ymax by default so the thin boundary
layer is visible against the 100-theta-tall domain.
"""
import argparse
import os
import sys

import h5py
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "tools"))
from compare_fields import load_field  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("h5")
    ap.add_argument("--out", default="blasius2d.png")
    ap.add_argument("--ymax", type=float, default=40.0)
    a = ap.parse_args()

    with h5py.File(a.h5, "r") as f:
        u = load_field(f, "un").mean(axis=0)   # (ny, nx), z-averaged
        v = load_field(f, "vn").mean(axis=0)
        p = load_field(f, "pn").mean(axis=0)
        xn = f["x"][...]
        yn = f["y"][...]

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    jmax = np.searchsorted(yn, a.ymax)
    fields = [("u", u, "RdBu_r"), ("v", v, "RdBu_r"), ("p", p, "RdBu_r")]
    fig, ax = plt.subplots(3, 1, figsize=(11, 8), sharex=True)
    for k, (name, fld, cmap) in enumerate(fields):
        fk = fld[:jmax, :]
        lim = np.max(np.abs(fk))
        if name == "u":            # u is one-signed; use its own range
            vmin, vmax = fk.min(), fk.max()
        else:
            vmin, vmax = -lim, lim
        pc = ax[k].pcolormesh(xn, yn[:jmax + 1], fk, cmap=cmap,
                              vmin=vmin, vmax=vmax, shading="auto")
        fig.colorbar(pc, ax=ax[k], pad=0.01, label=name)
        ax[k].set_ylabel("y")
        ax[k].set_ylim(0, a.ymax)
    ax[-1].set_xlabel("x")
    ax[0].set_title("Blasius boundary layer: u, v, p (z-averaged)")
    fig.tight_layout()
    fig.savefig(a.out, dpi=150)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
