#!/usr/bin/env python3
"""Streamwise-wall-normal (x-y) cross-section of a 2:1 wall-band refined channel
block-table field, for one or more runs side by side, all variables.

Reassembles the block-table field onto the finest-level uniform grid (coarse
cells replicated 2x per level jump; values treated cell-centred for viz) and
draws the x-y plane at mid-span (z = Lz/2). One figure: rows = variables
(u, v, w, p), columns = runs. The colour scale per variable is shared across the
runs (fixed to the first run's robust range) so any interface band / blow-up
stands out. The 2:1 interface y-locations (refined-region boundaries) are drawn
as dashed lines.

Usage:
  python3 tools/slice_channel.py OUT.png FIELD1.h5[:LABEL] [FIELD2.h5[:LABEL] ...]
"""
from __future__ import annotations
import sys
import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

VARS = [("un", "u"), ("vn", "v"), ("wn", "w"), ("pn", "p")]


def reassemble_xy(path, var):
    """Return (img[ny_fine, nx_fine], yfine_edges, interface_y) for the mid-z slice."""
    with h5py.File(path, "r") as f:
        nb = int(f.attrs["block_nb_x"])
        nx, ny, nz = int(f.attrs["nx"]), int(f.attrs["ny"]), int(f.attrs["nz"])
        blocks = f["blocks"][...]
        D = f[var][...]
        ybase = f["y"][...]
    maxlev = int(blocks[:, 3].max())
    s = 2 ** maxlev
    nxf, nyf, nzf = nx * s, ny * s, nz * s
    kz = nzf // 2  # mid-span fine index
    img = np.full((nyf, nxf), np.nan)
    refined_jcells = set()  # fine-y indices that belong to a level>0 block
    for bid, (ox, oy, oz, lev) in enumerate(blocks):
        f0 = 2 ** (maxlev - lev)            # fine cells per this-level cell
        z0, z1 = oz * f0, (oz + nb) * f0
        if not (z0 <= kz < z1):
            continue
        kc = (kz - z0) // f0                 # local k(z) cell in this block
        blk = D[bid]                         # h5py axes = [k(z), j(y), i(x)]
        for i in range(nb):                  # i -> x
            xi0 = (ox + i) * f0
            for j in range(nb):              # j -> y
                yj0 = (oy + j) * f0
                img[yj0:yj0 + f0, xi0:xi0 + f0] = blk[kc, j, i]
        if lev > 0:
            for j in range(nb):
                yj0 = (oy + j) * f0
                refined_jcells.update(range(yj0, yj0 + f0))
    # interface fine-y positions = boundary between refined and unrefined rows
    yf = np.linspace(ybase[0], ybase[-1], nyf + 1)
    iface = []
    rj = sorted(refined_jcells)
    if rj:
        edges = [rj[0]] + [rj[k] for k in range(1, len(rj)) if rj[k] != rj[k - 1] + 1] \
                        + [rj[k] + 1 for k in range(len(rj) - 1) if rj[k + 1] != rj[k] + 1] + [rj[-1] + 1]
        for e in sorted(set(edges)):
            iface.append(yf[e])
    return img, yf, iface


def main():
    out = sys.argv[1]
    specs = sys.argv[2:]
    runs = []
    for sp in specs:
        if ":" in sp and not sp[1:3] == ":\\":
            path, label = sp.rsplit(":", 1)
        else:
            path, label = sp, sp
        runs.append((path, label))

    nrun = len(runs)
    fig, axes = plt.subplots(len(VARS), nrun, figsize=(5.2 * nrun, 3.0 * len(VARS)),
                             squeeze=False)
    for vi, (vk, vn) in enumerate(VARS):
        imgs = [reassemble_xy(p, vk) for p, _ in runs]
        ref = imgs[0][0]
        finite = ref[np.isfinite(ref)]
        if vn == "u":
            vmin, vmax = np.nanpercentile(finite, 1), np.nanpercentile(finite, 99)
        else:
            a = np.nanpercentile(np.abs(finite), 99)
            vmin, vmax = -a, a
        cmap = "viridis" if vn == "u" else "RdBu_r"
        for ri, (img, yf, iface) in enumerate(imgs):
            ax = axes[vi][ri]
            im = ax.imshow(img, origin="lower", aspect="auto", cmap=cmap,
                           vmin=vmin, vmax=vmax,
                           extent=[0, img.shape[1], yf[0], yf[-1]])
            for yI in iface:
                ax.axhline(yI, color="k", ls="--", lw=0.6, alpha=0.7)
            if vi == 0:
                ax.set_title(runs[ri][1], fontsize=11)
            if ri == 0:
                ax.set_ylabel(f"{vn}\n y", fontsize=11)
            ax.set_xticks([])
            plt.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    fig.suptitle("Refined channel x-y slice (z = Lz/2), 2:1 wall bands (dashed = interface)",
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    fig.savefig(out, dpi=110)
    print("wrote", out)


if __name__ == "__main__":
    main()
