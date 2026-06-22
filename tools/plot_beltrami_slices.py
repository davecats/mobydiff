#!/usr/bin/env python3
"""Cross-section figures of the Beltrami flow over time: single-block vs refined.

Reads the multi-timestep field output of two Beltrami runs (a uniform single-
block run and a 2:1 block-refined run) and draws a z = Lz/2 cross-section of each
solver variable (u, v, w, p) separately, as a grid of rows = {single block, 2:1
refined} and columns = time. One figure per variable. The block structure is
overlaid (thin grey = every block; bold lines = the coarse/fine interfaces, i.e.
the boundary of the refined region), so any interface artifact is localised. A
diverging colour map centred at zero is shared per variable and fixed to the
single-block (physical) scale, so a refined-case blow-up saturates conspicuously.

Usage:
  python3 tools/plot_beltrami_slices.py UNIFORM_DIR REFINED_DIR [--out-prefix beltrami] [--ntimes 5]
"""

from __future__ import annotations

import argparse
import glob
import os

import h5py
import numpy as np

VARS = [("un", "u"), ("vn", "v"), ("wn", "w"), ("pn", "p")]


def snapshots(d):
    fs = glob.glob(os.path.join(d, "beltrami_*.h5"))
    return sorted(fs, key=lambda f: int(f.rsplit("_", 1)[1].split(".")[0]))


def blocks_on_slice(h5, dset, ztarget):
    """For blocks crossing z=ztarget: (x0, x1, y0, y1, lev, slab[j,i]) of dset."""
    lx = float(h5.attrs["lx"]); nx = int(h5.attrs["nx"])
    nbx = int(h5.attrs["block_nb_x"]); nby = int(h5.attrs["block_nb_y"]); nbz = int(h5.attrs["block_nb_z"])
    blocks = h5["blocks"][...]
    D = h5[dset]
    out = []
    for bid, (ox, oy, oz, lev) in enumerate(blocks):
        h = lx/(nx*2**lev)
        z0 = oz*h; z1 = (oz + nbz)*h
        if not (z0 <= ztarget < z1):
            continue
        k = min(max(int((ztarget - z0)/h), 0), nbz - 1)
        out.append((ox*h, (ox + nbx)*h, oy*h, (oy + nby)*h, int(lev), D[bid][k]))
    return out


def draw_blocks(ax, patches):
    """Thin outline for every block; bold outline around the refined (level>0) region."""
    from matplotlib.patches import Rectangle
    fine = []
    for x0, x1, y0, y1, lev, _ in patches:
        ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False,
                               edgecolor="0.4", lw=0.25, alpha=0.5))
        if lev > 0:
            fine.append((x0, x1, y0, y1))
    if fine:
        fx0 = min(f[0] for f in fine); fx1 = max(f[1] for f in fine)
        fy0 = min(f[2] for f in fine); fy1 = max(f[3] for f in fine)
        ax.add_patch(Rectangle((fx0, fy0), fx1 - fx0, fy1 - fy0, fill=False,
                               edgecolor="k", lw=1.6))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("uniform_dir"); p.add_argument("refined_dir")
    p.add_argument("--out-prefix", default="beltrami")
    p.add_argument("--ntimes", type=int, default=5)
    args = p.parse_args()

    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows = [("single block", snapshots(args.uniform_dir)),
            ("2:1 refined", snapshots(args.refined_dir))]
    for label, fs in rows:
        if not fs:
            raise SystemExit(f"no beltrami_*.h5 in {label} dir")
    n = min(args.ntimes, *(len(fs) for _, fs in rows))
    idx = np.linspace(0, min(len(fs) for _, fs in rows) - 1, n).round().astype(int)

    with h5py.File(rows[0][1][idx[0]], "r") as h:
        lx = float(h.attrs["lx"]); zt = lx/2.0

    for dset, vname in VARS:
        # symmetric colour scale fixed to the single-block (physical) field range
        vmax = 0.0
        for ti in idx:
            with h5py.File(rows[0][1][ti], "r") as h:
                for *_, slab in blocks_on_slice(h, dset, zt):
                    vmax = max(vmax, np.abs(slab).max())
        vmax = vmax or 1.0

        fig, axes = plt.subplots(2, n, figsize=(3.0*n, 6.2), squeeze=False, constrained_layout=True)
        im = None
        for r, (label, fs) in enumerate(rows):
            for c, ti in enumerate(idx):
                ax = axes[r][c]
                with h5py.File(fs[ti], "r") as h:
                    t = float(h.attrs["t_current"])
                    patches = blocks_on_slice(h, dset, zt)
                    for x0, x1, y0, y1, lev, slab in patches:
                        im = ax.imshow(slab, origin="lower", extent=[x0, x1, y0, y1],
                                       vmin=-vmax, vmax=vmax, cmap="RdBu_r", aspect="equal")
                    draw_blocks(ax, patches)
                ax.set_xlim(0, lx); ax.set_ylim(0, lx); ax.set_xticks([]); ax.set_yticks([])
                if r == 0:
                    ax.set_title(f"t={t:.0f}", fontsize=10)
                if c == 0:
                    ax.set_ylabel(label, fontsize=11)
        fig.colorbar(im, ax=axes, shrink=0.8, label=f"{vname}  at  z = Lz/2  (scale ±{vmax:.2f})")
        fig.suptitle(f"Beltrami: {vname} cross-section (z = Lz/2) over time "
                     f"-- bold line = 2:1 interface", fontsize=12)
        out = f"{args.out_prefix}_{vname}.png"
        fig.savefig(out, dpi=140)
        plt.close(fig)
        print(f"wrote {out}  (±{vmax:.3f})")


if __name__ == "__main__":
    main()
