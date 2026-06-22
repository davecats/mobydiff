#!/usr/bin/env python3
"""z = Lz/2 slice of the discrete divergence (MOBY_RKDIV dump, in the pn slot).

Shows where the 2:1 interface manufactures spurious divergence from the exact
divergence-free field -- which faces/cells (flat face vs edge vs corner) carry
the bogus Poisson RHS the corrector then 'fixes'. Per-block imshow, interface bold.

Usage: python3 tools/plot_div_slice.py UNIFORM_DIR REFINED_DIR [--out div.png]
"""
from __future__ import annotations
import argparse, os
import h5py, numpy as np


def slabs(h5, ztarget):
    lx = float(h5.attrs["lx"]); nx = int(h5.attrs["nx"])
    nbx = int(h5.attrs["block_nb_x"]); nby = int(h5.attrs["block_nb_y"]); nbz = int(h5.attrs["block_nb_z"])
    blocks = h5["blocks"][...]; D = h5["pn"]
    out = []
    for bid, (ox, oy, oz, lev) in enumerate(blocks):
        h = lx/(nx*2**lev); z0 = oz*h; z1 = (oz+nbz)*h
        if not (z0 <= ztarget < z1):
            continue
        k = min(max(int((ztarget-z0)/h), 0), nbz-1)
        out.append((ox*h, (ox+nbx)*h, oy*h, (oy+nby)*h, int(lev), D[bid][k]))
    return out


def draw(ax, patches):
    from matplotlib.patches import Rectangle
    fine = []
    for x0, x1, y0, y1, lev, _ in patches:
        ax.add_patch(Rectangle((x0, y0), x1-x0, y1-y0, fill=False, edgecolor="0.4", lw=0.25, alpha=0.5))
        if lev > 0:
            fine.append((x0, x1, y0, y1))
    if fine:
        ax.add_patch(Rectangle((min(f[0] for f in fine), min(f[2] for f in fine)),
                               max(f[1] for f in fine)-min(f[0] for f in fine),
                               max(f[3] for f in fine)-min(f[2] for f in fine),
                               fill=False, edgecolor="k", lw=1.6))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("uniform_dir"); p.add_argument("refined_dir")
    p.add_argument("--out", default="div_slice.png")
    a = p.parse_args()
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    rows = [("single block", a.uniform_dir), ("2:1 refined", a.refined_dir)]
    with h5py.File(os.path.join(a.uniform_dir, "beltrami_900004.h5")) as h:
        lx = float(h.attrs["lx"]); zt = lx/2.0
    vmax = 0.0
    data = {}
    for r, (_, d) in enumerate(rows):
        with h5py.File(os.path.join(d, "beltrami_900004.h5")) as h:
            data[r] = slabs(h, zt)
        for *_, s in data[r]:
            vmax = max(vmax, np.abs(s).max())
    fig, axes = plt.subplots(1, 2, figsize=(13, 6), constrained_layout=True)
    im = None
    for r, (label, _) in enumerate(rows):
        ax = axes[r]
        for x0, x1, y0, y1, lev, s in data[r]:
            im = ax.imshow(s, origin="lower", extent=[x0, x1, y0, y1], vmin=-vmax, vmax=vmax,
                           cmap="RdBu_r", aspect="equal")
        draw(ax, data[r])
        ax.set_xlim(0, lx); ax.set_ylim(0, lx); ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f"{label}   max|div|={max(np.abs(s).max() for *_,s in data[r]):.2e}")
    fig.colorbar(im, ax=axes, shrink=0.8, label=f"discrete divergence of EXACT field (z=Lz/2, ±{vmax:.2f})")
    fig.suptitle("Spurious divergence the 2:1 projection sees as its Poisson RHS", fontsize=13)
    fig.savefig(a.out, dpi=140)
    print(f"wrote {a.out}  (±{vmax:.3f})")


if __name__ == "__main__":
    main()
