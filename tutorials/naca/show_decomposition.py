#!/usr/bin/env python3
"""AS-BUILT decomposition report of a block-table coefficient (or field)
file: per-level leaf/cell counts + a block-outline figure of the nest.

  show_decomposition.py ibm_coeff_b11.h5 [--out decomposition_b11.png]
"""
import argparse

import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("h5")
    ap.add_argument("--out", default="decomposition_b11.png")
    ap.add_argument("--nose", type=float, nargs=2, default=[50.0, 48.0])
    a = ap.parse_args()

    with h5py.File(a.h5, "r") as f:
        blocks = f["blocks"][...]
        nb = int(f.attrs.get("block_nb", f.attrs.get("block_nb_x", 8)))
        nx = int(f.attrs["nx"]); ny = int(f.attrs["ny"]); nz = int(f.attrs["nz"])
        lx = float(f.attrs["lx"]); lz = float(f.attrs["lz"])
        mask = np.asarray(f.attrs.get("refine_dims", [1, 1, 1]), dtype=np.int64)
    lev = blocks[:, 3]
    lmax = int(lev.max())
    print(f"{a.h5}: {blocks.shape[0]} leaves, levels 0..{lmax}, "
          f"grid {nx}x{ny}x{nz}, refine_dims mask {mask.tolist()}")
    print(f"{'lvl':>3s} {'Delta [c]':>11s} {'blocks':>8s} {'Mcells':>8s}")
    tot = 0
    for k in range(lmax + 1):
        n = int((lev == k).sum())
        cells = n*nb**3
        tot += cells
        print(f"{k:3d} {lx/(nx*2**k):11.3e} {n:8d} {cells/1e6:8.3f}")
    print(f"TOTAL {blocks.shape[0]} blocks = {tot/1e6:.1f} M cells "
          f"(~{tot/1e6*0.42:.1f} GB device)")

    fig, axs = plt.subplots(1, 3, figsize=(18, 6))
    views = [(0, lx, 0, lz, "full domain"),
             (a.nose[0]-8, a.nose[0]+16, a.nose[1]-9, a.nose[1]+9, "nest"),
             (a.nose[0]-0.7, a.nose[0]+1.9, a.nose[1]-1.0, a.nose[1]+1.0, "body")]
    cmap = plt.cm.viridis
    for ax, (X0, X1, Z0, Z1, ttl) in zip(axs, views):
        patches, cols = [], []
        for ox, oy, oz, k in blocks:
            hx = lx/(nx*2**(int(k)*int(mask[0])))
            hz = lz/(nz*2**(int(k)*int(mask[2])))
            bx, bz, w, h = ox*hx, oz*hz, nb*hx, nb*hz
            if bx > X1 or bx + w < X0 or bz > Z1 or bz + h < Z0:
                continue
            patches.append(Rectangle((bx, bz), w, h))
            cols.append(int(k))
        pc = PatchCollection(patches, fc="none",
                             ec=[cmap(c/lmax) for c in cols], lw=0.35)
        ax.add_collection(pc)
        ax.plot([a.nose[0], a.nose[0]+1], [a.nose[1], a.nose[1]], "r-", lw=1.5)
        ax.set_xlim(X0, X1); ax.set_ylim(Z0, Z1)
        ax.set_aspect("equal"); ax.set_title(f"{ttl} ({len(patches)} blocks drawn)")
        ax.set_xlabel("x/c"); ax.set_ylabel("z/c")
    fig.suptitle(f"as-built decomposition: {blocks.shape[0]} leaves, "
                 f"{tot/1e6:.1f} M cells, levels 0..{lmax} (colour = level)")
    fig.tight_layout()
    fig.savefig(a.out, dpi=130)
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
