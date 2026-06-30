#!/usr/bin/env python3
"""Cross-section of all variables (u,v,w,p) for a refined-patch field vs a
base-grid control. Both are reassembled onto the finest lattice and drawn with
the TRUE (stretched) node coordinates, so coarse cells render as replicated 2x2
blocks and the refined patch shows per-cell detail -- the 2:1 interface is
visible directly. The patch outline (where a level>0 leaf exists in the plane) is
overlaid.

Usage:
  python3 tools/plot_patch_slice.py PATCH.h5 BASE.h5 --out slice.png \
      [--axis z|x|y] [--index FINE_INDEX]   # default: z through the patch centre
"""
from __future__ import annotations
import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare_fields import load_field, block_geometry  # noqa: E402
import h5py  # noqa: E402

VARS = [("u", "un"), ("v", "vn"), ("w", "wn"), ("p", "pn")]
AX = {"x": 2, "y": 1, "z": 0}      # array axis (z,y,x) for a slice NORMAL to axis
LABEL = {"x": ("y", "z"), "y": ("x", "z"), "z": ("x", "y")}


def subdivide(line, n):
    for _ in range(n):
        fine = np.empty(2 * len(line) - 1)
        fine[0::2] = line
        fine[1::2] = 0.5 * (line[:-1] + line[1:])
        line = fine
    return line


def fine_nodes(h5, lmax):
    return {d: subdivide(h5[d][...], lmax) for d in "xyz"}


def patch_span(h5, nb, lmax, axis):
    """Physical bounding box of refined leaves projected onto the plotted plane,
    as ((a0,a1),(b0,b1)) in the two in-plane axes."""
    blocks = h5["blocks"][...]
    xb, yb, zb = h5["x"][...], h5["y"][...], h5["z"][...]
    nodes = {"x": xb, "y": yb, "z": zb}
    nbm = {"x": nb[0], "y": nb[1], "z": nb[2]}
    a, b = LABEL[axis]
    span = {}
    for d, ax in ((a, 0), (b, 1)):
        f = 2 ** lmax
        cells = nbm[d]
        lo, hi = [], []
        for ox, oy, oz, lev in blocks:
            if lev == 0:
                continue
            o = {"x": ox, "y": oy, "z": oz}[d]   # level-lev cell origin
            ff = 2 ** (lmax - int(lev))
            # convert to finest-cell index range, then to node coords
            i0 = o * ff
            i1 = o * ff + cells * ff
            lo.append(i0); hi.append(i1)
        if lo:
            fn = subdivide(nodes[d], lmax)
            span[ax] = (fn[min(lo)], fn[max(hi)])
    return span.get(0), span.get(1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("patch")
    ap.add_argument("base")
    ap.add_argument("--out", required=True)
    ap.add_argument("--axis", default="z", choices=["x", "y", "z"])
    ap.add_argument("--index", type=int, default=None,
                    help="finest-lattice index of the slice (default: patch centre)")
    a = ap.parse_args()

    hp = h5py.File(a.patch, "r")
    hb = h5py.File(a.base, "r")
    blocks, nb, lmax, shape = block_geometry(hp)   # shape = (nz,ny,nx) finest
    _, _, lmax_b, _ = block_geometry(hb)
    up_b = 2 ** (lmax - lmax_b)                     # upsample base onto patch lattice
    fn = fine_nodes(hp, lmax)
    axis = a.axis
    nrm = AX[axis]

    # slice index: default = centre of the refined patch along the normal axis
    if a.index is not None:
        idx = a.index
    else:
        ref = blocks[blocks[:, 3] > 0]
        comp = {"x": 0, "y": 1, "z": 2}[axis]
        f = 2 ** lmax
        o = ref[:, comp]
        lev = ref[:, 3]
        ff = (2 ** (lmax - lev))
        centres = o * ff + (nb[comp] * ff) // 2
        idx = int(np.median(centres))
    a_lab, b_lab = LABEL[axis]
    an, bn = fn[a_lab], fn[b_lab]           # in-plane node lines
    sa, sb = patch_span(hp, nb, lmax, axis)

    # Add the LES eddy viscosity panel when both files carry it.
    vars_to_plot = list(VARS)
    if "nut" in hp and "nut" in hb:
        vars_to_plot.append(("nut", "nut"))

    fig, axes = plt.subplots(len(vars_to_plot), 2, figsize=(11, 4 * len(vars_to_plot)),
                             constrained_layout=True)
    for r, (name, dset) in enumerate(vars_to_plot):
        fp = load_field(hp, dset)
        fb = load_field(hb, dset)
        if up_b > 1:                               # base onto the patch lattice
            fb = fb.repeat(up_b, 0).repeat(up_b, 1).repeat(up_b, 2)
        sl = [slice(None)] * 3
        sl[nrm] = idx
        # remaining (z,y,x) axes minus the normal one are ordered [b, a]
        # (b = vertical, a = horizontal) for every choice of normal -- exactly
        # what pcolormesh(an, bn, C) expects with C[j(b), i(a)].
        pa = fp[tuple(sl)]
        ba = fb[tuple(sl)]
        if name == "p":            # pressure datum is arbitrary -> centre each panel
            pa = pa - pa.mean()
            ba = ba - ba.mean()
        vmin = min(pa.min(), ba.min())
        vmax = max(pa.max(), ba.max())
        for c, (arr, title) in enumerate(((ba, "without refinement"),
                                          (pa, "with refinement"))):
            ax = axes[r, c]
            lo, hi = vmin, vmax
            if name == "p":        # base p carries a huge null-mode -> scale each panel
                lim = float(np.percentile(np.abs(arr), 99)) or 1.0
                lo, hi = -lim, lim
            m = ax.pcolormesh(an, bn, arr, vmin=lo, vmax=hi,
                              cmap="RdBu_r", shading="flat")
            if sa and sb:
                ax.add_patch(Rectangle((sa[0], sb[0]), sa[1] - sa[0], sb[1] - sb[0],
                             fill=False, ec="k", lw=1.2, ls="--"))
            ax.set_title(f"{name}  ({title})")
            ax.set_xlabel(a_lab); ax.set_ylabel(b_lab)
            ax.set_aspect("equal")
            if c == 1:
                fig.colorbar(m, ax=axes[r, :].tolist(), shrink=0.85)
    fig.suptitle(f"{axis}-normal slice at finest index {idx}  "
                 f"({os.path.basename(a.patch)} vs {os.path.basename(a.base)})")
    fig.savefig(a.out, dpi=110)
    print(f"wrote {a.out}  ({axis}-normal, index {idx})")


if __name__ == "__main__":
    main()
