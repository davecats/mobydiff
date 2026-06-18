#!/usr/bin/env python3
"""Plot an x-y cross-section of u, v, w and p from a mobydiff field file.

Handles both layouts: the block-table format (per-block datasets, possibly
2:1 refined) and the legacy global-3D format. Each block is drawn on its own
refinement-level grid, so refined regions show their finer cells and the block
decomposition is visible (thin outlines mark block boundaries).

Values are plotted at the cells in which they are stored (u/v/w are
face-staggered, p is cell-centred); for a qualitative cross-section this is
fine. The plane is taken at the cell row nearest the requested z.

Usage:
  python3 tools/plot_field_section.py FIELD.h5 [--z Z] [--out out.png]
"""

from __future__ import annotations

import argparse
import os

import h5py
import numpy as np

VARS = [("u", "un"), ("v", "vn"), ("w", "wn"), ("p", "pn")]


def subdivide(line, times):
    """Midpoint subdivision `times` times (the solver's per-level node lines)."""
    for _ in range(times):
        fine = np.empty(2*(line.size - 1) + 1)
        fine[0::2] = line
        fine[1::2] = 0.5*(line[:-1] + line[1:])
        line = fine
    return line


class Field:
    """An x-y cross-section reader for one field file."""

    def __init__(self, path):
        self.h5 = h5py.File(path, "r")
        self.x0 = self.h5["x"][...]            # level-0 node lines
        self.y0 = self.h5["y"][...]
        self.z0 = self.h5["z"][...]
        self.block = self.h5["un"].ndim == 4
        if self.block:
            self.blocks = self.h5["blocks"][...]
            self.nb = int(self.h5.attrs["block_nb_x"])
            maxlev = int(self.blocks[:, 3].max())
            self.linesX = {l: subdivide(self.x0, l) for l in range(maxlev + 1)}
            self.linesY = {l: subdivide(self.y0, l) for l in range(maxlev + 1)}
            self.linesZ = {l: subdivide(self.z0, l) for l in range(maxlev + 1)}

    def lz(self):
        return float(self.h5.attrs["lz"])

    def patches(self, dset, z):
        """List of (x_edges, y_edges, plane[y, x]) covering the x-y plane at z."""
        data = self.h5[dset]
        if not self.block:
            zc = 0.5*(self.z0[:-1] + self.z0[1:])
            k = int(np.argmin(np.abs(zc - z)))
            return [(self.x0, self.y0, data[k, :, :])]
        out = []
        nb = self.nb
        for bid, (ox, oy, oz, lev) in enumerate(self.blocks):
            zl = self.linesZ[lev]
            if not (zl[oz] <= z < zl[oz + nb]):
                continue
            kk = int(np.clip(np.searchsorted(zl[oz:oz+nb+1], z) - 1, 0, nb - 1))
            out.append((self.linesX[lev][ox:ox+nb+1],
                        self.linesY[lev][oy:oy+nb+1],
                        data[bid][kk, :, :]))
        return out

    def close(self):
        self.h5.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("field", help="field HDF5 file")
    parser.add_argument("--z", type=float, default=None,
                        help="cross-section height z (default: mid-domain lz/2)")
    parser.add_argument("--out", default=None, help="output PNG (default: <field>_xy.png)")
    parser.add_argument("--cmap", default="RdBu_r")
    parser.add_argument("--edges", action="store_true",
                        help="draw every cell edge (default: only block outlines)")
    args = parser.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    f = Field(args.field)
    z = args.z if args.z is not None else 0.5*f.lz()
    out = args.out or (os.path.splitext(args.field)[0] + "_xy.png")

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    for ax, (label, dset) in zip(axes.ravel(), VARS):
        patches = f.patches(dset, z)
        vals = np.concatenate([p[2].ravel() for p in patches])
        if label == "p":
            vmax = np.abs(vals - vals.mean()).max() or 1.0
            vmin, vmax, off = -vmax, vmax, vals.mean()
        else:
            a = np.abs(vals).max() or 1.0
            vmin, vmax, off = -a, a, 0.0
        mesh = None
        for xe, ye, plane in patches:
            mesh = ax.pcolormesh(xe, ye, plane - off, cmap=args.cmap,
                                 vmin=vmin, vmax=vmax, shading="flat",
                                 edgecolors=("face" if not args.edges else "k"),
                                 linewidth=0.05)
            ax.add_patch(Rectangle((xe[0], ye[0]), xe[-1]-xe[0], ye[-1]-ye[0],
                                   fill=False, ec="0.4", lw=0.3))
        ax.set_title(f"{label}" + (" - mean" if label == "p" else ""))
        ax.set_xlabel("x"); ax.set_ylabel("y")
        if mesh is not None:
            fig.colorbar(mesh, ax=ax, shrink=0.85)
    fig.suptitle(f"{os.path.basename(args.field)}   x-y section at z = {z:.4f}")
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    f.close()
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
