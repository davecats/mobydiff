#!/usr/bin/env python3
"""Reduce the campaign's multi-GB outputs to the small arrays the figures use.

The run itself produces ~470 GB (snapshots of 5.6 GB each, solver statistics
files of 52 MB). Every figure in `asset/figures/` is built from a few hundred
kB of that: one z-plane of one snapshot, one z-averaged cross-section per
grid, and the block table. This script extracts exactly those, so
`plot_pipe.py` runs from the shipped assets and the large files can be deleted
after a campaign -- rerun it after re-running the case to refresh them.

    ./make_caches.py            # from the directory holding the run outputs

Writes: grid_residual.npz (Fig. 6b,c), pipe_prod_plane.npz (Fig. 7),
        pipe_prod_blocks.npz (Fig. 8).
"""

from __future__ import annotations

import os
import sys

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RUN = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..")
NAMES = ["c0", "c1", "c2", "c3", "c4", "mbc", "isof"]

# (tag, solver statistics file, nx) -- the resolution study of Fig. 6
RESID = [("coarse", "p_coarse2_stats.h5", 160),
         ("fine", "p_fine_stats.h5", 256),
         ("production", "p_pr_stats.h5", 256)]
SNAP = "p_pr_stat_198000.h5"          # the plane drawn in Fig. 7
CASE = "pipe_prod.h5"


def path(name):
    return os.path.join(RUN, name)


def grid_residual():
    """The z-averaged c0 mean over the cross-section, per grid.

    Fig. 6 measures the azimuthal spectrum of the residual about the
    axisymmetric mean; all it needs is this 2D field and its coordinates.
    """
    out = {}
    for tag, f, nx in RESID:
        if not os.path.exists(path(f)):
            print(f"   {tag}: {f} absent, skipped")
            continue
        with h5py.File(path(f), "r") as h:
            prof = h["profile"][...]
            x, y = h["xcoord"][...], h["ycoord"][...]
        sv = prof[:, 0].reshape(len(x), len(y))
        out[f"{tag}/sv"] = sv.astype(np.float64)
        out[f"{tag}/x"] = x.astype(np.float64)
        out[f"{tag}/y"] = y.astype(np.float64)
        out[f"{tag}/nx"] = np.array(nx)
        print(f"   {tag}: {sv.shape}")
    np.savez_compressed(os.path.join(HERE, "grid_residual.npz"), **out)


def plane():
    """One z-plane of one production snapshot: the instantaneous figure.

    Only the blocks that intersect the plane are read -- reassembling every
    field to draw one plane would move ~200x the data the figure needs.
    """
    def load_plane(h5, name, kg):
        blocks = h5["blocks"][...]
        nbx = int(h5.attrs["block_nb_x"]); nby = int(h5.attrs["block_nb_y"])
        nbz = int(h5.attrs["block_nb_z"])
        out = np.full((int(h5.attrs["ny"]), int(h5.attrs["nx"])), np.nan)
        d = h5[name]
        for bid, (ox, oy, oz, _l) in enumerate(blocks):
            if oz <= kg < oz + nbz:
                out[oy:oy + nby, ox:ox + nbx] = d[bid, kg - oz]
        return out

    with h5py.File(path(SNAP), "r") as h:
        node = {d: h[d][...] for d in "xyz"}
        kg = len(node["z"]) // 2
        out = {n: load_plane(h, n, kg).astype(np.float32)
               for n in ("un", "vn", "wn", "pn") + tuple(NAMES)}
        out["t_current"] = np.array(float(h.attrs["t_current"]))
    out["x"], out["y"] = node["x"], node["y"]
    out["z_plane"] = np.array(0.5 * (node["z"][kg] + node["z"][kg + 1]))
    out["source"] = np.array(SNAP)
    np.savez_compressed(os.path.join(HERE, "pipe_prod_plane.npz"), **out)
    print(f"   plane k={kg} of {SNAP}")


def blocks():
    """The leaf table and node lines: the block figure."""
    with h5py.File(path(CASE), "r") as h:
        out = {"blocks": h["blocks"][...],
               "block_nb": np.array(int(h.attrs["block_nb"]))}
        for d in "xyz":
            out[d] = h[f"{d}_nodes"][...]
    np.savez_compressed(os.path.join(HERE, "pipe_prod_blocks.npz"), **out)
    print(f"   {len(out['blocks'])} leaves of {out['block_nb']}^3")


def main():
    print(f"reading run outputs from {RUN}")
    grid_residual()
    plane()
    blocks()
    for f in ("grid_residual.npz", "pipe_prod_plane.npz", "pipe_prod_blocks.npz"):
        print(f"   {f}: {os.path.getsize(os.path.join(HERE, f)) / 1e6:.2f} MB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
