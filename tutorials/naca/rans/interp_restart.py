#!/usr/bin/env python3
"""Interpolate a block-table snapshot onto the leaf layout of another
(deeper-refined) case file, producing a restart the solver accepts.

  interp_restart.py <source_snapshot.h5> <target_case.h5> <out_restart.h5>

Piecewise-constant injection at the target cell centres: for each target
cell, the containing source leaf is looked up finest-level-first and its
cell value copied. Good enough as an initial condition — the projection
re-establishes discrete divergence-freedom within the first steps, and
the fresh solid classification re-damps penalized cells. refine_dims =
xz assumed (y is one never-refined block layer shared 1:1), so lookup is
2D in (x, z) with the full y pencil copied per hit. Grid, domain and nb
must match between source and target; attrs (t, step, bcs) are copied
from the SOURCE snapshot, the blocks table from the TARGET case.
"""
import sys

import h5py
import numpy as np

FIELDS = ["un", "vn", "wn", "pn", "k", "omega", "nut"]


def main():
    src_path, tgt_path, out_path = sys.argv[1:4]
    with h5py.File(src_path, "r") as f:
        sblocks = f["blocks"][...]
        nb = int(f.attrs["block_nb_x"])
        nx = int(f.attrs["nx"]); nz = int(f.attrs["nz"])
        lx = float(f.attrs["lx"]); lz = float(f.attrs["lz"])
        sdata = {n: f[n][...] for n in FIELDS if n in f}
        sattrs = dict(f.attrs)
        xline = f["x"][...]; yline = f["y"][...]; zline = f["z"][...]
    with h5py.File(tgt_path, "r") as f:
        tblocks = f["blocks"][...]

    # source lookup: (lev, ox, oz) -> bid; y never refined (oy == 0)
    lut = {}
    for bid, (ox, oy, oz, lev) in enumerate(sblocks):
        lut[(int(lev), int(ox), int(oz))] = bid
    slevs = sorted({int(b[3]) for b in sblocks}, reverse=True)

    nt = tblocks.shape[0]
    out = {n: np.empty((nt, nb, nb, nb)) for n in sdata}
    miss = 0
    for tb, (ox, oy, oz, lev) in enumerate(tblocks):
        hx = lx/(nx*2**int(lev)); hz = lz/(nz*2**int(lev))
        xc = (int(ox) + 0.5 + np.arange(nb))*hx
        zc = (int(oz) + 0.5 + np.arange(nb))*hz
        for kz in range(nb):
            for ix in range(nb):
                hit = None
                for sl in slevs:
                    shx = lx/(nx*2**sl); shz = lz/(nz*2**sl)
                    gi = int(xc[ix]/shx); gk = int(zc[kz]/shz)
                    key = (sl, (gi//nb)*nb, (gk//nb)*nb)
                    bid = lut.get(key)
                    if bid is not None:
                        hit = (bid, gk % nb, gi % nb)
                        break
                if hit is None:
                    miss += 1
                    for n in out:
                        out[n][tb, kz, :, ix] = 0.0
                    continue
                bid, kk, ii = hit
                for n in out:
                    out[n][tb, kz, :, ix] = sdata[n][bid, kk, :, ii]
    print(f"{nt} target leaves interpolated, {miss} uncovered cell "
          f"columns (expect 0: the leaf sets cover the same domain)")

    with h5py.File(out_path, "w") as f:
        f.create_dataset("blocks", data=tblocks)
        for n in out:
            f.create_dataset(n, data=out[n])
        f.create_dataset("x", data=xline)
        f.create_dataset("y", data=yline)
        f.create_dataset("z", data=zline)
        for k, v in sattrs.items():
            f.attrs[k] = v
        f.attrs["n_blocks"] = np.int32(nt)
    print(f"wrote {out_path} (t = {sattrs['t_current']:.6f}, "
          f"step = {int(sattrs['step'])})")


if __name__ == "__main__":
    main()
