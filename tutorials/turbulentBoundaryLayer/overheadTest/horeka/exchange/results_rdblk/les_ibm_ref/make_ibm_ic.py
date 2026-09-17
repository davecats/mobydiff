#!/usr/bin/env python3
"""Build the IBM-channel restart IC by mapping the developed KMM180 channel into
the fluid gap of the extended (uniform-y) domain.

The KMM180 restart is a wall-resolved channel of half-height 1 (ly=2.0, walls at
y=0 and y=2). The IBM channel has the SAME fluid gap (width 2.0) but embedded in a
taller domain (ly=2.5) with the walls at y_lo=0.259375 and y_hi=2.259375. So the
mapping is a pure y-shift (+y_lo) of the source field onto the IBM grid, with the
solid regions (y<y_lo, y>y_hi) zeroed (the IBM masks them anyway; the first
projection absorbs the residual divergence, so a transient is discarded).

Layout: the block-table format (un/vn/wn/pn = (n_blocks, NB, NB, NB)).
  * single level (default): copy a solver-minted cold-start field (cs_1.h5 -- carries
    the correct grid/BC/ibm_enabled attrs, blocks table and node lines) and overwrite
    the velocity/pressure block rows.  ->  IC.h5
  * refined (--leaves ibm_coeff_blocks.h5): build the multi-level leaf rows by
    scattering coarse + fine (subdivided) global fields per leaf level, reading the
    leaf table from the block-table coefficient file (the solver cross-checks it).
    Attrs/x/y/z copied from cs_1.h5.  ->  IC_refine.h5

Run with the geometry venv (h5py + numpy):
    /home/davide/ibmc/bin/python make_ibm_ic.py                          # IC.h5
    /home/davide/ibmc/bin/python make_ibm_ic.py --leaves ibm_coeff_blocks.h5 --out IC_refine.h5
"""
import argparse
import os
import shutil
import sys

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "tools"))
import make_channel_restart as mkr  # noqa: E402  (staggered helpers, validated)

SRC = os.path.join(ROOT, "tutorials/channel_kmm180/channel_kmm180_restart.h5")
COLD = os.path.join(HERE, "cs_1.h5")     # solver-minted template (attrs + layout)
OUT = os.path.join(HERE, "IC.h5")

LX, LY, LZ = 12.566370614359172, 2.5, 6.283185307179586
Y_LO = 8.3 * (LY / 80)          # 0.259375  (fluid gap = [Y_LO, Y_LO+2.0])
GAP = 2.0
NB = 8
NAMES = ("un", "vn", "wn", "pn")
VAR = {"un": 0, "vn": 1, "wn": 2, "pn": 3}
PERIODIC = (True, False, True)
LENGTHS = (LX, LY, LZ)


def global_field(src, nodes, var):
    """y-shifted KMM180 staggered field interpolated onto `nodes`, solid zeroed."""
    xs, ys, zs = src["nodes"]
    sp = mkr.staggered_positions((xs, ys, zs), var)
    sp = [sp[0], sp[1] + Y_LO, sp[2]]            # shift source up by Y_LO
    dp = mkr.staggered_positions(nodes, var)
    g = mkr.interp_field(src["fields"][["un", "vn", "wn", "pn"][var]],
                         sp, dp, PERIODIC, LENGTHS)
    ypos = dp[1]
    g[:, (ypos < Y_LO) | (ypos > Y_LO + GAP), :] = 0.0   # zero solid (axis 1 = y)
    return g


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--leaves", default=None,
                    help="block-table coef file -> build the multi-level refined IC")
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args()

    src = mkr.load_source(SRC)
    with h5py.File(COLD, "r") as f:
        xt, yt, zt = f["x"][...], f["y"][...], f["z"][...]
        cold_attrs = dict(f.attrs)
        ref_shape = f["un"].shape
    base_nodes = (xt, yt, zt)

    if a.leaves is None:                          # single level: overwrite cs_1.h5
        with h5py.File(COLD, "r") as f:
            blocks = f["blocks"][...]
        assert (blocks[:, 3] == 0).all()
        glob = {nm: global_field(src, base_nodes, VAR[nm]) for nm in NAMES}
        rows = {nm: np.zeros(ref_shape) for nm in NAMES}
        for bid, (ox, oy, oz, lev) in enumerate(blocks):
            for nm in NAMES:
                rows[nm][bid] = glob[nm][oz:oz + NB, oy:oy + NB, ox:ox + NB]
        shutil.copyfile(COLD, a.out)
        with h5py.File(a.out, "r+") as f:
            for nm in NAMES:
                f[nm][...] = rows[nm]
            f.attrs["step"] = np.int32(0)
            f.attrs["t_current"] = 0.0
    else:                                         # refined: leaves from the coef file
        with h5py.File(a.leaves, "r") as f:
            blocks = f["blocks"][...]             # (n_leaves, 4) origin in level-l cells
        fine_nodes = tuple(mkr.subdivide(n) for n in base_nodes)
        grids = {0: {nm: global_field(src, base_nodes, VAR[nm]) for nm in NAMES},
                 1: {nm: global_field(src, fine_nodes, VAR[nm]) for nm in NAMES}}
        nleaf = len(blocks)
        rows = {nm: np.zeros((nleaf, NB, NB, NB)) for nm in NAMES}
        for bid, (ox, oy, oz, lev) in enumerate(blocks):
            for nm in NAMES:
                rows[nm][bid] = grids[lev][nm][oz:oz + NB, oy:oy + NB, ox:ox + NB]
        attrs = dict(cold_attrs)
        attrs.update({"block_nb_x": np.int32(NB), "block_nb_y": np.int32(NB),
                      "block_nb_z": np.int32(NB), "n_blocks": np.int32(nleaf),
                      "step": np.int32(0), "t_current": 0.0})
        with h5py.File(a.out, "w") as f:
            for k, v in attrs.items():
                f.attrs[k] = v
            f.create_dataset("blocks", data=blocks.astype(np.int32))
            for nm in NAMES:
                f.create_dataset(nm, data=rows[nm])
            f.create_dataset("x", data=xt)
            f.create_dataset("y", data=yt)
            f.create_dataset("z", data=zt)
        glob = grids[0]

    umax = max(np.abs(glob[n]).max() for n in ("un", "vn", "wn"))
    print(f"wrote {a.out}: gap=[{Y_LO:.6f},{Y_LO+GAP:.6f}], "
          f"{'refined ' + str(len(blocks)) + ' leaves' if a.leaves else 'single level'}, "
          f"max|vel|={umax:.4f}")


if __name__ == "__main__":
    main()
