#!/usr/bin/env python3
"""Zero-force twin of a mobygeom block-table file for the uniform-flow gate.

The real 3-level file removes buried leaves and carries the penalization
coefficients, so uniform flow cannot survive it. This twin keeps the SAME
touch-driven refinement (the interfaces under test) but zeroes what breaks
a constant field:

  - block_buried_l{l} := 0 everywhere, and the blocks table is rebuilt with
    mobygeom's build_leaf_table_py with burial removal OFF, so neither the
    Python table nor the solver's own builder (which re-reads the zeroed
    masks) removes a leaf -> no closed faces in the flow;
  - coef_blocks := 0 for every leaf -> the IBM exerts no force.

Usage: make_uniform_twin.py <real.h5> <grid.h5> <twin.h5>
"""
import os
import sys

import h5py
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "tools"))
from mobygeom import build_leaf_table_py  # noqa: E402


def main() -> int:
    real, grid, twin = sys.argv[1], sys.argv[2], sys.argv[3]
    with h5py.File(grid, "r") as g:
        periodic = [bool(int(p)) for p in g.attrs["periodic"]]
    with h5py.File(real, "r") as src, h5py.File(twin, "w") as dst:
        for k, v in src.attrs.items():
            dst.attrs[k] = v
        nb = int(src.attrs["block_nb"])
        levels = int(src.attrs["block_levels"])
        gnbt = tuple(int(src.attrs[a]) // nb for a in ("nx", "ny", "nz"))
        # Per-direction refinement mask (absent = xyz octree); xz-quadtree
        # files scale only x,z per level (docs/next_session_refine2d.md).
        mask = np.asarray(src.attrs.get("refine_dims", [1, 1, 1]), dtype=np.int64)
        touch, mwlo, mwdim = [], [], []
        for l in range(levels):
            if f"mask_win_lo_l{l}" in src.attrs:
                lo = np.asarray(src.attrs[f"mask_win_lo_l{l}"], dtype=np.int64)
                dims = np.asarray(src.attrs[f"mask_win_dims_l{l}"], dtype=np.int64)
                dst.attrs[f"mask_win_lo_l{l}"] = src.attrs[f"mask_win_lo_l{l}"]
                dst.attrs[f"mask_win_dims_l{l}"] = src.attrs[f"mask_win_dims_l{l}"]
                if f"lev_win_lo_l{l}" in src.attrs:
                    dst.attrs[f"lev_win_lo_l{l}"] = src.attrs[f"lev_win_lo_l{l}"]
                    dst.attrs[f"lev_win_dims_l{l}"] = src.attrs[f"lev_win_dims_l{l}"]
            else:
                lo = np.zeros(3, dtype=np.int64)
                dims = np.array([g * 2**(l * int(m)) for g, m in zip(gnbt, mask)],
                                dtype=np.int64)
            shape_zyx = tuple(int(d) for d in dims[::-1])
            t = src[f"block_touch_l{l}"][...].reshape(shape_zyx).transpose(2, 1, 0)
            touch.append(t.astype(bool))
            mwlo.append(lo)
            mwdim.append(dims)
            dst.create_dataset(f"block_touch_l{l}", data=src[f"block_touch_l{l}"][...])
            dst.create_dataset(f"block_buried_l{l}",
                               data=np.zeros_like(src[f"block_buried_l{l}"][...]))
        # buried=None: keep every leaf (mirrors the solver reading the
        # zeroed masks); refine_box/lines unused on the body path.
        lev, crd, _, _ = build_leaf_table_py(gnbt, levels, periodic, touch, None, None,
                                             None, nb, mask=mask,
                                             touch_wlo=mwlo, touch_wdim=mwdim)
        n = lev.shape[0]
        blocks = np.empty((n, 4), dtype=np.int32)
        blocks[:, 0:3] = (crd * nb).astype(np.int32)
        blocks[:, 3] = lev.astype(np.int32)
        dst.create_dataset("blocks", data=blocks)
        dst.create_dataset("coef_blocks", shape=(n, nb + 2, nb + 2, nb + 2, 3),
                           dtype=np.float64, chunks=(1, nb + 2, nb + 2, nb + 2, 3))
        # The CELL-CENTRED companion (increment S3): a case with [scalar]
        # requires coef_p_blocks, and the scalar analogue of "the body
        # exerts no force" is "the body penalises no scalar" -- zero here
        # too, so a uniform scalar survives the twin exactly as a uniform
        # velocity does.
        dst.create_dataset("coef_p_blocks", shape=(n, nb + 2, nb + 2, nb + 2),
                           dtype=np.float64, chunks=(1, nb + 2, nb + 2, nb + 2))
        # created zero-filled; nothing to write
        n_real = src["blocks"].shape[0]
        print(f"twin: {n} leaves (real file: {n_real}; "
              f"{n - n_real} buried leaves restored), coef = 0")
    return 0


if __name__ == "__main__":
    sys.exit(main())
