#!/usr/bin/env python3
"""Compare field datasets in two HDF5 output files.

Handles both the legacy global 3D layout and the block-table layout
(per-variable datasets of shape (nBlocksGlobal, nbz, nby, nbx) plus the
`blocks` table). Block-table files are reassembled onto the finest-level
lattice, replicating coarser cells, so single-level files compare exactly
against global-layout references at the same resolution.
"""

from __future__ import annotations

import argparse
import sys

import h5py
import numpy as np

FIELDS = ["un", "vn", "wn", "pn"]


def is_block_format(h5: h5py.File) -> bool:
    return h5["un"].ndim == 4


def block_geometry(h5: h5py.File):
    blocks = h5["blocks"][...]
    nb = (int(h5.attrs["block_nb_x"]), int(h5.attrs["block_nb_y"]), int(h5.attrs["block_nb_z"]))
    lmax = int(blocks[:, 3].max())
    shape = tuple(int(h5.attrs[a]) * 2**lmax for a in ("nz", "ny", "nx"))
    return blocks, nb, lmax, shape


def load_field(h5: h5py.File, name: str) -> np.ndarray:
    if not is_block_format(h5):
        return h5[name][...]
    blocks, nb, lmax, shape = block_geometry(h5)
    arr = np.zeros(shape, dtype=np.float64)
    data = h5[name]
    for bid, (ox, oy, oz, lev) in enumerate(blocks):
        f = 2 ** (lmax - int(lev))
        row = data[bid]
        if f > 1:
            row = row.repeat(f, axis=0).repeat(f, axis=1).repeat(f, axis=2)
        arr[oz * f:oz * f + nb[2] * f, oy * f:oy * f + nb[1] * f, ox * f:ox * f + nb[0] * f] = row
    return arr


def surviving_block_mask(h5: h5py.File) -> np.ndarray:
    """Boolean finest-lattice mask of cells covered by the file's block table."""
    if is_block_format(h5):
        blocks, nb, lmax, shape = block_geometry(h5)
    else:
        blocks = h5["blocks"][...]
        nb = (int(h5.attrs["block_nb_x"]), int(h5.attrs["block_nb_y"]), int(h5.attrs["block_nb_z"]))
        lmax = int(blocks[:, 3].max())
        shape = tuple(int(h5.attrs[a]) * 2**lmax for a in ("nz", "ny", "nx"))
    mask = np.zeros(shape, dtype=bool)
    for ox, oy, oz, lev in blocks:
        f = 2 ** (lmax - int(lev))
        mask[oz * f:oz * f + nb[2] * f, oy * f:oy * f + nb[1] * f, ox * f:ox * f + nb[0] * f] = True
    return mask


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("reference")
    parser.add_argument("candidate")
    parser.add_argument("datasets", nargs="*", default=FIELDS)
    parser.add_argument("--tolerance", type=float, default=None)
    parser.add_argument("--mask-surviving", action="store_true",
                        help="compare only cells covered by the candidate's block table")
    parser.add_argument("--export-global", default=None,
                        help="write the candidate reassembled onto the finest lattice to this HDF5 file")
    args = parser.parse_args()

    failed = False
    with h5py.File(args.reference, "r") as ref, h5py.File(args.candidate, "r") as cand:
        mask = surviving_block_mask(cand) if args.mask_surviving else None
        if mask is not None:
            covered = int(mask.sum())
            print(f"masked to {covered} of {mask.size} cells "
                  f"({mask.size - covered} in removed blocks)")
        out = h5py.File(args.export_global, "w") if args.export_global else None
        for name in args.datasets:
            a = load_field(ref, name)
            b = load_field(cand, name)
            if out is not None:
                out.create_dataset(name, data=b)
            diff = a - b
            if mask is not None:
                diff = diff[mask]
            max_abs = float(np.max(np.abs(diff)))
            l2 = float(np.sqrt(np.mean(diff * diff)))
            print(f"{name:2s} max_abs={max_abs:.16e} l2={l2:.16e}")
            if args.tolerance is not None and max_abs > args.tolerance:
                failed = True
        if out is not None:
            out.close()

    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
