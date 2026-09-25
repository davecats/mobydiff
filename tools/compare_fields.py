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


def field_datasets(h5: h5py.File) -> list[str]:
    """Every FIELD dataset of a file: same rank as `un` (4 in the block-table
    layout, 3 in the legacy global one), which excludes the `blocks` table and
    the x/y/z node lines without naming them."""
    ndim = h5["un"].ndim
    return [k for k, v in h5.items()
            if isinstance(v, h5py.Dataset) and v.ndim == ndim]


def common_datasets(ref: h5py.File, cand: h5py.File) -> list[str]:
    """The datasets present in BOTH files, canonical variables first and the
    rest (scalars, nut, k, omega, gamma, rethetat, fd) alphabetically -- so a
    run with passive scalars is compared on every dataset it wrote without
    having to list them."""
    both = set(field_datasets(ref)) & set(field_datasets(cand))
    head = [f for f in FIELDS if f in both]
    return head + sorted(both - set(head))


def is_block_format(h5: h5py.File) -> bool:
    return h5["un"].ndim == 4


def block_geometry(h5: h5py.File):
    """Blocks table, block size, per-direction replication factors (z,y,x
    order) at each level, and the finest-lattice shape. The refine_dims
    attribute (absent = xyz octree) marks the xz-quadtree variant, where y
    never refines and block y origins are global cells already."""
    blocks = h5["blocks"][...]
    nb = (int(h5.attrs["block_nb_x"]), int(h5.attrs["block_nb_y"]), int(h5.attrs["block_nb_z"]))
    lmax = int(blocks[:, 3].max())
    mask = np.asarray(h5.attrs.get("refine_dims", [1, 1, 1]), dtype=np.int64)  # x,y,z
    mzyx = mask[::-1]
    shape = tuple(int(h5.attrs[a]) * 2**(lmax * int(m)) for a, m in zip(("nz", "ny", "nx"), mzyx))
    return blocks, nb, lmax, mzyx, shape


def load_field(h5: h5py.File, name: str) -> np.ndarray:
    if not is_block_format(h5):
        return h5[name][...]
    blocks, nb, lmax, mzyx, shape = block_geometry(h5)
    arr = np.zeros(shape, dtype=np.float64)
    data = h5[name]
    for bid, (ox, oy, oz, lev) in enumerate(blocks):
        fz, fy, fx = (2 ** ((lmax - int(lev)) * int(m)) for m in mzyx)
        row = data[bid]
        if max(fx, fy, fz) > 1:
            row = row.repeat(fz, axis=0).repeat(fy, axis=1).repeat(fx, axis=2)
        arr[oz * fz:oz * fz + nb[2] * fz,
            oy * fy:oy * fy + nb[1] * fy,
            ox * fx:ox * fx + nb[0] * fx] = row
    return arr


def surviving_block_mask(h5: h5py.File) -> np.ndarray:
    """Boolean finest-lattice mask of cells covered by the file's block table."""
    blocks, nb, lmax, mzyx, shape = block_geometry(h5)
    mask = np.zeros(shape, dtype=bool)
    for ox, oy, oz, lev in blocks:
        fz, fy, fx = (2 ** ((lmax - int(lev)) * int(m)) for m in mzyx)
        mask[oz * fz:oz * fz + nb[2] * fz,
             oy * fy:oy * fy + nb[1] * fy,
             ox * fx:ox * fx + nb[0] * fx] = True
    return mask


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("reference")
    parser.add_argument("candidate")
    parser.add_argument("datasets", nargs="*",
                        help="datasets to compare (default: those present in both files)")
    parser.add_argument("--tolerance", type=float, default=None)
    parser.add_argument("--mask-surviving", action="store_true",
                        help="compare only cells covered by the candidate's block table")
    parser.add_argument("--export-global", default=None,
                        help="write the candidate reassembled onto the finest lattice to this HDF5 file")
    args = parser.parse_args()

    failed = False
    with h5py.File(args.reference, "r") as ref, h5py.File(args.candidate, "r") as cand:
        datasets = args.datasets or common_datasets(ref, cand)
        if not args.datasets:
            print("datasets: " + " ".join(datasets))
        mask = surviving_block_mask(cand) if args.mask_surviving else None
        if mask is not None:
            covered = int(mask.sum())
            print(f"masked to {covered} of {mask.size} cells "
                  f"({mask.size - covered} in removed blocks)")
        out = h5py.File(args.export_global, "w") if args.export_global else None
        for name in datasets:
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
