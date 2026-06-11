#!/usr/bin/env python3
"""Compare field datasets in two HDF5 output files."""

from __future__ import annotations

import argparse
import sys

import h5py
import numpy as np


def surviving_block_mask(h5: h5py.File) -> np.ndarray:
    """Boolean (nz,ny,nx) mask of cells covered by the file's block table.

    Cells of removed (solid-buried) blocks are not written by the solver, so
    comparisons against runs without removal must be restricted to this mask.
    """
    nx = int(h5.attrs["nx"])
    ny = int(h5.attrs["ny"])
    nz = int(h5.attrs["nz"])
    nbx = int(h5.attrs["block_nb_x"])
    nby = int(h5.attrs["block_nb_y"])
    nbz = int(h5.attrs["block_nb_z"])
    blocks = h5["blocks"][...]
    mask = np.zeros((nz, ny, nx), dtype=bool)
    for ox, oy, oz, _level in blocks:
        mask[oz:oz + nbz, oy:oy + nby, ox:ox + nbx] = True
    return mask


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("reference")
    parser.add_argument("candidate")
    parser.add_argument("datasets", nargs="*", default=["un", "vn", "wn", "pn"])
    parser.add_argument("--tolerance", type=float, default=None)
    parser.add_argument("--mask-surviving", action="store_true",
                        help="compare only cells covered by the candidate's block table")
    args = parser.parse_args()

    failed = False
    with h5py.File(args.reference, "r") as ref, h5py.File(args.candidate, "r") as cand:
        mask = surviving_block_mask(cand) if args.mask_surviving else None
        if mask is not None:
            covered = int(mask.sum())
            print(f"masked to {covered} of {mask.size} cells "
                  f"({mask.size - covered} in removed blocks)")
        for name in args.datasets:
            diff = ref[name][...] - cand[name][...]
            if mask is not None:
                diff = diff[mask]
            max_abs = float(np.max(np.abs(diff)))
            l2 = float(np.sqrt(np.mean(diff * diff)))
            print(f"{name:2s} max_abs={max_abs:.16e} l2={l2:.16e}")
            if args.tolerance is not None and max_abs > args.tolerance:
                failed = True

    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
