#!/usr/bin/env python3
"""Chunked comparison of two BLOCK-TABLE field snapshots with identical
leaf tables (per-block rows compared directly, a few thousand at a time).

tools/compare_fields.py reassembles block-table files onto the FINEST
lattice, which for deep-refinement airfoil cases is tens of GB and gets
the process OOM-killed; this comparator never materializes more than a
chunk. Use compare_fields.py when the two files may have different block
layouts or a legacy global reference.

usage: compare_snapshots.py a.h5 b.h5 [--tolerance T] dataset...
exit 0 when every dataset's max |diff| <= tolerance (default: report only).
"""
import argparse
import sys

import h5py
import numpy as np

CHUNK = 2000


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("reference")
    ap.add_argument("candidate")
    ap.add_argument("datasets", nargs="+")
    ap.add_argument("--tolerance", type=float, default=None)
    args = ap.parse_args()

    ok = True
    with h5py.File(args.reference, "r") as a, h5py.File(args.candidate, "r") as b:
        if not np.array_equal(a["blocks"][...], b["blocks"][...]):
            print("blocks tables differ -- use tools/compare_fields.py")
            sys.exit(1)
        for name in args.datasets:
            da, db = a[name], b[name]
            if da.shape != db.shape:
                print(f"{name}: shape {da.shape} vs {db.shape}")
                ok = False
                continue
            max_abs = 0.0
            for lo in range(0, da.shape[0], CHUNK):
                hi = min(lo + CHUNK, da.shape[0])
                max_abs = max(max_abs, float(np.max(np.abs(
                    da[lo:hi][...] - db[lo:hi][...]))))
            print(f"{name} max_abs={max_abs:.16e}")
            if args.tolerance is not None and max_abs > args.tolerance:
                ok = False
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
