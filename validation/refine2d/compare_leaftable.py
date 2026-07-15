#!/usr/bin/env python3
"""R2D-0 gate (b): compare the solver leaf table (leaftable_test stdout)
row-by-row against the mobygeom block-table `blocks` dataset.

usage: compare_leaftable.py <leaftable_test.out> <block_table.h5> <dims>
"""
import sys

import h5py
import numpy as np


def main() -> int:
    out_path, h5_path, dims = sys.argv[1], sys.argv[2], sys.argv[3]

    rows = []
    seen = False
    with open(out_path) as fh:
        for line in fh:
            if line.strip().startswith("leaftable"):
                seen = True
                continue
            if seen and line.strip():
                rows.append([int(v) for v in line.split()])
    ftab = np.array(rows, dtype=np.int64)[:, 1:]  # drop the id column

    with h5py.File(h5_path, "r") as h5:
        ptab = h5["blocks"][...].astype(np.int64)
        attr = np.asarray(h5.attrs.get("refine_dims", [1, 1, 1]), dtype=np.int64)

    want = {"xyz": [1, 1, 1], "xz": [1, 0, 1]}[dims]
    if list(attr) != want:
        print(f"FAIL: refine_dims attribute {list(attr)} != expected {want}")
        return 1
    if ftab.shape != ptab.shape:
        print(f"FAIL: {ftab.shape[0]} solver leaves != {ptab.shape[0]} mobygeom leaves")
        return 1
    bad = np.nonzero((ftab != ptab).any(axis=1))[0]
    if bad.size:
        print(f"FAIL: {bad.size} mismatching rows; first: id {bad[0]} "
              f"solver {ftab[bad[0]]} mobygeom {ptab[bad[0]]}")
        return 1
    levels = np.bincount(ptab[:, 3])
    print(f"PASS {dims}: {ftab.shape[0]} leaves identical row-by-row "
          f"(levels {levels.tolist()})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
