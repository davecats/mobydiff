#!/usr/bin/env python3
"""Exact HDF5 equality: same datasets (names, shapes, bit-identical values)
and same root attributes. Exit 0 on identical, 1 otherwise (h5diff is not
installed on every machine)."""
import sys

import h5py
import numpy as np


def main():
    a_path, b_path = sys.argv[1], sys.argv[2]
    ok = True
    with h5py.File(a_path, "r") as a, h5py.File(b_path, "r") as b:
        a_sets = sorted(a.keys())
        b_sets = sorted(b.keys())
        if a_sets != b_sets:
            print(f"dataset lists differ: {a_sets} vs {b_sets}")
            ok = False
        for name in a_sets:
            if name not in b:
                continue
            da, db = a[name][...], b[name][...]
            if da.shape != db.shape:
                print(f"{name}: shape {da.shape} vs {db.shape}")
                ok = False
            elif not np.array_equal(da, db):
                print(f"{name}: values differ (max abs diff "
                      f"{np.max(np.abs(da.astype(np.float64) - db.astype(np.float64)))})")
                ok = False
        a_attrs = dict(a.attrs)
        b_attrs = dict(b.attrs)
        if sorted(a_attrs) != sorted(b_attrs):
            print(f"attribute lists differ: {sorted(a_attrs)} vs {sorted(b_attrs)}")
            ok = False
        for key in a_attrs:
            if key in b_attrs and not np.array_equal(a_attrs[key], b_attrs[key]):
                print(f"attr {key}: {a_attrs[key]} vs {b_attrs[key]}")
                ok = False
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
