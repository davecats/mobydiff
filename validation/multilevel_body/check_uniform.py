#!/usr/bin/env python3
"""Uniform-flow gate on a multi-level block-table field file.

Checks every STORED value (raw per-block rows -- no reassembly, so coarse
and fine cells are each checked once) against the constants, and asserts
the leaf table actually spans all expected levels (the interfaces under
test exist).

Usage: check_uniform.py <field.h5> --u0 U --v0 V --w0 W --levels 3
"""
import argparse
import sys

import h5py
import numpy as np


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("h5")
    ap.add_argument("--u0", type=float, required=True)
    ap.add_argument("--v0", type=float, required=True)
    ap.add_argument("--w0", type=float, required=True)
    ap.add_argument("--levels", type=int, default=3)
    a = ap.parse_args()

    with h5py.File(a.h5, "r") as f:
        levels = f["blocks"][:, 3]
        hist = np.bincount(levels, minlength=a.levels)
        dev = {n: float(np.max(np.abs(f[n][...] - c)))
               for n, c in (("un", a.u0), ("vn", a.v0), ("wn", a.w0))}
        p = f["pn"][...]
        p_spread = float(p.max() - p.min())

    print(f"{a.h5}: {int(levels.size)} leaves, level histogram {hist.tolist()}")
    for n, d in dev.items():
        print(f"max|{n} - const| = {d:.3e}")
    print(f"pn spread = {p_spread:.3e}")

    ok = all(d == 0.0 for d in dev.values())
    if np.count_nonzero(hist) != a.levels:
        print(f"FAIL: expected leaves on all {a.levels} levels")
        ok = False
    print("uniform gate:", "PASS (exact)" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
