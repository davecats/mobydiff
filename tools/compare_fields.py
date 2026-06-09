#!/usr/bin/env python3
"""Compare field datasets in two HDF5 output files."""

from __future__ import annotations

import argparse
import sys

import h5py
import numpy as np


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("reference")
    parser.add_argument("candidate")
    parser.add_argument("datasets", nargs="*", default=["un", "vn", "wn", "pn"])
    parser.add_argument("--tolerance", type=float, default=None)
    args = parser.parse_args()

    failed = False
    with h5py.File(args.reference, "r") as ref, h5py.File(args.candidate, "r") as cand:
        for name in args.datasets:
            diff = ref[name][...] - cand[name][...]
            max_abs = float(np.max(np.abs(diff)))
            l2 = float(np.sqrt(np.mean(diff * diff)))
            print(f"{name:2s} max_abs={max_abs:.16e} l2={l2:.16e}")
            if args.tolerance is not None and max_abs > args.tolerance:
                failed = True

    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
