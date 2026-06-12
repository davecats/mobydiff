#!/usr/bin/env python3
"""Assert that the interface-decay gate contracted.

Reads the step-20 and step-200 snapshots of tutorials/interface_decay
(block-table layout) and fails unless max|u|,|v|,|w|,|p| all decreased
and are finite. Step 20 is the baseline rather than step 0 so the
initial projection of the (divergent) white noise does not count as
growth.
"""

import argparse
import sys

import h5py
import numpy as np

FIELDS = ["un", "vn", "wn", "pn"]


def field_max(path):
    with h5py.File(path, "r") as h5:
        return {name: float(np.abs(h5[name][...]).max()) for name in FIELDS}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("rundir", help="directory containing decay_20.h5 and decay_200.h5")
    parser.add_argument("--prefix", default="decay")
    parser.add_argument("--first", type=int, default=20)
    parser.add_argument("--last", type=int, default=200)
    args = parser.parse_args()

    early = field_max(f"{args.rundir}/{args.prefix}_{args.first}.h5")
    late = field_max(f"{args.rundir}/{args.prefix}_{args.last}.h5")

    failed = False
    for name in FIELDS:
        ok = np.isfinite(late[name]) and late[name] < early[name]
        status = "ok " if ok else "FAIL"
        print(f"{status} {name}: step {args.first} max {early[name]:.6e} -> "
              f"step {args.last} max {late[name]:.6e}")
        failed = failed or not ok

    if failed:
        print("interface decay gate FAILED: noise did not contract")
        return 1
    print("interface decay gate passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
