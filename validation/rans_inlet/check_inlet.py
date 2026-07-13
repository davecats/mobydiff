#!/usr/bin/env python3
"""Gate (b) metrics for the RANS scalar inlet (A3 INCREMENT 1).

The inlet ghost holds the freestream values via the midpoint identity
(ghost = 2 v - interior), so at steady inflow the first interior column
relaxes to the freestream: check the mid-channel k/omega there against
k_inf = 1.5 (tu/100 U_inf)^2 and omega_inf = k_inf Re / nut_ratio
(smoke tolerance -- convection-dominated, one cell of adjustment).

Usage: check_inlet.py <field.h5> --tu 5.0 --nut-ratio 10.0 --re 100.0
"""
import argparse
import os
import sys

import h5py
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "tools"))
from compare_fields import load_field  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("h5")
    ap.add_argument("--tu", type=float, required=True)
    ap.add_argument("--nut-ratio", type=float, required=True)
    ap.add_argument("--re", type=float, required=True)
    ap.add_argument("--u-inf", type=float, default=1.0)
    ap.add_argument("--rtol", type=float, default=0.10)
    a = ap.parse_args()

    k_inf = 1.5 * (a.tu / 100.0 * a.u_inf) ** 2
    omg_inf = k_inf * a.re / a.nut_ratio

    with h5py.File(a.h5, "r") as f:
        k = load_field(f, "k")
        omg = load_field(f, "omega")

    ny = k.shape[1]
    jm = ny // 2  # mid-channel row: profile peak, u = U_inf, du/dy = 0
    k1 = float(k[:, jm, 0].mean())
    o1 = float(omg[:, jm, 0].mean())
    dk = abs(k1 / k_inf - 1.0)
    do = abs(o1 / omg_inf - 1.0)
    print(f"k_inf = {k_inf:.6e}   first-column mid-channel k = {k1:.6e}  (dev {dk:.2%})")
    print(f"omg_inf = {omg_inf:.6e} first-column mid-channel omega = {o1:.6e}  (dev {do:.2%})")
    ok = dk < a.rtol and do < a.rtol
    print("inlet gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
