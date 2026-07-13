#!/usr/bin/env python3
"""A3 INCREMENT 2 gate: NACA 0012 full-turbulent SST sanity vs XFOIL/
literature — a BAND, not a tight match (first-order penalization
D_eff ~ D + h, 12c Dirichlet blockage, y+_1 ~ 2-4 resolved walls).

  check_naca.py forces_aoa0.txt forces_aoa4.txt forces_aoa8.txt

Averages the last --tail fraction of each history, prints C_L/C_D per
alpha, the lift slope per degree, and gates:
  |C_L(0)| < 0.02, slope in [0.06, 0.13]/deg (2pi rad = 0.1097/deg),
  C_D(0) in [0.008, 0.035].
"""
import argparse
import re
import sys

import numpy as np


def tail_mean(path: str, frac: float) -> tuple[float, float, float]:
    data = np.loadtxt(path, skiprows=1)
    n = data.shape[0]
    t = data[int((1.0 - frac) * n):, :]
    return float(t[:, 2].mean()), float(t[:, 3].mean()), float(data[-1, 1])


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("forces", nargs="+", help="forces_aoa<deg>.txt files")
    ap.add_argument("--tail", type=float, default=0.2)
    a = ap.parse_args()

    rows = []
    for path in a.forces:
        m = re.search(r"aoa(\d+)", path)
        alpha = float(m.group(1)) if m else float("nan")
        cl, cd, t_end = tail_mean(path, a.tail)
        rows.append((alpha, cl, cd))
        print(f"alpha = {alpha:4.1f}: C_L = {cl:+.4f}  C_D = {cd:.4f}  (t_end = {t_end:.1f})")

    rows.sort()
    alphas = np.array([r[0] for r in rows])
    cls = np.array([r[1] for r in rows])
    ok = True
    if 0.0 in alphas:
        cl0 = cls[alphas == 0.0][0]
        cd0 = [r[2] for r in rows if r[0] == 0.0][0]
        if abs(cl0) >= 0.02:
            print(f"FAIL: |C_L(0)| = {abs(cl0):.4f} >= 0.02"); ok = False
        if not (0.008 <= cd0 <= 0.035):
            print(f"FAIL: C_D(0) = {cd0:.4f} outside [0.008, 0.035]"); ok = False
    if alphas.size >= 2:
        slope = np.polyfit(alphas, cls, 1)[0]
        print(f"lift slope = {slope:.4f}/deg  (2pi rad = 0.1097/deg)")
        if not (0.06 <= slope <= 0.13):
            print("FAIL: slope outside [0.06, 0.13]/deg"); ok = False
    print("naca0012 gate:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
