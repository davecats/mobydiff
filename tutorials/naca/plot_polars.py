#!/usr/bin/env python3
"""C_L(alpha) and C_L-C_D polars from the sweep force histories, with
optional XFOIL and OpenFOAM overlays.

  plot_polars.py [--tail 0.2] [--out polars.png]
                 [--xfoil polar.txt]  (XFOIL PACC polar dump: the table
                                       with columns alpha CL CD CDp CM ...)
                 [--openfoam of.csv]  (CSV/whitespace: alpha, CL, CD
                                       [, anything]; '#' comments)

Reads every forces_aoa<tag>.txt in this directory (tag m2 -> alpha = -2).
The mean over the LAST --tail fraction of each history is the converged
coefficient; the printed rms over the tail flags unconverged angles.
Writes the mobydiff polar to polar_mobydiff.dat as well.
"""
import argparse
import glob
import os
import re

import numpy as np


def read_forces(path, tail):
    d = np.loadtxt(path, skiprows=1)
    tl = d[int((1.0 - tail)*d.shape[0]):]
    return (float(tl[:, 2].mean()), float(tl[:, 3].mean()),
            float(tl[:, 2].std()), float(tl[:, 3].std()))


def read_xfoil_polar(path):
    """XFOIL PACC dump: skip the header down to the '----' rule, then
    columns alpha CL CD CDp CM ..."""
    rows = []
    started = False
    with open(path) as f:
        for line in f:
            if started:
                parts = line.split()
                if len(parts) >= 3:
                    try:
                        rows.append([float(parts[0]), float(parts[1]), float(parts[2])])
                    except ValueError:
                        pass
            elif set(line.strip()) <= {"-", " "} and "-" in line:
                started = True
    if not rows:
        raise SystemExit(f"{path}: no XFOIL polar table found")
    return np.array(rows)


def read_generic_polar(path):
    """alpha, CL, CD columns; '#' comments; comma or whitespace separated."""
    rows = []
    with open(path) as f:
        for line in f:
            line = line.split("#")[0].strip().replace(",", " ")
            parts = line.split()
            if len(parts) >= 3:
                try:
                    rows.append([float(parts[0]), float(parts[1]), float(parts[2])])
                except ValueError:
                    pass
    if not rows:
        raise SystemExit(f"{path}: no alpha/CL/CD rows found")
    return np.array(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tail", type=float, default=0.2)
    ap.add_argument("--out", default="polars.png")
    ap.add_argument("--xfoil", default=None)
    ap.add_argument("--openfoam", default=None)
    a = ap.parse_args()

    here = os.path.dirname(os.path.abspath(__file__))
    rows = []
    for p in sorted(glob.glob(os.path.join(here, "forces_aoa*.txt"))):
        m = re.search(r"forces_aoa(m?\d+)\.txt$", p)
        if not m:
            continue
        alpha = float(m.group(1).replace("m", "-"))
        cl, cd, scl, scd = read_forces(p, a.tail)
        rows.append((alpha, cl, cd, scl, scd))
        conv = "" if scl < 5e-3 else "   <-- tail rms high, check convergence"
        print(f"alpha {alpha:+5.1f}: C_L = {cl:+.4f} (rms {scl:.1e})  "
              f"C_D = {cd:.4f} (rms {scd:.1e}){conv}")
    if not rows:
        raise SystemExit("no forces_aoa*.txt found -- run the sweep first")
    rows.sort()
    pol = np.array(rows)
    if pol.shape[0] >= 2:
        slope = np.polyfit(pol[:, 0], pol[:, 1], 1)[0]
        print(f"lift slope (fit over all angles) = {slope:.4f}/deg "
              f"= {100.0*slope/0.10966:.0f}% of 2pi/rad")
    np.savetxt(os.path.join(here, "polar_mobydiff.dat"), pol[:, :3],
               header="alpha CL CD (tail-%.2f means)" % a.tail, fmt="%+.6e")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))
    ax1.errorbar(pol[:, 0], pol[:, 1], yerr=pol[:, 3], marker="o",
                 label="mobydiff (L6-xz IBM)", zorder=3)
    ax2.plot(pol[:, 2], pol[:, 1], marker="o", label="mobydiff (L6-xz IBM)", zorder=3)
    if a.xfoil:
        xf = read_xfoil_polar(a.xfoil)
        ax1.plot(xf[:, 0], xf[:, 1], "s--", ms=4, label="XFOIL")
        ax2.plot(xf[:, 2], xf[:, 1], "s--", ms=4, label="XFOIL")
    if a.openfoam:
        of = read_generic_polar(a.openfoam)
        ax1.plot(of[:, 0], of[:, 1], "^:", ms=5, label="OpenFOAM")
        ax2.plot(of[:, 2], of[:, 1], "^:", ms=5, label="OpenFOAM")
    ax1.set_xlabel(r"$\alpha$ [deg]"); ax1.set_ylabel(r"$C_L$")
    ax1.grid(alpha=0.3); ax1.legend()
    ax2.set_xlabel(r"$C_D$"); ax2.set_ylabel(r"$C_L$")
    ax2.grid(alpha=0.3); ax2.legend()
    fig.suptitle("NACA 0012, Re = 4e5, fully turbulent (SST)")
    fig.tight_layout()
    fig.savefig(a.out, dpi=150)
    print(f"wrote {a.out} and polar_mobydiff.dat")


if __name__ == "__main__":
    main()
