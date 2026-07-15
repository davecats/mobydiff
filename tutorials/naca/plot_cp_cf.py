#!/usr/bin/env python3
"""Overlay mobydiff surface Cp/Cf (surface_cp_cf.py npz) with XFOIL and
OpenFOAM references.

  plot_cp_cf.py cpcf_naca_aoa4_100000.npz [--out cpcf_aoa4.png]
      [--xfoil-cp cpx.txt]   XFOIL CPWR dump: columns x [y] Cp, one loop
      [--xfoil-cf cfx.txt]   XFOIL "VPLO -> CF" dump or 2-col x cf
      [--of-cp of_cp.csv]    generic 2+ col: x/c, Cp  ('#' comments, csv ok)
      [--of-cf of_cf.csv]    generic 2+ col: x/c, Cf

External files are plotted as given (single loops are fine -- they wrap
both sides); mobydiff curves are drawn per side.
"""
import argparse

import numpy as np


def read_cols(path, ncol_min=2):
    rows = []
    with open(path) as f:
        for line in f:
            line = line.split("#")[0].strip().replace(",", " ")
            parts = line.split()
            if len(parts) >= ncol_min:
                try:
                    rows.append([float(v) for v in parts])
                except ValueError:
                    continue
    if not rows:
        raise SystemExit(f"{path}: no numeric rows")
    n = min(len(r) for r in rows)
    return np.array([r[:n] for r in rows])


def xfoil_xy(path):
    """XFOIL CPWR/CF dumps: 2 cols (x, val) or 3 cols (x, y, val)."""
    d = read_cols(path)
    return (d[:, 0], d[:, -1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("npz")
    ap.add_argument("--out", default="cpcf_compare.png")
    ap.add_argument("--xfoil-cp", default=None)
    ap.add_argument("--xfoil-cf", default=None)
    ap.add_argument("--of-cp", default=None)
    ap.add_argument("--of-cf", default=None)
    a = ap.parse_args()

    d = np.load(a.npz)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.6))
    for side, lab, c in ((1.0, "mobydiff upper", "C0"), (0.0, "mobydiff lower", "C1")):
        sel = d["up"] == side
        o = np.argsort(d["xoc"][sel])
        ax1.plot(d["xoc"][sel][o], d["cp"][sel][o], c, label=lab)
        ax2.plot(d["xoc"][sel][o], d["cf"][sel][o], c, label=lab)
    if a.xfoil_cp:
        x, v = xfoil_xy(a.xfoil_cp)
        ax1.plot(x, v, "k--", lw=1, label="XFOIL")
    if a.xfoil_cf:
        x, v = xfoil_xy(a.xfoil_cf)
        ax2.plot(x, v, "k--", lw=1, label="XFOIL")
    if a.of_cp:
        x, v = xfoil_xy(a.of_cp)
        ax1.plot(x, v, "g:", lw=1.5, label="OpenFOAM")
    if a.of_cf:
        x, v = xfoil_xy(a.of_cf)
        ax2.plot(x, v, "g:", lw=1.5, label="OpenFOAM")
    ax1.invert_yaxis()
    ax1.set_xlabel("x/c"); ax1.set_ylabel(r"$C_p$")
    ax2.set_xlabel("x/c"); ax2.set_ylabel(r"$C_f$")
    ax2.axhline(0.0, color="0.7", lw=0.8)
    for ax in (ax1, ax2):
        ax.grid(alpha=0.3); ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(a.out, dpi=150)
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
