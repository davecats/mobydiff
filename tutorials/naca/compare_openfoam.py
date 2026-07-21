#!/usr/bin/env python3
"""Compare the C10 mobydiff aoa5 surface data and polar point against an
OpenFOAM kOmegaSST reference case (simpleFoam, body-fixed mesh, inlet
velocity rotated by alpha — the same convention as the C10 setup).

  compare_openfoam.py <of_case_dir> [--npz cpcf_c10_aoa5_150000.npz]
                      [--time latest] [--out cpcf_vs_openfoam.png]

OpenFOAM surface sampling (.raw): p_{suction,pressure}_side.raw with
columns x y z p (rho = 1, U = 1, p_inf = 0 -> Cp = 2 p) and
wallShearStress_*.raw with x y z tau_x tau_y tau_z (the traction on the
wall; attached TE-ward flow has tau_x < 0 -> Cf = 2 |tau| sign(-tau_x)).
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def load_raw(path):
    return np.loadtxt(path, comments="#")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("of_dir")
    ap.add_argument("--npz", default="cpcf_c10_aoa5_150000.npz")
    ap.add_argument("--time", default=None)
    ap.add_argument("--out", default="cpcf_vs_openfoam.png")
    a = ap.parse_args()

    samp = os.path.join(a.of_dir, "postProcessing", "surface_sampling")
    t = a.time or sorted(os.listdir(samp), key=float)[-1]
    of = {}
    for side in ("suction", "pressure"):
        p = load_raw(os.path.join(samp, t, f"p_{side}_side.raw"))
        tau = load_raw(os.path.join(samp, t, f"wallShearStress_{side}_side.raw"))
        o = np.argsort(p[:, 0])
        of[side] = {"x": p[o, 0], "cp": 2.0*p[o, 3]}
        o = np.argsort(tau[:, 0])
        mag = np.linalg.norm(tau[o, 3:6], axis=1)
        of[side]["xt"] = tau[o, 0]
        of[side]["cf"] = 2.0*mag*np.sign(-tau[o, 3])
    print(f"OpenFOAM sampling time {t}")

    d = np.load(a.npz)
    mb = {}
    for side, up in (("suction", 1.0), ("pressure", 0.0)):
        sel = (d["up"] == up)
        o = np.argsort(d["xoc"][sel])
        mb[side] = {"x": d["xoc"][sel][o], "cp": d["cp"][sel][o],
                    "cf": d["cf"][sel][o]}

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 4.6))
    cols = {"suction": "tab:blue", "pressure": "tab:orange"}
    for side in ("suction", "pressure"):
        ax1.plot(mb[side]["x"], mb[side]["cp"], "-", color=cols[side],
                 label=f"mobydiff {side}")
        ax1.plot(of[side]["x"], of[side]["cp"], "--", color=cols[side],
                 lw=1.2, label=f"OpenFOAM {side}")
        ax2.plot(mb[side]["x"], mb[side]["cf"], "-", color=cols[side],
                 label=f"mobydiff {side}")
        ax2.plot(of[side]["xt"], of[side]["cf"], "--", color=cols[side],
                 lw=1.2, label=f"OpenFOAM {side}")
    ax1.invert_yaxis()
    ax1.set_xlabel("x/c"); ax1.set_ylabel(r"$C_p$")
    ax2.set_xlabel("x/c"); ax2.set_ylabel(r"$C_f$")
    ax2.set_ylim(-0.005, 0.03)
    for ax in (ax1, ax2):
        ax.grid(alpha=0.3); ax.legend(fontsize=8)
    fig.suptitle(r"NACA 0012, Re = 4e5, $\alpha$ = 5$^\circ$: "
                 "mobydiff C10 IBM (SST, tu 5%) vs OpenFOAM kOmegaSST")
    fig.tight_layout()
    fig.savefig(a.out, dpi=150)
    print(f"wrote {a.out}")

    # scalar summary: suction peak, TE separation onset, integrated Cf
    for name, s in (("mobydiff", mb), ("OpenFOAM", of)):
        pk = np.nanmin(s["suction"]["cp"])
        xk = "x" if "cf" not in s["suction"] else ("x" if name == "mobydiff" else "xt")
        cf_s = s["suction"]["cf"]
        xs = s["suction"][xk if name == "OpenFOAM" else "x"]
        m = np.isfinite(cf_s)
        neg = xs[m][cf_s[m] < 0]
        xsep = neg[neg > 0.5].min() if (neg > 0.5).any() else np.nan
        print(f"{name:9s}: Cp_min = {pk:+.3f}   suction-side Cf<0 from "
              f"x/c = {xsep:.3f}" if np.isfinite(xsep) else
              f"{name:9s}: Cp_min = {pk:+.3f}   no TE separation")


if __name__ == "__main__":
    main()
