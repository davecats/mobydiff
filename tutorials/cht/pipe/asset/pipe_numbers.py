#!/usr/bin/env python3
"""Every table in pipe_report.md, computed from the accumulators.

Regenerating the report after a longer statistics leg is then one command,
and no number in the text is typed by hand.

    ./pipe_numbers.py > pipe_numbers.txt
"""
from __future__ import annotations
import sys
import numpy as np
import os
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))      # pipe_stats.py, beside the inis
from pipe_stats import profiles
from compare_neuhauser import ref_case


def asset(name):
    return os.path.join(HERE, name)

NU, U_TAU, DELTA_V, PR = 1/5300., 0.068379, 2.759318e-3, 0.71
R, LZ = 0.5, 12.5
A_WALL = 2*np.pi*R*LZ
NAMES = ["c0", "c1", "c2", "c3", "c4", "mbc", "isof"]
# The production grid, whose accumulators ship in this directory. The three
# coarser grids of the resolution study are reported in the README from the
# numbers this same script printed for them; their accumulators are not kept.
GRID = {"production": dict(stats="pipe_prod_statsD.npz",
                           snaps="pipe_prod_snaps.npz",
                           heat="pipe_prod_heat.txt",
                           dxp=1.84, dzp=5.06, n="256x256x896")}


def wall_value(r, t, interface=0.5):
    m = (r < interface) & np.isfinite(t)
    a, b = np.polyfit(r[m][-6:], t[m][-6:], 1)
    return a*interface + b


def main():
    REF = {nm: ref_case(nm, PR) for nm in NAMES}
    rr = REF["c0"]["r"]; mf = rr < 0.5
    ttau_ref = 1.00330/U_TAU
    tw_ref = {nm: wall_value(rr, REF[nm]["t"]) for nm in NAMES}

    print("== TABLE: interface heat balance (time mean, % of the discrete generation)")
    for tag, g in GRID.items():
        hf = g["heat"]
        h = np.vstack([np.loadtxt(asset(f)) for f in hf]) \
            if isinstance(hf, list) else np.loadtxt(asset(hf))
        tot = np.array([[row[2+3*i+2] for i in range(7)] for row in h])
        g["Q"] = -tot.mean(axis=0)[0]
        g["qw"] = g["Q"]/A_WALL
        gen = g["Q"]
        print(f"   {tag:7} t = {h[0,1]:.0f}..{h[-1,1]:.0f}  " +
              "  ".join(f"{n} {v:7.3f}" for n, v in zip(NAMES, -tot.mean(axis=0))))
    print()

    for tag, g in GRID.items():
        g["P"] = profiles(dict(np.load(asset(g["stats"]), allow_pickle=True)))
        try:
            g["S"] = profiles(dict(np.load(asset(g["snaps"]), allow_pickle=True)))
        except FileNotFoundError:
            g["S"] = None
        g["ttau"] = g["qw"]/U_TAU

    print("== TABLE: velocity statistics (snapshots)")
    for tag, g in GRID.items():
        S = g["S"]
        if S is None:
            continue
        f = (S["count"][:, 0] > 0) & (S["r"] < 0.5)
        r_, w = S["r"][f], S["uz"][f, 0]
        dev = lambda k, rk: 100*(S[k][f, 0] - np.interp(r_, rr, REF["c0"][rk])) / \
            np.interp(r_, rr, REF["c0"][rk])
        core = r_ < 0.1
        band = (r_ > 0.35) & (r_ < 0.47)
        tau = -NU*np.gradient(w, r_) + S["urz_cov"][f, 0]
        ex = U_TAU**2*r_/R
        m = r_ < 0.48
        print(f"   {tag:7} <uz> core {dev('uz','uz')[core].mean():+5.1f}%  "
              f"uz' core {dev('uz_rms','uz_rms')[core].mean():+5.1f}%  "
              f"ur' {dev('ur_rms','ur_rms')[band].mean():+5.1f}%  "
              f"ut' {dev('up_rms','ut_rms')[band].mean():+5.1f}%  "
              f"<ur uz> {dev('urz_cov','urz')[band].mean():+5.1f}%  "
              f"| stress law max {100*np.abs(tau[m]-ex[m]).max()/U_TAU**2:.1f}%"
              f" mean {100*np.abs(tau[m]-ex[m]).mean()/U_TAU**2:.2f}%")
    print()

    print("== TABLE: thermal profiles, c0 (plane statistics, balance-normalised)")
    for tag, g in GRID.items():
        P = g["P"]; f = (P["count"][:, 0] > 0) & (P["r"] < 0.5)
        r_ = P["r"][f]
        tm = (P["c0_mean"][f, 0] - wall_value(r_, P["c0_mean"][f, 0]))/g["ttau"]
        qm = np.interp(r_, rr, REF["c0"]["t"] - tw_ref["c0"])/ttau_ref
        tr = P["c0_rms"][f, 0]/g["ttau"]
        qr = np.interp(r_, rr, REF["c0"]["t_rms"])/ttau_ref
        axis = r_ < 0.05; wall = (0.5 - r_) / DELTA_V < 50
        print(f"   {tag:7} q_w {g['qw']:.5f}  theta_tau {g['ttau']:.4f}  "
              f"<theta>+ axis {100*(tm-qm)[axis].mean()/qm[axis].mean():+5.1f}%  "
              f"theta'+ axis {100*(tr-qr)[axis].mean()/qr[axis].mean():+5.1f}%  "
              f"theta'+ y+<50 {100*np.mean((tr-qr)[wall]/qr[wall]):+5.1f}%")
    print()

    # THE COMPARISON RADIUS MATTERS. The outermost fluid bin is a CUT CELL on
    # either grid, where what is stored is a penalization blend over a cell
    # straddling the wall -- comparing that against a body-fitted DNS is not
    # like-for-like, and it is also where the profiles are steepest. Quoted
    # just outside them (y+ ~ 7), every case agrees to 0.2-1 %; the cut-cell
    # row is kept below to show what that caveat is worth.
    print("== TABLE: interface signature, theta'/theta'_c0")
    for r_if in (0.48, 0.49, 0.4977):
        print(f"   at r = {r_if:.4f}  (y+ = {(0.5-r_if)/DELTA_V:.2f})"
              + ("   <- quoted" if r_if == 0.48 else
                 "   <- cut cells" if r_if > 0.497 else ""))
        ref_if = {nm: np.interp(r_if, rr, REF[nm]["t_rms"]) for nm in NAMES}
        for tag, g in GRID.items():
            Pg = g["P"]; fg = (Pg["count"][:, 0] > 0) & (Pg["r"] < 0.5)
            i = np.flatnonzero(fg)[np.argmin(np.abs(Pg["r"][fg] - r_if))]
            b = Pg["c0_rms"][i, 0]
            print(f"      {tag:7} " + "  ".join(f"{n} {Pg[f'{n}_rms'][i,0]/b:5.3f}"
                                               for n in NAMES))
        print("      ref     " + "  ".join(f"{n} {ref_if[n]/ref_if['c0']:5.3f}"
                                           for n in NAMES))
    print()

    print("== TABLE: through the solid, theta'/theta'(R), c0")
    for tag, g in GRID.items():
        P = g["P"]; sm = (P["count"][:, 1] > 0) & (P["r"] > 0.5) & (P["r"] < 0.5+0.095)
        i0 = np.flatnonzero(sm)[0]
        xs = (P["r"][sm]-0.5)/0.1
        vals = P["c0_rms"][sm, 1]/P["c0_rms"][i0, 1]
        sel = [np.argmin(np.abs(xs-t)) for t in (0.1, 0.2, 0.3, 0.5, 0.7, 0.9)]
        print(f"   {tag:7} " + "  ".join(f"{xs[i]:.2f}:{vals[i]:.3f}" for i in sel))
    rs = REF["c0"]["r"]; ms = (rs > 0.5) & (rs < 0.6)
    iface = np.interp(0.5001, rs, REF["c0"]["t_rms"])
    xs = (rs[ms]-0.5)/0.1; vals = REF["c0"]["t_rms"][ms]/iface
    sel = [np.argmin(np.abs(xs-t)) for t in (0.1, 0.2, 0.3, 0.5, 0.7, 0.9)]
    print("   ref     " + "  ".join(f"{xs[i]:.2f}:{vals[i]:.3f}" for i in sel))
    return 0


if __name__ == "__main__":
    sys.exit(main())
