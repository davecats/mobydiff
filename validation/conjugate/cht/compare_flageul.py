#!/usr/bin/env python3
"""Detailed comparison against Flageul et al.'s RAW data (not the digitised plot).

    ./compare_flageul.py [--stats flageul_stats.h5] [--out flageul_detail.png]

Their repository (repo.ijs.si/CFLAG/incompact3d, fetched by fetch_flageul.py)
publishes far more than the one curve of figure 5: mean and rms velocity, both
turbulent heat fluxes, the temperature variance in the FLUID and in the SOLID,
and a sweep of NINE (G, alpha) cases. That turns a single-curve comparison into
a multi-quantity one, and it removes the digitisation error entirely.

HOW WRONG THE DIGITISATION WAS, measured: rms 0.15, bias +0.13, and the peak
read as 6.208 where the truth is 5.941 (+4.5 %). It was systematically HIGH,
so every "ours is low" statement in the earlier campaigns was overstated by
about that much. The .dat files stay in the tree for provenance, but nothing
should be compared against them now.
"""

from __future__ import annotations

import argparse
import os

import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

NS = 7
S, SS, US, CLO, JLO, CHI, JHI = range(NS)
RE, PR = 149.0, 0.71
YLO, YHI = 1.0, 3.0
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "flageul_data")

# our sweep, and their nine cases (kappa_s, C_s) -> K
OURS = [("k1", 1.0, 1.0), ("k2", 1.0, 100.0), ("k3", 1.0, 1.0e4),
        ("a1", 0.1, 0.1), ("a2", 10.0, 10.0), ("a3", 100.0, 100.0)]
THEIRS = ["g1a1", "g05a05", "g05a1", "g05a2", "g1a05", "g1a2",
          "g2a05", "g2a1", "g2a2"]

C_US, C_REF = "#2a78d6", "#eb6834"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#8b8a85"
SURFACE, GRID = "#fcfcfb", "#e4e3df"


def their_params(name):
    g, y = name[1:].split("a")
    val = lambda s: float(s) / 10 ** (len(s) - 1) if s.startswith("0") else float(s)
    G, Y = val(g), val(y)
    ks, als = 1.0 / Y, 1.0 / G
    return ks, als, ks / als, (ks * ks / als) ** 0.5


def load_ours(path, isc=0):
    with h5py.File(path, "r") as f:
        P, y = f["profile"][...], f["coord"][...]
    m = (y > YLO) & (y < YHI)
    i0 = int(np.argmax(m))
    i1 = int(len(m) - 1 - np.argmax(m[::-1]))
    tt = 0.5 * (abs(P[i0, NS * isc + JLO]) + abs(P[i1, NS * isc + JHI]))
    mean = P[:, NS * isc + S]
    var = np.maximum(P[:, NS * isc + SS] - mean ** 2, 0.0) / tt ** 2
    low = m & (y < 0.5 * (YLO + YHI))
    sol = y < YLO
    return dict(yp=(y[low] - YLO) * RE, var=var[low],
                th=(np.abs(mean - mean[i0]) / tt)[low],
                ys=(y[sol] - YLO) * RE, vs=var[sol],
                ths=((mean[sol] - mean[sol][0]) / tt), P=P, low=low, tt=tt, y=y)


def stat(o_y, o, r_y, r, lo=None, hi=None):
    lo = max(o_y[0], r_y[r_y > 0].min() if lo is None else lo)
    hi = min(o_y[-1], r_y.max() if hi is None else hi)
    g = np.logspace(np.log10(lo), np.log10(hi), 300)
    d = np.interp(g, o_y, o) - np.interp(g, r_y, r)
    return np.sqrt((d ** 2).mean()), np.abs(d).max(), d.mean(), g, d


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stats", default=os.path.join(HERE, "flageul_stats.h5"))
    ap.add_argument("--vel", default=os.path.join(HERE, "Fstat_vel.h5"))
    ap.add_argument("--out", default=os.path.join(HERE, "flageul_detail.png"))
    a = ap.parse_args()

    d = load_ours(a.stats)
    ex = np.loadtxt(os.path.join(DATA, "g1a1_fluct1.dat"))
    m1 = np.loadtxt(os.path.join(DATA, "g1a1_moy1.dat"))
    fb = np.loadtxt(os.path.join(DATA, "g1a1_fluctb.dat"))
    with h5py.File(a.vel, "r") as f:
        Pv = f["profile"][...]
    lw = d["low"]
    U = Pv[lw, 0]
    uu = Pv[lw, 3] - Pv[lw, 0] ** 2
    vv = Pv[lw, 4] - Pv[lw, 1] ** 2
    ww = Pv[lw, 5] - Pv[lw, 2] ** 2
    uv = Pv[lw, 6] - Pv[lw, 0] * Pv[lw, 1]
    ut = (d["P"][lw, US] - Pv[lw, 0] * d["P"][lw, S]) / d["tt"]

    rows = [("U+", U, m1[:, 1], m1[:, 0]), ("theta+", d["th"], m1[:, 2], m1[:, 0]),
            ("u'2", uu, ex[:, 1], ex[:, 0]), ("v'2", vv, ex[:, 2], ex[:, 0]),
            ("w'2", ww, ex[:, 3], ex[:, 0]), ("-u'v'", -uv, ex[:, 4], ex[:, 0]),
            ("u'T'", ut, ex[:, 5], ex[:, 0]), ("T'2", d["var"], ex[:, 7], ex[:, 0])]
    print(f"  OURS vs Flageul g1a1 (RAW data), 0.5 < y+ < 137")
    print(f"  {'quantity':<10}{'rms':>9}{'max':>9}{'bias':>9}   % of that quantity's peak")
    for nm, o, r, ry in rows:
        rms, mx, bias, _, _ = stat(d["yp"], o, ry, r)
        print(f"  {nm:<10}{rms:9.4f}{mx:9.4f}{bias:+9.4f}   {100*rms/np.abs(r).max():6.1f} %")

    # ---- the figure -------------------------------------------------------
    fig, axg = plt.subplots(2, 3, figsize=(15.0, 8.0), facecolor=SURFACE)
    ax = axg.ravel()
    for x in ax:
        x.set_facecolor(SURFACE)
        for sp in ("top", "right"):
            x.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            x.spines[sp].set_color(GRID)
        x.tick_params(colors=INK2, labelsize=8.5)
        x.grid(True, color=GRID, lw=0.6, zorder=0)
        x.set_axisbelow(True)
        x.set_xlabel("$y^+$", color=INK2, fontsize=9.5)

    def panel(i, title, series, ylab, logx=True, logy=False):
        A = ax[i]
        for (lab, o, r, ry, c) in series:
            A.plot(ry, r, "-", color=C_REF, lw=3.2, alpha=0.35, zorder=3)
            A.plot(d["yp"], o, "-", color=c, lw=1.8, zorder=4)
            j = int(0.62 * len(o))
            A.annotate(lab, xy=(d["yp"][j], o[j]), xytext=(5, 5),
                       textcoords="offset points", color=c, fontsize=8.5,
                       bbox=dict(boxstyle="round,pad=0.12", fc=SURFACE,
                                 ec="none", alpha=0.8))
        if logx:
            A.set_xscale("log")
        if logy:
            A.set_yscale("log")
        A.set_ylabel(ylab, color=INK2, fontsize=9.5)
        A.set_title(title, color=INK, fontsize=10.5, loc="left", pad=6)

    panel(0, "(a)  mean profiles",
          [("$U^+$", U, m1[:, 1], m1[:, 0], C_US),
           ("$\\theta^+$", d["th"], m1[:, 2], m1[:, 0], "#1baf7a")],
          "$U^+$,  $\\theta^+$")
    panel(1, "(b)  temperature variance",
          [("$\\overline{T'^2}$", d["var"], ex[:, 7], ex[:, 0], C_US)],
          "$\\overline{T'^2}/T_\\tau^2$")
    panel(2, "(c)  velocity stresses",
          [("$\\overline{u'^2}$", uu, ex[:, 1], ex[:, 0], C_US),
           ("$\\overline{w'^2}$", ww, ex[:, 3], ex[:, 0], "#1baf7a"),
           ("$\\overline{v'^2}$", vv, ex[:, 2], ex[:, 0], "#4a3aa7"),
           ("$-\\overline{u'v'}$", -uv, ex[:, 4], ex[:, 0], "#104281")],
          "wall units")
    panel(3, "(d)  streamwise heat flux",
          [("$\\overline{u'T'}$", ut, ex[:, 5], ex[:, 0], C_US)],
          "$\\overline{u'T'}/(u_\\tau T_\\tau)$")

    # (e) the SOLID side -- the one real discrepancy
    E = ax[4]
    E.semilogy(-fb[:, 0], fb[:, 1], "-", color=C_REF, lw=3.2, alpha=0.45, zorder=3)
    E.semilogy(-d["ys"], d["vs"], "-", color=C_US, lw=1.8, zorder=4)
    E.set_xlabel("depth into the solid,  $-y^+$", color=INK2, fontsize=9.5)
    E.set_ylabel("$\\overline{T'^2}/T_\\tau^2$", color=INK2, fontsize=9.5)
    E.set_title("(e)  variance INSIDE the solid  — the discrepancy",
                color=INK, fontsize=10.5, loc="left", pad=6)
    E.annotate("ours decays more slowly and TURNS UP at the\n"
               "insulated outer face; theirs does not.\n"
               f"outer face {d['vs'][0]:.3f} vs {fb[0,1]:.4f} — a factor "
               f"{d['vs'][0]/fb[0,1]:.0f}.\nThe solid MEAN is exact "
               "(slopes agree to 0.05 %),\nso this is not a charging transient.",
               xy=(0.97, 0.55), xycoords="axes fraction", ha="right", va="top",
               color=MUTED, fontsize=8)

    # (f) the SWEEP: interface variance vs effusivity, ours AND their nine cases
    F = ax[5]
    yq = 0.5
    tk, tv = [], []
    for c in THEIRS:
        p = os.path.join(DATA, f"{c}_fluct1.dat")
        if not os.path.exists(p):
            continue
        e = np.loadtxt(p)
        ks, als, cs, K = their_params(c)
        tk.append(K)
        tv.append(float(np.interp(yq, e[:, 0], e[:, 7])))
    ok = [(K, float(np.interp(yq, load_ours(a.stats, i)["yp"],
                              load_ours(a.stats, i)["var"])))
          for i, (nm, ks, cs) in enumerate(OURS) for K in [np.sqrt(ks * cs)]]
    F.loglog(tk, tv, "D", color=C_REF, ms=9, mfc="none", mew=2.0, zorder=5)
    F.loglog(*zip(*sorted(ok)), "-o", color=C_US, lw=2.0, ms=7, zorder=4)
    F.set_xlabel("effusivity  $K=\\sqrt{\\kappa_s C_s}$", color=INK2, fontsize=9.5)
    F.set_ylabel(f"$\\overline{{T'^2}}$ at $y^+={yq}$", color=INK2, fontsize=9.5)
    F.set_title("(f)  the SWEEP: their nine cases vs our six",
                color=INK, fontsize=10.5, loc="left", pad=6)
    # THEIR OWN DATA DOES NOT COLLAPSE ON K. Two of their cases share
    # K = 1.414 (g05a05 and g2a1) and differ by 14 %; the two at K = 0.707 by
    # 6.5 %. So "same effusivity => same interface response" is not exact even
    # in the reference, which is what our own sweep has been reporting as a
    # failing gate. Mark the pairs so the figure says it.
    pair = {}
    for K, v in zip(tk, tv):
        pair.setdefault(round(K, 3), []).append(v)
    for K, vs_ in pair.items():
        if len(vs_) > 1:
            F.plot([K, K], [min(vs_), max(vs_)], "-", color=C_REF, lw=2.5,
                   alpha=0.5, zorder=4)
            F.annotate(f"same $K$,\n{100*(max(vs_)/min(vs_)-1):.0f} % apart",
                       xy=(K, max(vs_)), xytext=(6, 6), textcoords="offset points",
                       color=C_REF, fontsize=7.5, ha="left")
    F.annotate("diamonds: Flageul's nine $(G,\\alpha)$ cases; line: this work.\n"
               "THEIR data does not collapse on $K$ either — the vertical bars\n"
               "are their own same-$K$ pairs.",
               xy=(0.03, 0.05), xycoords="axes fraction", ha="left", va="bottom",
               color=MUTED, fontsize=8)

    fig.suptitle("Conjugate heat transfer, $Re_\\tau=149$, $Pr=0.71$ — against "
                 "Flageul et al.'s RAW data (thick pale orange), not the digitised plot",
                 color=INK, fontsize=12, x=0.008, ha="left", y=0.985)
    fig.tight_layout(rect=(0, 0.005, 1, 0.955))
    fig.savefig(a.out, dpi=160, facecolor=SURFACE)
    print(f"\n  {a.out} written")


if __name__ == "__main__":
    main()
