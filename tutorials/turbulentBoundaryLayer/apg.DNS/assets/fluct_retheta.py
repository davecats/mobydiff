#!/usr/bin/env python3
"""Turbulent fluctuations at matched Re_theta: ZPG DNS, ZPG RANS, FPG RANS.

RANS stresses by Boussinesq: -u'v' = nu_t dU/dy, u_i'u_i' = 2/3 k (the
normal stresses of an eddy-viscosity model are isotropic up to the small
strain terms). Wall units with the local u_tau of each case.

    fluct_retheta.py        -> assets/fluct_retheta.png (+ the peak table on stdout)

ZPG DNS: the committed ../smooth/assets/mobydiff/xyz_3200_384_136/data.nc.
RANS (local fields, regenerated from the committed inis):
../smooth.RANS/production_62365.h5, ../fpg.RANS/production_170012.h5.
"""
import os
import sys

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
TBL = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, HERE)
import dns_pg  # noqa: E402  (also puts pg_bl / rans_bl on the path)
import pg_bl  # noqa: E402
import rans_bl as rb  # noqa: E402

INK, INK2, MUTED, GRID, SURF = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#fcfcfb"
C = {"ZPG DNS": "#2a78d6", "ZPG RANS": "#eb6834", "FPG RANS (beta -0.3)": "#1baf7a"}
LS = {"ZPG DNS": "-", "ZPG RANS": "--", "FPG RANS (beta -0.3)": "-."}
RETH = (450.0, 650.0)
OUT = HERE


def dns():
    return dns_pg.window(os.path.join(TBL, "smooth", "assets", "mobydiff", "xyz_3200_384_136", "data.nc"))


def rans(path):
    f = rb.load(path)
    f["uv"] = -f["nut"] * np.gradient(f["U"], f["y"], axis=1)
    return f


def style(ax, xlabel, ylabel):
    ax.set_facecolor(SURF)
    ax.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("#c3c2b7")
    ax.tick_params(colors=INK2, labelsize=8)
    ax.set_xlabel(xlabel, color=INK2, fontsize=9)
    ax.set_ylabel(ylabel, color=INK2, fontsize=9)


def main():
    cases = {"ZPG DNS": dns(),
             "ZPG RANS": rans(os.path.join(TBL, "smooth.RANS", "production_62365.h5")),
             "FPG RANS (beta -0.3)": rans(os.path.join(TBL, "fpg.RANS", "production_170012.h5"))}
    q = {k: pg_bl.pg_integrals(f) for k, f in cases.items()}

    fig, axs = plt.subplots(2, 4, figsize=(16, 7.6), facecolor=SURF)
    table = []
    for row, rt in enumerate(RETH):
        for name, f in cases.items():
            qq = q[name]
            ok = qq["x"] > 60.0                        # past the trip / k-source region
            i = int(np.flatnonzero(ok)[np.argmin(abs(qq["Re_theta"][ok] - rt))])
            ut = np.sqrt(qq["tauw"][i])
            nu = 1.0 / f["re"]
            yp = f["y"] * ut / nu
            d99p = qq["delta99"][i] * ut / nu
            col, ls = C[name], LS[name]
            k = f["k"][i] / ut ** 2
            uv = -f["uv"][i] / ut ** 2
            axs[row, 0].semilogx(yp, k, color=col, ls=ls, lw=1.7, label=name)
            axs[row, 1].semilogx(yp, uv, color=col, ls=ls, lw=1.7)
            axs[row, 3].plot(f["y"] / qq["delta99"][i], uv, color=col, ls=ls, lw=1.7)
            if name == "ZPG DNS":
                for comp, lw in (("uu", 1.7), ("vv", 1.2), ("ww", 1.2)):
                    r = np.sqrt(np.maximum(f[comp][i], 0)) / ut
                    axs[row, 2].semilogx(yp, r, color=col, ls=ls, lw=lw)
                    j = int(np.argmax(r))
                    axs[row, 2].annotate(f"{comp[0]}'", (yp[j], r[j]), xytext=(4, 4),
                                         textcoords="offset points", color=INK2, fontsize=8)
            else:
                axs[row, 2].semilogx(yp, np.sqrt(2.0 / 3.0 * np.maximum(k, 0)), color=col, ls=ls, lw=1.7)
            for ax in axs[row, :3]:
                ax.axvline(d99p, color=col, ls=ls, lw=0.6, alpha=0.6)
            jk = int(np.argmax(k))
            table.append((rt, name, qq["x"][i], qq["Re_theta"][i], qq["beta"][i], qq["cf"][i],
                          qq["H"][i], d99p, k[jk], yp[jk], uv.max(), yp[int(np.argmax(uv))]))
        labels = (r"$k^+$", r"$-\overline{u'v'}^+$", r"$u'^+_{rms}$ (DNS: u', v', w'; RANS: $\sqrt{2k/3}$)",
                  r"$-\overline{u'v'}^+$")
        for c, ax in enumerate(axs[row]):
            if c < 3:
                ax.set_xlim(0.5, 2000)
                style(ax, r"$y^+$", labels[c])
            else:
                ax.set_xlim(0, 1.3)
                style(ax, r"$y/\delta_{99}$", labels[c])
            ax.set_title(rf"$Re_\theta$ = {rt:.0f}", color=INK, fontsize=9, loc="left")
    h, l = axs[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=3, frameon=False, fontsize=10, labelcolor=INK)
    fig.text(0.5, 0.005, r"thin vertical lines: $\delta_{99}^+$ of each case; local $u_\tau$ of each case",
             ha="center", color=INK2, fontsize=8)
    fig.tight_layout(rect=(0, 0.02, 1, 0.94))
    fig.savefig(os.path.join(OUT, "fluct_retheta.png"), dpi=130)

    hdr = f"{'Re_th':>5} {'case':22} {'x':>5} {'Re_th':>6} {'beta':>6} {'c_f':>8} {'H':>5} {'d99+':>5} {'k+max':>6} {'@y+':>5} {'uv+max':>6} {'@y+':>5}"
    print(hdr)
    for r in table:
        print(f"{r[0]:5.0f} {r[1]:22} {r[2]:5.0f} {r[3]:6.0f} {r[4]:+6.3f} {r[5]:8.5f} {r[6]:5.3f} {r[7]:5.0f} "
              f"{r[8]:6.2f} {r[9]:5.1f} {r[10]:6.3f} {r[11]:5.1f}")


if __name__ == "__main__":
    main()
