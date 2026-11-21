#!/usr/bin/env python3
"""The APG DNS window against the RANS: streamwise development and profiles.

    dns_vs_rans.py                                   (committed window -> assets/dns_vs_rans_*)
    dns_vs_rans.py STATS_LATE STATS_EARLY OUT_PREFIX (any window of the local statistics)

Writes OUT_PREFIX_x.png (beta, c_f, H, Re_theta, U_e, delta_99 against x) and
OUT_PREFIX_prof.png (U+, -u'v'+ and k+ against y+ at x = 300 and 500), plus
OUT_PREFIX_x.txt (the station table).

Series:
  DNS            assets/mobydiff/xyz_3200_384_272/data.nc (t 16 500-43 500), or
                 the window LATE - EARLY of the cumulative statistics
  RANS beta 0.5  ../apg.RANS/production_170182.h5 (table x 1, the RANS target)
  RANS same p    rans_check/check_169456.h5 (table x 1.087 = the DNS's top
                 pressure; coarse grid 800 x 192, c_f ~1.3 % below production)
  ZPG DNS        ../smooth/assets/mobydiff/xyz_3200_384_136/data.nc (grey)
All read with pg_bl.pg_integrals (Griffin edge on the local P0). The RANS
shear stress is nu_t dU/dy, its k the transported k. The RANS fields are
local (regenerated from ../apg.RANS and rans_check/check.ini).
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

INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#898781", "#e1e0d9"
C_DNS, C_RANS, C_SAME, C_ZPG = "#2a78d6", "#eb6834", "#1baf7a", "#a8a79f"


def rans_field(path):
    f = rb.load(path)
    dUdy = np.gradient(f["U"], f["y"], axis=1)
    f["uv"] = -f["nut"] * dUdy
    return f


def style(ax, xlabel, ylabel):
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
    if len(sys.argv) > 1:
        late, early, out = sys.argv[1:4]
        dns = dns_pg.window(late, early)
    else:
        dns = dns_pg.window(os.path.join(HERE, "mobydiff", "xyz_3200_384_272", "data.nc"))
        out = os.path.join(HERE, "dns_vs_rans")
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    series = [
        ("DNS APG, t %.0f-%.0f" % (dns["t0"], dns["t1"]), C_DNS, "-", 1.8, dns),
        ("RANS, beta_RANS = 0.5 (table x 1)", C_RANS, "--", 1.6,
         rans_field(os.path.join(TBL, "apg.RANS", "production_170182.h5"))),
        ("RANS, same p_top as DNS (x 1.087, coarse)", C_SAME, "-.", 1.6,
         rans_field(os.path.join(TBL, "apg.DNS", "rans_check", "check_169456.h5"))),
    ]
    zpg = dns_pg.window(os.path.join(TBL, "smooth", "assets", "mobydiff", "xyz_3200_384_136", "data.nc"))
    q = {lab: pg_bl.pg_integrals(f) for lab, _, _, _, f in series}
    qz = pg_bl.pg_integrals(zpg)

    # ---- streamwise development
    panels = [("beta", r"$\beta$", (-0.1, 0.75)), ("cf", r"$c_f\times10^3$", (2.5, 7.0)),
              ("H", r"$H$", (1.35, 1.75)), ("Re_theta", r"$Re_\theta$", (0, 1300)),
              ("Ue", r"$U_e$", (0.78, 1.01)), ("delta99", r"$\delta_{99}$", (0, 24))]
    fig, axs = plt.subplots(2, 3, figsize=(13, 7), facecolor="#fcfcfb")
    for ax, (key, lab, lim) in zip(axs.flat, panels):
        sc = 1e3 if key == "cf" else 1.0
        ax.set_facecolor("#fcfcfb")
        ax.plot(qz["x"], qz[key] * sc, color=C_ZPG, lw=1.2, ls=":", label="ZPG DNS")
        for name, col, ls, lw, _ in series:
            ax.plot(q[name]["x"], q[name][key] * sc, color=col, ls=ls, lw=lw, label=name)
        if key == "beta":
            ax.axhline(0.5, color=MUTED, lw=0.8)
            ax.axvspan(200, 550, color="#f0efec", zorder=0)
            ax.text(375, 0.70, "read-out range", ha="center", color=INK2, fontsize=8)
        ax.set_xlim(0, 640)
        ax.set_ylim(*lim)
        style(ax, "x", lab)
    h, l = axs[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=4, frameon=False, fontsize=9, labelcolor=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(out + "_x.png", dpi=130)

    # ---- profiles in wall units
    stations = (300.0, 500.0)
    fig, axs = plt.subplots(2, 3, figsize=(13, 7.5), facecolor="#fcfcfb")
    for row, xs in enumerate(stations):
        for name, col, ls, lw, f in series:
            qq = q[name]
            i = int(np.argmin(abs(qq["x"] - xs)))
            ut = np.sqrt(qq["tauw"][i])
            nu = 1.0 / f["re"]
            yp = f["y"] * ut / nu
            axs[row, 0].semilogx(yp, f["U"][i] / ut, color=col, ls=ls, lw=lw, label=name)
            axs[row, 1].semilogx(yp, -f["uv"][i] / ut ** 2, color=col, ls=ls, lw=lw)
            axs[row, 2].semilogx(yp, f["k"][i] / ut ** 2, color=col, ls=ls, lw=lw)
            d99p = qq["delta99"][i] * ut / nu
            for ax in axs[row]:
                ax.axvline(d99p, color=col, lw=0.6, ls=ls, alpha=0.7)
        yy = np.logspace(1.3, 3, 50)
        axs[row, 0].plot(yy, np.log(yy) / 0.41 + 5.0, color=MUTED, lw=0.8)
        for ax, lab in zip(axs[row], (r"$U^+$", r"$-\overline{u'v'}^+$", r"$k^+$")):
            ax.set_facecolor("#fcfcfb")
            ax.set_xlim(0.5, 3000)
            style(ax, r"$y^+$", lab)
            ax.set_title(f"x = {xs:.0f}", color=INK, fontsize=9, loc="left")
    axs[0, 0].set_ylim(0, 30)
    axs[1, 0].set_ylim(0, 30)
    h, l = axs[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=3, frameon=False, fontsize=9, labelcolor=INK)
    fig.text(0.5, 0.005, "thin vertical lines: delta_99 of each series; grey line: 1/0.41 ln y+ + 5",
             ha="center", color=INK2, fontsize=8)
    fig.tight_layout(rect=(0, 0.02, 1, 0.94))
    fig.savefig(out + "_prof.png", dpi=130)

    # ---- station table
    names = [s[0] for s in series]
    with open(out + "_x.txt", "w") as fo:
        fo.write("# x | " + " | ".join(f"{n}: beta cf H Re_theta" for n in ["ZPG DNS"] + names) + "\n")
        for xs in range(100, 651, 50):
            row = [f"{xs:4d}"]
            for qq in [qz] + [q[n] for n in names]:
                i = int(np.argmin(abs(qq["x"] - xs)))
                row.append(f"{qq['beta'][i]:+.3f} {qq['cf'][i]:.5f} {qq['H'][i]:.3f} {qq['Re_theta'][i]:6.0f}")
            fo.write(" | ".join(row) + "\n")
    print(open(out + "_x.txt").read())


if __name__ == "__main__":
    main()
