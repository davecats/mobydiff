#!/usr/bin/env python3
"""beta(x) of the APG DNS against the RANS, with the two half windows as the
sampling check -> assets/beta_dns_rans.png.

    beta_dns_rans.py

DNS: the committed window assets/mobydiff/xyz_3200_384_272/data.nc (t 16 500
-43 500). The half windows (t 16 500-29 500 and 29 500-43 500) need the local
cumulative statistics run/stats/stats_apg_{825000,1475000,2175000}.h5; when
those are present their beta is (re)written to assets/beta_halves.txt, else
read from it. RANS (local fields, regenerated from the committed inis):
../apg.RANS/production_170182.h5 (beta_RANS = 0.5, table x 1) and
rans_check/check_169456.h5 (the DNS's top pressure, table x 1.087, coarse).
"""
import os
import re
import sys

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
CASE = os.path.dirname(HERE)
TBL = os.path.dirname(CASE)
sys.path.insert(0, HERE)
import dns_pg  # noqa: E402  (also puts pg_bl / rans_bl on the path)
import pg_bl  # noqa: E402
import rans_bl as rb  # noqa: E402

NC = os.path.join(HERE, "mobydiff", "xyz_3200_384_272", "data.nc")
STATS = os.path.join(CASE, "run", "stats", "stats_apg_%d.h5")
HALVES = os.path.join(HERE, "beta_halves.txt")
INK, INK2, MUTED, GRID, SURF = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#fcfcfb"


def halves():
    """x, beta of the two half windows and their time ranges."""
    steps = (825000, 1475000, 2175000)
    if all(os.path.exists(STATS % s) for s in steps):
        out, cols, ts = None, [], []
        for a, b in zip(steps[:-1], steps[1:]):
            w = dns_pg.window(STATS % b, STATS % a)
            q = pg_bl.pg_integrals(w)
            out = q["x"] if out is None else out
            cols.append(q["beta"])
            ts.append((w["t0"], w["t1"]))
        hdr = "x  beta(t %.0f-%.0f)  beta(t %.0f-%.0f)" % (*ts[0], *ts[1])
        np.savetxt(HALVES, np.column_stack([out, *cols]), fmt="%.6e", header=hdr)
    d = np.loadtxt(HALVES)
    with open(HALVES) as f:
        hdr = f.readline()
    t = [float(v) for pair in re.findall(r"t (\d+)-(\d+)", hdr) for v in pair]
    return d[:, 0], d[:, 1], d[:, 2], t


def main():
    dns = dns_pg.window(NC)
    full = pg_bl.pg_integrals(dns)
    xh, b1, b2, th = halves()
    rans = pg_bl.pg_integrals(rb.load(os.path.join(TBL, "apg.RANS", "production_170182.h5")))
    same = pg_bl.pg_integrals(rb.load(os.path.join(CASE, "rans_check", "check_169456.h5")))

    fig, ax = plt.subplots(figsize=(10, 5.2), facecolor=SURF)
    ax.set_facecolor(SURF)
    ax.axvspan(200, 550, color="#f0efec", zorder=0)
    ax.text(375, 0.715, "read-out range x = 200-550", ha="center", color=INK2, fontsize=8)
    ax.axhspan(0.47, 0.53, color="#e8f0fb", zorder=0)
    ax.axhline(0.5, color=MUTED, lw=0.8)
    ax.plot(xh, b1, color="#2a78d6", lw=0.9, alpha=0.45,
            label="DNS half windows (t %.0f-%.0f, %.0f-%.0f)" % tuple(th))
    ax.plot(xh, b2, color="#2a78d6", lw=0.9, alpha=0.45)
    ax.plot(full["x"], full["beta"], color="#2a78d6", lw=2.0,
            label=f"DNS APG, t {dns['t0']:.0f}-{dns['t1']:.0f}")
    ax.plot(rans["x"], rans["beta"], color="#eb6834", lw=1.6, ls="--", label="RANS, beta_RANS = 0.5 (table x 1)")
    ax.plot(same["x"], same["beta"], color="#1baf7a", lw=1.6, ls="-.", label="RANS, same p_top as DNS (x 1.087, coarse)")
    ax.set_xlim(0, 640)
    ax.set_ylim(-0.05, 0.75)
    ax.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("#c3c2b7")
    ax.tick_params(colors=INK2, labelsize=9)
    ax.set_xlabel("x", color=INK2)
    ax.set_ylabel(r"$\beta = (\delta^*/\tau_w)\, dP_e/dx$", color=INK2)
    ax.legend(loc="lower right", frameon=False, fontsize=9, labelcolor=INK)
    ax.text(212, 0.537, "target 0.5 ± 0.03", color=INK2, fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "beta_dns_rans.png"), dpi=140)

    win = lambda x: (x >= 200) & (x <= 550)  # noqa: E731
    rows = [("DNS full", full["x"], full["beta"]), ("DNS half 1", xh, b1), ("DNS half 2", xh, b2),
            ("RANS 0.5", rans["x"], rans["beta"]), ("RANS same p", same["x"], same["beta"])]
    for lab, x, b in rows:
        bw = b[win(x)]
        print(f"{lab:12} mean {bw.mean():+.3f} min {bw.min():+.3f} max {bw.max():+.3f}  " +
              " ".join(f"{xs}:{b[np.argmin(abs(x - xs))]:+.3f}" for xs in (200, 300, 400, 500, 550)))
    print("max |half - full| on 200-550: %.4f" % max(abs(b1 - full["beta"])[win(xh)].max(),
                                                     abs(b2 - full["beta"])[win(xh)].max()))


if __name__ == "__main__":
    main()
