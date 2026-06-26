#!/usr/bin/env python3
"""Mean streamwise velocity in wall coordinates (U+ vs y+, semilog-x) from a
single channel snapshot, both walls folded, one or more runs.

For this forced channel the wall stress balances the body force, tau_w = forcing*h,
so u_tau = sqrt(forcing_x * h) and y+ = y_wall * Re * u_tau, U+ = U/u_tau. The
2:1 wall-band interface is drawn as a vertical line at its y+.

Usage: python3 tools/channel_loglaw.py OUT.png FIELD1.h5[:LABEL] [FIELD2.h5[:LABEL] ...]
"""
from __future__ import annotations
import sys
import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def level_line(base, lev):
    line = base.copy()
    for _ in range(lev):
        mid = 0.5 * (line[:-1] + line[1:])
        new = np.empty(2 * len(line) - 1)
        new[0::2] = line; new[1::2] = mid
        line = new
    return line


def mean_U(path):
    with h5py.File(path, "r") as f:
        nb = int(f.attrs["block_nb_x"]); yb = f["y"][...]
        Re = float(f.attrs["re"]); fx = float(f.attrs["forcing_x"]); Ly = float(f.attrs["ly"])
        blocks = f["blocks"][...]; U = f["un"][...]
    rows = {}
    for bid, (ox, oy, oz, lev) in enumerate(blocks):
        for jj in range(nb):
            rows.setdefault((int(lev), int(oy) + jj), []).append(U[bid, :, jj, :].ravel())
    lines = {L: level_line(yb, L) for L in set(k[0] for k in rows)}
    recs = []
    for (lev, gj), planes in rows.items():
        yc = 0.5 * (lines[lev][gj] + lines[lev][gj + 1])
        recs.append((yc, np.concatenate(planes).mean()))
    recs.sort()
    y = np.array([r[0] for r in recs]); Um = np.array([r[1] for r in recs])
    h = Ly / 2.0
    utau = np.sqrt(fx * h)
    return y, Um, Re, utau, h


def main():
    out = sys.argv[1]
    runs = [(sp.rsplit(":", 1) if ":" in sp[2:] else (sp, sp)) for sp in sys.argv[2:]]
    fig, ax = plt.subplots(figsize=(8, 6))
    colors = ["tab:blue", "tab:red", "tab:green"]
    yint_plus = None
    for ri, (path, label) in enumerate(runs):
        y, Um, Re, utau, h = mean_U(path)
        c = colors[ri % len(colors)]
        # lower wall: wall distance = y ; upper wall: 2h - y
        for wall, yd, mk in [("lower", y, "o"), ("upper", 2 * h - y, "s")]:
            m = yd > 0
            yp = yd[m] * Re * utau
            Up = Um[m] / utau
            order = np.argsort(yp)
            ax.semilogx(yp[order], Up[order], "-" + mk, ms=3, lw=1.0, color=c,
                        alpha=0.9 if wall == "lower" else 0.5,
                        label=f"{label} ({wall} wall)")
        yint_plus = 0.643 * Re * utau  # interface wall distance in + units
    yp_ref = np.logspace(0, np.log10(180), 50)
    ax.semilogx(yp_ref, yp_ref, "k:", lw=1, label="U+ = y+ (sublayer)")
    ax.semilogx(yp_ref, (1 / 0.41) * np.log(yp_ref) + 5.2, "k--", lw=1,
                label="U+ = ln(y+)/0.41 + 5.2 (log law)")
    if yint_plus:
        ax.axvline(yint_plus, color="grey", lw=1.2, ls="-.", label=f"2:1 interface (y+={yint_plus:.0f})")
    ax.set_xlabel("y+ (wall distance)"); ax.set_ylabel("U+")
    ax.set_xlim(0.8, 200); ax.set_ylim(0, 20)
    ax.grid(alpha=0.3, which="both"); ax.legend(fontsize=8, loc="upper left")
    ax.set_title("Mean velocity, wall coordinates (step 250, single snapshot, both walls)")
    fig.tight_layout(); fig.savefig(out, dpi=120); print("wrote", out)


if __name__ == "__main__":
    main()
