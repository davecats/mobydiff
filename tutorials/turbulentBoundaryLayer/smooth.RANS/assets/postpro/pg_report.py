#!/usr/bin/env python3
"""Report of the pressure-gradient boundary layers against the ZPG one.

    python3 pg_report.py OUT.png TABLE_DIR [ZPG=FIELD.h5] LABEL=FIELD.h5:PTOP.txt:BETA [...]

Each case is LABEL=FIELD:PTOP:BETA (the converged field, its top-pressure
table and the target beta). The reference ZPG RANS drawn in grey is
ZPG=FIELD (default ../../coarse_38117.h5, the coarse grid; pass
../../production_62365.h5 for the production grid); the beta ratio r is
always taken from ../../production_62365.h5. Writes the
figure (beta, p_top/p_wall, V_top, c_f, H, Re_theta) and, per case,
TABLE_DIR/<LABEL>_stations.txt with x, Re_theta, H, c_f, beta, the estimated
DNS beta (r beta, r from beta_ratio.py) and K.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rans_bl as rb  # noqa: E402
import pg_bl  # noqa: E402
import top_optimise as to  # noqa: E402
import beta_ratio  # noqa: E402

TUT = os.path.join(HERE, "..", "..")
COLORS = {"ZPG": "#7a7a76"}
PALETTE = ["#eb6834", "#2a78d6", "#1baf7a", "#6250d6"]


def main():
    out, tdir = sys.argv[1:3]
    cases = []
    zpg = os.path.join(TUT, "coarse_38117.h5")
    for arg in sys.argv[3:]:
        lab, spec = arg.split("=", 1)
        if lab == "ZPG":
            zpg = spec
            continue
        fn, ptop, beta = spec.split(":")
        cases.append((lab, fn, ptop, float(beta)))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    zf = rb.load(zpg)
    zq = pg_bl.pg_integrals(zf)
    rt = beta_ratio.ratio(os.path.join(TUT, "production_62365.h5"))
    d = rb.dns()
    fig, ax = plt.subplots(2, 3, figsize=(17, 9.5))
    ax = ax.ravel()
    k = (zq["x"] > 5) & (zq["x"] <= 640)
    ax[0].plot(zq["x"][k], zq["beta"][k], color=COLORS["ZPG"], lw=1.6, label="ZPG (RANS)")
    ax[3].plot(zq["x"][k], zq["cf"][k], color=COLORS["ZPG"], lw=1.6, label="ZPG (RANS)")
    ax[3].plot(d["x"][d["x"] <= 550], d["cf"][d["x"] <= 550], "k", lw=2.2, label="ZPG (DNS)")
    ax[4].plot(zq["x"][k], zq["H"][k], color=COLORS["ZPG"], lw=1.6)
    ax[5].plot(zq["x"][k], zq["Re_theta"][k], color=COLORS["ZPG"], lw=1.6)
    for c, (lab, fn, ptop, beta) in zip(PALETTE, cases):
        f = rb.load(fn)
        q = pg_bl.pg_integrals(f)
        x = q["x"]
        k = (x > 5) & (x <= 640)
        ax[0].plot(x[k], q["beta"][k], color=c, lw=2, label=f"{lab} (RANS)")
        rr = np.interp(x, rt["x"], rt["r"])
        s = (x >= 200) & (x <= 550)
        ax[0].plot(x[s], (rr * q["beta"])[s], color=c, lw=1.2, ls="--", label=f"{lab}: r·β, DNS estimate")
        ax[0].axhline(beta, color=c, lw=0.8, ls=":")
        pt = np.loadtxt(ptop)
        ax[1].plot(pt[:, 0], pt[:, 1] - pt[0, 1], color=c, lw=1.2, ls=":", label=f"{lab}: p_top − p_top(0)")
        ax[1].plot(x[k], (q["pw"] - q["pw"][np.argmin(abs(x - 30))])[k], color=c, lw=2, label=f"{lab}: p_wall − p_wall(30)")
        xv, vt = to.top_plane_v(fn)
        ax[2].plot(xv, vt, color=c, lw=2, label=lab)
        ax[3].plot(x[k], q["cf"][k], color=c, lw=2, label=lab)
        ax[4].plot(x[k], q["H"][k], color=c, lw=2)
        ax[5].plot(x[k], q["Re_theta"][k], color=c, lw=2)
        with open(os.path.join(tdir, f"{lab}_stations.txt"), "w") as fh:
            fh.write(f"# {lab}: {os.path.relpath(fn, tdir)} (step {f['step']}, t {f['t']:.0f}), target beta {beta:+.2f}\n")
            fh.write("# Griffin edge on the local free-stream P0, delta*/theta to delta_99 (pg_bl.py); beta_DNS_est = r beta with r from "
                     "beta_ratio_zpg.txt\n")
            fh.write("#   x   Re_theta    H      c_f      U_e     beta   beta_DNS_est  K*1e6   V_top\n")
            for xs in range(150, 651, 50):
                i = np.argmin(abs(x - min(xs, 640)))
                fh.write(f"{x[i]:6.0f} {q['Re_theta'][i]:8.1f} {q['H'][i]:6.3f} {q['cf'][i]:.5f} "
                         f"{q['Ue'][i]:.4f} {q['beta'][i]:+.3f} {rr[i] * q['beta'][i]:+.3f}     "
                         f"{q['K'][i] * 1e6:+.3f} {vt[i]:+.5f}\n")
            ss = (x >= 200) & (x <= 550)
            fh.write(f"# x 200-550: beta mean {q['beta'][ss].mean():+.3f} (min {q['beta'][ss].min():+.3f}, "
                     f"max {q['beta'][ss].max():+.3f}); DNS estimate mean {(rr * q['beta'])[ss].mean():+.3f}; "
                     f"min c_f (x 20-640) {q['cf'][(x > 20) & (x < 640)].min():.5f}\n")
    ax[0].set(xlabel="x", ylabel="β = δ*/τ_w dP_e/dx", xlim=(0, 650), title="(a) Clauser parameter")
    ax[1].set(xlabel="x", ylabel="Δp", xlim=(0, 650), title="(b) top (control) and wall pressure")
    ax[2].set(xlabel="x", ylabel="V_top", xlim=(0, 650), title="(c) resulting top velocity")
    ax[3].set(xlabel="x", ylabel="c_f", xlim=(0, 650), ylim=(0.002, 0.009), title="(d) skin friction")
    ax[4].set(xlabel="x", ylabel="H", xlim=(0, 650), ylim=(1.3, 1.9), title="(e) shape factor")
    ax[5].set(xlabel="x", ylabel="Re_θ", xlim=(0, 650), title="(f) momentum-thickness Reynolds number")
    for a in ax:
        a.grid(alpha=0.25, lw=0.6)
    ax[0].legend(fontsize=8, frameon=False)
    ax[1].legend(fontsize=8, frameon=False)
    ax[3].legend(fontsize=8, frameon=False)
    ax[0].axvspan(200, 550, color="#ddd", alpha=0.4, lw=0)
    fig.tight_layout()
    fig.savefig(out, dpi=100)
    print("wrote", out)


if __name__ == "__main__":
    main()
