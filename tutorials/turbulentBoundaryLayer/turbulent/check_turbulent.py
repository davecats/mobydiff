#!/usr/bin/env python3
"""Turbulent ZPG boundary-layer diagnostics from a snapshot (z-averaged).

  check_turbulent.py <field.h5> [--plot out.png] [--nu-re 450]

Reports, vs streamwise x: momentum-thickness Reynolds number Re_theta(x),
the shape factor H, and the skin friction c_f(x) = 2 nu (dU/dy)_wall / U_e^2,
compared with the turbulent correlation c_f = 0.024 Re_theta^-0.25 (and the
laminar Blasius c_f = 0.664/sqrt(Re_x-ish) upstream of the trip). At a chosen
station it plots the mean U+ vs y+ against the log law U+ = ln(y+)/0.41 + 5.

NOTE: a single snapshot is only z-averaged, not time-averaged, so the
turbulent region is instantaneous -- statistics converge only over many
snapshots (a production run). This is a transition/were-we-turbulent check,
not a converged-statistics tool.
"""
import argparse
import os
import sys

import h5py
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "tools"))
from compare_fields import load_field  # noqa: E402

KAPPA, B = 0.41, 5.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("h5")
    ap.add_argument("--plot", default=None)
    ap.add_argument("--station", type=float, default=0.7, help="x/lx for the U+ profile")
    a = ap.parse_args()

    with h5py.File(a.h5, "r") as f:
        u = load_field(f, "un").mean(axis=0)   # (ny, nx) z-averaged
        w = load_field(f, "wn")
        xn, yn = f["x"][...], f["y"][...]
        re = float(f.attrs["re"])
        t = float(f.attrs.get("t_current", np.nan))
    nu = 1.0 / re
    yc = 0.5 * (yn[:-1] + yn[1:])
    dy = np.diff(yn)
    nx = u.shape[1]

    # w rms(x): a turbulence indicator (0 in the laminar/2D region)
    wrms = np.sqrt((w ** 2).mean(axis=(0, 1)))

    xs, reth, Hs, cf = [], [], [], []
    for i in range(nx):
        up = u[:, i]
        ue = up[-1]
        if ue <= 0:
            continue
        th = np.sum((up / ue) * (1 - up / ue) * dy)
        ds = np.sum((1 - up / ue) * dy)
        tauw = nu * up[0] / yc[0]                 # du/dy at the wall (first cell)
        xs.append(xn[i]); reth.append(ue * th / nu)
        Hs.append(ds / th if th > 0 else np.nan)
        cf.append(2 * tauw / ue ** 2)
    xs, reth, Hs, cf = map(np.array, (xs, reth, Hs, cf))
    cf_corr = 0.024 * np.maximum(reth, 1) ** -0.25

    print(f"snapshot t = {t:.1f}   Re_delta*,0 = {re:.0f}")
    print(f"{'x':>7} {'Re_theta':>9} {'H':>6} {'c_f':>9} {'c_f corr':>9} {'w_rms':>9}")
    for frac in (0.1, 0.3, 0.5, 0.7, 0.9):
        i = int(frac * (len(xs) - 1))
        print(f"{xs[i]:7.1f} {reth[i]:9.1f} {Hs[i]:6.3f} {cf[i]:9.3e} "
              f"{cf_corr[i]:9.3e} {wrms[int(frac*(nx-1))]:9.2e}")

    # turbulent onset: first x where w_rms exceeds 1% of U_inf
    turb = np.where(wrms > 0.01)[0]
    if len(turb):
        print(f"turbulent onset (w_rms>0.01): x ~ {xn[turb[0]]:.1f}, "
              f"fully developed by x ~ {xn[turb[len(turb)//2]]:.1f}")
    else:
        print("no turbulence detected (w_rms < 0.01 everywhere)")

    if a.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 3, figsize=(15, 4))
        ax[0].plot(xs, cf, label="measured (z-avg)")
        ax[0].plot(xs, cf_corr, "k--", label=r"$0.024\,Re_\theta^{-1/4}$")
        ax[0].set_xlabel("x"); ax[0].set_ylabel(r"$c_f$"); ax[0].set_ylim(0, 0.01)
        ax[0].legend(); ax[0].set_title("skin friction")
        ax[1].plot(xn[:len(wrms)], wrms)
        ax[1].set_xlabel("x"); ax[1].set_ylabel(r"$w_{rms}$")
        ax[1].set_title("spanwise fluctuation (transition marker)")
        # U+ vs y+ at the station
        i = int(a.station * (nx - 1))
        up = u[:, i]; ue = up[-1]
        tauw = nu * up[0] / yc[0]; utau = np.sqrt(abs(tauw))
        yp = yc * utau / nu; Up = up / utau
        ax[2].semilogx(yp, Up, ".", ms=4, label=f"x={xn[i]:.0f}")
        yl = np.logspace(0, np.log10(yp.max()), 50)
        ax[2].semilogx(yl, np.log(yl) / KAPPA + B, "k--", label="log law")
        ax[2].semilogx(yl[yl < 12], yl[yl < 12], "k:", label=r"$U^+=y^+$")
        ax[2].set_xlabel(r"$y^+$"); ax[2].set_ylabel(r"$U^+$")
        ax[2].set_xlim(1, None); ax[2].legend(); ax[2].set_title("mean profile")
        fig.tight_layout(); fig.savefig(a.plot, dpi=150)
        print("wrote", a.plot)


if __name__ == "__main__":
    main()
