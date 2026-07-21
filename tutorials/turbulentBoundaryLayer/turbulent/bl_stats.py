#!/usr/bin/env python3
"""Post-process boundaryLayer statistics (bl_stats.h5): (x,y) profiles
averaged in span + time.

  bl_stats.py bl_stats.h5 [--plot out.png] [--retheta 450]

Six-panel figure:
  top:    c_f(Re_theta) vs 0.024 Re_theta^-1/4;  shape factor H(Re_theta)
  middle: Alfredsson diagnostic plot u'/U vs U/U_inf (with the linear fit
          u'/U = 0.286 - 0.255 U/U_inf) and the mean U+(y+), both at
          Re_theta = --retheta
  bottom: Reynolds stresses at Re_theta = --retheta; Clauser pressure-gradient
          parameter beta(Re_theta) = (delta*/tau_w) dp_e/dx = -(delta*/tau_w)
          U_e dU_e/dx (~0 for a ZPG layer)

Stats layout (see boundarylayer_stats.f90): profile[nx*ny, nstat] flattened
y-fastest; nstat = [u,v,w,uu,vv,ww,uv,uw,vw,p]. These are the time+span
means, so unlike a single snapshot they converge to true statistics.
"""
import argparse
import sys

import h5py
import numpy as np

U, V, W, UU, VV, WW, UV, UW, VW, P = range(10)
KAPPA, B = 0.41, 5.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("h5")
    ap.add_argument("--plot", default=None)
    ap.add_argument("--retheta", type=float, default=450.0,
                    help="Re_theta of the station for the profile panels")
    a = ap.parse_args()

    with h5py.File(a.h5, "r") as f:
        nx, ny = int(f.attrs["nx"]), int(f.attrs["ny"])
        re = float(f.attrs["re"])
        prof = f["profile"][...].reshape(nx, ny, -1)   # (x, y, stat)
        x = f["xcoord"][...]
        y = f["ycoord"][...]
        t = float(f.attrs["t_current"])
    nu = 1.0 / re
    # cell widths in y (for integrals); y are cell centres
    yf = np.empty(ny + 1)
    yf[1:-1] = 0.5 * (y[:-1] + y[1:])
    yf[0] = 0.0
    yf[-1] = y[-1] + (y[-1] - yf[-2])
    dy = np.diff(yf)

    Um = prof[:, :, U]           # (x, y) mean streamwise
    Ue = Um[:, -1]
    theta = np.sum((Um / Ue[:, None]) * (1 - Um / Ue[:, None]) * dy, axis=1)
    dstar = np.sum((1 - Um / Ue[:, None]) * dy, axis=1)
    H = dstar / theta
    reth = Ue * theta / nu
    tauw = nu * Um[:, 0] / y[0]
    cf = 2 * tauw / Ue ** 2
    cf_corr = 0.024 * np.maximum(reth, 1) ** -0.25

    print(f"stats t = {t:.1f}  Re_delta*,0 = {re:.0f}  (nx={nx}, ny={ny})")
    print(f"{'x':>7} {'Re_theta':>9} {'H':>6} {'c_f':>10} {'c_f corr':>10}")
    for frac in (0.2, 0.4, 0.6, 0.8):
        i = int(frac * (nx - 1))
        print(f"{x[i]:7.1f} {reth[i]:9.1f} {H[i]:6.3f} {cf[i]:10.3e} {cf_corr[i]:10.3e}")

    # Clauser pressure-gradient parameter beta = (delta*/tau_w) dp_e/dx, with
    # dp_e/dx = -U_e dU_e/dx (rho = 1). U_e(x) is noisy at short averaging, so
    # smooth it before differentiating; ~0 for a ZPG layer.
    Ue_s = np.convolve(Ue, np.ones(21) / 21, mode="same")
    dUedx = np.gradient(Ue_s, x)
    beta = -(dstar / np.maximum(tauw, 1e-30)) * Ue * dUedx

    # station nearest the requested Re_theta (only where the BL is turbulent)
    i0 = np.searchsorted(reth, 250)                       # skip transition
    i = i0 + int(np.argmin(np.abs(reth[i0:] - a.retheta)))
    utau = np.sqrt(abs(tauw[i]))
    yp = y * utau / nu
    Up = Um[i] / utau
    print(f"profile station: x={x[i]:.0f}, Re_theta={reth[i]:.0f} (target {a.retheta:.0f})")

    if a.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(3, 2, figsize=(13, 13))
        rt = f"$Re_\\theta={reth[i]:.0f}$"

        m = reth > 250          # turbulent range (skip transition)

        # (0,0) c_f(Re_theta): zoom on the turbulent band
        ax[0, 0].plot(reth[m], cf[m], label="measured")
        ax[0, 0].plot(reth[m], cf_corr[m], "k--", label=r"$0.024\,Re_\theta^{-1/4}$")
        ax[0, 0].set_xlabel(r"$Re_\theta$"); ax[0, 0].set_ylabel(r"$c_f$")
        ax[0, 0].set_ylim(0.003, 0.006); ax[0, 0].legend(); ax[0, 0].set_title("skin friction")

        # (0,1) H(Re_theta): zoom on the turbulent plateau
        ax[0, 1].plot(reth[m], H[m])
        ax[0, 1].axhline(1.4, ls=":", c="grey", label="H→1.4 (high-Re TBL)")
        ax[0, 1].set_xlabel(r"$Re_\theta$"); ax[0, 1].set_ylabel("H")
        ax[0, 1].set_ylim(1.35, 1.75); ax[0, 1].legend(); ax[0, 1].set_title("shape factor")

        # (1,0) Alfredsson diagnostic plot: u'/U vs U/U_inf, with the linear fit
        Uinf = Ue[i]
        urms = np.sqrt(np.maximum(prof[i, :, UU] - prof[i, :, U] ** 2, 0.0))
        good = Um[i] > 0.15 * Uinf                       # drop the near-wall U->0 divergence
        ax[1, 0].plot(Um[i][good] / Uinf, urms[good] / Um[i][good], ".", ms=4, label="DNS")
        xl = np.linspace(0.15, 1.0, 50)
        ax[1, 0].plot(xl, 0.286 - 0.255 * xl, "k--",
                      label=r"$0.286-0.255\,U/U_\infty$")
        ax[1, 0].set_xlabel(r"$U/U_\infty$"); ax[1, 0].set_ylabel(r"$u'/U$")
        ax[1, 0].set_xlim(0.15, 1.0); ax[1, 0].set_ylim(0, 0.3); ax[1, 0].legend()
        ax[1, 0].set_title(f"diagnostic plot ({rt})")

        # (1,1) mean profile U+(y+)
        ax[1, 1].semilogx(yp, Up, ".", ms=4, label="DNS")
        yl = np.logspace(0, np.log10(yp.max()), 50)
        ax[1, 1].semilogx(yl, np.log(yl) / KAPPA + B, "k--", label="log law")
        ax[1, 1].semilogx(yl[yl < 12], yl[yl < 12], "k:", label=r"$U^+=y^+$")
        ax[1, 1].set_xlabel(r"$y^+$"); ax[1, 1].set_ylabel(r"$U^+$")
        ax[1, 1].set_xlim(1, None); ax[1, 1].set_ylim(0, 25); ax[1, 1].legend()
        ax[1, 1].set_title(f"mean profile ({rt})")

        # (2,0) Reynolds stresses (fluctuations, wall units)
        uu = prof[i, :, UU] - prof[i, :, U] ** 2
        vv = prof[i, :, VV] - prof[i, :, V] ** 2
        ww = prof[i, :, WW] - prof[i, :, W] ** 2
        uv = prof[i, :, UV] - prof[i, :, U] * prof[i, :, V]
        ax[2, 0].plot(yp, np.sqrt(np.maximum(uu, 0)) / utau, label=r"$u'_{rms}$")
        ax[2, 0].plot(yp, np.sqrt(np.maximum(vv, 0)) / utau, label=r"$v'_{rms}$")
        ax[2, 0].plot(yp, np.sqrt(np.maximum(ww, 0)) / utau, label=r"$w'_{rms}$")
        ax[2, 0].plot(yp, -uv / utau ** 2, label=r"$-\overline{u'v'}$")
        ax[2, 0].set_xlabel(r"$y^+$"); ax[2, 0].set_xlim(0, min(yp.max(), 400))
        ax[2, 0].set_ylim(0, 3.0); ax[2, 0].legend()
        ax[2, 0].set_title(f"Reynolds stresses ({rt})")

        # (2,1) Clauser beta(Re_theta); exclude the outlet region where the
        # edge-velocity gradient (hence beta) is corrupted by the outflow BC
        mb = m & (x < 0.9 * x[-1])
        ax[2, 1].plot(reth[mb], beta[mb])
        ax[2, 1].axhline(0.0, ls="--", c="k", label=r"$\beta=0$ (ZPG)")
        ax[2, 1].set_xlabel(r"$Re_\theta$"); ax[2, 1].set_ylabel(r"$\beta$")
        ax[2, 1].set_ylim(-0.3, 0.3); ax[2, 1].legend()
        ax[2, 1].set_title("Clauser pressure-gradient parameter")

        fig.suptitle(f"boundaryLayer statistics — t={t:.0f}, Re_δ*,0={re:.0f}", y=1.0)
        fig.tight_layout(); fig.savefig(a.plot, dpi=140)
        print("wrote", a.plot)


if __name__ == "__main__":
    main()
