#!/usr/bin/env python3
"""Post-process boundaryLayer statistics (bl_stats.h5): (x,y) profiles
averaged in span + time.

  bl_stats.py bl_stats.h5 [--plot out.png] [--station 0.6]

Reports, vs streamwise x: Re_theta(x), shape factor H(x), skin friction
c_f(x) = 2 nu (dU/dy)_wall / U_e^2 vs the turbulent correlation
0.024 Re_theta^-1/4; and at a chosen station the mean U+(y+) vs the log law
plus the Reynolds stresses u'v', u'_rms, v'_rms, w'_rms.

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
    ap.add_argument("--station", type=float, default=0.6)
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

    if a.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 3, figsize=(15, 4))
        ax[0].plot(reth, cf, label="measured")
        ax[0].plot(reth, cf_corr, "k--", label=r"$0.024\,Re_\theta^{-1/4}$")
        ax[0].set_xlabel(r"$Re_\theta$"); ax[0].set_ylabel(r"$c_f$")
        ax[0].set_ylim(0, 0.008); ax[0].legend(); ax[0].set_title("skin friction")

        i = int(a.station * (nx - 1))
        utau = np.sqrt(abs(tauw[i]))
        yp = y * utau / nu
        ax[1].semilogx(yp, Um[i] / utau, ".", ms=4, label=f"x={x[i]:.0f}")
        yl = np.logspace(0, np.log10(yp.max()), 50)
        ax[1].semilogx(yl, np.log(yl) / KAPPA + B, "k--", label="log law")
        ax[1].semilogx(yl[yl < 12], yl[yl < 12], "k:", label=r"$U^+=y^+$")
        ax[1].set_xlabel(r"$y^+$"); ax[1].set_ylabel(r"$U^+$")
        ax[1].set_xlim(1, None); ax[1].legend(); ax[1].set_title("mean profile")

        # Reynolds stresses (fluctuation = <ab> - <a><b>), in wall units
        uu = prof[i, :, UU] - prof[i, :, U] ** 2
        vv = prof[i, :, VV] - prof[i, :, V] ** 2
        ww = prof[i, :, WW] - prof[i, :, W] ** 2
        uv = prof[i, :, UV] - prof[i, :, U] * prof[i, :, V]
        ax[2].plot(yp, np.sqrt(np.maximum(uu, 0)) / utau, label=r"$u'_{rms}$")
        ax[2].plot(yp, np.sqrt(np.maximum(vv, 0)) / utau, label=r"$v'_{rms}$")
        ax[2].plot(yp, np.sqrt(np.maximum(ww, 0)) / utau, label=r"$w'_{rms}$")
        ax[2].plot(yp, -uv / utau ** 2, label=r"$-u'v'$")
        ax[2].set_xlabel(r"$y^+$"); ax[2].set_xlim(0, min(yp.max(), 300))
        ax[2].legend(); ax[2].set_title("Reynolds stresses")
        fig.tight_layout(); fig.savefig(a.plot, dpi=150)
        print("wrote", a.plot)


if __name__ == "__main__":
    main()
