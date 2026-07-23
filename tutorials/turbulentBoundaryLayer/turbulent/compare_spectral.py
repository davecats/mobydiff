#!/usr/bin/env python3
"""Compare the converged boundaryLayer DNS statistics against an in-house
pseudo-spectral ZPG turbulent-boundary-layer DNS (tbl_uncontrolled.mat).

  compare_spectral.py <stats.h5> [--mat tbl_uncontrolled.mat]
                      [--retheta 677] [--out spectral_compare.png]

The spectral set stores mean U+(y+) profiles at 2400 streamwise stations
spanning Re_theta ~ 380-2480 (plus u_tau, Re_theta, Re_tau) — no Reynolds
stresses, but a wide Re_theta range. So the comparison is three-fold:
  (a) the mean U+(y+) profile at the matched Re_theta,
  (b) c_f(Re_theta)   over the whole turbulent range,
  (c) H(Re_theta)     over the whole turbulent range,
the last two being trends the single-station Schlatter & Örlü profile cannot
provide. c_f = 2/(U_e+)^2 and H = int(1-U/Ue)dy / int (U/Ue)(1-U/Ue)dy are
taken straight from the spectral U+ profiles (both independent of l+).
"""
import argparse

import h5py
import numpy as np
import scipy.io as sio

U, V, W, UU, VV, WW, UV = 0, 1, 2, 3, 4, 5, 6


def spectral_cf_H(mat):
    """c_f(Re_theta) and H(Re_theta) from the spectral U+ profiles."""
    ret = mat["re_theta"].ravel()
    Up = mat["u_mean"]           # (nstat, ny) = U+  = U/u_tau
    yp = mat["y_plus"]           # (nstat, ny)
    Ue = Up[:, -1]               # edge U+
    cf = 2.0 / Ue ** 2
    f = Up / Ue[:, None]
    H = np.array([np.trapz(1 - f[i], yp[i]) / np.trapz(f[i] * (1 - f[i]), yp[i])
                  for i in range(len(ret))])
    order = np.argsort(ret)
    return ret[order], cf[order], H[order]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stats")
    ap.add_argument("--mat", default="tbl_uncontrolled.mat")
    ap.add_argument("--retheta", type=float, default=677.0)
    ap.add_argument("--out", default="spectral_compare.png")
    a = ap.parse_args()

    mat = sio.loadmat(a.mat)
    ret_s = mat["re_theta"].ravel()
    Up_s, yp_s = mat["u_mean"], mat["y_plus"]
    js = int(np.argmin(np.abs(ret_s - a.retheta)))
    reth_s = ret_s[js]
    ret_cf, cf_s, H_s = spectral_cf_H(mat)

    with h5py.File(a.stats, "r") as f:
        nx, ny = int(f.attrs["nx"]), int(f.attrs["ny"])
        re = float(f.attrs["re"]); t = float(f.attrs["t_current"])
        prof = f["profile"][...].reshape(nx, ny, -1)
        y = f["ycoord"][...]
    nu = 1.0 / re
    dy = np.diff(np.concatenate(([0.0], 0.5 * (y[:-1] + y[1:]), [y[-1]])))
    Um = prof[:, :, U]; Ue = Um[:, -1]
    theta = np.sum((Um / Ue[:, None]) * (1 - Um / Ue[:, None]) * dy, axis=1)
    dstar = np.sum((1 - Um / Ue[:, None]) * dy, axis=1)
    H = dstar / theta
    reth = Ue * theta / nu
    tauw = nu * Um[:, 0] / y[0]
    cf = 2 * tauw / Ue ** 2

    i0 = np.searchsorted(reth, 300)
    i = i0 + int(np.argmin(np.abs(reth[i0:] - a.retheta)))
    utau = np.sqrt(abs(tauw[i]))
    yp = y * utau / nu
    Up = Um[i] / utau

    print(f"                     DNS (this work)     spectral DNS (ours)")
    print(f"  Re_theta        {reth[i]:14.1f}     {reth_s:12.1f}")
    print(f"  c_f             {cf[i]:14.5f}     {2/Up_s[js,-1]**2:12.5f}"
          f"   ({100*(cf[i]-2/Up_s[js,-1]**2)/(2/Up_s[js,-1]**2):+.1f}%)")
    Hs = np.interp(reth[i], ret_cf, H_s)
    print(f"  H               {H[i]:14.3f}     {Hs:12.3f}   ({100*(H[i]-Hs)/Hs:+.1f}%)")
    print(f"  U_e+            {Up[-1]:14.2f}     {Up_s[js,-1]:12.2f}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(17, 5))

    ax[0].semilogx(yp_s[js], Up_s[js], "k-", lw=1.5,
                   label=f"spectral DNS (Re$_\\theta$={reth_s:.0f})")
    ax[0].semilogx(yp, Up, "o", ms=3, mfc="none",
                   label=f"this work (Re$_\\theta$={reth[i]:.0f})")
    ax[0].set_xlabel(r"$y^+$"); ax[0].set_ylabel(r"$U^+$")
    ax[0].set_xlim(1, 400); ax[0].set_ylim(0, 22); ax[0].legend(loc="upper left")
    ax[0].set_title("mean velocity")

    m = reth > 300
    ax[1].plot(ret_cf, cf_s, "k-", lw=1.5, label="spectral DNS")
    ax[1].plot(reth[m], cf[m], "-", c="C0", lw=1.5, label="this work")
    ax[1].plot(reth[m], 0.024 * reth[m] ** -0.25, "--", c="grey",
               label=r"$0.024\,Re_\theta^{-1/4}$")
    ax[1].set_xlabel(r"$Re_\theta$"); ax[1].set_ylabel(r"$c_f$")
    ax[1].set_xlim(300, reth[m].max() * 1.02); ax[1].set_ylim(0.003, 0.006)
    ax[1].legend(); ax[1].set_title("skin friction")

    ax[2].plot(ret_cf, H_s, "k-", lw=1.5, label="spectral DNS")
    ax[2].plot(reth[m], H[m], "-", c="C0", lw=1.5, label="this work")
    ax[2].set_xlabel(r"$Re_\theta$"); ax[2].set_ylabel("H")
    ax[2].set_xlim(300, reth[m].max() * 1.02); ax[2].set_ylim(1.35, 1.75)
    ax[2].legend(); ax[2].set_title("shape factor")

    fig.suptitle(f"boundaryLayer DNS vs in-house spectral DNS — t={t:.0f}")
    fig.tight_layout(); fig.savefig(a.out, dpi=150)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
