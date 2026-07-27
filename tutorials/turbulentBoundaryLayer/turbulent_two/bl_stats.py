#!/usr/bin/env python3
"""Post-process boundaryLayer statistics (bl_stats.h5): (x,y) profiles
averaged in span + time.

  bl_stats.py bl_stats.h5 [--plot out.png] [--retheta 450] [--ref passivewall.hdf5]

Six-panel figure:
  top:    c_f(Re_theta) vs 0.024 Re_theta^-1/4;  shape factor H(Re_theta)
  middle: Alfredsson diagnostic plot u'/U vs U/U_inf (with the linear fit
          u'/U = 0.286 - 0.255 U/U_inf) and the mean U+(y+), both at
          Re_theta = --retheta
  bottom: Reynolds stresses at Re_theta = --retheta; Clauser pressure-gradient
          parameter beta(Re_theta) = (delta*/tau_w) dp_e/dx = -(delta*/tau_w)
          U_e dU_e/dx (~0 for a ZPG layer)

Reference overlay (--ref, default passivewall.hdf5): the SIMSON spectral ZPG-TBL
DNS (Schmitt/KIT; Re_delta*,0 = 450, same as this case) carries mean profiles
AND Reynolds stresses over the full streamwise development, so it is overlaid on
every panel it can populate: c_f, H, U+, the diagnostic plot and the Reynolds
stresses. (A legacy .mat with mean profiles only -- tbl_uncontrolled.mat -- is
also accepted; it then populates only c_f/H/U+.) The Clauser-beta panel has no
reference counterpart.

Stats layout (see boundarylayer_stats.f90): profile[nx*ny, nstat] flattened
y-fastest; nstat = [u,v,w,uu,vv,ww,uv,uw,vw,p]. These are the time+span
means, so unlike a single snapshot they converge to true statistics.
"""
import argparse
import os
import sys

import h5py
import numpy as np

U, V, W, UU, VV, WW, UV, UW, VW, P = range(10)
KAPPA, B = 0.41, 5.0


def load_reference(path):
    """Load a spectral ZPG-TBL reference. Returns a dict with Re_theta-sorted
    trends (reth, cf, H) and a station lookup profile(target_reth) giving
    y+, U+ and -- when available -- the Reynolds stresses (urms+, vrms+,
    wrms+, -uv+). Handles the SIMSON passivewall.hdf5 (has stresses, in outer
    units) and the legacy tbl_uncontrolled.mat (mean U+ only). None if absent."""
    if not path or not os.path.exists(path):
        return None

    if path.endswith((".hdf5", ".h5")):
        f = h5py.File(path, "r")
        nu = 1.0 / float(f["parameters"].attrs["re_d1_0"])
        x, y = f["mesh"]["x"][...], f["mesh"]["y"][...]
        u = f["mean"]["u"][...]                         # (nx, ny), outer units
        cov = f["covariance"]
        uu, vv, ww, uv = (cov["uu"][...], cov["vv"][...],
                          cov["ww"][...], cov["uv"][...])
        Ue = u[:, -1]
        theta = np.trapz((u / Ue[:, None]) * (1 - u / Ue[:, None]), y, axis=1)
        dstar = np.trapz(1 - u / Ue[:, None], y, axis=1)
        H = dstar / theta
        reth = Ue * theta / nu
        h1, h2 = y[1], y[2]                              # 2nd-order wall gradient
        dudy = (h2 / (h1 * (h2 - h1))) * u[:, 1] - (h1 / (h2 * (h2 - h1))) * u[:, 2]
        utau = np.sqrt(np.abs(nu * dudy))
        cf = 2 * (nu * dudy) / Ue ** 2
        valid = np.zeros(len(x), bool)                   # drop the SIMSON outflow fringe
        valid[: np.argmax(reth) + 1] = True

        def profile(target):
            idx = np.where(valid)[0]
            j = idx[np.argmin(np.abs(reth[idx] - target))]
            ut = utau[j]
            return dict(reth=reth[j], yp=y * ut / nu, Up=u[j] / ut,
                        urms=np.sqrt(np.maximum(uu[j], 0)) / ut,
                        vrms=np.sqrt(np.maximum(vv[j], 0)) / ut,
                        wrms=np.sqrt(np.maximum(ww[j], 0)) / ut,
                        muv=-uv[j] / ut ** 2, Uratio=u[j] / Ue[j])
        o = np.argsort(reth[valid])
        return dict(name="spectral DNS (SIMSON)", has_stress=True, profile=profile,
                    reth=reth[valid][o], cf=cf[valid][o], H=H[valid][o])

    import scipy.io as sio                               # legacy .mat, mean only
    m = sio.loadmat(path)
    ret = m["re_theta"].ravel()
    Up, yp = m["u_mean"], m["y_plus"]
    Ue = Up[:, -1]
    cf = 2.0 / Ue ** 2
    ff = Up / Ue[:, None]
    H = np.array([np.trapz(1 - ff[i], yp[i]) / np.trapz(ff[i] * (1 - ff[i]), yp[i])
                  for i in range(len(ret))])
    o = np.argsort(ret)

    def profile(target):
        j = int(np.argmin(np.abs(ret - target)))
        return dict(reth=ret[j], yp=yp[j], Up=Up[j], Uratio=Up[j] / Ue[j])
    return dict(name="spectral DNS", has_stress=False, profile=profile,
                reth=ret[o], cf=cf[o], H=H[o])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("h5")
    ap.add_argument("--plot", default=None)
    ap.add_argument("--retheta", type=float, default=450.0,
                    help="Re_theta of the station for the profile panels")
    ap.add_argument("--ref", default="passivewall.hdf5",
                    help="spectral reference to overlay (passivewall.hdf5 with "
                         "stresses, or a mean-only .mat); '' disables")
    a = ap.parse_args()
    ref = load_reference(a.ref) if a.ref else None

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
    rp = ref["profile"](reth[i]) if ref else None         # matched-Re_theta reference profile

    if a.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(3, 2, figsize=(13, 13))
        rt = f"$Re_\\theta={reth[i]:.0f}$"
        RC = "C3"                                         # reference colour

        m = reth > 250          # turbulent range (skip transition)
        if ref is not None:
            sr = (ref["reth"] >= 250) & (ref["reth"] <= reth[m].max())
        rlab = ref["name"] if ref else None
        rlab_rt = (f"{ref['name']} (Re$_\\theta$={rp['reth']:.0f})" if ref else None)

        # (0,0) c_f(Re_theta): zoom on the turbulent band
        ax[0, 0].plot(reth[m], cf[m], label="measured")
        if ref is not None:
            ax[0, 0].plot(ref["reth"][sr], ref["cf"][sr], "-", c=RC, lw=1.3, label=rlab)
        ax[0, 0].plot(reth[m], cf_corr[m], "k--", label=r"$0.024\,Re_\theta^{-1/4}$")
        ax[0, 0].set_xlabel(r"$Re_\theta$"); ax[0, 0].set_ylabel(r"$c_f$")
        ax[0, 0].set_ylim(0.003, 0.006); ax[0, 0].legend(); ax[0, 0].set_title("skin friction")

        # (0,1) H(Re_theta): zoom on the turbulent plateau
        ax[0, 1].plot(reth[m], H[m], label="measured")
        if ref is not None:
            ax[0, 1].plot(ref["reth"][sr], ref["H"][sr], "-", c=RC, lw=1.3, label=rlab)
        ax[0, 1].axhline(1.4, ls=":", c="grey", label="H→1.4 (high-Re TBL)")
        ax[0, 1].set_xlabel(r"$Re_\theta$"); ax[0, 1].set_ylabel("H")
        ax[0, 1].set_ylim(1.35, 1.75); ax[0, 1].legend(); ax[0, 1].set_title("shape factor")

        # (1,0) Alfredsson diagnostic plot: u'/U vs U/U_inf, with the linear fit
        Uinf = Ue[i]
        urms = np.sqrt(np.maximum(prof[i, :, UU] - prof[i, :, U] ** 2, 0.0))
        good = Um[i] > 0.15 * Uinf                       # drop the near-wall U->0 divergence
        ax[1, 0].plot(Um[i][good] / Uinf, urms[good] / Um[i][good], ".", ms=4, label="DNS")
        if rp is not None and "urms" in rp:
            gr = rp["Uratio"] > 0.15
            ax[1, 0].plot(rp["Uratio"][gr], rp["urms"][gr] / (rp["Up"][gr]),
                          "-", c=RC, lw=1.3, label=rlab)
        xl = np.linspace(0.15, 1.0, 50)
        ax[1, 0].plot(xl, 0.286 - 0.255 * xl, "k--",
                      label=r"$0.286-0.255\,U/U_\infty$")
        ax[1, 0].set_xlabel(r"$U/U_\infty$"); ax[1, 0].set_ylabel(r"$u'/U$")
        ax[1, 0].set_xlim(0.15, 1.0); ax[1, 0].set_ylim(0, 0.3); ax[1, 0].legend()
        ax[1, 0].set_title(f"diagnostic plot ({rt})")

        # (1,1) mean profile U+(y+)
        ax[1, 1].semilogx(yp, Up, ".", ms=4, label="DNS")
        if rp is not None:
            ax[1, 1].semilogx(rp["yp"], rp["Up"], "-", c=RC, lw=1.3, label=rlab_rt)
        yl = np.logspace(0, np.log10(yp.max()), 50)
        ax[1, 1].semilogx(yl, np.log(yl) / KAPPA + B, "k--", label="log law")
        ax[1, 1].semilogx(yl[yl < 12], yl[yl < 12], "k:", label=r"$U^+=y^+$")
        ax[1, 1].set_xlabel(r"$y^+$"); ax[1, 1].set_ylabel(r"$U^+$")
        ax[1, 1].set_xlim(1, None); ax[1, 1].set_ylim(0, 25); ax[1, 1].legend()
        ax[1, 1].set_title(f"mean profile ({rt})")

        # (2,0) Reynolds stresses (fluctuations, wall units); reference as lines
        uu = prof[i, :, UU] - prof[i, :, U] ** 2
        vv = prof[i, :, VV] - prof[i, :, V] ** 2
        ww = prof[i, :, WW] - prof[i, :, W] ** 2
        uv = prof[i, :, UV] - prof[i, :, U] * prof[i, :, V]
        for arr, op, c, lbl in [(uu, 1, "C0", r"$u'_{rms}$"), (vv, 1, "C1", r"$v'_{rms}$"),
                                (ww, 1, "C2", r"$w'_{rms}$")]:
            ax[2, 0].plot(yp, np.sqrt(np.maximum(arr, 0)) / utau, c=c, label=lbl)
        ax[2, 0].plot(yp, -uv / utau ** 2, c="C3", label=r"$-\overline{u'v'}$")
        if rp is not None and "urms" in rp:
            for key, c in [("urms", "C0"), ("vrms", "C1"), ("wrms", "C2"), ("muv", "C3")]:
                ax[2, 0].plot(rp["yp"], rp[key], "--", c=c, lw=1.1)
            ax[2, 0].plot([], [], "k--", lw=1.1, label=f"{ref['name']}")
        ax[2, 0].set_xlabel(r"$y^+$"); ax[2, 0].set_xlim(0, min(yp.max(), 400))
        ax[2, 0].set_ylim(0, 3.0); ax[2, 0].legend(ncol=2)
        ax[2, 0].set_title(f"Reynolds stresses ({rt}; dashed = {ref['name'] if ref else 'ref'})")

        # (2,1) Clauser beta(Re_theta); exclude the outlet region where the
        # edge-velocity gradient (hence beta) is corrupted by the outflow BC
        mb = m & (x < 0.9 * x[-1])
        ax[2, 1].plot(reth[mb], beta[mb])
        ax[2, 1].axhline(0.0, ls="--", c="k", label=r"$\beta=0$ (ZPG)")
        ax[2, 1].set_xlabel(r"$Re_\theta$"); ax[2, 1].set_ylabel(r"$\beta$")
        ax[2, 1].set_ylim(-0.001, 0.001); ax[2, 1].legend()
        ax[2, 1].set_title("Clauser pressure-gradient parameter")

        fig.suptitle(f"boundaryLayer statistics — t={t:.0f}, Re_δ*,0={re:.0f}", y=1.0)
        fig.tight_layout(); fig.savefig(a.plot, dpi=140)
        print("wrote", a.plot)


if __name__ == "__main__":
    main()
