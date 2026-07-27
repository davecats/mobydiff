#!/usr/bin/env python3
"""Compare the converged boundaryLayer DNS against the SIMSON spectral
ZPG-TBL reference (passivewall.hdf5, Schmitt/KIT; Re_delta*,0 = 450, same as
this case, trip at x=10). Unlike the earlier .mat reference this one carries
Reynolds stresses and the full streamwise development from the inlet, so the
comparison spans mean flow AND fluctuations.

  compare_passivewall.py <stats.h5> [--ref passivewall.hdf5]
                         [--retheta 677] [--out passivewall_compare.png]

Four panels: c_f(Re_theta), H(Re_theta), mean U+(y+), and the Reynolds
stresses u'/v'/w'_rms+ and -u'v'+ at the matched Re_theta.
"""
import argparse

import h5py
import numpy as np

U, V, W, UU, VV, WW, UV = 0, 1, 2, 3, 4, 5, 6


def load_ref(path):
    f = h5py.File(path, "r")
    nu = 1.0 / float(f["parameters"].attrs["re_d1_0"])
    x, y = f["mesh"]["x"][...], f["mesh"]["y"][...]
    u = f["mean"]["u"][...]
    cov = f["covariance"]
    uu, vv, ww, uv = cov["uu"][...], cov["vv"][...], cov["ww"][...], cov["uv"][...]
    Ue = u[:, -1]
    theta = np.trapz((u / Ue[:, None]) * (1 - u / Ue[:, None]), y, axis=1)
    dstar = np.trapz(1 - u / Ue[:, None], y, axis=1)
    H = dstar / theta
    reth = Ue * theta / nu
    h1, h2 = y[1], y[2]
    dudy = (h2 / (h1 * (h2 - h1))) * u[:, 1] - (h1 / (h2 * (h2 - h1))) * u[:, 2]
    utau = np.sqrt(np.abs(nu * dudy))
    cf = 2 * (nu * dudy) / Ue ** 2
    valid = np.zeros(len(x), bool); valid[: np.argmax(reth) + 1] = True
    return dict(nu=nu, y=y, u=u, uu=uu, vv=vv, ww=ww, uv=uv,
                Ue=Ue, utau=utau, cf=cf, H=H, reth=reth, valid=valid)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stats")
    ap.add_argument("--ref", default="passivewall.hdf5")
    ap.add_argument("--retheta", type=float, default=677.0)
    ap.add_argument("--out", default="passivewall_compare.png")
    a = ap.parse_args()

    R = load_ref(a.ref)
    idx = np.where(R["valid"])[0]
    jr = idx[np.argmin(np.abs(R["reth"][idx] - a.retheta))]
    ur = R["utau"][jr]
    yp_r = R["y"] * ur / R["nu"]
    Up_r = R["u"][jr] / ur
    urms_r = np.sqrt(np.maximum(R["uu"][jr], 0)) / ur
    vrms_r = np.sqrt(np.maximum(R["vv"][jr], 0)) / ur
    wrms_r = np.sqrt(np.maximum(R["ww"][jr], 0)) / ur
    uv_r = -R["uv"][jr] / ur ** 2
    o = np.argsort(R["reth"][R["valid"]])
    reth_rt, cf_rt, H_rt = R["reth"][R["valid"]][o], R["cf"][R["valid"]][o], R["H"][R["valid"]][o]

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
    H = dstar / theta; reth = Ue * theta / nu
    tauw = nu * Um[:, 0] / y[0]; cf = 2 * tauw / Ue ** 2

    i0 = np.searchsorted(reth, 300)
    i = i0 + int(np.argmin(np.abs(reth[i0:] - a.retheta)))
    utau = np.sqrt(abs(tauw[i])); yp = y * utau / nu; Up = Um[i] / utau
    urms = np.sqrt(np.maximum(prof[i, :, UU] - prof[i, :, U] ** 2, 0)) / utau
    vrms = np.sqrt(np.maximum(prof[i, :, VV] - prof[i, :, V] ** 2, 0)) / utau
    wrms = np.sqrt(np.maximum(prof[i, :, WW] - prof[i, :, W] ** 2, 0)) / utau
    uv = (prof[i, :, UV] - prof[i, :, U] * prof[i, :, V]) / utau ** 2

    print(f"                     DNS (this work)   SIMSON spectral")
    print(f"  Re_theta        {reth[i]:14.1f}   {R['reth'][jr]:12.1f}")
    print(f"  c_f             {cf[i]:14.5f}   {R['cf'][jr]:12.5f}   ({100*(cf[i]-R['cf'][jr])/R['cf'][jr]:+.1f}%)")
    print(f"  H               {H[i]:14.3f}   {R['H'][jr]:12.3f}   ({100*(H[i]-R['H'][jr])/R['H'][jr]:+.1f}%)")
    print(f"  u_tau           {utau:14.4f}   {ur:12.4f}")
    print(f"  u'_rms peak     {urms.max():14.3f}   {urms_r.max():12.3f}")
    print(f"  v'_rms peak     {vrms.max():14.3f}   {vrms_r.max():12.3f}")
    print(f"  w'_rms peak     {wrms.max():14.3f}   {wrms_r.max():12.3f}")
    print(f"  -u'v' peak      {(-uv).max():14.3f}   {uv_r.max():12.3f}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 2, figsize=(13, 10))

    m = reth > 150          # include the transition/overshoot
    sr = (reth_rt >= 150) & (reth_rt <= reth[m].max())
    ax[0, 0].plot(reth[m], cf[m], label="this work")
    ax[0, 0].plot(reth_rt[sr], cf_rt[sr], "-", c="C3", lw=1.3, label="SIMSON")
    ax[0, 0].plot(reth[m], 0.024 * reth[m] ** -0.25, "k--", label=r"$0.024\,Re_\theta^{-1/4}$")
    ax[0, 0].set_xlabel(r"$Re_\theta$"); ax[0, 0].set_ylabel(r"$c_f$")
    ax[0, 0].set_ylim(0.002, 0.019); ax[0, 0].legend(); ax[0, 0].set_title("skin friction (with transition)")

    ax[0, 1].plot(reth[m], H[m], label="this work")
    ax[0, 1].plot(reth_rt[sr], H_rt[sr], "-", c="C3", lw=1.3, label="SIMSON")
    ax[0, 1].set_xlabel(r"$Re_\theta$"); ax[0, 1].set_ylabel("H")
    ax[0, 1].set_ylim(1.35, 2.7); ax[0, 1].legend(); ax[0, 1].set_title("shape factor (with transition)")

    ax[1, 0].semilogx(yp_r, Up_r, "-", c="C3", lw=1.5, label=f"SIMSON (Re$_\\theta$={R['reth'][jr]:.0f})")
    ax[1, 0].semilogx(yp, Up, "o", ms=3, mfc="none", label=f"this work (Re$_\\theta$={reth[i]:.0f})")
    ax[1, 0].set_xlabel(r"$y^+$"); ax[1, 0].set_ylabel(r"$U^+$")
    ax[1, 0].set_xlim(1, 400); ax[1, 0].set_ylim(0, 22); ax[1, 0].legend(loc="upper left")
    ax[1, 0].set_title("mean velocity")

    for r, d, c, lbl in [(urms_r, urms, "C0", r"$u'_{rms}$"), (wrms_r, wrms, "C2", r"$w'_{rms}$"),
                         (vrms_r, vrms, "C1", r"$v'_{rms}$"), (uv_r, -uv, "C3", r"$-\overline{u'v'}$")]:
        ax[1, 1].plot(yp_r, r, "-", c=c, lw=1.5)
        ax[1, 1].plot(yp, d, "o", c=c, ms=3, mfc="none", label=lbl)
    ax[1, 1].plot([], [], "k-", lw=1.5, label="SIMSON (lines)")
    ax[1, 1].set_xlabel(r"$y^+$"); ax[1, 1].set_xlim(0, 300); ax[1, 1].set_ylim(0, 3)
    ax[1, 1].legend(ncol=2); ax[1, 1].set_title("Reynolds stresses (symbols: DNS, lines: SIMSON)")

    fig.suptitle(f"boundaryLayer DNS vs SIMSON spectral (passivewall) — Re$_\\theta$≈{a.retheta:.0f}, t={t:.0f}")
    fig.tight_layout(); fig.savefig(a.out, dpi=150)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
