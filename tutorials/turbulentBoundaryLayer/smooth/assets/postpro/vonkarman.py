#!/usr/bin/env python3
"""von Karman momentum-integral check for a (nearly) ZPG boundary layer:

    c_f/2 = d(theta)/dx  +  (delta* + 2 theta)/Ue * d(Ue)/dx
            \___term1__/    \_________term2 (pressure gradient)_________/

The identity is exact for a 2D steady incompressible BL, so how well it closes
(residual = c_f/2 - term1 - term2, with c_f the DIRECT wall-gradient value) is a
convergence/consistency check independent of the wall shear. c_f estimated FROM
the balance, c_f_vk = 2(term1+term2), is compared with the direct c_f and SIMSON.

  vonkarman.py [stats.h5] [--ref passivewall.hdf5] [--out vonkarman.png]
"""
import argparse
import h5py
import numpy as np

nu = 1 / 450.


def smooth(a, w=41):
    k = np.ones(w) / w
    return np.convolve(a, k, mode="same")


def integrals_dns(fn):
    f = h5py.File(fn, "r")
    nx, ny = int(f.attrs["nx"]), int(f.attrs["ny"])
    p = f["profile"][...].reshape(nx, ny, -1); x = f["xcoord"][...]; y = f["ycoord"][...]
    dy = np.diff(np.concatenate(([0], 0.5 * (y[:-1] + y[1:]), [y[-1]])))
    U = p[:, :, 0]; Ue = U[:, -1]
    theta = np.sum((U / Ue[:, None]) * (1 - U / Ue[:, None]) * dy, axis=1)
    dstar = np.sum((1 - U / Ue[:, None]) * dy, axis=1)
    cf = 2 * (nu * U[:, 0] / y[0]) / Ue ** 2
    reth = Ue * theta / nu
    return x, theta, dstar, Ue, cf, reth


def integrals_simson(fn):
    f = h5py.File(fn, "r")
    x = f["mesh"]["x"][...]; y = f["mesh"]["y"][...]; u = f["mean"]["u"][...]
    Ue = u[:, -1]
    theta = np.trapz((u / Ue[:, None]) * (1 - u / Ue[:, None]), y, axis=1)
    dstar = np.trapz(1 - u / Ue[:, None], y, axis=1)
    h1, h2 = y[1], y[2]
    du = (h2 / (h1 * (h2 - h1))) * u[:, 1] - (h1 / (h2 * (h2 - h1))) * u[:, 2]
    cf = 2 * nu * du / Ue ** 2
    reth = Ue * theta / nu
    valid = np.zeros(len(x), bool); valid[: np.argmax(reth) + 1] = True
    return x, theta, dstar, Ue, cf, reth, valid


def vk_terms(x, theta, dstar, Ue):
    dthetadx = np.gradient(smooth(theta), x)          # term1
    dUedx = np.gradient(smooth(Ue), x)
    term2 = (dstar + 2 * theta) / Ue * dUedx          # pressure-gradient term
    return dthetadx, term2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stats", nargs="?", default="production_stats.h5")
    ap.add_argument("--ref", default="passivewall.hdf5")
    ap.add_argument("--out", default="vonkarman.png")
    a = ap.parse_args()

    x, th, ds, Ue, cf, reth = integrals_dns(a.stats)
    t1, t2 = vk_terms(x, th, ds, Ue)
    cf_vk = 2 * (t1 + t2)
    resid = cf / 2 - t1 - t2

    xs, ths, dss, Ues, cfs, reths, val = integrals_simson(a.ref)
    t1s, t2s = vk_terms(xs, ths, dss, Ues)
    cf_vks = 2 * (t1s + t2s)
    resids = cfs / 2 - t1s - t2s

    print(f"{'Re_th':>6} | {'term1 dth/dx':>12} {'term2 (PG)':>11} | "
          f"{'cf_vk':>9} {'cf_direct':>9} {'cf_SIMSON':>9} | {'resid/(cf/2)':>12}")
    for tgt in (450, 550, 677, 800):
        i = np.argmin(np.abs(reth - tgt))
        k = np.where(val)[0][np.argmin(np.abs(reths[val] - reth[i]))]
        print(f"{reth[i]:6.0f} | {t1[i]:12.3e} {t2[i]:11.3e} | "
              f"{cf_vk[i]:9.5f} {cf[i]:9.5f} {cfs[k]:9.5f} | {100*resid[i]/(cf[i]/2):11.1f}%")

    print("\n  SIMSON's own von Karman closure (spectral reference):")
    for tgt in (450, 550, 677, 800):
        k = np.where(val)[0][np.argmin(np.abs(reths[val] - tgt))]
        print(f"{reths[k]:6.0f} | term1={t1s[k]:.3e} term2={t2s[k]:.3e} | "
              f"cf_vk={cf_vks[k]:.5f} cf_direct={cfs[k]:.5f} | resid/(cf/2)={100*resids[k]/(cfs[k]/2):.1f}%")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    m = (reth > 400) & (reth < 900) & (x < 0.9 * x[-1])
    fig, ax = plt.subplots(1, 2, figsize=(14, 5.5))
    ax[0].plot(reth[m], cf[m] / 2, "C0", lw=1.6, label=r"$c_f/2$ (direct wall shear)")
    ax[0].plot(reth[m], t1[m], "C1", lw=1.3, label=r"term1 $d\theta/dx$")
    ax[0].plot(reth[m], t2[m], "C2", lw=1.3, label=r"term2 $(\delta^*{+}2\theta)/U_e\,dU_e/dx$")
    ax[0].plot(reth[m], (t1 + t2)[m], "C3--", lw=1.3, label=r"term1+term2 ($c_{f,vk}/2$)")
    ax[0].set_xlabel(r"$Re_\theta$"); ax[0].set_ylabel("momentum budget")
    ax[0].legend(fontsize=8); ax[0].set_title("von Karman terms (case i)")
    mvs = val & (reths > 400) & (reths < 900)
    ax[1].plot(reth[m], cf[m], "C0", lw=1.6, label="c_f direct (case i)")
    ax[1].plot(reth[m], cf_vk[m], "C3--", lw=1.4, label="c_f from von Karman (case i)")
    ax[1].plot(reths[mvs], cfs[mvs], "k", lw=1.6, label="c_f SIMSON")
    ax[1].set_xlabel(r"$Re_\theta$"); ax[1].set_ylabel(r"$c_f$")
    ax[1].set_ylim(0.0042, 0.0052); ax[1].legend(fontsize=8); ax[1].set_title("c_f: direct vs von Karman vs SIMSON")
    fig.tight_layout(); fig.savefig(a.out, dpi=150); print("\nwrote", a.out)


if __name__ == "__main__":
    main()
