#!/usr/bin/env python3
"""Compare the converged boundaryLayer DNS statistics against the
Schlatter & Örlü (2010) ZPG turbulent-boundary-layer reference DNS.

  compare_schlatter_orlu.py <stats.h5> [--ref ref_schlatter_orlu_Re670.prof]
                            [--out schlatter_orlu_compare.png]

The reference (KTH database, vel_0670_dns.prof) is at Re_theta = 677; the
comparison station is our DNS profile at the matched Re_theta. Overlays the
mean U+(y+) and the Reynolds stresses u'/v'/w'_rms+ and -u'v'+; prints the
integral quantities (c_f, H, u_tau).

Reference columns: y/d99, y+, U+, urms+, vrms+, wrms+, uv+, prms+, ...
"""
import argparse
import re
import sys

import h5py
import numpy as np

U, V, W, UU, VV, WW, UV = 0, 1, 2, 3, 4, 5, 6


def read_ref(path):
    meta = {}
    with open(path) as f:
        for line in f:
            if not line.startswith("%") or "=" not in line:
                continue
            val = float(line.split("=")[1])
            if r"Re_{\theta}" in line:
                meta["reth"] = val
            elif "H_{12}" in line:
                meta["H"] = val
            elif "c_f" in line:
                meta["cf"] = val
    data = np.loadtxt(path, comments="%")
    return data, meta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stats")
    ap.add_argument("--ref", default="ref_schlatter_orlu_Re670.prof")
    ap.add_argument("--out", default="schlatter_orlu_compare.png")
    a = ap.parse_args()

    ref, meta = read_ref(a.ref)
    yp_r, Up_r = ref[:, 1], ref[:, 2]
    urms_r, vrms_r, wrms_r, uv_r = ref[:, 3], ref[:, 4], ref[:, 5], ref[:, 6]
    reth_ref = meta.get("reth", 677.0)
    cf_ref, H_ref = meta.get("cf", np.nan), meta.get("H", np.nan)

    with h5py.File(a.stats, "r") as f:
        nx, ny = int(f.attrs["nx"]), int(f.attrs["ny"])
        re = float(f.attrs["re"]); t = float(f.attrs["t_current"])
        prof = f["profile"][...].reshape(nx, ny, -1)
        x = f["xcoord"][...]; y = f["ycoord"][...]
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
    i = i0 + int(np.argmin(np.abs(reth[i0:] - reth_ref)))
    utau = np.sqrt(abs(tauw[i]))
    yp = y * utau / nu
    Up = Um[i] / utau
    urms = np.sqrt(np.maximum(prof[i, :, UU] - prof[i, :, U] ** 2, 0)) / utau
    vrms = np.sqrt(np.maximum(prof[i, :, VV] - prof[i, :, V] ** 2, 0)) / utau
    wrms = np.sqrt(np.maximum(prof[i, :, WW] - prof[i, :, W] ** 2, 0)) / utau
    uv = (prof[i, :, UV] - prof[i, :, U] * prof[i, :, V]) / utau ** 2

    print(f"                     DNS (this work)     Schlatter & Örlü")
    print(f"  Re_theta        {reth[i]:14.1f}     {reth_ref:12.1f}")
    print(f"  c_f             {cf[i]:14.5f}     {cf_ref:12.5f}   ({100*(cf[i]-cf_ref)/cf_ref:+.1f}%)")
    print(f"  H               {H[i]:14.3f}     {H_ref:12.3f}   ({100*(H[i]-H_ref)/H_ref:+.1f}%)")
    print(f"  u_tau           {utau:14.4f}     {np.sqrt(cf_ref/2):12.4f}")
    print(f"  u'_rms peak     {urms.max():14.3f}     {urms_r.max():12.3f}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))

    ax[0].semilogx(yp_r, Up_r, "k-", lw=1.5, label="Schlatter & Örlü (Re$_\\theta$=677)")
    ax[0].semilogx(yp, Up, "o", ms=3, mfc="none", label=f"DNS (Re$_\\theta$={reth[i]:.0f})")
    ax[0].set_xlabel(r"$y^+$"); ax[0].set_ylabel(r"$U^+$")
    ax[0].set_xlim(1, 400); ax[0].set_ylim(0, 22); ax[0].legend(); ax[0].set_title("mean velocity")

    for r, d, c, lbl in [(urms_r, urms, "C0", r"$u'_{rms}$"),
                         (wrms_r, wrms, "C2", r"$w'_{rms}$"),
                         (vrms_r, vrms, "C1", r"$v'_{rms}$"),
                         (-uv_r, -uv, "C3", r"$-\overline{u'v'}$")]:
        ax[1].plot(yp_r, r, "-", c=c, lw=1.5)
        ax[1].plot(yp, d, "o", c=c, ms=3, mfc="none", label=lbl)
    ax[1].plot([], [], "k-", lw=1.5, label="S&Ö (lines)")
    ax[1].set_xlabel(r"$y^+$"); ax[1].set_xlim(0, 300); ax[1].set_ylim(0, 3)
    ax[1].legend(ncol=2); ax[1].set_title("Reynolds stresses (symbols: DNS, lines: S&Ö)")

    fig.suptitle(f"boundaryLayer DNS vs Schlatter & Örlü (2010) — Re$_\\theta$≈{reth_ref:.0f}, t={t:.0f}")
    fig.tight_layout(); fig.savefig(a.out, dpi=150)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
