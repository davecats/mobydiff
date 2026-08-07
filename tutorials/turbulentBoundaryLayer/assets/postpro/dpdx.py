#!/usr/bin/env python3
"""Streamwise pressure gradient of the boundaryLayer statistics.

  dpdx.py bl_stats.h5 [--out dpdx.png]

<p>_yz(x): the mean pressure averaged over the wall-normal (y) and spanwise
(z) directions (the stats are already span+time averaged, so this is the
extra y-average, dy-weighted over the domain height). Its streamwise
derivative d<p>_yz/dx is the mean streamwise pressure gradient -- the direct
ZPG diagnostic (should be ~0). Also shown against the edge-pressure gradient
implied by U_e(x) via Bernoulli, dp_e/dx = -U_e dU_e/dx.
"""
import argparse

import h5py
import numpy as np

U, P = 0, 9


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("h5")
    ap.add_argument("--out", default="dpdx.png")
    a = ap.parse_args()

    with h5py.File(a.h5, "r") as f:
        nx, ny = int(f.attrs["nx"]), int(f.attrs["ny"])
        re = float(f.attrs["re"]); t = float(f.attrs["t_current"])
        prof = f["profile"][...].reshape(nx, ny, -1)
        x = f["xcoord"][...]; y = f["ycoord"][...]
    nu = 1.0 / re
    yf = np.empty(ny + 1)
    yf[1:-1] = 0.5 * (y[:-1] + y[1:]); yf[0] = 0.0; yf[-1] = y[-1] + (y[-1] - yf[-2])
    dy = np.diff(yf); Ly = yf[-1] - yf[0]

    p = prof[:, :, P]                       # <p>_z(x, y), span+time mean
    p_yz = np.sum(p * dy, axis=1) / Ly      # wall-normal average -> <p>_yz(x)
    dpdx = np.gradient(p_yz, x)

    # edge-pressure gradient from the freestream velocity (Bernoulli)
    Ue = prof[:, -1, U]
    Ue_s = np.convolve(Ue, np.ones(21) / 21, mode="same")
    dpedx = -Ue_s * np.gradient(Ue_s, x)

    # Re_theta(x) for context / turbulent mask
    Um = prof[:, :, U]; Uet = Um[:, -1]
    theta = np.sum((Um / Uet[:, None]) * (1 - Um / Uet[:, None]) * dy, axis=1)
    reth = Uet * theta / nu
    m = (reth > 300) & (x < 0.92 * x[-1])   # turbulent, exclude outflow zone

    print(f"stats t = {t:.0f}")
    print(f"  mean |d<p>_yz/dx| over turbulent range: {np.mean(np.abs(dpdx[m])):.2e}")
    print(f"  peak (transition):                      {np.max(np.abs(dpdx[reth>150])):.2e}")
    print(f"  d<p>_yz/dx developed (Re_theta>500):    {np.mean(dpdx[reth>500]):.2e}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    ax[0].plot(x, p_yz - p_yz[m][0], lw=1.2)
    ax[0].axhline(0, ls=":", c="grey")
    ax[0].set_ylabel(r"$\langle p\rangle_{yz}(x) - $ ref")
    ax[0].set_title(f"wall-normal+spanwise mean pressure — t={t:.0f}")

    ax[1].plot(x[m], dpdx[m], lw=1.2, label=r"$d\langle p\rangle_{yz}/dx$")
    ax[1].plot(x[m], dpedx[m], lw=1.0, c="C1", alpha=0.8,
               label=r"$-U_e\,dU_e/dx$ (edge)")
    ax[1].axhline(0, ls="--", c="k")
    ax[1].set_xlabel("x"); ax[1].set_ylabel(r"$d\langle p\rangle_{yz}/dx$")
    ax[1].set_ylim(-2e-4, 2e-4); ax[1].legend()
    ax[1].set_title("mean streamwise pressure gradient (ZPG diagnostic)")

    # secondary Re_theta axis ticks on top
    axt = ax[0].twiny()
    axt.set_xlim(ax[0].get_xlim())
    xt = [x[np.argmin(np.abs(reth - r))] for r in (400, 600, 800)]
    axt.set_xticks(xt); axt.set_xticklabels([f"{r}" for r in (400, 600, 800)])
    axt.set_xlabel(r"$Re_\theta$")

    fig.tight_layout(); fig.savefig(a.out, dpi=150)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
