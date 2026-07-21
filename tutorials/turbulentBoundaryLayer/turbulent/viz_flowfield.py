#!/usr/bin/env python3
"""Visualize a turbulent-boundary-layer snapshot: the growing layer in a
side (x-y) plane, the near-wall streaks from above (x-z), and the
wall-normal ejections/sweeps.

  viz_flowfield.py <field.h5> [--out flowfield.png] [--ymax 25]
"""
import argparse
import os
import sys

import h5py
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "tools"))
from compare_fields import load_field  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("h5")
    ap.add_argument("--out", default="flowfield.png")
    ap.add_argument("--ymax", type=float, default=25.0)
    a = ap.parse_args()

    with h5py.File(a.h5, "r") as f:
        xn, yn, zn = f["x"][...], f["y"][...], f["z"][...]
        t, re = float(f.attrs["t_current"]), float(f.attrs["re"])
        u = load_field(f, "un")            # (nz, ny, nx)
        v = load_field(f, "vn")
    nu = 1.0 / re
    xc = 0.5 * (xn[:-1] + xn[1:]); yc = 0.5 * (yn[:-1] + yn[1:]); zc = 0.5 * (zn[:-1] + zn[1:])
    nz, ny, nx = u.shape

    # friction velocity at mid-domain (for y+): u_tau = sqrt(nu dU/dy|wall)
    Um = u.mean(axis=0)                      # (ny, nx) spanwise mean
    imid = nx // 2
    utau = np.sqrt(nu * Um[0, imid] / yc[0])
    lplus = nu / utau
    # near-wall plane at y+ ~ 12 (the streak layer)
    jwall = int(np.argmin(np.abs(yc - 12.0 * lplus)))
    yplus_plane = yc[jwall] / lplus

    jmax = np.searchsorted(yn, a.ymax)
    zmid = nz // 2

    # boundary-layer thickness delta99(x): where <u>_z first reaches 0.99 Ue
    Ue = Um[-1, :]
    d99 = np.array([yc[np.searchsorted(Um[:, i], 0.99 * Ue[i])] if np.any(Um[:, i] >= 0.99 * Ue[i])
                    else yc[-1] for i in range(nx)])

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(3, 1, figsize=(15, 10))

    # (A) side view: instantaneous u in the mid-span x-y plane
    pc = ax[0].pcolormesh(xc, yc[:jmax], u[zmid, :jmax, :], cmap="turbo",
                          vmin=0.0, vmax=1.1, shading="auto", rasterized=True)
    ax[0].plot(xc, d99, "w-", lw=0.8, alpha=0.8)
    fig.colorbar(pc, ax=ax[0], pad=0.01, label="u")
    ax[0].set_ylabel("y"); ax[0].set_ylim(0, a.ymax)
    ax[0].set_title(f"streamwise velocity u — mid-span x-y plane (white = δ99); "
                    f"t={t:.0f}, Re_δ*,0={re:.0f}")

    # (B) near-wall streaks: u' = u - <u>_z in an x-z plane at y+ ~ 12
    up = u[:, jwall, :] - Um[jwall, :]
    lim = np.percentile(np.abs(up), 99)
    pc = ax[1].pcolormesh(xc, zc, up, cmap="RdBu_r", vmin=-lim, vmax=lim,
                          shading="auto", rasterized=True)
    fig.colorbar(pc, ax=ax[1], pad=0.01, label="u'")
    ax[1].set_ylabel("z"); ax[1].set_aspect("equal")
    ax[1].set_title(f"near-wall streaks: u' in the x-z plane at y⁺≈{yplus_plane:.0f} "
                    f"(elongated low/high-speed streaks)")

    # (C) wall-normal velocity v (ejections/sweeps) in the mid-span x-y plane
    lim = np.percentile(np.abs(v[zmid, :jmax, :]), 99)
    pc = ax[2].pcolormesh(xc, yc[:jmax], v[zmid, :jmax, :], cmap="RdBu_r",
                          vmin=-lim, vmax=lim, shading="auto", rasterized=True)
    fig.colorbar(pc, ax=ax[2], pad=0.01, label="v")
    ax[2].set_xlabel("x"); ax[2].set_ylabel("y"); ax[2].set_ylim(0, a.ymax)
    ax[2].set_title("wall-normal velocity v — mid-span x-y plane (ejections v>0 / sweeps v<0)")

    fig.tight_layout()
    fig.savefig(a.out, dpi=130)
    print(f"u_tau≈{utau:.4f}, l+≈{lplus:.4f}, streak plane y+≈{yplus_plane:.1f}")
    print("wrote", a.out)


if __name__ == "__main__":
    main()
