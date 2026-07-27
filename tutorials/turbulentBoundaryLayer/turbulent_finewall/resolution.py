#!/usr/bin/env python3
"""Grid resolution in wall units vs Re_theta for the boundaryLayer DNS.

  resolution.py <stats.h5> <field.h5> [--out resolution.png]

Uses the converged mean profile (stats) for u_tau(x), Re_theta(x) and
delta99(x), and a field file for the x/y/z node lines (dx, dy(y), dz). Plots
Dx+, Dz+, Dy+_min (wall cell) and Dy+_max (the largest dy cell that lies
INSIDE the boundary-layer thickness delta99) as functions of Re_theta.
"""
import argparse
import sys

import h5py
import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stats")
    ap.add_argument("field")
    ap.add_argument("--out", default="resolution.png")
    a = ap.parse_args()

    with h5py.File(a.stats, "r") as f:
        nx, ny = int(f.attrs["nx"]), int(f.attrs["ny"])
        re = float(f.attrs["re"])
        prof = f["profile"][...].reshape(nx, ny, -1)
        xc = f["xcoord"][...]
        yc = f["ycoord"][...]
    with h5py.File(a.field, "r") as g:
        xnode = g["x"][...]; ynode = g["y"][...]; znode = g["z"][...]
    nu = 1.0 / re

    dx = xnode[1] - xnode[0]                 # uniform
    dz = znode[1] - znode[0]                 # uniform
    dyc = np.diff(ynode)                     # (ny,) cell widths, fixed in x
    yface = ynode                            # (ny+1,)

    Um = prof[:, :, 0]                        # mean streamwise (x, y)
    Ue = Um[:, -1]
    dy = np.diff(np.concatenate(([0.0], 0.5 * (yc[:-1] + yc[1:]), [yc[-1]])))
    theta = np.sum((Um / Ue[:, None]) * (1 - Um / Ue[:, None]) * dy, axis=1)
    reth = Ue * theta / nu
    tauw = nu * Um[:, 0] / yc[0]
    utau = np.sqrt(np.abs(tauw))
    lplus = nu / utau                        # viscous length (x)

    # delta99(x): first y where <U> reaches 0.99 Ue
    d99 = np.array([yc[np.searchsorted(Um[i], 0.99 * Ue[i])]
                    if np.any(Um[i] >= 0.99 * Ue[i]) else yc[-1] for i in range(nx)])

    dxp = dx / lplus
    dzp = dz / lplus
    dyp_min = dyc[0] / lplus                  # wall cell (fixed dy, varying l+)
    # largest dy cell whose face is within delta99
    dyp_max = np.empty(nx)
    for i in range(nx):
        inside = yface[1:] <= d99[i]
        dyp_max[i] = (dyc[inside].max() if inside.any() else dyc[0]) / lplus[i]

    m = reth > 250                            # turbulent range
    print(f"resolution (turbulent range Re_theta {reth[m].min():.0f}-{reth[m].max():.0f}):")
    print(f"  Dx+     = {dxp[m].min():.1f} - {dxp[m].max():.1f}")
    print(f"  Dz+     = {dzp[m].min():.1f} - {dzp[m].max():.1f}")
    print(f"  Dy+_min = {dyp_min[m].min():.2f} - {dyp_min[m].max():.2f}")
    print(f"  Dy+_max = {dyp_max[m].min():.1f} - {dyp_max[m].max():.1f}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8, 5.5))
    ax.plot(reth[m], dxp[m], label=r"$\Delta x^+$")
    ax.plot(reth[m], dzp[m], label=r"$\Delta z^+$")
    ax.plot(reth[m], dyp_max[m], label=r"$\Delta y^+_{max}$ (within $\delta_{99}$)")
    ax.plot(reth[m], dyp_min[m], label=r"$\Delta y^+_{min}$ (wall)")
    ax.set_xlabel(r"$Re_\theta$"); ax.set_ylabel(r"grid spacing in wall units")
    nxg, nzg = len(xnode) - 1, len(znode) - 1
    ax.set_title(f"DNS resolution vs $Re_\\theta$  ({nxg}×{ny}×{nzg}, $Re_{{\\delta^*,0}}$={re:.0f})")
    ax.legend(); ax.grid(alpha=0.3); ax.set_ylim(0, None)
    fig.tight_layout(); fig.savefig(a.out, dpi=150)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
