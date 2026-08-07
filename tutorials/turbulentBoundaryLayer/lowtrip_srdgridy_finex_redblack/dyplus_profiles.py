#!/usr/bin/env python3
"""Wall-normal grid spacing in wall units, dy+ vs y+, over the WHOLE domain
height at 10 streamwise stations (labelled by Re_theta).

  dyplus_profiles.py <stats.h5> <field.h5> [--out dyplus_profiles.png]

The y node line is fixed in x; each station just rescales it by the local
u_tau(x): y+ = y*u_tau/nu, dy+ = dy*u_tau/nu. Shows the near-wall clustering,
the resolved band inside the BL, and the freestream coarsening (blayer grid).
"""
import argparse
import h5py
import numpy as np

nu = 1 / 450.


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stats")
    ap.add_argument("field")
    ap.add_argument("--out", default="dyplus_profiles.png")
    ap.add_argument("--n", type=int, default=10)
    a = ap.parse_args()

    with h5py.File(a.stats, "r") as f:
        nx, ny = int(f.attrs["nx"]), int(f.attrs["ny"])
        p = f["profile"][...].reshape(nx, ny, -1); xc = f["xcoord"][...]; yc = f["ycoord"][...]
    with h5py.File(a.field, "r") as g:
        ynode = g["y"][...]                      # faces (ny+1)
    U = p[:, :, 0]; Ue = U[:, -1]
    utau = np.sqrt(np.abs(nu * U[:, 0] / yc[0]))
    dy = np.diff(ynode)                          # (ny,) cell widths
    ycen = 0.5 * (ynode[:-1] + ynode[1:])
    # Re_theta(x) for station labels
    dyc = np.diff(np.concatenate(([0], 0.5 * (yc[:-1] + yc[1:]), [yc[-1]])))
    theta = np.sum((U / Ue[:, None]) * (1 - U / Ue[:, None]) * dyc, axis=1)
    reth = Ue * theta / nu
    delta99 = np.array([yc[np.searchsorted(U[i], 0.99 * Ue[i])] if np.any(U[i] >= 0.99 * Ue[i])
                        else yc[-1] for i in range(nx)])

    # 10 stations evenly in x over the developed plate (skip trip + outlet)
    i0 = np.searchsorted(xc, 40); i1 = np.searchsorted(xc, 0.92 * xc[-1])
    idx = np.linspace(i0, i1 - 1, a.n).astype(int)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    cmap = plt.cm.viridis(np.linspace(0, 1, len(idx)))
    fig, ax = plt.subplots(figsize=(9, 6))
    for c, i in zip(cmap, idx):
        yp = ycen * utau[i] / nu
        dyp = dy * utau[i] / nu
        ax.loglog(yp, dyp, "-", c=c, lw=1.2, label=f"Re$_\\theta$={reth[i]:.0f}")
        d99p = delta99[i] * utau[i] / nu
        ax.plot(d99p, np.interp(delta99[i], ycen, dyp), "o", c=c, ms=5, mfc="none")  # BL edge marker
    ax.axhline(4, ls=":", c="grey"); ax.text(1.2, 4.3, r"$\Delta y^+=4$", color="grey", fontsize=9)
    ax.set_xlabel(r"$y^+$"); ax.set_ylabel(r"$\Delta y^+$")
    ax.set_title(r"wall-normal spacing $\Delta y^+(y^+)$, full domain height (o = $\delta_{99}$)")
    ax.legend(fontsize=8, ncol=2, title="station"); ax.grid(alpha=0.3, which="both")
    fig.tight_layout(); fig.savefig(a.out, dpi=150)
    print("wrote", a.out)
    print(f"  wall dy+ range: {dy[0]*utau[idx].min()/nu:.2f}-{dy[0]*utau[idx].max()/nu:.2f}; "
          f"domain-top dy+ up to {dy[-1]*utau[idx].max()/nu:.0f}")


if __name__ == "__main__":
    main()
