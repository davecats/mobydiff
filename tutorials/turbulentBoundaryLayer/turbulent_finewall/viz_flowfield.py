#!/usr/bin/env python3
"""Visualize a turbulent-boundary-layer snapshot: side (x-y) views of the
streamwise u, wall-normal v, spanwise w and pressure p, plus the near-wall
streaks from above (x-z).

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
        w = load_field(f, "wn")
        p = load_field(f, "pn")
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

    # Spanwise-mean pressure <p>_z (the large-scale structure, incl. the ZPG
    # offset and the velocity-neutral pn-drift), and the turbulent
    # fluctuation p' = p - <p>_z (the low-pressure cores of the near-wall
    # vortices). Centre <p>_z on its own mean so the diverging map reads.
    pmean_z = p.mean(axis=0)                        # (ny, nx)
    p = p - pmean_z[np.newaxis, :, :]              # p' (nz, ny, nx)
    pmean_z = pmean_z - pmean_z.mean()

    fig, ax = plt.subplots(6, 1, figsize=(15, 19))

    def side(axis, field, cmap, label, title, vmin=None, vmax=None):
        # field is a full (nz,ny,nx) volume (mid-span slice) or a 2D (ny,nx)
        # spanwise-averaged field.
        f2 = field[zmid, :jmax, :] if field.ndim == 3 else field[:jmax, :]
        if vmin is None:
            lim = np.percentile(np.abs(f2), 99)
            vmin, vmax = -lim, lim
        pc = axis.pcolormesh(xc, yc[:jmax], f2, cmap=cmap, vmin=vmin, vmax=vmax,
                             shading="auto", rasterized=True)
        fig.colorbar(pc, ax=axis, pad=0.01, label=label)
        axis.set_ylabel("y"); axis.set_ylim(0, a.ymax); axis.set_title(title)

    # (A) streamwise velocity u (side) with delta99
    side(ax[0], u, "turbo", "u",
         f"streamwise velocity u — mid-span x-y plane (white = δ99); "
         f"t={t:.0f}, Re_δ*,0={re:.0f}", vmin=0.0, vmax=1.1)
    ax[0].plot(xc, d99, "w-", lw=0.8, alpha=0.8)

    # (B) wall-normal velocity v (side): ejections / sweeps
    side(ax[1], v, "RdBu_r", "v",
         "wall-normal velocity v — mid-span x-y plane (ejections v>0 / sweeps v<0)")

    # (C) spanwise velocity w (side): the spanwise turbulent motions
    side(ax[2], w, "RdBu_r", "w",
         "spanwise velocity w — mid-span x-y plane")

    # (D) turbulent pressure fluctuation p' (side): low-pressure eddy cores
    side(ax[3], p, "RdBu_r", "p'",
         "turbulent pressure fluctuation p' = p − ⟨p⟩_z — mid-span x-y plane")

    # (E) spanwise-mean pressure <p>_z (side): the large-scale structure /
    # velocity-neutral pn-drift that p' removes
    side(ax[4], pmean_z, "RdBu_r", "⟨p⟩_z",
         "spanwise-mean pressure ⟨p⟩_z − mean (the large-scale / pn-drift mode)")

    # (F) near-wall streaks: u' = u - <u>_z in an x-z plane at y+ ~ 12
    up = u[:, jwall, :] - Um[jwall, :]
    lim = np.percentile(np.abs(up), 99)
    pc = ax[5].pcolormesh(xc, zc, up, cmap="RdBu_r", vmin=-lim, vmax=lim,
                          shading="auto", rasterized=True)
    fig.colorbar(pc, ax=ax[5], pad=0.01, label="u'")
    ax[5].set_ylabel("z"); ax[5].set_xlabel("x"); ax[5].set_aspect("equal")
    ax[5].set_title(f"near-wall streaks: u' in the x-z plane at y⁺≈{yplus_plane:.0f} "
                    f"(elongated low/high-speed streaks)")

    fig.tight_layout()
    fig.savefig(a.out, dpi=130)
    print(f"u_tau≈{utau:.4f}, l+≈{lplus:.4f}, streak plane y+≈{yplus_plane:.1f}")
    print("wrote", a.out)


if __name__ == "__main__":
    main()
