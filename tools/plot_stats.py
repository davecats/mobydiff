#!/usr/bin/env python3
"""Read and plot the channel-flow statistics written by mobydiff.

The solver accumulates wall-normal profiles into an HDF5 file (default name
``channel_stats.h5``). This little tool reads that file and shows the
quantities you usually want from a turbulent channel:

  * the mean streamwise velocity profile        U(y)   and   U+(y+)
  * the Reynolds stresses                        u'_rms, v'_rms, w'_rms, -<u'v'>
  * the total stress and its viscous / turbulent breakdown

It also prints the friction velocity u_tau and friction Reynolds number Re_tau.

Usage:
    python3 tools/plot_stats.py channel_stats.h5
    python3 tools/plot_stats.py channel_stats.h5 --save stats.png

The file layout (for the curious):
    dataset "coord"   -> y at the cell centres            shape (Ny,)
    dataset "count"   -> area-weighted sample count        shape (Ny,)
    dataset "profile" -> the AVERAGED statistics           shape (Ny, 14)
    attributes: "re", "forcing_x", "t_current", ...
The 14 profile columns are (already divided by "count", i.e. real averages):
    0:<u>  1:<v>  2:<w>  3:<uu>  4:<vv>  5:<ww>  6:<uv>  7:<uw>  8:<vw> ...
"""

import argparse
import glob
import os

import h5py
import numpy as np
import matplotlib.pyplot as plt

# Column indices inside the "profile" dataset.
U, V, W, UU, VV, WW, UV = 0, 1, 2, 3, 4, 5, 6


def read_stats(path):
    """Read one stats file (plus any refinement-level companions).

    Refined runs split the profile over several files: ``channel_stats.h5``
    for the coarse level and ``channel_stats_l1.h5``, ``_l2.h5`` ... for the
    finer levels. We read them all and stitch the rows together, sorted by y.
    """
    base, ext = os.path.splitext(path)
    files = [path] + sorted(glob.glob(base + "_l*" + ext))

    ys, profs = [], []
    attrs = {}
    for f in files:
        if not os.path.exists(f):
            continue
        with h5py.File(f, "r") as h:
            coord = h["coord"][...]
            prof = h["profile"][...]
            count = h["count"][...]
            keep = count > 0.0          # rows that actually collected samples
            ys.append(coord[keep])
            profs.append(prof[keep])
            if not attrs:
                attrs = {k: h.attrs[k] for k in h.attrs}

    y = np.concatenate(ys)
    prof = np.concatenate(profs, axis=0)
    order = np.argsort(y)               # merge levels into one monotone profile
    return y[order], prof[order], attrs


def main():
    ap = argparse.ArgumentParser(description="Plot mobydiff channel statistics.")
    ap.add_argument("file", nargs="?", default="channel_stats.h5",
                    help="stats HDF5 file (default: channel_stats.h5)")
    ap.add_argument("--save", metavar="PNG", default=None,
                    help="save the figure to this file instead of showing it")
    args = ap.parse_args()

    y, P, attrs = read_stats(args.file)

    re = float(attrs.get("re", 1.0))
    fx = float(attrs.get("forcing_x", 0.0))
    nu = 1.0 / re                       # kinematic viscosity (mu, since rho = 1)

    # --- mean profile and fluctuations ------------------------------------
    Umean = P[:, U]
    # rms of the fluctuations: sqrt(<u u> - <u><u>), clipped at 0 for round-off.
    urms = np.sqrt(np.maximum(P[:, UU] - P[:, U] ** 2, 0.0))
    vrms = np.sqrt(np.maximum(P[:, VV] - P[:, V] ** 2, 0.0))
    wrms = np.sqrt(np.maximum(P[:, WW] - P[:, W] ** 2, 0.0))
    # Reynolds shear stress -<u'v'>
    uv = -(P[:, UV] - P[:, U] * P[:, V])

    # --- stress budget ----------------------------------------------------
    # viscous stress mu*dU/dy and Reynolds stress -<u'v'>; their sum should be
    # a straight line across the channel (fully developed balance).
    tau_visc = nu * np.gradient(Umean, y)
    tau_reyn = uv
    tau_tot = tau_visc + tau_reyn

    # --- friction velocity ------------------------------------------------
    h = 0.5 * (y.min() + y.max())       # channel half-height
    # (1) from the driving force balance  tau_wall = fx * h
    utau_force = np.sqrt(abs(fx * h)) if fx != 0.0 else 0.0
    # (2) from the wall velocity gradient tau_wall = nu * dU/dy at the wall
    utau_wall = np.sqrt(abs(nu * (Umean[1] - Umean[0]) / (y[1] - y[0])))
    utau = utau_force if utau_force > 0.0 else utau_wall
    retau = utau / nu

    print(f"file            : {args.file}")
    print(f"time            : {float(attrs.get('t_current', 0.0)):.4g}")
    print(f"Re (bulk)       : {re:.4g}")
    print(f"u_tau (forcing) : {utau_force:.5g}")
    print(f"u_tau (wall dU) : {utau_wall:.5g}")
    print(f"Re_tau          : {retau:.5g}")

    # wall units for the mean profile (distance from the nearest wall)
    ywall = np.minimum(y - y.min(), y.max() - y)
    yplus = ywall * retau
    uplus = Umean / utau if utau > 0.0 else Umean

    # --- plot -------------------------------------------------------------
    fig, ax = plt.subplots(2, 2, figsize=(11, 8))

    ax[0, 0].plot(y, Umean, "-o", ms=3)
    ax[0, 0].set(xlabel="y", ylabel="U", title="Mean streamwise velocity")

    ax[0, 1].semilogx(yplus, uplus, "-o", ms=3, label="DNS")
    yv = np.logspace(0, np.log10(max(yplus.max(), 2)), 50)
    ax[0, 1].plot(yv[yv < 12], yv[yv < 12], "k--", lw=1, label="U+ = y+")
    ax[0, 1].plot(yv[yv > 8], np.log(yv[yv > 8]) / 0.41 + 5.0, "k:", lw=1,
                  label="log law")
    ax[0, 1].set(xlabel="y+", ylabel="U+", title="Mean profile in wall units")
    ax[0, 1].legend()

    ax[1, 0].plot(y, urms, label="u'_rms")
    ax[1, 0].plot(y, vrms, label="v'_rms")
    ax[1, 0].plot(y, wrms, label="w'_rms")
    ax[1, 0].plot(y, uv, label="-<u'v'>")
    ax[1, 0].set(xlabel="y", ylabel="stress", title="Reynolds stresses")
    ax[1, 0].legend()

    ax[1, 1].plot(y, tau_visc, label="viscous  nu dU/dy")
    ax[1, 1].plot(y, tau_reyn, label="Reynolds -<u'v'>")
    ax[1, 1].plot(y, tau_tot, "k-", lw=2, label="total")
    ax[1, 1].set(xlabel="y", ylabel="stress",
                 title="Total stress and breakdown")
    ax[1, 1].legend()

    fig.tight_layout()
    if args.save:
        fig.savefig(args.save, dpi=130)
        print(f"saved figure    : {args.save}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
