#!/usr/bin/env python3
"""Steady-state convergence diagnostics for the Blasius run.

  check_convergence.py conv_*.h5 [--plot convergence.png]

For each consecutive snapshot pair reports the drift RATE
max|u(t2)-u(t1)| / (t2-t1) and the rms drift rate, and per snapshot the
momentum-thickness error vs Blasius at x/lx = 0.5 and the max interior
divergence. A converged run has the drift rate falling to a floor and the
gate metrics flat; a still-developing run has them trending.
"""
import argparse
import os
import sys

import h5py
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "tools"))
from compare_fields import load_field  # noqa: E402
from compare_blasius import solve_blasius  # noqa: E402

BETA = 0.664114672


def interior_maxdiv(u, v, w, xn, yn, zn):
    dx, dy, dz = np.diff(xn), np.diff(yn), np.diff(zn)
    nz, ny, nx = u.shape
    div = ((np.roll(u, -1, 2) - u) / dx[None, None, :]
           + (np.roll(v, -1, 1) - v) / dy[None, :, None]
           + (np.roll(w, -1, 0) - w) / dz[:, None, None])
    return float(np.max(np.abs(div[:, :ny - 1, :nx - 1])))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("h5", nargs="+")
    ap.add_argument("--plot", default=None)
    a = ap.parse_args()
    files = sorted(a.h5, key=lambda s: int(s.split("_")[-1].split(".")[0]))

    eta_b, f_b, fp_b = solve_blasius()

    prev_u, prev_t = None, None
    rows = []
    print(f"{'file':22s} {'t':>8} {'drift/dt(max)':>14} {'drift/dt(rms)':>14} "
          f"{'theta_err%':>11} {'maxdiv':>10}")
    for fn in files:
        with h5py.File(fn, "r") as f:
            u = load_field(f, "un"); v = load_field(f, "vn"); w = load_field(f, "wn")
            xn, yn, zn = f["x"][...], f["y"][...], f["z"][...]
            re = float(f.attrs["re"]); t = float(f.attrs["t_current"])
        nu = 1.0 / re
        um = u.mean(axis=0)
        # theta error at x/lx = 0.5
        nx = um.shape[1]; i = nx // 2
        dy = np.diff(yn)
        up = um[:, i]; ue = up[-1]
        theta = np.sum((up / ue) * (1 - up / ue) * dy)
        xv = (nu ** -1) / BETA ** 2
        th_th = BETA * np.sqrt(nu * (xn[i] + xv) / ue)
        th_err = 100.0 * (theta - th_th) / th_th
        mdiv = interior_maxdiv(u, v, w, xn, yn, zn)
        if prev_u is not None:
            dt = t - prev_t
            dmax = float(np.max(np.abs(u - prev_u))) / dt
            drms = float(np.sqrt(np.mean((u - prev_u) ** 2))) / dt
        else:
            dmax = drms = np.nan
        rows.append((t, dmax, drms, th_err, mdiv))
        print(f"{os.path.basename(fn):22s} {t:8.1f} {dmax:14.3e} {drms:14.3e} "
              f"{th_err:11.3f} {mdiv:10.2e}")
        prev_u, prev_t = u, t

    if a.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        arr = np.array(rows)
        fig, ax = plt.subplots(1, 2, figsize=(11, 4))
        ax[0].semilogy(arr[1:, 0], arr[1:, 1], "o-", label="max|du|/dt")
        ax[0].semilogy(arr[1:, 0], arr[1:, 2], "s-", label="rms|du|/dt")
        ax[0].set_xlabel("t"); ax[0].set_ylabel("drift rate"); ax[0].legend()
        ax[0].set_title("field drift rate -> floor = converged")
        ax[1].plot(arr[:, 0], arr[:, 3], "o-")
        ax[1].set_xlabel("t"); ax[1].set_ylabel(r"$\theta$ error % (x/lx=0.5)")
        ax[1].set_title("gate metric stationary = converged")
        fig.tight_layout(); fig.savefig(a.plot, dpi=150)
        print("wrote", a.plot)


if __name__ == "__main__":
    main()
