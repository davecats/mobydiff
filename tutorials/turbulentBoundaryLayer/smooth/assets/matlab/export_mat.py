#!/usr/bin/env python3
"""Export the shipped boundary-layer statistics to ONE MATLAB file for
external validation (e.g. a RANS transition model):

    python3 export_mat.py [--stats ../../production_stats.h5]
                          [--since ../../production_stats_t0.h5]
                          [--xmax 550] [--out tbl_dns_smooth.mat]

The window is the difference of the two cumulative statistics files (raw sums
and counts), as in make_mobydiff_nc.py. Everything stored is span + time
averaged; 2D arrays are (ny, nx) -- rows y, columns x, MATLAB convention -- so
that contourf(x, y, U) works directly. The nondimensionalisation and every
definition are documented in read_tbl_dns.m, which is the reader to hand out.
"""
import argparse
import os

import h5py
import numpy as np
import scipy.io as sio

U, V, W, UU, VV, WW, UV, UW, VW, P = range(10)
BLASIUS_DSTAR = 1.7207876573        # delta* = 1.72079 sqrt(nu x / U_inf)

HERE = os.path.dirname(os.path.abspath(__file__))
TUT = os.path.join(HERE, "..", "..")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stats", default=os.path.join(TUT, "production_stats.h5"))
    ap.add_argument("--since", default=os.path.join(TUT, "production_stats_t0.h5"))
    ap.add_argument("--ini", default=os.path.join(TUT, "production_stats.ini"))
    ap.add_argument("--xmax", type=float, default=550.0,
                    help="drop x > xmax, the outlet adjustment zone (README 'The outlet zone')")
    ap.add_argument("--out", default=os.path.join(HERE, "tbl_dns_smooth.mat"))
    a = ap.parse_args()

    with h5py.File(a.stats, "r") as f:
        nx, ny = int(f.attrs["nx"]), int(f.attrs["ny"])
        re = float(f.attrs["re"])
        raw, cnt = f["raw_sum"][...], f["count"][...]
        x, y = f["xcoord"][...], f["ycoord"][...]
        t1 = float(f.attrs["t_current"])
    with h5py.File(a.since, "r") as f:
        raw = raw - f["raw_sum"][...]
        cnt = cnt - f["count"][...]
        t0 = float(f.attrs["t_current"])
    assert np.all(cnt > 0)
    prof = (raw / cnt[:, None]).reshape(nx, ny, -1)
    nu = 1.0 / re

    # Means and Reynolds stresses (covariances of the fluctuations).
    Um, Vm, Wm, Pm = (prof[:, :, k] for k in (U, V, W, P))
    uu = prof[:, :, UU] - Um * Um
    vv = prof[:, :, VV] - Vm * Vm
    ww = prof[:, :, WW] - Wm * Wm
    uv = prof[:, :, UV] - Um * Vm
    k = 0.5 * (uu + vv + ww)

    # Integral thicknesses: midpoint rule on the cell-centre profile, the
    # wall at y = 0 and the top at y = Ly bounding the cells.
    yf = np.r_[0.0, 0.5 * (y[1:] + y[:-1]), 2 * y[-1] - 0.5 * (y[-1] + y[-2])]
    dy = np.diff(yf)
    Ue = Um[:, -1]
    r = Um / Ue[:, None]
    dstar = np.sum((1 - r) * dy, axis=1)
    theta = np.sum(r * (1 - r) * dy, axis=1)
    j99 = np.argmax(r >= 0.99, axis=1)               # first cell above 0.99 U_e
    jm = j99 - 1                                     # linear interpolation below it
    i = np.arange(nx)
    d99 = y[jm] + (0.99 - r[i, jm]) * (y[j99] - y[jm]) / (r[i, j99] - r[i, jm])

    # Wall shear: second-order one-sided gradient through U(0) = 0 and the
    # first two cell centres.
    y0, y1 = y[0], y[1]
    dUdy_w = (Um[:, 0] * y1**2 - Um[:, 1] * y0**2) / (y0 * y1 * (y1 - y0))
    tauw = nu * dUdy_w
    utau = np.sqrt(np.abs(tauw))

    # Eddy viscosity from the DNS by the Boussinesq hypothesis on the one
    # shear component that matters in a thin layer: nu_t = -<u'v'> / S_xy*2,
    # 2 S_xy = dU/dy + dV/dx. Undefined where the shear vanishes (outside the
    # layer): NaN above delta_99.
    dUdy = np.gradient(np.c_[np.zeros(nx), Um], np.r_[0.0, y], axis=1, edge_order=2)[:, 1:]
    dVdx = np.gradient(Vm, x, axis=0, edge_order=2)
    nut = -uv / (dUdy + dVdx)
    nut[y[None, :] > d99[:, None]] = np.nan

    # Clauser pressure-gradient parameter (zero for a ZPG layer), U_e smoothed
    # over 21 cells (edge-padded) before differentiating.
    Ue_s = np.convolve(np.pad(Ue, 10, mode="edge"), np.ones(21) / 21, mode="valid")
    beta = -(dstar / tauw) * Ue * np.gradient(Ue_s, x)

    # Re_x: the inflow is the Blasius profile at Re_delta* = 450, i.e. at a
    # virtual leading edge x_v = (Re_delta*,in / 1.72079)^2 nu / U_inf upstream.
    rex_in = (re / BLASIUS_DSTAR) ** 2
    rex = rex_in + x * re

    keep = x <= a.xmax
    T = lambda f2: np.ascontiguousarray(f2[keep].T)          # (ny, nx_kept)
    C = lambda f1: f1[keep]
    out = dict(
        x=C(x), y=y, nu=nu, Re_dstar_in=re, Re_x_in=rex_in,
        U=T(Um), V=T(Vm), P=T(Pm),
        uu=T(uu), vv=T(vv), ww=T(ww), uv=T(uv), k=T(k),
        dUdy=T(dUdy), nut_over_nu=T(nut / nu),
        Ue=C(Ue), delta99=C(d99), dstar=C(dstar), theta=C(theta), H=C(dstar / theta),
        tauw=C(tauw), utau=C(utau), cf=C(2 * tauw / Ue**2), p_wall=C(Pm[:, 0]),
        Re_x=C(rex), Re_theta=C(Ue * theta / nu), Re_dstar=C(Ue * dstar / nu),
        Re_tau=C(utau * d99 / nu), beta_clauser=C(beta),
        urms_max=C(np.sqrt(np.max(uu, axis=1))),
        time_window=np.array([t0, t1]), n_samples=float(cnt[0]),
        input_ini=open(a.ini).read(),
        source=("mobydiff DNS, smooth-wall ZPG turbulent boundary layer, CaNS grid "
                "3200x384x136 on 650x100x26 delta*_in; span+time averages over "
                f"t = {t0:.1f}..{t1:.1f} ({cnt[0]:.0f} samples); x <= {a.xmax:g}"),
    )
    sio.savemat(a.out, out, do_compression=True, oned_as="column")
    print(f"wrote {a.out}: {keep.sum()} x {ny} points, window {t0:.0f}..{t1:.0f}, "
          f"Re_theta {out['Re_theta'][0]:.0f}..{out['Re_theta'][-1]:.0f}")


if __name__ == "__main__":
    main()
