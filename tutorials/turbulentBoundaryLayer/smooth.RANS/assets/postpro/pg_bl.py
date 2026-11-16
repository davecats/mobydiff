#!/usr/bin/env python3
"""Boundary-layer integrals of a RANS field WITH a pressure gradient.

    python3 pg_bl.py FIELD.h5 [--every 50]

Under a pressure gradient the potential flow above the layer varies with y
(dU/dy = dV/dx), so "U_e = U at the top row" (the ZPG convention of
rans_bl.py and of the DNS export) misplaces the edge. Here the edge follows
Griffin, Fu & Moin (Phys. Rev. Fluids 6, 2021): the local inviscid velocity

    U_I(x, y) = sqrt(2 (P0_e - p) - V^2),

delta_99 is the first height where U/U_I = 0.99, U_e = U_I(delta_99), and
delta*, theta integrate (1 - U/U_I) and U/U_I (1 - U/U_I) from the wall to
delta_99 (the pressure-gradient convention; integrating to the top would add
the free stream's U_I uncertainty times its height).

The reference total pressure P0_e(x) is that of the free stream JUST ABOVE
the layer, not of the top row as in Griffin et al.: the free stream here is
not iso-energetic. The inlet fixes U (Blasius) but not p, so P0 = p + U^2/2
varies across the inlet streamlines, by ~0.5 % of U^2 with no gradient and
up to 1.5 % under the top-pressure tables of apg/fpg.RANS (the DNS shows the
same weak outer shear: U is 0.3 % below U_top just above the layer). With
the top-row P0, U/U_I levels off at ~0.991 above an FPG layer and delta_99
lands anywhere on that plateau (20 instead of 7.4 at x = 200). The edge
height y_e is the first point, going up, where |dP0/dy| < EDGE_DP0DY (the
free-stream level is <= 5e-4, the layer's >= 2e-3) and U/U_I,top > 0.9;
P0_e = P0(y_e). For the DNS export (no pressure, no gradient) the same test
reads |U dU/dy| < EDGE_DP0DY and U_e = U(y_e).

    P_e = p(delta_99), dP_e/dx smoothed over SMOOTH_L,
    beta = delta*/tau_w dP_e/dx (Clauser), K = nu/U_e^2 dU_e/dx.
"""
import argparse
import os
import sys

import numpy as np
from scipy.signal import savgol_filter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rans_bl as rb  # noqa: E402

SMOOTH_L = 30.0          # length of the Savitzky-Golay window for d/dx, in delta*_0
EDGE_DP0DY = 1.0e-3      # |dP0/dy| below which a point is free stream, U_inf^2/delta*_0


def smooth_ddx(x, q, length=SMOOTH_L):
    """d/dx of q on a uniform x line with a quadratic Savitzky-Golay filter."""
    dx = x[1] - x[0]
    win = int(round(length / dx)) | 1
    return savgol_filter(q, win, 2, deriv=1, delta=dx, mode="interp")


def edge_index(y, ptot, ratio):
    """Per x row: the first j where |d ptot/dy| < EDGE_DP0DY and ratio > 0.9
    (the top row if none). ptot, ratio are (nx, ny)."""
    g = np.abs(np.gradient(ptot, y, axis=1))
    ok = (g < EDGE_DP0DY) & (ratio > 0.9)
    j = np.argmax(ok, axis=1)
    j[~ok.any(axis=1)] = len(y) - 1
    return j


def edge_total_pressure(y, U, P0, p, V):
    """P0 of the free stream just above the layer (module docstring)."""
    UIt = np.sqrt(np.maximum(2.0 * (P0[:, -1:] - p) - V ** 2, 1e-12))
    j = edge_index(y, P0, U / UIt)
    return P0[np.arange(len(j)), j]


def pg_integrals(f):
    """1D quantities of a loaded field (rans_bl.load): dict of x-arrays."""
    x, y, U, V, p, re = f["x"], f["y"], f["U"], f["V"], f["pn"], f["re"]
    nu = 1.0 / re
    P0 = edge_total_pressure(y, U, p + 0.5 * (U ** 2 + V ** 2), p, V)
    UI = np.sqrt(np.maximum(2.0 * (P0[:, None] - p) - V ** 2, 1e-12))
    r = U / UI
    nx = len(x)
    d99, Ue, Pe, dstar, theta = (np.empty(nx) for _ in range(5))
    for i in range(nx):
        j = int(np.argmax(r[i] >= 0.99))
        if j == 0:
            j = 1
        w = (0.99 - r[i, j - 1]) / (r[i, j] - r[i, j - 1])
        d99[i] = y[j - 1] + w * (y[j] - y[j - 1])
        Ue[i] = UI[i, j - 1] + w * (UI[i, j] - UI[i, j - 1])
        Pe[i] = p[i, j - 1] + w * (p[i, j] - p[i, j - 1])
        # trapezoid from the wall (U = 0) through the cell centres to delta_99
        yy = np.r_[0.0, y[:j], d99[i]]
        rr = np.r_[0.0, r[i, :j], 0.99]
        dstar[i] = np.trapz(1.0 - rr, yy)
        theta[i] = np.trapz(rr * (1.0 - rr), yy)
    y0, y1 = y[0], y[1]
    tauw = nu * (U[:, 0] * y1 ** 2 - U[:, 1] * y0 ** 2) / (y0 * y1 * (y1 - y0))
    dPe = smooth_ddx(x, Pe)
    dUe = smooth_ddx(x, Ue)
    return dict(x=x, Ue=Ue, Pe=Pe, P0=P0, delta99=d99, dstar=dstar, theta=theta,
                H=dstar / theta, tauw=tauw, cf=2.0 * tauw / Ue ** 2,
                Re_theta=Ue * theta / nu, Re_tau=np.sqrt(np.abs(tauw)) * d99 / nu,
                dPedx=dPe, dUedx=dUe, beta=dstar / tauw * dPe, K=nu / Ue ** 2 * dUe,
                pw=p[:, 0], Vtop=V[:, -1])


def dns_truncated(d):
    """The DNS export's ZPG field with the same edge and truncation: no
    pressure in the export and none needed with no gradient, so the edge test
    is on U^2/2 and U_e = U(y_e) (module docstring); delta*, theta from the
    wall to delta_99. Returns x, Ue, dstar, theta, cf."""
    x, y, U = d["x"], d["y"], d["U"]           # U is (ny, nx)
    nx = len(x)
    je = edge_index(y, 0.5 * U.T ** 2, (U / U[-1]).T)
    Ue = U[je, np.arange(nx)]
    dstar, theta = np.empty(nx), np.empty(nx)
    for i in range(nx):
        r = U[:, i] / Ue[i]
        j = max(int(np.argmax(r >= 0.99)), 1)
        w = (0.99 - r[j - 1]) / (r[j] - r[j - 1])
        d99 = y[j - 1] + w * (y[j] - y[j - 1])
        yy = np.r_[0.0, y[:j], d99]
        rr = np.r_[0.0, r[:j], 0.99]
        dstar[i] = np.trapz(1.0 - rr, yy)
        theta[i] = np.trapz(rr * (1.0 - rr), yy)
    return dict(x=x, Ue=Ue, dstar=dstar, theta=theta, cf=d["cf"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("field")
    ap.add_argument("--every", type=float, default=50.0)
    a = ap.parse_args()
    f = rb.load(a.field)
    q = pg_integrals(f)
    print(f"== {a.field} (step {f['step']}, t {f['t']:.0f})")
    print("    x |   Ue     Re_th   H      cf       beta    K*1e6  delta99")
    for xs in np.arange(a.every, 600.0, a.every):
        i = np.argmin(abs(q["x"] - xs))
        print(f" {xs:4.0f} | {q['Ue'][i]:.4f} {q['Re_theta'][i]:6.0f} {q['H'][i]:.3f} {q['cf'][i]:.5f} "
              f"{q['beta'][i]:+.3f} {q['K'][i] * 1e6:+.3f} {q['delta99'][i]:6.2f}")


if __name__ == "__main__":
    main()
