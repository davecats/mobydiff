#!/usr/bin/env python3
"""Compare a boundary-layer snapshot against the Blasius similarity solution.

  compare_blasius.py <field.h5> [--theta 1.0] [--re-theta from-file]
                     [--stations 0.15 0.3 0.5 0.7 0.85] [--plot out.png]

Lengths are assumed nondimensionalized with the inlet momentum thickness
(theta_in = --theta), velocities with U_inf, so Re_theta,in = [flow] re of
the run. The virtual origin follows from Blasius growth:

    theta(x) = beta*sqrt(nu*(x + x_v)/U),  x_v = Re_theta,in*theta_in/beta^2

with beta = 2 f''(0) = 0.664115. At each station the script measures the
momentum thickness, displacement thickness/shape factor and the u and v
profiles, and compares them with the Blasius solution at the same
(virtual-origin-shifted) x. Gates: theta within 2%, H within 2%, u-profile
deviation < 1% of U_e, v-profile deviation < 15% of the local entrainment
scale (the p = 0 top holds the entrainment ~9% below Blasius aloft — the
finite-height bias). Default stations stop at x/lx = 0.7: the last ~15% of
the domain is the outlet influence zone (theta bends by a few % there).
The independent Blasius ODE solve below double-checks the solver's inlet
table.
"""
import argparse
import os
import sys

import h5py
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "tools"))
from compare_fields import load_field  # noqa: E402


def solve_blasius(eta_max=15.0, n=15000):
    """RK4 + secant shooting for f''' = -f f''/2; returns eta, f, f'."""
    h = eta_max / n

    def integrate(fpp0):
        y = np.array([0.0, 0.0, fpp0])
        f = np.empty(n + 1)
        fp = np.empty(n + 1)
        f[0], fp[0] = 0.0, 0.0
        rhs = lambda y: np.array([y[1], y[2], -0.5 * y[0] * y[2]])
        for i in range(n):
            k1 = rhs(y)
            k2 = rhs(y + 0.5 * h * k1)
            k3 = rhs(y + 0.5 * h * k2)
            k4 = rhs(y + h * k3)
            y = y + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
            f[i + 1], fp[i + 1] = y[0], y[1]
        return f, fp, y[1] - 1.0

    a, b = 0.33, 0.34
    _, _, ea = integrate(a)
    _, _, eb = integrate(b)
    for _ in range(50):
        c = b - eb * (b - a) / (eb - ea)
        f, fp, ec = integrate(c)
        if abs(ec) < 1e-13:
            break
        a, ea, b, eb = b, eb, c, ec
    eta = np.linspace(0.0, eta_max, n + 1)
    return eta, f, fp


def blasius_eval(eta, eta_b, f_b, fp_b):
    """f and f' at eta, following the outer asymptote f = eta - disp beyond
    the table (np.interp would clamp f and corrupt the entrainment v)."""
    disp = eta_b[-1] - f_b[-1]
    f = np.where(eta >= eta_b[-1], eta - disp, np.interp(eta, eta_b, f_b))
    fp = np.where(eta >= eta_b[-1], 1.0, np.interp(eta, eta_b, fp_b))
    return f, fp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("h5")
    ap.add_argument("--theta", type=float, default=1.0, help="inlet momentum thickness")
    ap.add_argument("--stations", type=float, nargs="+", default=[0.15, 0.3, 0.5, 0.7])
    ap.add_argument("--plot", default=None)
    a = ap.parse_args()

    with h5py.File(a.h5, "r") as f:
        u = load_field(f, "un")   # (nz, ny, nx); u[..,i] lives at x_node[i]
        v = load_field(f, "vn")   # v[.., j, .] lives at y_node[j], x cell centre
        xn = f["x"][...]
        yn = f["y"][...]
        re = float(f.attrs["re"])
    nu = 1.0 / re
    uinf = 1.0

    eta_b, f_b, fp_b = solve_blasius()
    beta = np.trapz(fp_b * (1.0 - fp_b), eta_b)          # momentum-thickness constant
    disp = np.trapz(1.0 - fp_b, eta_b)                   # displacement constant
    print(f"Blasius solve: beta = {beta:.6f} (0.664115), delta* const = {disp:.6f} (1.720788), "
          f"H = {disp/beta:.4f} (2.5911)")

    re_theta_in = uinf * a.theta / nu
    x_v = re_theta_in * a.theta / beta**2
    yc = 0.5 * (yn[:-1] + yn[1:])
    dy = np.diff(yn)
    nx = u.shape[2]
    lx = xn[-1]

    print(f"Re_theta,in = {re_theta_in:.1f}  virtual origin x_v = {x_v:.2f}")
    print(f"{'x/lx':>6} {'x':>8} {'Ue':>7} {'theta':>8} {'th_err%':>8} {'H':>7} "
          f"{'H_err%':>7} {'du_max':>9} {'dv_max':>9} {'v_top':>8}")

    worst = dict(th=0.0, H=0.0, du=0.0, dv=0.0)
    plots = []
    for frac in a.stations:
        i = min(int(round(frac * (nx - 1))), nx - 1)
        x = xn[i]
        up = u[:, :, i].mean(axis=0)                      # u(y) at x_node[i], z-avg
        ue = up[-1]
        # v at the same x: v index j at y_node[j], x cell-centred -> centre nearest
        iv = min(i, v.shape[2] - 1)
        vp = v[:, :, iv].mean(axis=0)

        theta = np.sum((up / ue) * (1.0 - up / ue) * dy)
        dstar = np.sum((1.0 - up / ue) * dy)
        H = dstar / theta
        x_tot = x + x_v
        th_th = beta * np.sqrt(nu * x_tot / ue)
        th_err = 100.0 * (theta - th_th) / th_th
        H_err = 100.0 * (H - disp / beta) / (disp / beta)

        # profiles vs Blasius at the same virtual-origin-shifted station
        eta_u = yc * np.sqrt(ue / (nu * x_tot))
        _, fp_u = blasius_eval(eta_u, eta_b, f_b, fp_b)
        u_th = ue * fp_u
        du = np.max(np.abs(up - u_th)) / ue

        # v rows sit on the LOW y faces; the top face lives in solver halos
        # and is not written
        eta_v = yn[:-1] * np.sqrt(ue / (nu * x_tot))
        vscale = 0.5 * np.sqrt(nu * ue / x_tot)
        f_v, fp_v = blasius_eval(eta_v, eta_b, f_b, fp_b)
        v_th = vscale * (eta_v * fp_v - f_v)
        v_edge = vscale * disp                            # entrainment velocity
        # gate v inside/just above the layer; far aloft the p = 0 top pulls
        # the entrainment a few % below Blasius (reported as v_top, info only)
        band = eta_v <= 8.0
        dv = np.max(np.abs(vp[band] - v_th[band])) / v_edge
        v_top = (vp[-1] - v_th[-1]) / v_edge

        worst["th"] = max(worst["th"], abs(th_err))
        worst["H"] = max(worst["H"], abs(H_err))
        worst["du"] = max(worst["du"], du)
        worst["dv"] = max(worst["dv"], dv)
        plots.append((x, x_tot, ue, up, eta_u, eta_v, vp / vscale))
        print(f"{frac:6.2f} {x:8.1f} {ue:7.4f} {theta:8.4f} {th_err:8.2f} {H:7.4f} "
              f"{H_err:7.2f} {du:9.2e} {dv:9.2e} {v_top:+8.3f}")

    ok = worst["th"] < 2.0 and worst["H"] < 2.0 and worst["du"] < 1e-2 and worst["dv"] < 1.5e-1
    print(f"worst: theta {worst['th']:.2f}%  H {worst['H']:.2f}%  "
          f"du/Ue {worst['du']:.2e}  dv/v_edge {worst['dv']:.2e}")
    print("blasius gate:", "PASS" if ok else "FAIL")

    if a.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 3, figsize=(14, 4))
        for x, x_tot, ue, up, eta_u, eta_v, vsim in plots:
            ax[0].plot(up / ue, eta_u, ".", ms=3, label=f"x={x:.0f}")
        ax[0].plot(fp_b, eta_b, "k-", lw=1, label="Blasius")
        ax[0].set_xlabel("u/Ue"); ax[0].set_ylabel(r"$\eta$"); ax[0].set_ylim(0, 8)
        ax[0].legend(fontsize=7)
        # v similarity form: v / (0.5 sqrt(nu Ue / x)) = eta f' - f
        for x, x_tot, ue, up, eta_u, eta_v, vsim in plots:
            ax[1].plot(vsim, eta_v, ".", ms=3, label=f"x={x:.0f}")
        ax[1].plot(eta_b*fp_b - f_b, eta_b, "k-", lw=1, label="Blasius")
        ax[1].set_xlabel(r"$v\,/\,\frac{1}{2}\sqrt{\nu U_e/x}$")
        ax[1].set_ylabel(r"$\eta$"); ax[1].set_ylim(0, 8); ax[1].set_xlim(0, 2)
        ax[1].legend(fontsize=7)
        xs = np.array([p[0] for p in plots])
        th_meas = [np.sum((p[3]/p[2])*(1-p[3]/p[2])*dy) for p in plots]
        xfine = np.linspace(0, lx, 200)
        ax[2].plot(xs, th_meas, "o", label="measured")
        ax[2].plot(xfine, beta*np.sqrt(nu*(xfine + x_v)), "k-", lw=1, label="Blasius")
        ax[2].set_xlabel("x"); ax[2].set_ylabel(r"$\theta$"); ax[2].legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(a.plot, dpi=150)
        print("wrote", a.plot)

    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
