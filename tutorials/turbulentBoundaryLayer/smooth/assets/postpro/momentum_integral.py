#!/usr/bin/env python3
"""von Karman momentum-integral balance of a boundary-layer statistics file.

    momentum_integral.py <bl_stats.h5> [<second bl_stats.h5>] [--re 450] [--ymax 60]
                         [--bands 450 500 550 600 650 700 730 745 748 750]

For each x band prints the growth of the momentum thickness, d(theta)/dx,
against what the integral balance of a flat-plate layer allows,

    c_f/2 - (H + 2) (theta/U_e) dU_e/dx ,

and the mean wall and freestream pressure. In a zero-pressure-gradient layer
the two columns agree (to the sampling scatter and the neglected normal-stress
terms) and p_wall = p_e. Where they part, the boundary is acting on the layer:
this is how the outlet zone of the case was measured (README, "The outlet
zone"). With two files the bands are printed side by side, then theta, c_f and
H of both at a few stations.

The statistics file is the solver's (boundarylayer_stats.f90): profiles on the
(x, y) plane, y fastest, columns U V W UU VV WW UV UW VW P at cell centres.
"""
import argparse
import h5py
import numpy as np

U, P = 0, 9
trapezoid = getattr(np, "trapezoid", None) or np.trapz    # numpy < 2


def load(name, re, ymax):
    with h5py.File(name, "r") as f:
        nx, ny = int(f.attrs["nx"]), int(f.attrs["ny"])
        prof = f["profile"][...].reshape(nx, ny, -1)
        x, y = f["xcoord"][...], f["ycoord"][...]
    top = np.searchsorted(y, ymax)                   # freestream row, above the layer
    ue = prof[:, top, U]
    yw = np.concatenate(([0.0], y[:top]))            # the wall, then the cell centres
    un = np.concatenate((np.zeros((nx, 1)), prof[:, :top, U]), axis=1)/ue[:, None]
    dstar = trapezoid(1.0 - un, yw, axis=1)
    theta = trapezoid(un*(1.0 - un), yw, axis=1)
    cf = 2.0*prof[:, 0, U]/y[0]/re/ue**2
    return dict(x=x, ue=ue, theta=theta, H=dstar/theta, cf=cf, pw=prof[:, 0, P], pe=prof[:, top, P])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("stats", nargs="+")
    ap.add_argument("--re", type=float, default=450.0)
    ap.add_argument("--ymax", type=float, default=60.0)
    ap.add_argument("--bands", type=float, nargs="+",
                    default=[450, 500, 550, 600, 650, 680, 700, 715, 730, 740, 745, 748, 749.5, 750])
    a = ap.parse_args()
    runs = [load(name, a.re, a.ymax) for name in a.stats[:2]]

    print("x band            " + "".join("| d theta/dx  balance    p_wall     p_e       " for _ in runs))
    for lo, hi in zip(a.bands[:-1], a.bands[1:]):
        line = f"{lo:7.1f} ..{hi:7.1f}  "
        for r in runs:
            m = np.where((r["x"] >= lo) & (r["x"] < hi))[0]
            i0, i1 = m[0], m[-1]
            dx = r["x"][i1] - r["x"][i0]
            dth = (r["theta"][i1] - r["theta"][i0])/dx
            due = (r["ue"][i1] - r["ue"][i0])/dx
            bal = r["cf"][m].mean()/2 - (r["H"][m].mean() + 2)*r["theta"][m].mean()/r["ue"][m].mean()*due
            line += f"| {dth: .2e}  {bal: .2e}  {r['pw'][m].mean(): .2e}  {r['pe'][m].mean(): .2e}  "
        print(line)

    print("\nstation      " + "".join("| theta    c_f       H       " for _ in runs) + ("| theta, c_f of the second against the first" if len(runs) == 2 else ""))
    for xs in (300, 400, 500, 600, 650, 700, 720, 740, runs[0]["x"][-1]):
        k = int(np.argmin(abs(runs[0]["x"] - xs)))
        line = f"x {runs[0]['x'][k]:7.2f}    "
        for r in runs:
            line += f"| {r['theta'][k]:.4f}  {r['cf'][k]:.6f}  {r['H'][k]:.4f}  "
        if len(runs) == 2:
            line += f"| {runs[1]['theta'][k]/runs[0]['theta'][k] - 1:+.2%}  {runs[1]['cf'][k]/runs[0]['cf'][k] - 1:+.2%}"
        print(line)


if __name__ == "__main__":
    main()
