#!/usr/bin/env python3
"""Smooth against rough wall: the momentum balance to the outlet, the roughness
Reynolds number along the plate and the Hama roughness function.

    rough_compare.py --smooth <stats.h5> [<earlier stats.h5>] --rough <stats.h5> [<earlier>]
                     [--re 450] [--ymax 60] [--k 0.38348] [--crest 0.767]
                     [--stations 300 385 450 550 600]

A statistics file is the solver's cumulative one (raw sums + counts); giving
an EARLIER copy of the same file makes the window the difference of the two
(what `submit_rough_tbl.sh` keeps as stats_<side>_<step>.h5), so the
transient after the restart is excluded.

Printed:
 1. per x band, for both sides: d theta/dx against the von Karman balance
    c_f/2 - (H + 2)(theta/U_e) dU_e/dx, with c_f from the WALL CELL on the
    smooth side and from the BALANCE itself on the rough side (the wall-cell
    velocity sits inside the roughness there), and the mean wall / freestream
    pressure -- the outlet-zone measurement of the README, now with a body
    crossing the outlet plane;
 2. along x: theta, c_f, u_tau and k+ = k u_tau Re of the rough side, next to
    the smooth side's;
 3. at the stations: the Hama roughness function Delta U+ = U+_smooth(y+) -
    U+_rough(y+) averaged over the log region ABOVE THE CRESTS (y > crest,
    30 < y+ < 0.2 delta99+), each side in its own u_tau, the rough y measured
    from the mean plane (--k); and the two c_f.

The rough-wall c_f is 2 d theta/dx + (H + 2) theta/U_e dU_e/dx (the balance
solved for c_f), smoothed over +-10 cells; it is the total drag per unit
area, form drag included, which is what the roughness function is defined
with. The mean profile over the roughness is the plane average INCLUDING the
penalized solid (superficial average); above the crests it is the fluid
average, which is all the comparison uses.
"""
import argparse
import h5py
import numpy as np

U, P = 0, 9
trapezoid = getattr(np, "trapezoid", None) or np.trapz    # numpy < 2
KAPPA = 0.41


def window(names):
    """Mean profiles of the window between two cumulative files (or one)."""
    with h5py.File(names[0], "r") as f:
        nx, ny = int(f.attrs["nx"]), int(f.attrs["ny"])
        raw = f["raw_sum"][...]
        cnt = f["count"][...]
        x, y = f["xcoord"][...], f["ycoord"][...]
        t1 = float(f.attrs.get("t_current", np.nan))
    t0 = np.nan
    if len(names) > 1:
        with h5py.File(names[1], "r") as f:
            raw = raw - f["raw_sum"][...]
            cnt = cnt - f["count"][...]
            t0 = float(f.attrs.get("t_current", np.nan))
    prof = (raw/np.maximum(cnt, 1)[:, None]).reshape(nx, ny, -1)
    return x, y, prof, (t0, t1)


def integral(x, y, prof, re, ymax):
    top = np.searchsorted(y, ymax)
    ue = prof[:, top, U]
    yw = np.concatenate(([0.0], y[:top]))
    un = np.concatenate((np.zeros((len(x), 1)), prof[:, :top, U]), axis=1)/ue[:, None]
    dstar = trapezoid(1.0 - un, yw, axis=1)
    theta = trapezoid(un*(1.0 - un), yw, axis=1)
    cf_wall = 2.0*prof[:, 0, U]/y[0]/re/ue**2
    d99 = np.array([yw[np.searchsorted(un[i], 0.99)] for i in range(len(x))])
    return dict(x=x, y=y, ue=ue, theta=theta, H=dstar/theta, cf_wall=cf_wall, d99=d99,
                pw=prof[:, 0, P], pe=prof[:, top, P], prof=prof)


def smooth_series(v, half):
    k = np.ones(2*half + 1)/(2*half + 1)
    return np.convolve(np.pad(v, half, mode="edge"), k, mode="valid")


def cf_balance(r, half=10):
    """c_f from the integral balance, smoothed over +-half cells."""
    x = r["x"]
    th = smooth_series(r["theta"], half)
    ue = smooth_series(r["ue"], half)
    H = smooth_series(r["H"], half)
    dth = np.gradient(th, x)
    due = np.gradient(ue, x)
    return 2.0*dth + 2.0*(H + 2.0)*th/ue*due


def bands(runs, edges):
    print("x band          " + "".join(f"| {name:6s}: dth/dx  balance   p_wall    p_e      " for name, _ in runs))
    for lo, hi in zip(edges[:-1], edges[1:]):
        line = f"{lo:7.1f}..{hi:7.1f} "
        for name, r in runs:
            m = np.where((r["x"] >= lo) & (r["x"] < hi))[0]
            i0, i1 = m[0], m[-1]
            dx = r["x"][i1] - r["x"][i0]
            dth = (r["theta"][i1] - r["theta"][i0])/dx
            due = (r["ue"][i1] - r["ue"][i0])/dx
            cf = r["cf"][m].mean()
            bal = cf/2 - (r["H"][m].mean() + 2)*r["theta"][m].mean()/r["ue"][m].mean()*due
            line += f"| {dth: .2e} {bal: .2e} {r['pw'][m].mean(): .2e} {r['pe'][m].mean(): .2e} "
        print(line)


def delta_u_plus(rs, rr, k, crest, re, xs):
    """Hama roughness function at station xs."""
    i = int(np.argmin(abs(rs["x"] - xs)))
    out = {}
    for name, r, y0 in (("smooth", rs, 0.0), ("rough", rr, k)):
        ut = np.sqrt(r["cf"][i]/2)*r["ue"][i]
        yp = (r["y"] - y0)*ut*re
        up = r["prof"][i, :, U]/ut
        out[name] = (ut, yp, up, r["d99"][i]*ut*re)
    uts, yps, ups, d99p = out["smooth"]
    utr, ypr, upr, _ = out["rough"]
    lo = max(30.0, (crest - k)*utr*re*1.2)        # above the crests, in the log region
    hi = 0.2*d99p
    yy = np.geomspace(lo, hi, 40)
    du = np.interp(yy, yps, ups) - np.interp(yy, ypr, upr)
    return dict(x=rs["x"][i], ut_s=uts, ut_r=utr, kplus=k*utr*re, dU=du.mean(), dU_sd=du.std(),
                lo=lo, hi=hi, cf_s=rs["cf"][i], cf_r=rr["cf"][i], d99p=d99p)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--smooth", nargs="+", required=True)
    ap.add_argument("--rough", nargs="+")
    ap.add_argument("--re", type=float, default=450.0)
    ap.add_argument("--ymax", type=float, default=60.0)
    ap.add_argument("--k", type=float, default=0.38348082595870, help="semi-amplitude = mean-plane height")
    ap.add_argument("--crest", type=float, default=2*0.38348082595870)
    ap.add_argument("--stations", type=float, nargs="+", default=[300, 385, 450, 500, 550, 600])
    ap.add_argument("--bands", type=float, nargs="+",
                    default=[100, 120, 160, 200, 250, 300, 385, 450, 500, 550, 600, 630, 640, 645, 648, 649.5, 650])
    a = ap.parse_args()

    runs = []
    for name, files in (("smooth", a.smooth), ("rough", a.rough)):
        if not files:
            continue
        x, y, prof, (t0, t1) = window(files)
        r = integral(x, y, prof, a.re, a.ymax)
        r["cf"] = r["cf_wall"] if name == "smooth" else cf_balance(r)
        r["window"] = (t0, t1)
        runs.append((name, r))
        print(f"{name}: window t = {t0:.1f} .. {t1:.1f}" if np.isfinite(t0) else f"{name}: cumulative to t = {t1:.1f}")

    print("\n== momentum balance per x band (rough c_f = the balance, smoothed) ==")
    bands(runs, a.bands)

    print("\n== along the plate ==")
    print("x        " + "".join(f"| {name}: theta   H      c_f      u_tau   k+   " for name, _ in runs))
    for xs in (100, 150, 200, 250, 300, 385, 450, 500, 550, 600, 640, runs[0][1]["x"][-1]):
        i = int(np.argmin(abs(runs[0][1]["x"] - xs)))
        line = f"{runs[0][1]['x'][i]:7.2f} "
        for name, r in runs:
            ut = np.sqrt(max(r["cf"][i], 0.0)/2)*r["ue"][i]
            line += f"| {r['theta'][i]:.4f} {r['H'][i]:.3f} {r['cf'][i]:.5f} {ut:.4f} {a.k*ut*a.re:5.1f} "
        print(line)

    if len(runs) == 2:
        rs, rr = runs[0][1], runs[1][1]
        print("\n== Hama roughness function above the crests (each side its own u_tau; rough y from the mean plane) ==")
        print("x        k+     u_tau s/r        c_f s/r            dU+ (mean +- sd over y+ in [lo, hi])")
        for xs in a.stations:
            d = delta_u_plus(rs, rr, a.k, a.crest, a.re, xs)
            print(f"{d['x']:7.2f} {d['kplus']:5.1f}  {d['ut_s']:.4f}/{d['ut_r']:.4f}  {d['cf_s']:.5f}/{d['cf_r']:.5f}"
                  f"  {d['dU']:+.2f} +- {d['dU_sd']:.2f}  [{d['lo']:.0f}, {d['hi']:.0f}]")
        print("reference: MacDonald, Chan, Chung, Hutchins & Ooi 2016, full-span channel k+ 10 / lambda+ 113: dU+ = 3.72")


if __name__ == "__main__":
    main()
