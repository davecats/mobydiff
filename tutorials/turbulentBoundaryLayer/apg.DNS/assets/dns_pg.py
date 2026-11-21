#!/usr/bin/env python3
"""Boundary-layer integrals and beta of a DNS statistics window, and the
one-step correction of the top-pressure table.

    dns_pg.py STATS_LATE.h5 [STATS_EARLY.h5] [--every 50]
    dns_pg.py assets/mobydiff/xyz_3200_384_272/data.nc     (the committed window)
    dns_pg.py STATS_LATE.h5 STATS_EARLY.h5 --correct OUTDIR [--target 0.5]
              [--value 1.087] [--kappa 1.3] [--relax 0.9]

The statistics file of the boundaryLayer case is CUMULATIVE (raw sums and
sample counts per (x, y) cell, z averaged); the window between two copies is
their difference (one file: everything since the run started). Its means U, V
and P (accumulators 0, 1, 9; make_mobydiff_nc.py) go through the RANS
read-out ../../smooth.RANS/assets/postpro/pg_bl.py unchanged: the Griffin
edge on the local free-stream total pressure, delta_99, delta* and theta to
delta_99, tau_w from the two first cell centres, beta = delta*/tau_w dP_e/dx.

--correct writes OUTDIR/ptop.txt: the table that the RANS response model of
top_optimise.py predicts to bring beta to the target on x = 200 .. 550 (the
DNS read-out range; the RANS fit ran to 640). The current coefficients are
those of ../../apg.RANS/opt/it14 times --value (the ini's y_max_p_value: the
table times the value is the applied top pressure, and p_top is linear in the
coefficients); d beta = kappa (delta*/tau_w) R da on THIS window's state, kappa
the secant gain of the last RANS iterations (0.9-1.4). The written table is
divided by --value again, so the ini keeps its y_max_p_value and only
ptop.txt is replaced.
"""
import argparse
import json
import os
import sys

import h5py
import numpy as np
from scipy.signal import savgol_filter

HERE = os.path.dirname(os.path.abspath(__file__))
TBL = os.path.normpath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(TBL, "smooth.RANS", "assets", "postpro"))
import pg_bl  # noqa: E402
import top_optimise as to  # noqa: E402

U_, V_, P_ = 0, 1, 9                    # accumulator columns (make_mobydiff_nc.py)
OPT = os.path.join(TBL, "apg.RANS", "opt")


def window(late, early=None):
    """The mean field of the window as pg_bl wants it (x, y and (nx, ny) U, V,
    pn) plus the Reynolds stresses uu, vv, ww, uv and k.

    LATE is either a cumulative statistics file (with EARLY: the window
    between the two copies) or the exported data.nc of a window
    (../../smooth/assets/postpro/make_mobydiff_nc.py, the committed form)."""
    if late.endswith(".nc"):
        import netCDF4
        with netCDF4.Dataset(late) as d:
            g = lambda k: np.asarray(d[k][...], dtype=float)  # noqa: E731
            w = dict(x=g("x"), y=g("y"), U=g("u_mean"), V=g("v_mean"), pn=g("p_mean"),
                     uu=g("uu_stress"), vv=g("vv_stress"), ww=g("ww_stress"), uv=g("uv_stress"),
                     re=float(d.re_dstar), t1=float(d.time_window_end if "time_window_end" in d.ncattrs() else d.time),
                     t0=float(d.time_window_start) if "time_window_start" in d.ncattrs() else None)
    else:
        with h5py.File(late, "r") as f:
            nx, ny = int(f.attrs["nx"]), int(f.attrs["ny"])
            s, c = f["raw_sum"][...], f["count"][...].astype(float)
            x, y, re = f["xcoord"][...], f["ycoord"][...], float(f.attrs["re"])
            t1 = float(f.attrs["t_current"])
        t0 = None
        if early:
            with h5py.File(early, "r") as f:
                s, c = s - f["raw_sum"][...], c - f["count"][...]
                t0 = float(f.attrs["t_current"])
        m = (s / c[:, None]).reshape(nx, ny, -1)
        U, V, W = m[..., U_], m[..., V_], m[..., 2]
        w = dict(x=x, y=y, U=U, V=V, pn=m[..., P_], re=re, t0=t0, t1=t1,
                 uu=m[..., 3] - U * U, vv=m[..., 4] - V * V, ww=m[..., 5] - W * W, uv=m[..., 6] - U * V)
    w["k"] = 0.5 * (w["uu"] + w["vv"] + w["ww"])
    return w


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("late")
    ap.add_argument("early", nargs="?")
    ap.add_argument("--every", type=float, default=50.0)
    ap.add_argument("--correct", metavar="OUTDIR")
    ap.add_argument("--target", type=float, default=0.5)
    ap.add_argument("--value", type=float, default=1.087)
    ap.add_argument("--kappa", type=float, default=1.3)
    ap.add_argument("--relax", type=float, default=0.9)
    a = ap.parse_args()

    w = window(a.late, a.early)
    q = pg_bl.pg_integrals(w)
    x = q["x"]
    print("window t %s .. %.0f" % ("start" if w["t0"] is None else "%.0f" % w["t0"], w["t1"]))
    print(f"{'x':>6} {'Re_theta':>8} {'H':>6} {'c_f':>8} {'U_e':>6} {'delta99':>7} {'beta':>7} {'K*1e6':>6}")
    for xs in np.arange(a.every, x[-1], a.every):
        i = np.argmin(abs(x - xs))
        print(f"{x[i]:6.0f} {q['Re_theta'][i]:8.1f} {q['H'][i]:6.3f} {q['cf'][i]:8.5f} {q['Ue'][i]:6.3f} "
              f"{q['delta99'][i]:7.2f} {q['beta'][i]:+7.3f} {q['K'][i] * 1e6:+6.2f}")
    win = (x >= 200) & (x <= 550)
    b = q["beta"][win]
    print(f"beta on x 200-550: mean {b.mean():+.3f}  min {b.min():+.3f}  max {b.max():+.3f}")
    if not a.correct:
        return

    cfg = json.load(open(os.path.join(OPT, "case.json")))
    p_ramp, bump = cfg["p_ramp"], cfg.get("bump")
    npoly = cfg["npoly"] + (1 if bump else 0)
    a_cur = a.value * to.load_coeffs(os.path.join(OPT, "it14", "coeffs.txt"), npoly)
    pr = to.PressureResponse()
    R = []
    for m in range(npoly):
        e = np.zeros(npoly)
        e[m] = 1.0
        R.append(np.interp(x, pr.xg, np.gradient(pr(to.top_pressure(pr.xg, e, p_ramp, bump)), pr.xg)))
    R = np.array(R).T
    n = int(round(33.0 / (x[1] - x[0]))) | 1      # the RANS model's smoothing length
    model = (savgol_filter(q["dstar"], n, 2) / savgol_filter(q["tauw"], n, 2))[:, None] * R
    bt = a.target * to.ramp(x, *cfg["x_ramp"])
    da = np.linalg.lstsq(a.kappa * model[win], (bt - q["beta"])[win], rcond=None)[0]
    a_new = a_cur + a.relax * da
    pred = q["beta"] + a.kappa * model @ (a_new - a_cur)
    os.makedirs(a.correct, exist_ok=True)
    xt = np.linspace(0.0, to.LX, 1301)
    np.savetxt(os.path.join(a.correct, "coeffs.txt"), a_new,
               header=f"applied coefficients (= table x {a.value}); dp_top/dx = S sum_m a_m xi^m [+ a_bump (S_bump - S)]")
    np.savetxt(os.path.join(a.correct, "ptop.txt"), np.c_[xt, to.top_pressure(xt, a_new, p_ramp, bump) / a.value],
               fmt="%.10e", header=f"x  p_top / {a.value}  (DNS correction from the window t {w['t0']}..{w['t1']:.0f}, "
               f"target {a.target:+.3f}, kappa {a.kappa}, relax {a.relax}; p_top(650) = 0)")
    print(f"a: {np.array2string(a_cur, precision=6)} -> {np.array2string(a_new, precision=6)}")
    print(f"predicted beta on x 200-550: mean {pred[win].mean():+.3f}  min {pred[win].min():+.3f}  "
          f"max {pred[win].max():+.3f}; wrote {a.correct}/ptop.txt (use with y_max_p_value = {a.value})")


if __name__ == "__main__":
    main()
