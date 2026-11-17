#!/usr/bin/env python3
"""The bubble and the transition front of the SD7003 runs from the reduced
snapshots (sd_reduce.py), one or several runs side by side.

    sd_bubble.py --run LABEL 'GLOB.red.npz' [--run ...] [--tmin 14] [--tmax T] [--out viz]

Time average of the span-averaged moments over the snapshots with t >=
tmin; the Reynolds stresses are those of the time-and-span average. Per
run, on the suction side:
  separation / reattachment  the sign changes of the mean tangential
        velocity 1.5 cells above the wall (x/c 0.05..0.98);
  transition  the first x/c where max_y(-u'v') reaches 0.001 U^2 (the
        criterion of Galbraith & Visbal 2010 and of the experiments
        they compare with);
  the amplification curve  max_y of u'u', w'w' against x/c: its level
        ahead of the bubble is the disturbance floor of the run, its
        slope the growth in the detached shear layer.
The wall is the upper side of sd7003.dat (leading edge at (4.5, 6.0),
chord 1, not rotated: the incidence is in the free stream).
"""
import argparse, glob, os
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--run", nargs=2, action="append", metavar=("LABEL", "GLOB"), required=True)
ap.add_argument("--tmin", type=float, default=14.0)
ap.add_argument("--tmax", type=float, default=1e30, help="end of the averaging window (half-window scatter)")
ap.add_argument("--out", default="viz")
ap.add_argument("--re", type=float, default=6.0e4)
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
XLE, YLE = 4.5, 6.0

d = np.loadtxt(os.path.join(os.path.dirname(os.path.abspath(__file__)), "sd7003.dat"), skiprows=1)
ile = np.argmin(d[:, 0])
xu, yu = d[:ile + 1][::-1, 0], d[:ile + 1][::-1, 1]           # upper side, leading to trailing edge
slope = np.gradient(yu, xu)

def average(pattern, tmin, tmax):
    files = sorted(glob.glob(pattern), key=lambda s: int(s.split("_")[-1].split(".")[0]))
    acc, n, ts = None, 0, []
    for fn in files:
        z = np.load(fn)
        if not tmin <= float(z["t"]) <= tmax:
            continue
        if acc is None:
            x, y, h = z["x"], z["y"], float(z["h"])
            acc = {k: np.zeros_like(z["u"]) for k in ("u", "v", "w", "p", "uu", "vv", "ww", "uv")}
        for k in acc:
            acc[k] += z[k]
        n += 1; ts.append(float(z["t"]))
    if not n:
        raise SystemExit("%s: no snapshot at t >= %g" % (pattern, tmin))
    m = {k: v / n for k, v in acc.items()}
    return x, y, h, m, n, ts

def bilinear(f, x, y, h, xp, yp):
    fi = (xp - x[0]) / h; fj = (yp - y[0]) / h
    i = np.clip(np.floor(fi).astype(int), 0, len(x) - 2); j = np.clip(np.floor(fj).astype(int), 0, len(y) - 2)
    s, t = fi - i, fj - j
    return (f[j, i] * (1 - s) * (1 - t) + f[j, i + 1] * s * (1 - t) + f[j + 1, i] * (1 - s) * t + f[j + 1, i + 1] * s * t)

res = {}
for label, pattern in a.run:
    x, y, h, m, n, ts = average(pattern, a.tmin, a.tmax)
    xc = (x - XLE)
    on = (xc > 0.02) & (xc < 0.995)
    ys = YLE + np.interp(xc[on], xu, yu); sl = np.interp(xc[on], xu, slope)
    tx, ty = 1 / np.hypot(1, sl), sl / np.hypot(1, sl)
    dist = 1.5 * h                                            # wall distance of the near-wall sample
    xs, yp = x[on] - dist * ty, ys + dist * tx
    ut = bilinear(m["u"], x, y, h, xs, yp) * tx + bilinear(m["v"], x, y, h, xs, yp) * ty
    cf = 2 * ut / dist / a.re
    # stresses of the time-and-span average, maxima over the fluid above the wall
    Y = y[:, None] * np.ones_like(x)[None, :]
    above = np.zeros_like(m["u"], dtype=bool)
    above[:, on] = Y[:, on] > (ys + 2 * h)[None, :]
    uv = np.where(above, -(m["uv"] - m["u"] * m["v"]), np.nan)
    uu = np.where(above, m["uu"] - m["u"] ** 2, np.nan)
    ww = np.where(above, m["ww"] - m["w"] ** 2, np.nan)
    with np.errstate(all="ignore"):
        uvm, uum, wwm = (np.nanmax(q[:, on], axis=0) for q in (uv, uu, ww))
    X = xc[on]
    win = (X > 0.05) & (X < 0.98)
    sgn = np.sign(ut[win]); k = np.nonzero(np.diff(sgn) != 0)[0]
    cross = X[win][k] + (X[win][k + 1] - X[win][k]) * ut[win][k] / (ut[win][k] - ut[win][k + 1])
    down = [c for c, i in zip(cross, k) if sgn[i] > 0]; up = [c for c, i in zip(cross, k) if sgn[i] < 0]
    xs_, xr_ = (down[0] if down else np.nan), (up[-1] if up else np.nan)
    it = np.nonzero((uvm > 1e-3) & (X > 0.05))[0]
    xt_ = X[it[0]] if len(it) else np.nan
    floor = np.exp(np.nanmean(np.log(np.maximum(uum[(X > 0.10) & (X < 0.20)], 1e-300))))
    floor_w = np.exp(np.nanmean(np.log(np.maximum(wwm[(X > 0.10) & (X < 0.20)], 1e-300))))
    print("%-8s %3d snapshots t %.2f..%.2f: separation %.3f, transition (-u'v' 1e-3) %.3f, reattachment %.3f, "
          "bubble %.3f c; floor x/c 0.1..0.2: u'u' %.2e, w'w' %.2e; max u'u' %.3f"
          % (label, n, ts[0], ts[-1], xs_, xt_, xr_, xr_ - xs_, floor, floor_w, np.nanmax(uum)))
    np.savetxt(os.path.join(a.out, "bubble_%s.csv" % label), np.c_[X, cf, uvm, uum, wwm], delimiter=",",
               header="x/c, c_f (u_t at 1.5 cells), max -u'v', max u'u', max w'w'", fmt="%.6e")
    res[label] = (X, cf, uvm, uum, wwm, xs_, xt_, xr_)

import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
COL = ["#2a78d6", "#d6452a", "#2a9d5c", "#8a4fd6", "#555555"]
fig, ax = plt.subplots(3, 1, figsize=(8, 9.5), sharex=True)
for c, (label, (X, cf, uvm, uum, wwm, xs_, xt_, xr_)) in zip(COL, res.items()):
    ax[0].plot(X, cf, color=c, lw=1.6, label="%s: sep %.3f, trans %.3f, reatt %.3f" % (label, xs_, xt_, xr_))
    ax[1].semilogy(X, np.maximum(uum, 1e-16), color=c, lw=1.6, label=label)
    ax[2].semilogy(X, np.maximum(wwm, 1e-16), color=c, lw=1.6, label=label)
    ax[1].semilogy(X, np.maximum(uvm, 1e-16), color=c, lw=1.0, ls="--")
ax[0].axhline(0, color="0.5", lw=0.8); ax[0].set_ylabel("$c_f$ (from $u_t$ 1.5 cells off the wall)"); ax[0].set_ylim(-0.012, 0.03)
ax[0].legend(frameon=False, fontsize=9, title="suction side, mean over t and span")
ax[1].axhline(1e-3, color="0.5", lw=0.8, ls=":"); ax[1].set_ylabel("max$_y$ $u'u'$ (solid), $-u'v'$ (dashed)")
ax[2].set_ylabel("max$_y$ $w'w'$"); ax[2].set_xlabel("x/c")
for q in ax:
    q.grid(alpha=0.25); q.set_xlim(0, 1)
fig.tight_layout(); fig.savefig(os.path.join(a.out, "bubble.png"), dpi=140)
print("wrote", os.path.join(a.out, "bubble.png"))
