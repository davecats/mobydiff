#!/usr/bin/env python3
"""pre-F1 vs F1 developed KMM180 statistics (20 t.u. windows, same restart)."""
import h5py, numpy as np, sys
U, V, W, UU, VV, WW, UV = 0, 1, 2, 3, 4, 5, 6
def load(fn):
    with h5py.File(fn) as f:
        y = f["coord"][...]; p = f["profile"][...]; c = f["count"][...]; re = float(f.attrs["re"]); t = float(f.attrs["t_current"]); step = int(f.attrs["step"])
    keep = c > 0; y, p = y[keep], p[keep]
    h = 0.5 * (y.max() + y.min()) if False else 1.0          # ly = 2, walls at 0 and 2 -> h = 1
    utau = 1.0                                               # forcing_x = 1, h = 1: tau_w = 1
    yw = np.minimum(y, 2.0 - y); yp = yw * re * utau
    Up = p[:, U] / utau
    urms = np.sqrt(np.maximum(p[:, UU] - p[:, U]**2, 0)); vrms = np.sqrt(np.maximum(p[:, VV] - p[:, V]**2, 0)); wrms = np.sqrt(np.maximum(p[:, WW] - p[:, W]**2, 0))
    uv = -(p[:, UV] - p[:, U] * p[:, V])
    return dict(y=y, yp=yp, Up=Up, u=urms, v=vrms, w=wrms, uv=uv, t=t, step=step, re=re)
a, b = load(sys.argv[1]), load(sys.argv[2])
print(f"pre-F1: step {a['step']} t {a['t']:.2f} | F1: step {b['step']} t {b['t']:.2f}  (Re_tau {a['re']:.0f})")
def fold(d, q):  # both walls folded by y+ (rows come in symmetric pairs)
    o = np.argsort(d["yp"]); return d["yp"][o], q[o]
print(f"{'quantity':14s} {'pre-F1':>10s} {'F1':>10s} {'diff':>10s}  where")
n = len(a["y"]); mid = n // 2
for name, q, where in [("U+ centreline", "Up", "y = 1"), ("u'+ peak", "u", "max"), ("v'+ peak", "v", "max"), ("w'+ peak", "w", "max"), ("-<u'v'>+ peak", "uv", "max")]:
    if where == "y = 1":
        va, vb = 0.5 * (a["Up"][mid - 1] + a["Up"][mid]), 0.5 * (b["Up"][mid - 1] + b["Up"][mid]); loc = "y+ 180"
    else:
        ia, ib = np.argmax(a[q]), np.argmax(b[q]); va, vb = a[q][ia], b[q][ib]; loc = f"y+ {a['yp'][ia]:.1f}/{b['yp'][ib]:.1f}"
    print(f"{name:14s} {va:10.4f} {vb:10.4f} {vb - va:10.2e}  {loc}")
Ub_a = np.trapz(a["Up"], a["y"]) / 2.0; Ub_b = np.trapz(b["Up"], b["y"]) / 2.0
print(f"{'bulk U+':14s} {Ub_a:10.4f} {Ub_b:10.4f} {Ub_b - Ub_a:10.2e}")
print("\nmax |F1 - pre| over the profile, and the profile's own scale:")
for name, q in [("U+", "Up"), ("u'+", "u"), ("v'+", "v"), ("w'+", "w"), ("-<u'v'>+", "uv")]:
    d = np.abs(b[q] - a[q]); print(f"  {name:10s} max diff {d.max():.3e} at y+ {a['yp'][np.argmax(d)]:.1f}   (scale {np.abs(a[q]).max():.3f})")
print("\nlog-law slice (both walls folded), U+ at y+ ~ 30, 60, 100, 180:")
for target in (30, 60, 100, 180):
    ia = np.argmin(np.abs(a["yp"] - target)); print(f"  y+ {a['yp'][ia]:6.1f}: pre {a['Up'][ia]:.3f}  F1 {b['Up'][ia]:.3f}  loglaw {np.log(a['yp'][ia])/0.41 + 5.2:.3f}")
