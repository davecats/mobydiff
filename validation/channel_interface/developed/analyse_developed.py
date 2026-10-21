#!/usr/bin/env python3
"""Developed band-refined channel: the statistics of one or more runs as a table.

Reads a run's channel_stats.h5 (+ _l1 beside it), combines the levels (rows
with count > 0), folds the two walls by y+ and prints U+ centreline, bulk U+,
the rms / shear-stress peaks against KMM / Moser Re_tau 180, and the
INTERFACE-BAND diagnostic: at the 2:1 interface (y = 0.5986 and its mirror)
the jump of each profile between the last fine row and the first coarse row,
against the variation over the two rows on either side -- a ratio near 1 is
a smooth profile, a ratio well above 1 an interface band.

    python3 analyse_developed.py runs/cheb12/stats/channel_stats.h5 [runs/redblack/stats/channel_stats.h5 ...]
"""
import glob, os, sys
import h5py, numpy as np

U, V, W, UU, VV, WW, UV = 0, 1, 2, 3, 4, 5, 6
KMM = {"U+ centreline": 18.3, "u'+ peak": 2.66, "v'+ peak": 0.84, "w'+ peak": 1.08, "-<u'v'>+ peak": 0.73}


def load(level0):
    ys, ps = [], []
    for fn in [level0] + sorted(glob.glob(level0.replace(".h5", "_l*.h5"))):
        with h5py.File(fn) as f:
            y, p, c = f["coord"][...], f["profile"][...], f["count"][...]
            re, t, step = float(f.attrs["re"]), float(f.attrs["t_current"]), int(f.attrs["step"])
        keep = c > 0
        ys.append(y[keep]); ps.append(p[keep])
    y = np.concatenate(ys); p = np.concatenate(ps); o = np.argsort(y); y, p = y[o], p[o]
    d = dict(y=y, yp=np.minimum(y, 2.0 - y) * re, Up=p[:, U], re=re, t=t, step=step)
    d["u"] = np.sqrt(np.maximum(p[:, UU] - p[:, U] ** 2, 0))
    d["v"] = np.sqrt(np.maximum(p[:, VV] - p[:, V] ** 2, 0))
    d["w"] = np.sqrt(np.maximum(p[:, WW] - p[:, W] ** 2, 0))
    d["uv"] = -(p[:, UV] - p[:, U] * p[:, V])
    return d


def band(d, q, yi):
    """jump across the interface at y = yi vs the variation two rows either side."""
    below = np.where(d["y"] < yi)[0][-2:]; above = np.where(d["y"] > yi)[0][:2]
    jump = abs(d[q][above[0]] - d[q][below[-1]])
    local = 0.5 * (abs(d[q][below[-1]] - d[q][below[0]]) + abs(d[q][above[1]] - d[q][above[0]]))
    return jump, local


runs = [(a, load(a)) for a in sys.argv[1:]]
w = max(len(os.path.dirname(os.path.dirname(a))) for a, _ in runs)
print(f"{'quantity':16s}" + "".join(f"{os.path.dirname(os.path.dirname(a)):>{w + 2}s}" for a, _ in runs) + f"{'KMM':>8s}")
for name, q in [("U+ centreline", "Up"), ("bulk U+", "Up"), ("u'+ peak", "u"), ("v'+ peak", "v"),
                ("w'+ peak", "w"), ("-<u'v'>+ peak", "uv")]:
    vals = []
    for _, d in runs:
        if name == "U+ centreline":
            n = len(d["y"]); vals.append(0.5 * (d["Up"][n // 2 - 1] + d["Up"][n // 2]))
        elif name == "bulk U+":
            vals.append(float(np.trapezoid(d["Up"], d["y"]) / 2.0))
        else:
            vals.append(float(d[q].max()))
    print(f"{name:16s}" + "".join(f"{v:>{w + 2}.4f}" for v in vals) + f"{KMM.get(name, float('nan')):>8.2f}")
for _, d in runs:
    print(f"  (t = {d['t']:.1f}, step {d['step']}, {len(d['y'])} rows)")
print("\ninterface band at y = 0.5986 / 1.4014: jump vs local variation (ratio)")
for a, d in runs:
    out = []
    for name, q in [("U+", "Up"), ("u'+", "u"), ("v'+", "v"), ("w'+", "w"), ("-<u'v'>+", "uv")]:
        r = []
        for yi in (0.598617012688, 1.401382987312):
            j, l = band(d, q, yi); r.append(j / max(l, 1e-300))
        out.append(f"{name} {max(r):.2f}")
    print(f"  {os.path.dirname(os.path.dirname(a)):>{w}s}: " + "  ".join(out))
