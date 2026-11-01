#!/usr/bin/env python3
"""Extrema of a block-table snapshot, leaf by leaf, for a restart that blows
up (the c/2048 calibration box, 2026-11-01: the injected c/512 field dies
within a few steps). For every field: global min / max, the number of
leaves holding a non-finite value, and the leaves with the largest |value|
(level, index, physical origin), so the blow-up can be placed: at the
body, at the level-7 box edge, at the outlet. Streams the file in batches
of leaves (an 82 GB file on a 64 GB node).
Usage: field_extrema.py FILE.h5 [--top 8] [--batch 4000] [--fields un,vn,wn,pn,k,omega,nut]"""
import argparse, time
import numpy as np, h5py

ap = argparse.ArgumentParser()
ap.add_argument("file"); ap.add_argument("--top", type=int, default=8)
ap.add_argument("--batch", type=int, default=4000)
ap.add_argument("--fields", default="un,vn,wn,pn,k,omega,nut")
a = ap.parse_args()
f = h5py.File(a.file, "r")
blocks = f["blocks"][...].astype(np.int64); n = len(blocks); lev = blocks[:, 3]
nb = int(f.attrs["block_nb_x"]) if "block_nb_x" in f.attrs else f["un"].shape[1]
N = np.array([f.attrs[k] for k in ("nx", "ny", "nz")], dtype=float)
L = np.array([f.attrs[k] for k in ("lx", "ly", "lz")], dtype=float)
h_leaf = L[None, :] / (N[None, :] / nb * 2.0 ** lev[:, None])     # leaf edge lengths
origin = blocks[:, :3] * h_leaf
print(f"{a.file}: {n} leaves, levels {np.bincount(lev)}, attrs step {f.attrs.get('step')} "
      f"t {f.attrs.get('t_current')} dt {f.attrs.get('dt')}")
fields = [x for x in a.fields.split(",") if x in f]
amax = {x: np.zeros(n) for x in fields}; gmin = {x: np.inf for x in fields}; gmax = {x: -np.inf for x in fields}
nonfin = {x: np.zeros(n, dtype=bool) for x in fields}
t0 = time.time()
for s in range(0, n, a.batch):
    e = min(s + a.batch, n)
    for x in fields:
        d = f[x][s:e]
        fin = np.isfinite(d)
        nonfin[x][s:e] = ~fin.all(axis=(1, 2, 3))
        dd = np.where(fin, d, 0.0)
        amax[x][s:e] = np.abs(dd).max(axis=(1, 2, 3))
        gmin[x] = min(gmin[x], dd.min()); gmax[x] = max(gmax[x], dd.max())
    if (s // a.batch) % 50 == 0:
        print(f"  {e}/{n} leaves, {time.time() - t0:.0f} s", flush=True)
for x in fields:
    print(f"\n{x}: min {gmin[x]:.4e} max {gmax[x]:.4e} (finite values); leaves with non-finite values: "
          f"{nonfin[x].sum()} (by level {np.bincount(lev[nonfin[x]], minlength=lev.max() + 1)})")
    order = np.argsort(-amax[x])[:a.top]
    for r in order:
        print(f"  leaf {r:8d} level {lev[r]} index {tuple(blocks[r, :3])} origin "
              f"({origin[r, 0]:.4f}, {origin[r, 1]:.4f}, {origin[r, 2]:.4f}) max|{x}| {amax[x][r]:.4e}")
    big = amax[x] > 10.0
    if big.any():
        print(f"  leaves with max|{x}| > 10: {big.sum()} (by level {np.bincount(lev[big], minlength=lev.max() + 1)}); "
              f"origin box x {origin[big, 0].min():.3f}..{origin[big, 0].max():.3f} "
              f"y {origin[big, 1].min():.3f}..{origin[big, 1].max():.3f} z {origin[big, 2].min():.3f}..{origin[big, 2].max():.3f}")
print(f"done in {time.time() - t0:.0f} s")
