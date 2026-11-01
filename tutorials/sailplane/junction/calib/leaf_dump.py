#!/usr/bin/env python3
"""The cells of one leaf of a block-table snapshot and the leaves around it:
where a blow-up starts (the c/2048 calibration box, 2026-11-01). Prints the
leaf's level, index and extent, the 8 x 8 x 8 values of one field as z-planes
(the datasets are stored (leaf, z, y, x)), and every leaf whose extent touches
this one with its level and max |field| -- a halo nobody writes shows as a
plane of garbage against sane neighbours, an interface error as a jump on
the face toward a leaf of another level.
Usage: leaf_dump.py FILE.h5 LEAF [--field un] [--planes 0,7] [--cells]"""
import argparse
import numpy as np, h5py

ap = argparse.ArgumentParser()
ap.add_argument("file"); ap.add_argument("leaf", type=int)
ap.add_argument("--field", default="un"); ap.add_argument("--planes", default="0,3,7")
ap.add_argument("--cells", action="store_true", help="print the full 8x8x8 block of every neighbour too")
a = ap.parse_args()
f = h5py.File(a.file, "r")
blocks = f["blocks"][...].astype(np.int64); lev = blocks[:, 3]
nb = f[a.field].shape[1]
N = np.array([f.attrs[k] for k in ("nx", "ny", "nz")], dtype=float)
L = np.array([f.attrs[k] for k in ("lx", "ly", "lz")], dtype=float)
h = L[None, :] / (N[None, :] * 2.0 ** lev[:, None])            # cell size per leaf
lo = blocks[:, :3] * h; hi = lo + nb * h
r = a.leaf
print(f"leaf {r}: level {lev[r]}, index {tuple(blocks[r, :3])}, h {h[r, 0]:.3e}, "
      f"extent x {lo[r, 0]:.4f}..{hi[r, 0]:.4f} y {lo[r, 1]:.4f}..{hi[r, 1]:.4f} z {lo[r, 2]:.4f}..{hi[r, 2]:.4f}")
d = f[a.field][r]
np.set_printoptions(linewidth=200, precision=3)
for k in [int(p) for p in a.planes.split(",")]:
    print(f"  {a.field} plane k = {k} (rows y, columns x):\n{d[k]}")
# neighbours: extents overlapping the leaf's extent widened by half a cell
eps = 0.5 * h[r]
touch = np.all((lo <= hi[r] + eps) & (hi >= lo[r] - eps), axis=1); touch[r] = False
nbrs = np.where(touch)[0]
print(f"  {len(nbrs)} touching leaves (levels {np.bincount(lev[nbrs], minlength=lev.max() + 1)}):")
for q in nbrs:
    side = []
    for dname, dd in zip("xyz", range(3)):
        if hi[q, dd] <= lo[r, dd] + eps[dd]: side.append(f"{dname}-")
        elif lo[q, dd] >= hi[r, dd] - eps[dd]: side.append(f"{dname}+")
    dq = f[a.field][q]
    print(f"    leaf {q:8d} level {lev[q]} side {' '.join(side) or 'overlap'} max|{a.field}| {np.abs(dq).max():.3e}")
    if a.cells:
        for k in [int(p) for p in a.planes.split(",")]:
            print(f"      plane k = {k}:\n{dq[k]}")
