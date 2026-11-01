#!/usr/bin/env python3
"""Interpolate a block-table snapshot onto the leaf layout of a deeper-refined
case file of the SAME grid (refine_dims = xyz): the initial field of the
c/2048 calibration case from the c/512 final field (2026-10-29). Piecewise-
constant injection: every target cell takes the value of the source cell
that contains its centre (the source leaf found by strip_loads.Field.locate
at the target leaf's origin; the target leaf is never coarser than the
source leaf under body refinement, which the script checks; a target leaf
the source table does not hold -- a buried leaf kept by the 2:1 coverage
rule of 2026-11-01 -- is filled from the nearest source leaf, cells
clamped). Written leaf
by leaf (the target holds 750 M cells), attrs copied from the source with
dt reset to --dt and the outlet plane `un_xmax` DROPPED (a solver-minted
template carries the plane of its own run; the solver rebuilds it).
Usage: interp_restart3d.py SOURCE.h5 SOURCE_CASE.h5 TARGET_CASE.h5 OUT.h5 [--dt 2e-5]"""
import sys, os, argparse, time
import numpy as np, h5py
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from strip_loads import Field

ap = argparse.ArgumentParser()
ap.add_argument("source"); ap.add_argument("source_case"); ap.add_argument("target_case"); ap.add_argument("out")
ap.add_argument("--dt", type=float, default=2e-5); ap.add_argument("--batch", type=int, default=4000)
a = ap.parse_args()
FIELDS = ["un", "vn", "wn", "pn", "k", "omega", "nut"]
S = Field(a.source, a.source_case); nb = S.nb
with h5py.File(a.target_case, "r") as f:
    tblocks = f["blocks"][...].astype(np.int64)
    for k in ("nx", "ny", "nz", "lx", "ly", "lz", "block_nb_x"):
        if k in f.attrs:
            assert f.attrs[k] == S.f.attrs[k], k
nt = len(tblocks); tlev = tblocks[:, 3]
n_at = np.array([S.n0 * 2 ** l for l in range(int(tlev.max()) + 1)])
# the source leaf holding each target leaf's first cell centre
h_t = S.L[None, :] / n_at[tlev]
origin_pt = (tblocks[:, :3] + 0.5) * h_t
srow, slev = S.locate(origin_pt)
# A target leaf over a region the source table does not hold is a buried
# leaf the source prepare dropped and the target prepare keeps next to a
# coarser leaf (the 2:1 coverage rule, 2026-11-01): fill it from the nearest
# source leaf with every cell index clamped to that leaf, which carries the
# neighbour's solid edge values (u ~ 0, the pressure continued) into it.
miss = np.nonzero(srow < 0)[0]
if len(miss):
    h_s = S.L[None, :] / S.n_at[S.blocks[:, 3]]
    lo = S.blocks[:, :3] * h_s; hi = lo + nb * h_s
    for i in miss:
        c = origin_pt[i] + (nb / 2 - 0.5) * h_t[i]              # the target leaf's centre
        d = np.maximum(lo - c, 0) + np.maximum(c - hi, 0)       # distance to each source leaf's box
        srow[i] = np.argmin((d * d).sum(axis=1)); slev[i] = S.blocks[srow[i], 3]
    print(f"{len(miss)} target leaves outside the source table (buried leaves the source dropped): "
          f"filled from the nearest source leaf, cells clamped")
clamp = np.zeros(nt, dtype=bool); clamp[miss] = True
assert (slev[~clamp] <= tlev[~clamp]).all(), "a target leaf coarser than its source leaf: injection undefined"
print(f"target {nt} leaves (levels {np.bincount(tlev)}); source levels under them {np.bincount(slev)}")
src = {n: S.f[n] for n in FIELDS if n in S.f}
out = h5py.File(a.out, "w")
for k, v in S.f.attrs.items():
    out.attrs[k] = v
out.attrs["dt"] = a.dt
out.create_dataset("blocks", data=tblocks.astype(S.f["blocks"].dtype))
for k in ("x", "y", "z"):
    out.create_dataset(k, data=S.f[k][...])
dsets = {n: out.create_dataset(n, shape=(nt, nb, nb, nb), dtype="f8", chunks=(1, nb, nb, nb)) for n in src}
idx = np.arange(nb)
t0 = time.time()
for s in range(0, nt, a.batch):
    e = min(s + a.batch, nt); r = np.arange(s, e)
    ratio = 2.0 ** (tlev[r] - slev[r])                       # target cells per source cell, per leaf
    # source cell index (0..nb-1) of each target cell, per direction: the
    # source cell holding the target cell's centre, floor((o_t + i + 1/2) / ratio) - o_s
    o_t = tblocks[r, :3]; o_s = S.blocks[srow[r], :3]
    rows = np.unique(srow[r])
    cache = {n: dict(zip(rows, [src[n][int(q)] for q in rows])) for n in src}
    for j, leaf in enumerate(r):
        q = ratio[j]
        ix, iy, iz = (np.floor((o_t[j, d] + idx + 0.5) / q).astype(np.int64) - o_s[j, d] for d in range(3))
        if clamp[leaf]:
            ix, iy, iz = (np.clip(v, 0, nb - 1) for v in (ix, iy, iz))
        assert ix.min() >= 0 and ix.max() < nb and iy.min() >= 0 and iy.max() < nb and iz.min() >= 0 and iz.max() < nb
        for n in src:
            blk = cache[n][srow[leaf]]                       # (z, y, x)
            dsets[n][leaf] = blk[np.ix_(iz, iy, ix)]
    if (s // a.batch) % 25 == 0:
        print(f"  {e}/{nt} leaves, {time.time() - t0:.0f} s", flush=True)
out.close()
print(f"wrote {a.out}: {nt} leaves, dt {a.dt}, un_xmax dropped, {time.time() - t0:.0f} s")
