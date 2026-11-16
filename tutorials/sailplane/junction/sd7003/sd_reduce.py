#!/usr/bin/env python3
"""Span-averaged moments of one SD7003 snapshot on the finest level.

    sd_reduce.py FIELD.h5 [FIELD.h5 ...]       ->  FIELD.red.npz each

Every leaf of the finest level whose block lies in the window (the suction
side, the bubble and the near wake) contributes the z-average of u, v, w, p
and of uu, vv, ww, uv over its cells to a 2D array on the finest lattice
(NaN where the finest level is absent or a block is buried). The velocities
are taken at their own staggered positions: half a cell is irrelevant for
the bubble length and the transition front read from these arrays
(sd_bubble.py averages them in time). A snapshot is 2.5 GB, its reduction
2 MB: the run directory keeps the reductions and the restart states only.
"""
import sys
import numpy as np, h5py

WIN = (4.40, 5.70, 5.90, 6.20)               # x0, x1, y0, y1 (chords)

for name in sys.argv[1:]:
    with h5py.File(name, "r") as f:
        blk = f["blocks"][...].astype(np.int64)
        nb = int(f.attrs["block_nb_x"])
        n0 = np.array([int(f.attrs[k]) for k in ("nx", "ny", "nz")])
        L = np.array([float(f.attrs[k]) for k in ("lx", "ly", "lz")])
        lmax = int(blk[:, 3].max())
        h = L[0] / (n0[0] * 2 ** lmax)
        i0, i1 = int(WIN[0] / h) // nb * nb, -(-int(WIN[1] / h) // nb) * nb
        j0, j1 = int(WIN[2] / h) // nb * nb, -(-int(WIN[3] / h) // nb) * nb
        nzf = n0[2] * 2 ** lmax
        rows = np.nonzero((blk[:, 3] == lmax) & (blk[:, 0] >= i0) & (blk[:, 0] < i1)
                          & (blk[:, 1] >= j0) & (blk[:, 1] < j1))[0]
        names = ("u", "v", "w", "p", "uu", "vv", "ww", "uv")
        S = {k: np.zeros((j1 - j0, i1 - i0)) for k in names}
        cnt = np.zeros((j1 - j0, i1 - i0))
        # hdf5 wants increasing indices; read in chunks of rows
        for a in range(0, len(rows), 4096):
            r = rows[a:a + 4096]
            q = {k: f[d][r] for k, d in (("u", "un"), ("v", "vn"), ("w", "wn"), ("p", "pn"))}   # (rows, z, y, x)
            q["uu"], q["vv"], q["ww"], q["uv"] = q["u"] ** 2, q["v"] ** 2, q["w"] ** 2, q["u"] * q["v"]
            for m, row in enumerate(r):
                i, j = blk[row, 0] - i0, blk[row, 1] - j0
                for k in names:
                    S[k][j:j + nb, i:i + nb] += q[k][m].sum(axis=0)
                cnt[j:j + nb, i:i + nb] += nb
        full = cnt == nzf
        out = {k: np.where(full, S[k] / np.maximum(cnt, 1), np.nan) for k in names}
        t = float(f.attrs["t_current"]) if "t_current" in f.attrs else np.nan
        step = int(f.attrs["step"]) if "step" in f.attrs else -1
    np.savez_compressed(name[:-3] + ".red.npz", x=(i0 + np.arange(i1 - i0) + 0.5) * h,
                        y=(j0 + np.arange(j1 - j0) + 0.5) * h, h=h, t=t, step=step, **out)
    print("%s: t %.4f, %d finest leaves in the window, %d of %d columns complete"
          % (name, t, len(rows), full.sum(), full.size))
