#!/usr/bin/env python3
"""Discrete divergence of a block-table snapshot per level, on the cells
whose high faces lie inside the same block (indices 0..nb-2 in each
direction; u/v/w are stored on the LOW faces, datasets (n, z, y, x)).
Solid cells (all six face velocities zero) are skipped. Prints max |div| h
(the dimensionless mass error per cell) and its rms, per level.
Usage: snapshot_divergence.py SNAP.h5 [...]"""
import sys, h5py, numpy as np
for p in sys.argv[1:]:
    f = h5py.File(p, 'r'); b = f['blocks'][...]; lev = b[:, 3]
    h0 = f.attrs['lx']/f.attrs['nx']; n = len(lev)
    mx = np.zeros(6); ss = np.zeros(6); cnt = np.zeros(6)
    U, V, Wd = f['un'], f['vn'], f['wn']
    for s in range(0, n, 10000):
        u = U[s:s+10000]; v = V[s:s+10000]; w = Wd[s:s+10000]; L = lev[s:s+10000]
        du = u[:, :-1, :-1, 1:] - u[:, :-1, :-1, :-1]
        dv = v[:, :-1, 1:, :-1] - v[:, :-1, :-1, :-1]
        dw = w[:, 1:, :-1, :-1] - w[:, :-1, :-1, :-1]
        fluid = (np.abs(u[:, :-1, :-1, 1:]) + np.abs(u[:, :-1, :-1, :-1]) + np.abs(v[:, :-1, 1:, :-1])
                 + np.abs(v[:, :-1, :-1, :-1]) + np.abs(w[:, 1:, :-1, :-1]) + np.abs(w[:, :-1, :-1, :-1])) > 0
        d = np.abs(du + dv + dw)   # = |div| h on a uniform isotropic level
        for l in range(6):
            sel = L == l
            if sel.any():
                dd = d[sel][fluid[sel]]
                if dd.size:
                    mx[l] = max(mx[l], dd.max()); ss[l] += (dd**2).sum(); cnt[l] += dd.size
    print(f"{p}: step {int(f.attrs['step'])} dt {float(f.attrs['dt']):.3e}")
    for l in range(6):
        if cnt[l]: print(f"   L{l}: fluid cells {int(cnt[l]):>10d}  max |div| h {mx[l]:.3e}  rms {np.sqrt(ss[l]/cnt[l]):.3e}")
