#!/usr/bin/env python3
"""Extrema of a block-table snapshot: dt / Courant / Peclet attributes, max
|u|, |v|, |w|, |p|, nut/nu, the cell count above a threshold and where the
largest values sit (chord units). The diagnostic of the calibration box's
cold start (2026-10-27).  Usage: snapshot_extrema.py SNAP.h5 [SNAP2.h5 ...]
"""
import sys, h5py, numpy as np

def extrema(path, thr=3.0):
    f = h5py.File(path, 'r')
    b = f['blocks'][...]; lev = b[:, 3]
    h0 = f.attrs['lx']/f.attrs['nx']; nu = 1/f.attrs['re']
    out = {'step': int(f.attrs['step']), 't': float(f.attrs['t_current']),
           'dt': float(f.attrs['dt']), 'cfl': f.attrs['cfl']}
    for name in ('un', 'vn', 'wn', 'pn', 'nut'):
        d = f[name]; n = d.shape[0]; best = (0.0, 0, (0, 0, 0)); cnt = 0
        for s in range(0, n, 20000):
            x = d[s:s+20000]; a = np.abs(x).reshape(x.shape[0], -1)
            m = a.max(1); i = int(m.argmax())
            if m[i] > best[0]:
                best = (float(m[i]), s+i, np.unravel_index(int(a[i].argmax()), x[i].shape))
            if name != 'nut': cnt += int((a > thr).sum())
        v, ib, (kz, jy, ix) = best; h = h0/2**lev[ib]
        pos = ((b[ib, 0]+ix+0.5)*h, (b[ib, 1]+jy+0.5)*h, (b[ib, 2]+kz+0.5)*h)
        if name == 'nut': v /= nu
        out[name] = (v, int(lev[ib]), pos, cnt)
    return out

for p in sys.argv[1:]:
    o = extrema(p)
    print(f"{p}: step {o['step']} t {o['t']:.4f} dt {o['dt']:.3e} (Courant, Peclet) {o['cfl']}")
    for k in ('un', 'vn', 'wn', 'pn', 'nut'):
        v, L, pos, cnt = o[k]
        extra = f"  cells>3: {cnt}" if k != 'nut' else " (nut/nu)"
        print(f"   {k}: max {v:10.3f} at L{L} ({pos[0]:.3f}, {pos[1]:.3f}, {pos[2]:.3f}){extra}")
