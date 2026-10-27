#!/usr/bin/env python3
"""Far-field (level-0) pressure of a block-table snapshot: mean pressure and
mean u in the first chord behind the inlet, in the chord before the outlet,
and the max |<p>| over 1 c slabs -- the amplitude of the global pressure
mode of an approximate projection (2026-10-27 diagnostics).
Usage: farfield_pressure.py SNAP.h5 [...]"""
import sys, h5py, numpy as np
for path in sys.argv[1:]:
    f = h5py.File(path, 'r'); b = f['blocks'][...]; lev = b[:, 3]
    h0 = f.attrs['lx']/f.attrs['nx']; idx = np.where(lev == 0)[0]
    pn, un = f['pn'], f['un']
    P = np.array([pn[i].mean() for i in idx]); U = np.array([un[i].mean() for i in idx])
    x = (b[idx, 0] + 4)*h0
    slabs = [(P[(x >= x0) & (x < x0+1)].mean()) for x0 in range(int(f.attrs['lx']))]
    print(f"{path.split('/')[-1]}: step {int(f.attrs['step'])} t {float(f.attrs['t_current']):.4f} dt {float(f.attrs['dt']):.2e}"
          f"  <p> inlet slab {slabs[0]:+.3f}  outlet slab {slabs[-1]:+.3f}  max|slab| {max(abs(s) for s in slabs):.3f}"
          f"  <u> inlet {U[x < 1].mean():.4f} outlet {U[x >= int(f.attrs['lx'])-1].mean():.4f}")
