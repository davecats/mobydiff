#!/usr/bin/env python3
"""The pocket gate: a rough wall whose trough sits ON the outlet plane.

The last column at the spanwise trough is a dimple open to the outlet: the
row-3 outlet face has a neighbour held near zero by its coefficient and a
solid face behind that. With the 2026-10-01 face equation alone that face
grew as inflow and blew up at t = 34 (docs/next_session_rough_outlet.md,
section 4). This reads the last snapshot and gates: every field finite, the
divergence small, and the outlet face in the dimple rows bounded (the rows
above it read the physical near-wall 0.01-0.03; the face itself settles at a
few 1e-2 of inflow).

    python3 check_pocket.py <snapshot> [--case pocket.case.h5] [--umax 0.5] [--div 1e-4]
"""
import argparse, sys
import h5py, numpy as np

ap = argparse.ArgumentParser()
ap.add_argument('snapshot'); ap.add_argument('--case', default='pocket.case.h5')
ap.add_argument('--umax', type=float, default=0.5, help='bound on |u| at the outlet plane in the dimple rows')
ap.add_argument('--div', type=float, default=1e-3, help='bound on the interior rms divergence (1.0e-4 at t = 50, still the cold-start transient at red-black 6; 7.9e-8 by t = 150; 1e31 blown up)')
a = ap.parse_args()

c = h5py.File(a.case, 'r'); blocks = c['blocks'][:]
f = h5py.File(a.snapshot, 'r')
ok = True
for name in ('un', 'vn', 'wn', 'pn', 'un_xmax'):
    if not np.isfinite(f[name][:]).all():
        print('NOT FINITE:', name); ok = False
if not ok:
    print('pocket gate: FAIL'); sys.exit(1)

# the wall block on the outlet plane; the trough is at mid-span
b0 = int(np.where((blocks[:, 0] == blocks[:, 0].max()) & (blocks[:, 1] == 0) & (blocks[:, 2] == 0))[0][0])
coef = c['coef_blocks'][b0]; nbx = coef.shape[0] - 2; nbz = coef.shape[2] - 2
k = nbz // 2
ux = f['un_xmax'][b0].reshape(nbz, -1)
dimple = ux[k - 1, 0:6]
print('outlet face at the trough, rows 1..6: ' + ' '.join('%9.2e' % v for v in dimple))
print('coefficient on the plane there:        ' + ' '.join('%9.2e' % v for v in coef[nbx + 1, 1:7, k, 0]))
print('coefficient one cell behind:           ' + ' '.join('%9.2e' % v for v in coef[nbx, 1:7, k, 0]))
umax = np.abs(dimple).max()
print('max|u| in the dimple rows = %.3e (bound %.2e)' % (umax, a.umax))
if umax > a.umax: ok = False

# interior divergence of the whole field (cells with six stored faces)
x = c['x_nodes'][:]; y = c['y_nodes'][:]; z = c['z_nodes'][:]
un = f['un'][:]; vn = f['vn'][:]; wn = f['wn'][:]
nb = un.shape[3], un.shape[2], un.shape[1]
rms = 0.0; n = 0
for b in range(un.shape[0]):
    ox, oy, oz = blocks[b, 0], blocks[b, 1], blocks[b, 2]
    dx = np.diff(x[ox:ox + nb[0] + 1]); dy = np.diff(y[oy:oy + nb[1] + 1]); dz = np.diff(z[oz:oz + nb[2] + 1])
    du = (un[b, :, :, 1:] - un[b, :, :, :-1]) / dx[None, None, :-1]
    dv = (vn[b, :, 1:, :] - vn[b, :, :-1, :]) / dy[None, :-1, None]
    dw = (wn[b, 1:, :, :] - wn[b, :-1, :, :]) / dz[:-1, None, None]
    d = du[:-1, :-1, :] + dv[:-1, :, :-1] + dw[:, :-1, :-1]
    rms += (d ** 2).sum(); n += d.size
rms = np.sqrt(rms / n)
print('interior rms divergence = %.3e (bound %.1e)' % (rms, a.div))
if rms > a.div: ok = False
print('pocket gate:', 'PASS' if ok else 'FAIL')
sys.exit(0 if ok else 1)
