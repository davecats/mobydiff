#!/usr/bin/env python3
"""Last column of a rough-TBL snapshot: outlet-face velocity (un_xmax), the dimple
faces (fluid face on the plane with a solid face one cell upstream), the stored
pressure of the last column vs an interior column."""
import sys, h5py, numpy as np
fn = sys.argv[1]
k_amp = 0.38348082595870; off = k_amp; nwx, nwz = 150, 6; lx, lz = 650.0, 26.0
f = h5py.File(fn, 'r')
blocks = f['blocks'][...]           # origin(3)+level per block
x, y, z = f['x'][...], f['y'][...], f['z'][...]
nbx, nby, nbz = f.attrs['block_nb_x'], f.attrs['block_nb_y'], f.attrs['block_nb_z']
nx, ny, nz = f.attrs['nx'], f.attrs['ny'], f.attrs['nz']
yc = 0.5*(y[:-1]+y[1:]); zc = 0.5*(z[:-1]+z[1:])
def h(xx, zz):
    return off + k_amp*np.sin(2*np.pi*nwx*xx/lx + np.pi/2)*np.sin(2*np.pi*nwz*zz/lz + np.pi/2)
# assemble the outlet plane (ny, nz) and the last columns
uf = np.full((ny, nz), np.nan); u_nx = np.full((ny, nz), np.nan); p_last = np.full((ny, nz), np.nan)
p_int = np.full((ny, nz), np.nan); p_mid = np.full((ny, nz), np.nan)
for b in range(blocks.shape[0]):
    ox, oy, oz = blocks[b, :3]
    if ox + nbx != nx: continue
    js = slice(oy, oy+nby); ks = slice(oz, oz+nbz)
    uf[js, ks] = f['un_xmax'][b].reshape(nbz, nby).T
    u_nx[js, ks] = f['un'][b, :, :, nbx-1].T       # face one cell upstream of the plane
    p_last[js, ks] = f['pn'][b, :, :, nbx-1].T
    p_int[js, ks] = f['pn'][b, :, :, nbx-50].T
    p_mid[js, ks] = f['pn'][b, :, :, nbx//2].T
# the plane faces: solid where the surface at x = Lx is above the face centre
hs_plane = h(lx, zc)[None, :]
solid_plane = yc[:, None] < hs_plane
# pocket faces: fluid on the plane, solid somewhere in the same (y,z) row within
# half a roughness wavelength upstream (the trough behind the plane is a dead end)
xc = 0.5*(x[:-1]+x[1:]); lam = lx/nwx; up = xc[xc > lx - lam/2]
hs_up = np.max(np.array([h(xi, zc) for xi in up]), axis=0)[None, :]
dimple = (~solid_plane) & (yc[:, None] < hs_up)
print(f"t = {f.attrs['t_current']:.1f}  step {f.attrs['step']}")
print(f"outlet plane: {uf.size} faces, solid on the plane {solid_plane.sum()}, pocket faces (fluid plane / solid within lambda/2 upstream) {dimple.sum()}")
print(f"u_f on solid plane faces: max|u| {np.abs(uf[solid_plane]).max():.2e}")
fl = ~solid_plane
print(f"u_f fluid faces: mean {uf[fl].mean():.4f}, min {uf[fl].min():.4f}, inflow faces {np.sum(uf[fl]<0)} ({100*np.mean(uf[fl]<0):.2f} %), mean of inflow {uf[fl][uf[fl]<0].mean() if np.any(uf[fl]<0) else 0:.4f}")
if dimple.any():
    d = uf[dimple]
    print(f"pocket faces: mean {d.mean():+.4f}, rms {d.std():.4f}, min {d.min():+.4f}, max {d.max():+.4f}, inflow {np.sum(d<0)} of {d.size}")
    for yl, yh in ((0, 0.1), (0.1, 0.2), (0.2, 0.4), (0.4, 0.6), (0.6, 0.8)):
        m = dimple & (yc[:, None] >= yl) & (yc[:, None] < yh)
        if m.any(): print(f"   y in [{yl},{yh}): n {m.sum():4d} mean {uf[m].mean():+.4f} min {uf[m].min():+.4f}")
# near-wall fluid faces not dimple, y < 2k
m = fl & ~dimple & (yc[:, None] < 2*k_amp)
print(f"other fluid faces below the crests (y < {2*k_amp:.3f}): n {m.sum()}, mean {uf[m].mean():+.4f}, min {uf[m].min():+.4f}")
def prms(p, name):
    pf = p[~solid_plane]
    zm = np.nanmean(np.where(solid_plane, np.nan, p), axis=1, keepdims=True)
    print(f"{name}: mean {pf.mean():+.3e}, rms about mean {pf.std():.3e}, rms about the span mean {np.nanstd(np.where(solid_plane, np.nan, p - zm)):.3e}")
prms(p_last, "p last column  ")
prms(p_int,  "p column nx-50 ")
prms(p_mid,  "p column nx/2  ")
