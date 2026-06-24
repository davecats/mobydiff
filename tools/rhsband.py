#!/usr/bin/env python3
"""Momentum-operator-truncation gate for the 2:1 interface.

Reads a MOBY_RHSDUMP field ("<prefix>_rhs.h5": un/vn/wn = the discrete momentum
RHS L_h(u) the predictor stored, captured on the pristine `initial = tgv3d`
field), evaluates the ANALYTIC RHS

    L(u)_i = -div(u u)_i + (1/Re) lap(u)_i

at each component's staggered point, and reports the error
discrete - analytic, split into interior / coarse-band / fine-band, per
velocity component. Pass several files (32/64/128) to also print the
convergence ORDER between consecutive resolutions:

    interior order ~ 2   (sanity, run on a uniform grid)
    fine-band  order      = the verdict (a broken interface stencil drops it).

The tgv3d field (k = 2*pi/Lx, on a cube):
    u = sin(kx)cos(ky)cos(kz)
    v = cos(kx)sin(ky)cos(kz)
    w = cos(kx)cos(ky)sin(kz)
lap(u_i) = -3 k^2 u_i, so the analytic diffusion is trivial; the advection is
assembled from the 9 analytic first derivatives (div form, no algebra by hand):
    div(uu)_i = sum_j [ (d_j u_i) u_j + u_i (d_j u_j) ].

ADV vs DIFF split: run a second dump at a HUGE Re (e.g. re=1e30, ire~0) so the
discrete RHS is advection only; pass --re 1e30 and the analytic diffusion
vanishes too -> that file/order isolates the advection operator. The normal-Re
file is the combined operator.

NORMAL vs TANGENTIAL: at a y-normal interface (slab_y) the normal component is
v; the per-component report (u/v/w) is the normal/tangential split. For the 3D
patch every orientation is present, so read all three components.

Usage:
    rhsband.py FILE_rhs.h5 [MORE_rhs.h5 ...] [--re RE] [--band B]
"""
import sys
import h5py
import numpy as np


def parse_args(argv):
    files, re_override, band = [], None, None
    i = 0
    while i < len(argv):
        a = argv[i]
        if a == "--re":
            re_override = float(argv[i + 1]); i += 2
        elif a == "--band":
            band = float(argv[i + 1]); i += 2
        else:
            files.append(a); i += 1
    return files, re_override, band


def analytic_L(var, x, y, z, k, ire):
    """Analytic momentum RHS for component `var` (1=u,2=v,3=w) at (x,y,z)."""
    kx, ky, kz = k * x, k * y, k * z
    sx, cx = np.sin(kx), np.cos(kx)
    sy, cy = np.sin(ky), np.cos(ky)
    sz, cz = np.sin(kz), np.cos(kz)
    # components
    u = sx * cy * cz
    v = cx * sy * cz
    w = cx * cy * sz
    comp = (u, v, w)
    # 9 first derivatives d_j u_i  (rows i = u,v,w ; cols j = x,y,z)
    dudx, dudy, dudz = k * cx * cy * cz, -k * sx * sy * cz, -k * sx * cy * sz
    dvdx, dvdy, dvdz = -k * sx * sy * cz, k * cx * cy * cz, -k * cx * sy * sz
    dwdx, dwdy, dwdz = -k * sx * cy * sz, -k * cx * sy * sz, k * cx * cy * cz
    d = ((dudx, dudy, dudz), (dvdx, dvdy, dvdz), (dwdx, dwdy, dwdz))
    divu = dudx + dvdy + dwdz                       # d_j u_j
    i = var - 1
    # div(u u)_i = sum_j (d_j u_i) u_j + u_i (d_j u_j) = (u.grad)u_i + u_i div u
    adv = (d[i][0] * comp[0] + d[i][1] * comp[1] + d[i][2] * comp[2]
           + comp[i] * divu)
    diff = -3.0 * k * k * comp[i]
    return -adv + ire * diff


def staggered_coords(var, ox, oy, oz, nb, h):
    """(z,y,x) meshgrids of the staggered position of component `var`."""
    p = np.arange(nb)
    # half-cell offset on the two non-staggered axes; face (no offset) on the
    # staggered axis (var 1->x, 2->y, 3->z).
    offx = 0.0 if var == 1 else 0.5
    offy = 0.0 if var == 2 else 0.5
    offz = 0.0 if var == 3 else 0.5
    x = (ox + p + offx) * h
    y = (oy + p + offy) * h
    z = (oz + p + offz) * h
    Z, Y, X = np.meshgrid(z, y, x, indexing="ij")
    return X, Y, Z


def band_masks(blocks, nx, nb, lx, band):
    """Per block: (level, |D|-style band membership mask (nb,nb,nb))."""
    fine = blocks[blocks[:, 3] == 1]
    hf = lx / (nx * 2)
    if fine.size:
        lo = np.array([fine[:, d].min() * hf for d in range(3)])
        hi = np.array([(fine[:, d].max() + nb) * hf for d in range(3)])
    else:
        lo = np.array([np.inf] * 3); hi = np.array([-np.inf] * 3)
    masks = []
    for ox, oy, oz, lev in blocks:
        h = lx / (nx * 2 ** lev)
        xc = (ox + np.arange(nb) + 0.5) * h
        yc = (oy + np.arange(nb) + 0.5) * h
        zc = (oz + np.arange(nb) + 0.5) * h
        inO = ((zc[:, None, None] >= lo[2] - band) & (zc[:, None, None] <= hi[2] + band)
             & (yc[None, :, None] >= lo[1] - band) & (yc[None, :, None] <= hi[1] + band)
             & (xc[None, None, :] >= lo[0] - band) & (xc[None, None, :] <= hi[0] + band))
        inI = ((zc[:, None, None] >= lo[2] + band) & (zc[:, None, None] <= hi[2] - band)
             & (yc[None, :, None] >= lo[1] + band) & (yc[None, :, None] <= hi[1] - band)
             & (xc[None, None, :] >= lo[0] + band) & (xc[None, None, :] <= hi[0] - band))
        masks.append((int(lev), np.broadcast_to(inO & ~inI, (nb, nb, nb))))
    return masks


def analyze(fn, re_override, band):
    h5 = h5py.File(fn, "r")
    lx = float(h5.attrs["lx"]); nx = int(h5.attrs["nx"])
    nb = int(h5.attrs["block_nb_x"]); re = float(h5.attrs["re"])
    blocks = h5["blocks"][...]
    data = {1: h5["un"][...], 2: h5["vn"][...], 3: h5["wn"][...]}  # (nblk,k,j,i)
    h5.close()
    if re_override is not None:
        re = re_override
    ire = 1.0 / re
    k = 2.0 * np.pi / lx
    hc = lx / nx
    if band is None:
        band = 1.01 * hc
    masks = band_masks(blocks, nx, nb, lx, band)

    # region -> component -> list of |error|
    regions = ("interior", "coarse_band", "fine_band")
    acc = {r: {v: [] for v in (1, 2, 3)} for r in regions}
    for bid, (ox, oy, oz, lev) in enumerate(blocks):
        h = lx / (nx * 2 ** lev)
        levb, inband = masks[bid]
        band_key = "fine_band" if lev == 1 else "coarse_band"
        for v in (1, 2, 3):
            X, Y, Z = staggered_coords(v, ox, oy, oz, nb, h)
            ana = analytic_L(v, X, Y, Z, k, ire)
            err = np.abs(data[v][bid] - ana)
            acc[band_key][v].append(err[inband])
            acc["interior"][v].append(err[~inband])

    res = {}
    for r in regions:
        res[r] = {}
        for v in (1, 2, 3):
            a = np.concatenate(acc[r][v]) if acc[r][v] else np.array([0.0])
            a = a[np.isfinite(a)]
            if a.size == 0:
                a = np.array([0.0])
            res[r][v] = (a.max(), np.sqrt(np.mean(a ** 2)), a.size)
    return nx, band, hc, res


def main():
    files, re_override, band = parse_args(sys.argv[1:])
    if not files:
        print(__doc__); sys.exit(1)
    vn = {1: "u", 2: "v", 3: "w"}
    rows = []
    for fn in files:
        nx, b, hc, res = analyze(fn, re_override, band)
        rows.append((fn, nx, res))
        print(f"\nfile={fn}  nx={nx}  band={b:.4f} (hc={hc:.4f})")
        for r in ("interior", "coarse_band", "fine_band"):
            parts = "  ".join(
                f"{vn[v]}: max={res[r][v][0]:.3e} rms={res[r][v][1]:.3e}"
                for v in (1, 2, 3))
            n = res["interior" if r == "interior" else r][1][2]
            print(f"  {r:12s} {parts}")
    # convergence order between consecutive (sorted by nx)
    rows.sort(key=lambda t: t[1])
    if len(rows) >= 2:
        print("\nconvergence order (rms, consecutive nx):")
        for r in ("interior", "coarse_band", "fine_band"):
            line = [f"  {r:12s}"]
            for a, b2 in zip(rows[:-1], rows[1:]):
                for v in (1, 2, 3):
                    ea, eb = a[2][r][v][1], b2[2][r][v][1]
                    if ea > 0 and eb > 0:
                        p = np.log(ea / eb) / np.log(b2[1] / a[1])
                        line.append(f"{vn[v]}{a[1]}->{b2[1]}:{p:.2f}")
                    else:
                        line.append(f"{vn[v]}{a[1]}->{b2[1]}:--")
            print("  ".join(line))


if __name__ == "__main__":
    main()
