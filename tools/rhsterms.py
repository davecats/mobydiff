#!/usr/bin/env python3
"""Term-by-term momentum-operator check at the 2:1 interface (MOBY_TERMDUMP).

Reads "<prefix>_adv.h5" (un/vn/wn = the x/y/z advection flux-divergence terms of
component `var`) and "<prefix>_dif.h5" (un/vn/wn = the x/y/z Laplacian terms,
times 1/Re), evaluates each term's ANALYTIC value for the tgv3d field at the
component's staggered point, and reports the error of EACH term separately,
split interior / fine-band, and (for the interface-facing row) per y-row of the
fine blocks. This pins a broken interface stencil to the exact term.

Divergence-form advection term in direction d for component i:
    adv_d = -[ (d_d u_d) * u_i  +  u_d * (d_d u_i) ]
Laplacian term in direction d:  dif_d = (1/Re) d2(u_i)/dx_d2 = -(1/Re) k^2 u_i.

Usage: rhsterms.py PREFIX_adv.h5 PREFIX_dif.h5 [--var 2] [--re RE] [--rows]
   or: rhsterms.py P32 P64 P128 ... (multiple prefixes for order; pass _adv/_dif)
"""
import sys
import h5py
import numpy as np

TERMS = ["adv_x", "adv_y", "adv_z", "dif_x", "dif_y", "dif_z"]
DSET = ["un", "vn", "wn"]


def parse(argv):
    files, var, re_o, rows = [], 2, None, False
    i = 0
    while i < len(argv):
        a = argv[i]
        if a == "--var": var = int(argv[i+1]); i += 2
        elif a == "--re": re_o = float(argv[i+1]); i += 2
        elif a == "--rows": rows = True; i += 1
        else: files.append(a); i += 1
    return files, var, re_o, rows


def derivs(x, y, z, k):
    sx, cx = np.sin(k*x), np.cos(k*x)
    sy, cy = np.sin(k*y), np.cos(k*y)
    sz, cz = np.sin(k*z), np.cos(k*z)
    comp = (sx*cy*cz, cx*sy*cz, cx*cy*sz)            # u, v, w
    d = ((k*cx*cy*cz, -k*sx*sy*cz, -k*sx*cy*sz),      # du/dx,du/dy,du/dz
         (-k*sx*sy*cz, k*cx*cy*cz, -k*cx*sy*sz),      # dv/...
         (-k*sx*cy*sz, -k*cx*sy*sz, k*cx*cy*cz))      # dw/...
    return comp, d


def analytic_terms(var, x, y, z, k, ire):
    comp, d = derivs(x, y, z, k)
    ui = comp[var-1]
    out = {}
    for dd in range(3):                              # advection x/y/z
        out[f"adv_{'xyz'[dd]}"] = -(d[dd][dd]*ui + comp[dd]*d[var-1][dd])
    diff = -ire * k*k * ui
    for dd in range(3):
        out[f"dif_{'xyz'[dd]}"] = diff
    return out


def stag(var, ox, oy, oz, nb, h):
    p = np.arange(nb)
    x = (ox + p + (0 if var == 1 else .5)) * h
    y = (oy + p + (0 if var == 2 else .5)) * h
    z = (oz + p + (0 if var == 3 else .5)) * h
    Z, Y, X = np.meshgrid(z, y, x, indexing="ij")
    return X, Y, Z


def load(prefix_adv, prefix_dif):
    a = h5py.File(prefix_adv, "r"); d = h5py.File(prefix_dif, "r")
    meta = (float(a.attrs["lx"]), int(a.attrs["nx"]), int(a.attrs["block_nb_x"]),
            float(a.attrs["re"]), a["blocks"][...])
    disc = {}
    for i, t in enumerate(["adv_x", "adv_y", "adv_z"]):
        disc[t] = a[DSET[i]][...]
    for i, t in enumerate(["dif_x", "dif_y", "dif_z"]):
        disc[t] = d[DSET[i]][...]
    a.close(); d.close()
    return meta, disc


def analyze(prefix_adv, prefix_dif, var, re_o, rows):
    (lx, nx, nb, re, blk), disc = load(prefix_adv, prefix_dif)
    if re_o is not None: re = re_o
    ire = 1.0/re; k = 2*np.pi/lx
    band = 1.01*lx/nx
    fine = blk[blk[:, 3] == 1]
    hf = lx/(nx*2)
    lo = np.array([fine[:, dd].min()*hf for dd in range(3)])
    hi = np.array([(fine[:, dd].max()+nb)*hf for dd in range(3)])
    ylo = fine[:, 1].min() if fine.size else None
    acc = {t: {"int": [], "fb": []} for t in TERMS}
    rowacc = {t: {} for t in TERMS}                  # row-> list (bottom fine blocks)
    for bid, (ox, oy, oz, lev) in enumerate(blk):
        if lev != 1:
            continue
        h = lx/(nx*2)
        X, Y, Z = stag(var, ox, oy, oz, nb, h)
        ana = analytic_terms(var, X, Y, Z, k, ire)
        xc = (ox+np.arange(nb)+.5)*h; yc = (oy+np.arange(nb)+.5)*h; zc = (oz+np.arange(nb)+.5)*h
        inO = ((zc[:, None, None] >= lo[2]-band) & (zc[:, None, None] <= hi[2]+band)
             & (yc[None, :, None] >= lo[1]-band) & (yc[None, :, None] <= hi[1]+band)
             & (xc[None, None, :] >= lo[0]-band) & (xc[None, None, :] <= hi[0]+band))
        inI = ((zc[:, None, None] >= lo[2]+band) & (zc[:, None, None] <= hi[2]-band)
             & (yc[None, :, None] >= lo[1]+band) & (yc[None, :, None] <= hi[1]-band)
             & (xc[None, None, :] >= lo[0]+band) & (xc[None, None, :] <= hi[0]-band))
        m = np.broadcast_to(inO & ~inI, (nb, nb, nb))
        for t in TERMS:
            e = np.abs(disc[t][bid] - ana[t])
            acc[t]["fb"].append(e[m]); acc[t]["int"].append(e[~m])
            if rows and oy == ylo:                    # orientation-B bottom blocks
                for jrow in range(nb):
                    rowacc[t].setdefault(jrow, []).append(e[:, jrow, :].ravel())
    res = {}
    for t in TERMS:
        fb = np.concatenate(acc[t]["fb"]) if acc[t]["fb"] else np.array([0.])
        it = np.concatenate(acc[t]["int"]) if acc[t]["int"] else np.array([0.])
        res[t] = (np.sqrt(np.mean(it**2)), np.sqrt(np.mean(fb**2)))
    rowres = {t: {jr: np.sqrt(np.mean(np.concatenate(v)**2))
                  for jr, v in rowacc[t].items()} for t in TERMS} if rows else None
    return nx, res, rowres


def main():
    files, var, re_o, rows = parse(sys.argv[1:])
    advs = [f for f in files if "_adv" in f]
    difs = [f for f in files if "_dif" in f]
    runs = []
    for fa, fd in zip(advs, difs):
        nx, res, rowres = analyze(fa, fd, var, re_o, rows)
        runs.append((nx, res, rowres))
        print(f"\n=== var={var} nx={nx}  (interior / fine-band rms error per term) ===")
        for t in TERMS:
            print(f"  {t}: interior={res[t][0]:.3e}   fine_band={res[t][1]:.3e}")
        if rows and rowres:
            print("  fine-band-bottom (orientation-B) per y-row j=1..nb:")
            for t in TERMS:
                nbrows = sorted(rowres[t])
                s = " ".join(f"{rowres[t][jr]:.2e}" for jr in nbrows)
                print(f"    {t}: {s}")
    runs.sort(key=lambda r: r[0])
    if len(runs) >= 2:
        print("\n=== fine-band convergence order per term ===")
        for t in TERMS:
            seg = []
            for a, b in zip(runs[:-1], runs[1:]):
                ea, eb = a[1][t][1], b[1][t][1]
                p = np.log(ea/eb)/np.log(b[0]/a[0]) if ea > 0 and eb > 0 else float('nan')
                seg.append(f"{a[0]}->{b[0]}:{p:.2f}")
            print(f"  {t}: " + "  ".join(seg))


if __name__ == "__main__":
    main()
