#!/usr/bin/env python3
"""Where does the 2:1 interface inject error within ONE RK substep?

Reads the MOBY_RKDUMP output (field before the predictor, after the predictor,
after the corrector -- _900001/_900002/_900003.h5) of a uniform single-block run
and a 2:1 block-refined run of the Beltrami flow, on a z = Lz/2 slice. For each
solver variable (u, v, w, p) it writes TWO figures, each a
{single block, refined} x {before, after predictor, after corrector} grid:

  <prefix>_<var>_field.png  -- the variable itself (single block vs refined)
  <prefix>_<var>_error.png  -- the error vs the analytic Beltrami (t=0 reference;
                               one substep barely advances the field, so any
                               structure is numerical -- the interface stands out)

The block grid is overlaid (thin grey) and the refined-region boundary -- the 2:1
interface -- is drawn bold. Comparing the single-block row (smooth) to the
refined row localises exactly which stage and which variable the interface
corrupts.

Usage:
  python3 tools/plot_rkdump.py UNIFORM_DIR REFINED_DIR [--out-prefix rkdump]
"""

from __future__ import annotations

import argparse
import os

import h5py
import numpy as np

VARS = [("un", "u"), ("vn", "v"), ("wn", "w"), ("pn", "p")]
STAGES = [("900001", "before predictor"), ("900002", "after predictor"),
          ("900003", "after corrector")]


def exact_slab(vname, xc, yc, zc, k0):
    """Analytic Beltrami (t=0) for one variable on the slice, shape (j, i).
    zc is the slab's actual cell-centred z (per block; up to h/2 from Lz/2)."""
    sx, cx = np.sin(k0*xc)[None, :], np.cos(k0*xc)[None, :]   # (1, i)
    sy, cy = np.sin(k0*yc)[:, None], np.cos(k0*yc)[:, None]   # (j, 1)
    sz, cz = np.sin(k0*zc), np.cos(k0*zc)                     # scalars
    if vname == "u":   # sin(kz)+cos(ky)
        return np.broadcast_to(sz + cy, (yc.size, xc.size))
    if vname == "v":   # sin(kx)+cos(kz)
        return np.broadcast_to(sx + cz, (yc.size, xc.size))
    if vname == "w":   # sin(ky)+cos(kx)
        return sy + cx
    return -0.5*((sz + cy)**2 + (sx + cz)**2 + (sy + cx)**2)  # p = -1/2|u|^2


def slabs(h5, dset, ztarget):
    """For blocks crossing z=ztarget: (x0,x1,y0,y1,lev, xc, yc, slab[j,i])."""
    lx = float(h5.attrs["lx"]); nx = int(h5.attrs["nx"])
    nbx = int(h5.attrs["block_nb_x"]); nby = int(h5.attrs["block_nb_y"]); nbz = int(h5.attrs["block_nb_z"])
    blocks = h5["blocks"][...]; D = h5[dset]
    out = []
    for bid, (ox, oy, oz, lev) in enumerate(blocks):
        h = lx/(nx*2**lev)
        z0 = oz*h; z1 = (oz + nbz)*h
        if not (z0 <= ztarget < z1):
            continue
        k = min(max(int((ztarget - z0)/h), 0), nbz - 1)
        xc = (ox + np.arange(nbx) + 0.5)*h
        yc = (oy + np.arange(nby) + 0.5)*h
        zc = (oz + k + 0.5)*h
        out.append((ox*h, (ox + nbx)*h, oy*h, (oy + nby)*h, int(lev), xc, yc, zc, D[bid][k]))
    return out


def draw_blocks(ax, patches):
    from matplotlib.patches import Rectangle
    fine = []
    for x0, x1, y0, y1, lev, *_ in patches:
        ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False, edgecolor="0.4", lw=0.25, alpha=0.5))
        if lev > 0:
            fine.append((x0, x1, y0, y1))
    if fine:
        fx0 = min(f[0] for f in fine); fx1 = max(f[1] for f in fine)
        fy0 = min(f[2] for f in fine); fy1 = max(f[3] for f in fine)
        ax.add_patch(Rectangle((fx0, fy0), fx1 - fx0, fy1 - fy0, fill=False, edgecolor="k", lw=1.6))


def panel_data(path, dset, vname, zt, k0, error):
    with h5py.File(path, "r") as h:
        patches = slabs(h, dset, zt)
    rendered = []
    for x0, x1, y0, y1, lev, xc, yc, zc, slab in patches:
        field = slab.astype(float)
        val = field - exact_slab(vname, xc, yc, zc, k0) if error else field
        rendered.append((x0, x1, y0, y1, lev, val))
    return rendered


def figure(rows, dirs, dset, vname, zt, k0, lx, error, out):
    import matplotlib.pyplot as plt
    data = {}  # (r, s) -> list of (x0,x1,y0,y1,lev,val)
    vmax = 0.0
    for r, d in enumerate(dirs):
        for s, (sid, _) in enumerate(STAGES):
            rd = panel_data(os.path.join(d, f"beltrami_{sid}.h5"), dset, vname, zt, k0, error)
            data[(r, s)] = rd
            scale_row = (r == 0) if error else True   # error: scale to single-block? use both
            for *_, val in rd:
                vmax = max(vmax, np.abs(val).max())
    vmax = vmax or 1.0
    fig, axes = plt.subplots(2, 3, figsize=(10.5, 6.6), squeeze=False, constrained_layout=True)
    im = None
    for r, label in enumerate(rows):
        for s, (sid, sname) in enumerate(STAGES):
            ax = axes[r][s]
            for x0, x1, y0, y1, lev, val in data[(r, s)]:
                im = ax.imshow(val, origin="lower", extent=[x0, x1, y0, y1],
                               vmin=-vmax, vmax=vmax, cmap="RdBu_r", aspect="equal")
            draw_blocks(ax, data[(r, s)])
            ax.set_xlim(0, lx); ax.set_ylim(0, lx); ax.set_xticks([]); ax.set_yticks([])
            if r == 0:
                ax.set_title(sname, fontsize=11)
            if s == 0:
                ax.set_ylabel(label, fontsize=11)
    kind = "error vs analytic" if error else "value"
    fig.colorbar(im, ax=axes, shrink=0.85, label=f"{vname} {kind}  (±{vmax:.2e})")
    fig.suptitle(f"Beltrami {vname}: {kind} through one RK substep "
                 f"(z=Lz/2, bold = 2:1 interface)", fontsize=12)
    fig.savefig(out, dpi=140); plt.close(fig)
    print(f"wrote {out}  (±{vmax:.3e})")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("uniform_dir"); p.add_argument("refined_dir")
    p.add_argument("--out-prefix", default="rkdump")
    args = p.parse_args()
    import matplotlib; matplotlib.use("Agg")

    with h5py.File(os.path.join(args.uniform_dir, "beltrami_900001.h5"), "r") as h:
        lx = float(h.attrs["lx"]); k0 = 2.0*np.pi/lx; zt = lx/2.0
    dirs = [args.uniform_dir, args.refined_dir]
    rows = ["single block", "2:1 refined"]
    for dset, vname in VARS:
        figure(rows, dirs, dset, vname, zt, k0, lx, False, f"{args.out_prefix}_{vname}_field.png")
        figure(rows, dirs, dset, vname, zt, k0, lx, True, f"{args.out_prefix}_{vname}_error.png")


if __name__ == "__main__":
    main()
