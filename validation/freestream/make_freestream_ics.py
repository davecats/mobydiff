#!/usr/bin/env python3
"""Build the freestream-gate restart ICs from solver-minted templates.

The restart metadata carries periodicity/BC rows, so a periodic-run field
cannot be restarted directly into an inflow/outflow ini: mint a 0-step
template WITH the target ini (run_gates.sh does it), then overwrite the
velocity/pressure block rows here.

  pois:   copy the developed pois_ref field rows into the pois_io template
          (same grid, same nb, same single-level block table)  -> IC_pois.h5
  vortex: superpose a Lamb-Oseen vortex on the uniform (1,0,0) template
          (staggered evaluation, quasi-2D)                     -> IC_vortex.h5

--mirror builds the x-mirror image of either (x -> lx - x, u -> -u): the
initial condition of the same flow run in -x with the outlet at x_min (the
low-side outlet legs). The template must then come from the mirrored ini.

The normal velocity ON the outlet face is part of the state. With the outlet
at x_min it is face 0 of un; with the outlet at x_max it is not in un and
goes into the un_xmax dataset (io.f90 write_outlet_planes), which both
generators write: an IC without it starts that face from the zero-gradient
value instead, and the two orientations would not start from the same field.
"""
import argparse
import shutil

import h5py
import numpy as np

NAMES = ("un", "vn", "wn", "pn")


def staggered_positions(nodes, var):
    cent = [0.5*(n[:-1] + n[1:]) for n in nodes]
    pos = list(cent)
    if var < 3:
        pos[var] = nodes[var][:-1]
    return pos


def block_scatter(f, name, global_field):
    """Overwrite dataset rows (nblocks, nbz, nby, nbx) from a global (z,y,x) array."""
    blocks = f["blocks"][...]
    data = f[name]
    nbz, nby, nbx = data.shape[1:]
    for bid, (ox, oy, oz, lev) in enumerate(blocks):
        assert lev == 0
        data[bid] = global_field[oz:oz + nbz, oy:oy + nby, ox:ox + nbx]


def write_xmax_plane(f, plane):
    """un_xmax: the outlet-face u(z, y) at x = lx, one row per block (the
    block's own (y, z) patch, y fastest; rows of blocks away from x_max are
    unused and left zero)."""
    blocks = f["blocks"][...]
    nbz, nby, nbx = f["un"].shape[1:]
    nx = int(f.attrs["nx"])
    rows = np.zeros((len(blocks), nby*nbz))
    for bid, (ox, oy, oz, lev) in enumerate(blocks):
        assert lev == 0
        if ox + nbx == nx:
            rows[bid] = plane[oz:oz + nbz, oy:oy + nby].reshape(-1)   # (z, y): y fastest
    if "un_xmax" in f:
        del f["un_xmax"]
    f.create_dataset("un_xmax", data=rows)


def block_gather(f, name):
    """The global (z,y,x) array of a single-level block dataset."""
    blocks = f["blocks"][...]
    data = f[name]
    nbz, nby, nbx = data.shape[1:]
    shape = [int(f.attrs[a]) for a in ("nz", "ny", "nx")]
    g = np.zeros(shape)
    for bid, (ox, oy, oz, lev) in enumerate(blocks):
        assert lev == 0
        g[oz:oz + nbz, oy:oy + nby, ox:ox + nbx] = data[bid]
    return g


def cmd_pois_mirror(args):
    """pois, mirrored: flow in -x, the outlet face at x = 0."""
    G = 0.08
    shutil.copy(args.template, args.out)
    with h5py.File(args.src, "r") as s, h5py.File(args.out, "r+") as d:
        src = {name: block_gather(s, name) for name in NAMES}
        # u faces: face i of the image is face nx - i of the (periodic) source
        u = -np.roll(src["un"][:, :, ::-1], 1, axis=2)
        p = src["pn"][:, :, ::-1]
        x = d["x"][...]
        shift = G*float(x[1] - x[0])/2.0 - float(np.mean(p[:, :, 0]))
        block_scatter(d, "un", u)
        block_scatter(d, "vn", src["vn"][:, :, ::-1])
        block_scatter(d, "wn", src["wn"][:, :, ::-1])
        block_scatter(d, "pn", p + shift)
    print(f"wrote {args.out} (mirrored rows from {args.src}, p shifted by {shift:+.3e})")


def cmd_pois(args):
    if args.mirror:
        return cmd_pois_mirror(args)
    G = 0.08                    # pois_ref forcing_x = the developed -dp/dx
    shutil.copy(args.template, args.out)
    with h5py.File(args.src, "r") as s, h5py.File(args.out, "r+") as d:
        for name in NAMES:
            d[name][...] = s[name][...]
        # The periodic run's p level is arbitrary (all-Neumann null mode); the
        # outlet pins the level, so pre-shift p to the pinned line
        # p(x) = -G*(x - lx) (0 at the outlet FACE) to avoid a spurious
        # level-adjustment transient.
        x = d["x"][...]
        h_last = float(x[-1] - x[-2])
        blocks = d["blocks"][...]
        p = d["pn"]
        # outlet column = blocks whose ox is max, local last x row
        ox_max = blocks[:, 0].max()
        vals = [p[bid][:, :, -1] for bid, blk in enumerate(blocks) if blk[0] == ox_max]
        shift = G*h_last/2.0 - float(np.mean(vals))
        p[...] = p[...] + shift
        # the outlet face: the periodic image of the source's face 0
        write_xmax_plane(d, block_gather(s, "un")[:, :, 0])
    print(f"wrote {args.out} (rows from {args.src}, p shifted by {shift:+.3e})")


def cmd_vortex(args):
    GAMMA, RC, XC, YC, UINF = 0.443, 0.15, 0.75, 1.0, 1.0
    shutil.copy(args.template, args.out)
    with h5py.File(args.out, "r+") as f:
        nodes = (f["x"][...], f["y"][...], f["z"][...])
        nz = len(nodes[2]) - 1
        for var, name in enumerate(NAMES):
            if name == "wn":
                continue
            px, py, _ = staggered_positions(nodes, var)
            X, Y = np.meshgrid(px, py, indexing="xy")   # (ny, nx)
            if args.mirror:                             # evaluate the image point
                X = nodes[0][-1] - X
            dx, dy = X - XC, Y - YC
            r2 = dx*dx + dy*dy
            r2 = np.maximum(r2, 1e-30)
            # u_theta/r, finite at the axis
            swirl = GAMMA/(2.0*np.pi*r2)*(1.0 - np.exp(-r2/RC**2))
            if name == "un":
                plane = UINF - swirl*dy
            elif name == "vn":
                plane = swirl*dx
            else:                                       # pn: leave the template value
                continue
            if args.mirror and name == "un":
                plane = -plane
            g = np.repeat(plane[np.newaxis, :, :], nz, axis=0)  # (z, y, x)
            block_scatter(f, name, g)
            if name == "un" and not args.mirror:
                # the outlet face x = lx, evaluated like every other u face
                dxo, dyo = nodes[0][-1] - XC, py - YC
                r2o = np.maximum(dxo*dxo + dyo*dyo, 1e-30)
                uo = UINF - GAMMA/(2.0*np.pi*r2o)*(1.0 - np.exp(-r2o/RC**2))*dyo
                write_xmax_plane(f, np.repeat(uo[np.newaxis, :], nz, axis=0))
    print(f"wrote {args.out} (Lamb-Oseen at ({XC},{YC}), rc={RC}, Gamma={GAMMA})")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("pois")
    p.add_argument("--src", default="pois_ref_final.h5")
    p.add_argument("--template", default="pois_template.h5")
    p.add_argument("--out", default="IC_pois.h5")
    p.add_argument("--mirror", action="store_true")
    p.set_defaults(func=cmd_pois)
    v = sub.add_parser("vortex")
    v.add_argument("--template", default="vortex_template.h5")
    v.add_argument("--out", default="IC_vortex.h5")
    v.add_argument("--mirror", action="store_true")
    v.set_defaults(func=cmd_vortex)
    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
