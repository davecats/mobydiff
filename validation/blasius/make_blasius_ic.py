#!/usr/bin/env python3
"""Fill a solver-minted restart template with the Blasius similarity field.

The impulsive uniform start seeds grid-scale ringing that takes very long
to flush; starting from the Blasius field (with the same virtual origin as
the inlet profile) removes the transient and turns the run into a sharp
gate: the solution must STAY on Blasius up to discretization error.

  1) mint a template (1 step, no restart):
       mpirun -n 4 ../../build_cpu/moby_solve template.ini
  2) ./make_blasius_ic.py            (template_1.h5 -> IC_blasius.h5)
  3) run blasius2d.ini (its [restart] points at IC_blasius.h5)

u/v rows follow the staggered storage: u[..,i] at x_node[i] (the inlet
face included), v[..,j] at y_node[j], cell centres elsewhere. p is left at
the template value (Blasius has uniform p; the outlet pins the level).

The two OUTLET faces (x_max and the top) are state too, and index nb+1 is
not in un/vn: they go into un_xmax / vn_ymax (io.f90 write_outlet_planes),
filled here with the same analytic field. The template carries those
datasets from its own one-step run; left as they are they would restart the
outlet faces from that unrelated field.
"""
import argparse
import os
import shutil
import sys

import h5py
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare_blasius import solve_blasius, blasius_eval  # noqa: E402

BETA = 0.664114672


def staggered_positions(nodes, var):
    pos = [0.5*(n[:-1] + n[1:]) for n in nodes]
    if var < 3:
        pos[var] = nodes[var][:-1]
    return pos


def block_scatter(f, name, global_field):
    blocks = f["blocks"][...]
    data = f[name]
    nbz, nby, nbx = data.shape[1:]
    for bid, (ox, oy, oz, lev) in enumerate(blocks):
        assert lev == 0
        data[bid] = global_field[oz:oz + nbz, oy:oy + nby, ox:ox + nbx]


def write_outlet_plane(f, name, plane, axis):
    """One row per block of the outlet-normal velocity on a HIGH domain face
    (axis 0 = x_max, plane (z, y); axis 1 = y_max, plane (z, x)); the first
    tangential index runs fastest, rows of blocks away from the face are
    unused."""
    blocks = f["blocks"][...]
    nbz, nby, nbx = f["un"].shape[1:]
    nglob = (int(f.attrs["nx"]), int(f.attrs["ny"]))[axis]
    nb_n, nb_t = ((nbx, nby), (nby, nbx))[axis]
    rows = np.zeros((len(blocks), nb_t*nbz))
    for bid, (ox, oy, oz, lev) in enumerate(blocks):
        assert lev == 0
        on, ot = ((ox, oy), (oy, ox))[axis]
        if on + nb_n == nglob:
            rows[bid] = plane[oz:oz + nbz, ot:ot + nb_t].reshape(-1)
    if name in f:
        del f[name]
    f.create_dataset(name, data=rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--template", default="template_1.h5")
    ap.add_argument("--out", default="IC_blasius.h5")
    ap.add_argument("--theta", type=float, default=1.0)
    args = ap.parse_args()

    eta_b, f_b, fp_b = solve_blasius()
    shutil.copy(args.template, args.out)
    with h5py.File(args.out, "r+") as f:
        nodes = (f["x"][...], f["y"][...], f["z"][...])
        re = float(f.attrs["re"])
        nu, uinf = 1.0/re, 1.0
        x_v = (uinf*args.theta/nu)*args.theta/BETA**2

        nz = len(nodes[2]) - 1
        for var, name in enumerate(("un", "vn")):
            px, py, _ = staggered_positions(nodes, var)
            X, Y = np.meshgrid(px, py, indexing="xy")       # (ny, nx)
            eta = Y*np.sqrt(uinf/(nu*(X + x_v)))
            f_e, fp_e = blasius_eval(eta, eta_b, f_b, fp_b)
            if name == "un":
                plane = uinf*fp_e
            else:
                vscale = 0.5*np.sqrt(nu*uinf/(X + x_v))
                plane = vscale*(eta*fp_e - f_e)
            g = np.repeat(plane[np.newaxis, :, :], nz, axis=0)
            block_scatter(f, name, g)
            # the outlet face of this component: u at x = lx, v at y = ly
            if name == "un":
                xo, yo = nodes[0][-1], py
            else:
                xo, yo = px, nodes[1][-1]
            eta = yo*np.sqrt(uinf/(nu*(xo + x_v)))
            f_e, fp_e = blasius_eval(eta, eta_b, f_b, fp_b)
            if name == "un":
                line = uinf*fp_e
            else:
                line = 0.5*np.sqrt(nu*uinf/(xo + x_v))*(eta*fp_e - f_e)
            write_outlet_plane(f, ("un_xmax", "vn_ymax")[var],
                               np.repeat(line[np.newaxis, :], nz, axis=0), var)
        f["wn"][...] = 0.0
        f["pn"][...] = 0.0
    print(f"wrote {args.out} (Blasius field, x_v = {x_v:.2f})")


if __name__ == "__main__":
    main()
