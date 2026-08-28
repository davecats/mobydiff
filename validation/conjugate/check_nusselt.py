#!/usr/bin/env python3
"""The conjugate interface-heat (Nusselt) diagnostic (increment C3, gate 2).

C1 left the conjugate flux columns of `scalar_stats.f90` SMOKE-gated only.
This puts numbers on them, three ways, each independent of the others:

  `sum`     an independent Python transcription of the SAME discrete sum, from
            the case file's phi and the snapshot's theta. Solver vs Python, so
            it catches a transcription slip in either -- the S4 idiom.

  `slab`    a CLOSED FORM. A slab with an insulated outer face and a
            volumetric source in the solid must, at steady state, deliver
            every watt it generates across the interface:
                H = C_s S y_w A .
            No reference run, no discretisation tolerance -- and because the
            source is weighted by the fluid fraction and that fraction is
            EXACT for a plane, the identity is exact too. (The idiom of the
            S5a ibmwf180 gate: build the case so its budget telescopes.)

  `budget`  the CONTROL-VOLUME cross-check, the way A2 validated C_L/C_D
            against a Gauss border flux -- but stated so it holds at any time
            rather than only at a steady state. Over the set of cells the
            solver classifies FLUID, the discrete flux form telescopes and
            leaves exactly the interface faces, so

                d/dt sum_fluid C theta dV = H  -  (flux out through a plane)

            with the plane flux (convective AND diffusive) evaluated in Python
            from the snapshot. Three snapshots give the centred dE/dt, so the
            statement is second order in dt rather than an assumption of
            steadiness -- which matters because a converged conjugate channel
            is expensive and a transient one is not.

            LANDMINE inherited from A2: the border flux reads the VELOCITY, so
            it is only as good as the projection's divergence. The driver runs
            it at two `niter` and reports max|div u| with the result.

    ./check_nusselt.py sum    field.h5 --case case.h5 --kappa K [--heat h.txt]
    ./check_nusselt.py slab   --heat h.txt --y-wall Y --kappa K --cap C
                              --source S --area A
    ./check_nusselt.py budget a.h5 b.h5 c.h5 --heat h.txt --cap C --wavy
                              [--plane Y]
"""

from __future__ import annotations

import argparse
import math
import os
import sys

import numpy as np
import h5py

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scalar"))
from scalar_tools import BlockGeometry                       # noqa: E402
from check_oblique import load_phi                           # noqa: E402
from check_conjugate import wavy_wall_height                 # noqa: E402

MIN_COSINE = 5.0e-2          # scalar.f90 CONJ_MIN_COSINE


def global_field(h5, name):
    """A single-level block-table dataset as one global (k, j, i) array.

    Field files carry the node lines as x/y/z and the block size as
    block_nb_{x,y,z}; CASE files (load_phi above) use x_nodes/... and
    block_nb. Keeping the two straight is not cosmetic -- C1's weight checker
    could not tell the two index orders apart because its slab is symmetric.
    """
    blocks = h5["blocks"][...]
    nb = [int(h5.attrs["block_nb_x"]), int(h5.attrs["block_nb_y"]),
          int(h5.attrs["block_nb_z"])]
    if int(blocks[:, 3].max()) != 0:
        raise SystemExit("single-level snapshots only")
    d = h5[name][...]
    n = [int(blocks[:, m].max()) + nb[m] for m in range(3)]
    out = np.zeros((n[2], n[1], n[0]))
    for bid in range(blocks.shape[0]):
        ox, oy, oz = (int(v) for v in blocks[bid, :3])
        out[oz:oz + nb[2], oy:oy + nb[1], ox:ox + nb[0]] = d[bid]
    return out


def cmd_sum(a):
    """The interface sum, transcribed independently from the case file."""
    phi, centres = load_phi(a.case)
    with h5py.File(a.field, "r") as h5:
        th = global_field(h5, a.name)
        node = [h5["x"][...], h5["y"][...], h5["z"][...]]
    dm = 1.0 / (a.re * a.pr)
    # cell sizes and the areas of the three face directions
    dl = [np.diff(n) for n in node]
    total = 0.0
    nfaces = 0
    for d in range(3):
        axis = 2 - d
        # the face between global cell m-1 and m along `axis`; phi carries one
        # ghost, so phi index m+1 is cell m.
        pl = np.take(phi, range(1, phi.shape[axis] - 2), axis=axis)
        pr = np.take(phi, range(2, phi.shape[axis] - 1), axis=axis)
        sl = tuple(slice(1, phi.shape[m] - 1) if m != axis else slice(None)
                   for m in range(3))
        pl, pr = pl[sl], pr[sl]
        tl = np.take(th, range(0, th.shape[axis] - 1), axis=axis)
        tr = np.take(th, range(1, th.shape[axis]), axis=axis)
        cut = (pl < 0.0) != (pr < 0.0)
        if not cut.any():
            continue
        hd = float(dl[d][0])            # uniform in the cut direction
        invd = 1.0 / hd
        w = np.where(np.abs(pl - pr) * invd < MIN_COSINE, 0.5, pl / (pl - pr))
        kl = np.where(pl < 0.0, a.kappa, 1.0)
        kr = np.where(pr < 0.0, a.kappa, 1.0)
        kface = dm / (w / kl + (1.0 - w) / kr + a.rc * dm * invd)
        # the flux in the +axis sense, positive INTO the fluid: -k dT/dn,
        # signed by which side is solid (the solver's own convention).
        flux = -kface * (tr - tl) * invd
        sgn = np.where(pl < 0.0, 1.0, -1.0)
        area = np.ones_like(flux)
        for t in range(3):
            if t == d:
                continue
            shape = [1, 1, 1]
            shape[2 - t] = len(dl[t])
            area = area * dl[t].reshape(shape)
        total += float((sgn * flux * area)[cut].sum())
        nfaces += int(cut.sum())

    print(f"   python: {nfaces} cut faces   interface heat = {total:.16e}")
    ok = True
    if a.heat:
        rows = read_heat(a.heat)
        _, solver = rows[-1]
        rel = abs(solver - total) / max(abs(total), 1e-300)
        print(f"   solver: {solver:.16e}   relative difference {rel:.3e}")
        ok = rel <= a.tolerance
        print("   PASS" if ok else f"   FAIL (tolerance {a.tolerance:g})")
    return 0 if ok else 1


def read_heat(path, column=2):
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            v = [float(x) for x in line.split()]
            if len(v) >= 5:
                rows.append((v[1], v[2 + column]))
    if not rows:
        raise SystemExit(f"{path}: no data rows")
    return rows


def cmd_slab(a):
    rows = read_heat(a.heat)
    t, solver = rows[-1]
    exact = a.cap * a.source * a.y_wall * a.area
    rel = abs(solver - exact) / abs(exact)
    print(f"   t = {t:.6f}   solver {solver:.16e}   closed form {exact:.16e}")
    print(f"   relative difference {rel:.3e}")
    if len(rows) >= 2:
        print(f"   (steadiness: the previous sample was {rows[-2][1]:.16e}, "
              f"relative change {abs(rows[-1][1]/rows[-2][1] - 1.0):.3e})")
    ok = rel <= a.tolerance
    print("   PASS" if ok else f"   FAIL (tolerance {a.tolerance:g})")
    return 0 if ok else 1


def case_marker(case, nblocks, shape):
    """The SOLVER's own material marker, straight out of the case file: a cell
    centre is solid iff |coef_p| > SOLID/Re, which is the test
    init_scalar_conjugate signs phi with. Reading it beats reproducing the
    geometry analytically -- the les_ibm wall planes, for instance, are only
    the requested planes to float32 (the T1 lesson), and a single misclassified
    cell would silently break a budget that is supposed to telescope.
    Tiles are ghost-inclusive and stored (i, j, k); the interior in the field
    files' (k, j, i) order is the transpose."""
    with h5py.File(case, "r") as h5:
        coef = h5["coef_p_blocks"][...]
        nb = int(h5.attrs["block_nb"])
        if coef.shape[0] != nblocks:
            raise SystemExit("case file and snapshot have different block tables")
    out = np.abs(coef[:, 1:nb + 1, 1:nb + 1, 1:nb + 1].transpose(0, 3, 2, 1)) > 1.0e20
    if out.shape[1:] != shape:
        raise SystemExit(f"marker shape {out.shape[1:]} vs field {shape}")
    return out


def fluid_energy(path, name, cap, wavy, plane, case=None):
    """sum C theta dV over the cells the SOLVER classifies fluid, up to
    `plane` -- with the solver's own capacity, read from its `vfrac`."""
    with h5py.File(path, "r") as h5:
        geo = BlockGeometry(h5)
        th = h5[name][...]
        vf = h5["vfrac"][...]
        t = float(h5.attrs["t_current"])
        marker = case_marker(case, geo.n_blocks, th.shape[1:]) if case else None
        e, bad = 0.0, 0
        for bid in range(geo.n_blocks):
            x, y, z, dV = geo.mesh(bid)
            xx = np.broadcast_to(x, th[bid].shape)
            yy = np.broadcast_to(y, th[bid].shape)
            if marker is not None:
                solid = marker[bid]
            else:
                solid = yy < wavy_wall_height(xx, lx=geo.leng[0]) if wavy else yy < 0.0
            keep = (~solid) & (yy < plane)
            f = vf[bid]
            cc = f + (1.0 - f) * cap
            e += float((cc * th[bid] * np.broadcast_to(dV, th[bid].shape))[keep].sum())
            # a contradiction between the analytic marker and the solver's
            # fraction would silently break the classification this rests on
            bad += int(np.sum(solid & (f == 1.0)) + np.sum((~solid) & (f == 0.0)))
    return e, t, bad


def plane_flux(path, name, plane, re, pr):
    """The convective + diffusive theta flux UPWARD through one global y face,
    from the snapshot alone -- the A2 border-flux idiom, on a plane that lies
    entirely in the fluid.

    The face is the global y node nearest `plane`, moved off a block boundary
    if it lands on one: there the two cells belong to different blocks and the
    interior-only datasets carry no halo to join them.
    """
    with h5py.File(path, "r") as h5:
        geo = BlockGeometry(h5)
        th = h5[name][...]
        vn = h5["vn"][...]
        ynode = geo.lines[1][0]
        nb = geo.nb[1]
        jf = int(np.argmin(np.abs(ynode - plane)))
        if jf % nb == 0:
            jf += 1
        dm = 1.0/(re*pr)
        flux = 0.0
        for bid in range(geo.n_blocks):
            oy = int(geo.blocks[bid, 1])
            j = jf - oy
            if not 1 <= j <= nb - 1:
                continue
            x, y, z, dV = geo.mesh(bid)
            (_, _), (yc, dyc), (_, _) = geo.block_axes(bid)
            da = np.broadcast_to(dV, th[bid].shape)[:, j, :]/dyc[j]   # dx dz
            tlo, thi = th[bid][:, j - 1, :], th[bid][:, j, :]
            dy = float(yc[j] - yc[j - 1])
            flux += float((vn[bid][:, j, :]*0.5*(tlo + thi)*da).sum())
            flux -= float((dm*(thi - tlo)/dy*da).sum())
    return flux, float(ynode[jf])


def cmd_budget(a):
    """dE/dt over the fluid cell set == interface heat - plane flux."""
    # The CV's upper boundary must be the SAME face the flux is evaluated on,
    # or the balance is short by a whole cell row -- and the plane is snapped
    # off block boundaries, so it is the FACE that decides, not the request.
    yface = 1e30
    fl = 0.0
    if a.plane is not None:
        fl, yface = plane_flux(a.second, a.name, a.plane, a.re, a.pr)
    e0, t0, b0 = fluid_energy(a.first, a.name, a.cap, a.wavy, yface, a.case)
    e1, t1, b1 = fluid_energy(a.second, a.name, a.cap, a.wavy, yface, a.case)
    e2, t2, b2 = fluid_energy(a.third, a.name, a.cap, a.wavy, yface, a.case)
    dedt = (e2 - e0)/(t2 - t0)
    rows = dict(read_heat(a.heat))
    key = min(rows, key=lambda t: abs(t - t1))
    heat = rows[key]
    rhs = heat - fl
    scale = max(abs(dedt), abs(heat), 1e-300)
    rel = abs(dedt - rhs)/scale
    print(f"   fluid cells: E({t0:.6f}) = {e0:.12e}   E({t2:.6f}) = {e2:.12e}")
    print(f"   dE/dt (centred) = {dedt:.12e}")
    if a.plane is not None:
        print(f"   interface heat = {heat:.12e}   plane flux (y = {yface:g}) "
              f"= {fl:.12e}")
    else:
        print(f"   interface heat = {heat:.12e}   (no plane: the whole fluid)")
    print(f"   residual = {dedt - rhs:.3e}   relative = {rel:.3e}"
          f"   (marker contradictions {b0 + b1 + b2})")
    if a.emit:
        with open(a.emit, "a") as fh:
            fh.write("%s %.10e %.10e %.10e\n" % (a.tag, t2 - t0, dedt - rhs, rel))
    ok = rel <= a.tolerance and (b0 + b1 + b2) == 0
    print("   PASS" if ok else f"   FAIL (tolerance {a.tolerance:g})")
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("sum")
    p.add_argument("field")
    p.add_argument("--case", required=True)
    p.add_argument("--kappa", type=float, required=True)
    p.add_argument("--rc", type=float, default=0.0)
    p.add_argument("--re", type=float, default=1.0)
    p.add_argument("--pr", type=float, default=1.0)
    p.add_argument("--heat", default=None)
    p.add_argument("--name", default="theta")
    p.add_argument("--tolerance", type=float, default=1e-13)
    p.set_defaults(func=cmd_sum)

    p = sub.add_parser("slab")
    p.add_argument("--heat", required=True)
    p.add_argument("--y-wall", type=float, required=True)
    p.add_argument("--kappa", type=float, default=1.0)
    p.add_argument("--cap", type=float, required=True)
    p.add_argument("--source", type=float, required=True)
    p.add_argument("--area", type=float, required=True)
    p.add_argument("--tolerance", type=float, default=1e-10)
    p.set_defaults(func=cmd_slab)

    p = sub.add_parser("budget")
    p.add_argument("first")
    p.add_argument("second")
    p.add_argument("third")
    p.add_argument("--heat", required=True)
    p.add_argument("--cap", type=float, required=True)
    p.add_argument("--wavy", action="store_true")
    p.add_argument("--case", default=None,
                   help="classify cells with the case file's own IBM marker")
    p.add_argument("--plane", type=float, default=None)
    p.add_argument("--re", type=float, default=100.0)
    p.add_argument("--pr", type=float, default=0.71)
    p.add_argument("--name", default="theta")
    p.add_argument("--emit", default=None)
    p.add_argument("--tag", default="run")
    p.add_argument("--tolerance", type=float, default=1e-6)
    p.set_defaults(func=cmd_budget)

    a = ap.parse_args()
    return a.func(a)


if __name__ == "__main__":
    sys.exit(main())
