#!/usr/bin/env python3
"""Checkers for the conjugate-interface gates (increment C1).

    ./check_conjugate.py slab     <field.h5> --y-wall Y --kappa K [--prev P]
    ./check_conjugate.py weight   <case.h5>  --y-wall Y
    ./check_conjugate.py conserve <a.h5> <b.h5> --capacity C
    ./check_conjugate.py limit    <conjugate.h5> <reference.h5> --y-wall Y

`slab`     the exact two-material steady solution (gate 1).
`weight`   the level-set weight w, straight out of the case file's own
           dwall_blocks / coef_p_blocks, against the analytic cut position --
           the one genuinely new ingredient of the increment, checked without
           going through the solver at all (gate 1).
`conserve` sum(C theta dV) between two snapshots (gate 3).
`limit`    a conjugate run against its kappa_s -> infinity / 0 twin (gate 2).

Geometry comes from scalar_tools.BlockGeometry, so cell centres and volumes
are the solver's to the last bit.
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import h5py

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scalar"))
from scalar_tools import BlockGeometry                       # noqa: E402

SOLID_THRESHOLD = 1.0e20


def wavy_wall_height(x, amp=0.025, nwave=1, phase=0.0, lx=1.0):
    """ibm.f90 wavy_wall_height, verbatim (the analytic-path geometry)."""
    return amp * 0.5 * (1.0 + np.sin(2.0 * np.pi * nwave * x / lx + phase)) + 1.0e-2


def slab_exact(y, y_wall, kappa, ly, rc=0.0, depth=0.0, kappa_outer=None,
               outer_init=0.0):
    """The piecewise-linear multi-material steady solution, T(0)=0, T(ly)=1.

    A contact resistance rc simply adds to the series and puts a JUMP q*rc in
    T at the interface; the profile stays piecewise linear, so it stays an
    exact fixed point of the discrete operator.

    With depth > 0 the solid carries the F2 second material band: a SHELL of
    thickness `depth` against the wall plane with conductivity `kappa`, and
    everything deeper with `kappa_outer`. Three layers in series, still
    piecewise linear, so still an exact fixed point -- and now the fixed
    point depends on where the band boundary is, which is exactly what the
    level-set weight on the shifted level set psi = phi + depth decides.

    kappa_outer = 0 is the INSULATOR the pipe campaign actually uses. The
    series resistance is then infinite: no heat crosses the band boundary,
    the shell and the fluid sit at the y = ly Dirichlet value, and the outer
    band is inert at whatever it started from (`outer_init`).
    """
    if depth <= 0.0:
        q = 1.0 / (y_wall / kappa + rc + (ly - y_wall))
        return np.where(y <= y_wall, q * y / kappa,
                        q * (y_wall / kappa + rc + (y - y_wall)))
    y_band = y_wall - depth
    if kappa_outer == 0.0:
        return np.where(y <= y_band, outer_init, 1.0)
    r_out = y_band / kappa_outer
    r_shell = depth / kappa
    q = 1.0 / (r_out + r_shell + rc + (ly - y_wall))
    return np.where(
        y <= y_band, q * y / kappa_outer,
        np.where(y <= y_wall,
                 q * (r_out + (y - y_band) / kappa),
                 q * (r_out + r_shell + rc + (y - y_wall))))


def read_scalar(path, name="theta"):
    h5 = h5py.File(path, "r")
    return h5, BlockGeometry(h5), h5[name][...]


def cmd_slab(a):
    h5, geo, th = read_scalar(a.field, a.name)
    ly = geo.leng[1]
    worst = -1.0
    worst_at = (float("nan"),) * 3
    for bid in range(geo.n_blocks):
        x, y, z, _ = geo.mesh(bid)
        yy = np.broadcast_to(y, th[bid].shape)
        ref = slab_exact(yy, a.y_wall, a.kappa, ly, a.contact,
                         a.depth, a.outer_kappa, a.outer_init)
        err = np.abs(th[bid] - ref)
        if err.max() > worst:
            worst = float(err.max())
            k, j, i = np.unravel_index(int(err.argmax()), err.shape)
            worst_at = (float(yy[k, j, i]), float(th[bid][k, j, i]), float(ref[k, j, i]))
    h5.close()

    drift = float("nan")
    if a.prev:
        p5 = h5py.File(a.prev, "r")
        drift = float(np.abs(p5[a.name][...] - th).max())
        p5.close()

    if a.depth > 0.0:
        print(f"   kappa_s = {a.kappa:g}  y_wall = {a.y_wall!r}"
              f"  band depth = {a.depth!r} (boundary at y = {a.y_wall - a.depth!r})"
              f"  kappa_outer = {a.outer_kappa:g}")
    else:
        q = 1.0 / (a.y_wall / a.kappa + a.contact + (ly - a.y_wall))
        print(f"   kappa_s = {a.kappa:g}  y_wall = {a.y_wall!r}  R_c = {a.contact:g}"
              f"  q_exact = {q:.16e}")
    print(f"   max|theta - exact| = {worst:.6e}   at y = {worst_at[0]:.6f} "
          f"({worst_at[1]:.16e} vs {worst_at[2]:.16e})")
    if a.prev:
        print(f"   transient residual (last write interval) = {drift:.3e}")
    # A NaN field must FAIL, not slip through a `<= tolerance` comparison.
    ok = np.isfinite(worst) and 0.0 <= worst <= a.tolerance
    print("   PASS" if ok else f"   FAIL (tolerance {a.tolerance:g})")
    return 0 if ok else 1


def cmd_weight(a):
    """w = phi_L/(phi_L - phi_R) from the CASE FILE, vs the analytic cut.

    Independent of the solver: it rebuilds phi from the two datasets the
    scheme reads -- dwall_blocks for the magnitude and coef_p_blocks for the
    sign -- and compares the level-set weight of every cut y-arm with the
    analytic cut fraction of the STL plane. Case-file tiles are
    ghost-inclusive, shape (nBlocks, nb+2, nb+2, nb+2) in (i, j, k) order --
    CASE-FILE tiles, unlike the FIELD datasets, which are (k, j, i); see
    check_oblique.py load_phi. This checker cannot tell the two apart (its
    slab is x/z symmetric) and only ever indexes axis 1, which is j either
    way, so the distinction does not reach the number below.
    """
    with h5py.File(a.case, "r") as h5:
        blocks = h5["blocks"][...]
        nb = int(h5.attrs["block_nb"])
        ynode = h5["y_nodes"][...]
        dwall = h5["dwall_blocks"][...]
        coef = h5["coef_p_blocks"][...]
    if int(blocks[:, 3].max()) != 0:
        raise SystemExit("weight: single-level case files only")

    # One ghost node at each end, so ghost-cell centres are addressable: the
    # cut arm may well sit ON a block boundary, which is exactly the case the
    # ghost-inclusive tiles exist for (strategy doc Section 8, arrangement 2).
    ext = np.concatenate(([2 * ynode[0] - ynode[1]], ynode,
                          [2 * ynode[-1] - ynode[-2]]))

    worst, n = 0.0, 0
    for bid in range(blocks.shape[0]):
        oy = int(blocks[bid, 1])
        yc = 0.5 * (ext[oy:oy + nb + 2] + ext[oy + 1:oy + nb + 3])   # j = 0..nb+1
        solid = np.abs(coef[bid]) > SOLID_THRESHOLD
        phi = np.where(solid, -np.maximum(dwall[bid], 1e-300), dwall[bid])
        for j in range(nb + 1):           # every y arm, ghost arms included
            cut = solid[:, j, :] != solid[:, j + 1, :]
            if not cut.any():
                continue
            pl, pr = phi[:, j, :][cut], phi[:, j + 1, :][cut]
            w = pl / (pl - pr)
            # Analytic: the fraction of the arm on the LOW cell's side.
            w_ref = (a.y_wall - yc[j]) / (yc[j + 1] - yc[j])
            worst = max(worst, float(np.abs(w - w_ref).max()))
            n += int(cut.sum())

    print(f"   cut y-arms = {n}   max|w - w_exact| = {worst:.3e}")
    ok = n > 0 and worst <= a.tolerance
    print("   PASS" if ok else f"   FAIL (tolerance {a.tolerance:g}, arms {n})")
    return 0 if ok else 1


def capacity_integral(path, name, capacity, geometry):
    """sum(C theta dV) with the SOLVER's own cell capacity.

    C3: the capacity of a cut cell is FLUID-FRACTION WEIGHTED,
    C = f + (1-f) C_s, so the invariant the scheme conserves is built on f
    -- and f is a discrete quantity (a plane reconstructed from the central
    differences of phi), not something a checker can rebuild from the
    analytic geometry to the last bit. An analytic capacity map instead of
    the solver's would leave sum (C_ref - C) dtheta dV, which is ~1e-3
    relative here, i.e. the gate would measure the map mismatch rather than
    conservation. So the solver writes f to the snapshot as `vfrac` and the
    checker reads it. The analytic marker is still used, as a CROSS-CHECK of
    the fraction (reported, not gated: they must agree away from the
    interface exactly and inside it to O(h)).
    """
    h5, geo, th = read_scalar(path, name)
    vf = h5["vfrac"][...] if "vfrac" in h5 else None
    if vf is None:
        raise SystemExit(f"{path}: no vfrac dataset -- conjugate runs write it")
    total = 0.0
    worst = 0.0
    for bid in range(geo.n_blocks):
        x, y, z, dV = geo.mesh(bid)
        if geometry == "wavy":
            solid = np.broadcast_to(y, th[bid].shape) < wavy_wall_height(
                np.broadcast_to(x, th[bid].shape), lx=geo.leng[0])
        else:
            solid = np.broadcast_to(y, th[bid].shape) < geometry
        f = vf[bid]
        cc = f + (1.0 - f) * capacity
        total += float((cc * th[bid] * np.broadcast_to(dV, th[bid].shape)).sum())
        worst = max(worst, float(np.abs(f - np.where(solid, 0.0, 1.0)).max()))
    h5.close()
    return total, worst


def cmd_conserve(a):
    geometry = "wavy" if a.wavy else a.y_wall
    ia, da = capacity_integral(a.first, a.name, a.capacity, geometry)
    ib, db = capacity_integral(a.second, a.name, a.capacity, geometry)
    scale = max(abs(ia), abs(ib), 1e-300)
    rel = abs(ib - ia) / scale
    print(f"   sum(C theta dV): {ia:.16e} -> {ib:.16e}")
    print(f"   drift = {ib - ia:.3e}   relative = {rel:.3e}")
    print(f"   (fluid fraction vs the analytic pointwise marker: max dev "
          f"{max(da, db):.3f} -- a cut cell, by construction)")
    ok = rel <= a.tolerance
    print("   PASS" if ok else f"   FAIL (tolerance {a.tolerance:g})")
    return 0 if ok else 1


def cmd_limit(a):
    """Fluid-side comparison of a conjugate run and its limiting twin."""
    ha, ga, ta = read_scalar(a.conjugate, a.name)
    hb, gb, tb = read_scalar(a.reference, a.name)
    worst = 0.0
    for bid in range(ga.n_blocks):
        x, y, z, _ = ga.mesh(bid)
        yy = np.broadcast_to(y, ta[bid].shape)
        fluid = yy > a.y_wall
        if not fluid.any():
            continue
        worst = max(worst, float(np.abs(ta[bid][fluid] - tb[bid][fluid]).max()))
    ha.close()
    hb.close()
    print(f"   max|theta_conjugate - theta_reference| over the FLUID = {worst:.6e}")
    ok = worst <= a.tolerance
    print("   PASS" if ok else f"   FAIL (tolerance {a.tolerance:g})")
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("slab")
    p.add_argument("field")
    p.add_argument("--y-wall", type=float, required=True)
    p.add_argument("--kappa", type=float, required=True)
    p.add_argument("--prev", default=None)
    p.add_argument("--contact", type=float, default=0.0)
    # The F2 second material band: shell depth, the material beyond it, and
    # the value an INSULATED outer band is frozen at (its solid_init).
    p.add_argument("--depth", type=float, default=0.0)
    p.add_argument("--outer-kappa", type=float, default=None)
    p.add_argument("--outer-init", type=float, default=0.0)
    p.add_argument("--name", default="theta")
    p.add_argument("--tolerance", type=float, default=1e-13)
    p.set_defaults(func=cmd_slab)

    p = sub.add_parser("weight")
    p.add_argument("case")
    p.add_argument("--y-wall", type=float, required=True)
    p.add_argument("--tolerance", type=float, default=1e-12)
    p.set_defaults(func=cmd_weight)

    p = sub.add_parser("conserve")
    p.add_argument("first")
    p.add_argument("second")
    p.add_argument("--capacity", type=float, required=True)
    p.add_argument("--wavy", action="store_true")
    p.add_argument("--y-wall", type=float, default=0.0)
    p.add_argument("--name", default="theta")
    p.add_argument("--tolerance", type=float, default=1e-13)
    p.set_defaults(func=cmd_conserve)

    p = sub.add_parser("limit")
    p.add_argument("conjugate")
    p.add_argument("reference")
    p.add_argument("--y-wall", type=float, required=True)
    p.add_argument("--name", default="theta")
    p.add_argument("--tolerance", type=float, default=1e-3)
    p.set_defaults(func=cmd_limit)

    a = ap.parse_args()
    return a.func(a)


if __name__ == "__main__":
    sys.exit(main())
