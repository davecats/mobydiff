#!/usr/bin/env python3
"""Gates for the THERMAL WALL FUNCTION AT AN IMMERSED WALL (S5a's open item).

S5a validated the Kader/Jayatilleke closure on DOMAIN walls only. Here the
wall is an immersed body, so the classified wall cell is a CUT cell: the
Dirichlet penalization pins it to `ibm_value` while the wall function
installs an eddy diffusivity on the same cell. `ibmwf180.ini` puts that cell
in the LOG branch of both wall functions (T3's only IBM wall-function case
sits at y+ ~ 2-3, the conduction branch, where nothing new happens).

Driving: isothermal walls at `ibm_value` and a constant volumetric `source`,
because ONE body value serves both walls. That makes the steady budget
closed-form, which is what the `budget` subcommand below exploits.

Subcommands
  wall    <ransgeom.h5> <snapshot.h5>
          which cells are wall cells, their k-based y+, whether the LOG
          branch fires (y+ > y+_T), and the momentum wall function's nu_t
          against an INDEPENDENT transcription of nu(y+ kappa/ln(E y+) - 1).
  budget  <snapshot.h5> <case.h5>
          the exact steady balance. Summed over all interior cells the
          convective and diffusive fluxes telescope to the domain boundary
          (periodic in x/z, zero-flux walls in y), so at steady state

              sum_cells coef_p (s_body - s) dV / Pr  =  -source * V_total

          EXACTLY -- and, restricted to the cells the penalization does not
          pin, the heat crossing into the body is source * V_fluid. Both are
          checked against the snapshot, and against the solver's own runtime
          heat file when it is given (--heat), which is what turned up the
          double count documented in README.md.
"""

from __future__ import annotations

import argparse
import sys

import h5py
import numpy as np

KAPPA = 0.41
ELOG = 9.8
CMU25 = 0.09 ** 0.25
SOLID_FACE_THRESHOLD = 1.0e20


def jayatilleke_p(prat: float) -> float:
    return 9.24 * (prat ** 0.75 - 1.0) * (1.0 + 0.28 * np.exp(-0.007 * prat))


def thermal_yplus(pr: float, prt: float, p: float) -> float:
    """Where Pr y+ meets Pr_t[ln(E y+)/kappa + P], beyond the minimum."""
    f = lambda y: pr * y - prt * (np.log(ELOG * y) / KAPPA + p)
    lo = prt / (KAPPA * pr)
    if f(lo) >= 0.0:
        return lo
    hi = lo
    while f(hi) < 0.0:
        hi *= 2.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) < 0.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def reassemble(h5: h5py.File, name: str) -> np.ndarray:
    """Block-table dataset -> the global (z,y,x) array (single level). The
    shape comes from the block table, not from nx/ny/nz: the ransgeom dump
    is self-contained and carries no grid attributes."""
    blocks = h5["blocks"][...]
    nb = (int(h5.attrs["block_nb_x"]), int(h5.attrs["block_nb_y"]),
          int(h5.attrs["block_nb_z"]))
    shape = (int(blocks[:, 2].max()) + nb[2],
             int(blocks[:, 1].max()) + nb[1],
             int(blocks[:, 0].max()) + nb[0])
    out = np.zeros(shape)
    data = h5[name]
    for bid, (ox, oy, oz, lev) in enumerate(blocks):
        out[oz:oz + nb[2], oy:oy + nb[1], ox:ox + nb[0]] = data[bid]
    return out


def case_coef_p(case: str) -> np.ndarray:
    """coef_p_blocks carries a one-cell ghost layer per leaf; strip it."""
    with h5py.File(case, "r") as f:
        blocks = f["blocks"][...]
        tiles = f["coef_p_blocks"]
        nbt = tiles.shape[1] - 2                    # ghost-inclusive extent
        nz = int(blocks[:, 2].max()) + nbt
        ny = int(blocks[:, 1].max()) + nbt
        nx = int(blocks[:, 0].max()) + nbt
        out = np.zeros((nz, ny, nx))
        for bid, (ox, oy, oz, lev) in enumerate(blocks):
            out[oz:oz + nbt, oy:oy + nbt, ox:ox + nbt] = tiles[bid][1:-1, 1:-1, 1:-1]
    return out


def cmd_wall(a) -> int:
    g = h5py.File(a.ransgeom, "r")
    f = h5py.File(a.snapshot, "r")
    wallcell = reassemble(g, "wallcell")
    yeff = reassemble(g, "yeff")
    dwall = reassemble(g, "dwall")
    k = reassemble(f, "k")
    nut = reassemble(f, "nut")
    nu = 1.0 / a.re

    p = jayatilleke_p(a.pr / a.prt)
    ypt = thermal_yplus(a.pr, a.prt, p)
    yplus = CMU25 * np.sqrt(np.maximum(k, 0.0)) * yeff / nu

    mask = wallcell == 1
    n = int(mask.sum())
    print(f"wall cells: {n} of {wallcell.size}   "
          f"(solid cells, wallcell == 2: {int((wallcell == 2).sum())})")
    if n == 0:
        print("FAIL: no IBM wall cell — the case is not exercising the gate")
        return 1
    yw = yplus[mask]
    print(f"Jayatilleke P({a.pr}/{a.prt}) = {p:.6f}   y+_T = {ypt:.4f}   "
          f"(momentum switch y+_lam = 11.5301)")
    print(f"wall-cell y+_k : min {yw.min():.4f}  max {yw.max():.4f}  "
          f"mean {yw.mean():.4f}")
    print(f"wall-cell dwall: min {dwall[mask].min():.6f}  "
          f"max {dwall[mask].max():.6f}")
    log_th = yw > ypt
    log_mom = yw > 11.5301
    print(f"THERMAL log branch fires on {int(log_th.sum())}/{n} wall cells; "
          f"MOMENTUM log branch on {int(log_mom.sum())}/{n}")

    # The momentum wall function's nu_t, transcribed independently.
    ref = np.where(log_mom, nu * (yw * KAPPA / np.log(ELOG * np.maximum(yw, 1e-30)) - 1.0), 0.0)
    got = nut[mask]
    rel = np.abs(got - ref) / np.maximum(np.abs(ref), 1e-300)
    print(f"wall-cell nu_t vs nu(y+ kappa/ln(E y+) - 1): max rel dev "
          f"{rel.max():.3e}   (nu_t {got.min():.6e} .. {got.max():.6e})")

    ok = bool(log_th.all()) and rel.max() <= a.tol
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


def cmd_budget(a) -> int:
    f = h5py.File(a.snapshot, "r")
    s = reassemble(f, a.scalar)
    cp = case_coef_p(a.case)
    if cp.shape != s.shape:
        raise SystemExit(f"coef_p {cp.shape} vs field {s.shape}")

    nx, ny, nz = s.shape[2], s.shape[1], s.shape[0]
    dv = (a.lx / nx) * (a.ly / ny) * (a.lz / nz)
    vol = a.lx * a.ly * a.lz

    pen = cp * (a.value - s) * dv / a.pr
    solid = np.abs(cp) > SOLID_FACE_THRESHOLD
    total_pen = pen.sum()
    unpinned_pen = pen[~solid].sum()
    v_solid = float(solid.sum()) * dv

    src_all = a.source * vol
    src_fluid = a.source * (vol - v_solid)

    print(f"cells {s.size}  ({int(solid.sum())} above the solid threshold)  "
          f"dV = {dv:.10e}")
    d1 = abs(total_pen + src_all) / src_all
    print(f"[1] IDENTITY over ALL interior cells:"
          f" sum coef_p (s_body - s) dV/Pr")
    print(f"    snapshot {total_pen:.10e}    -source*V_total"
          f" {-src_all:.10e}    rel dev {d1:.3e}")
    print(f"    (at steady state the whole source input leaves through the"
          f" penalization,")
    print(f"     INCLUDING what is deposited inside the body and never"
          f" crossed the interface)")
    print(f"[2] the same restricted to the UNPINNED cells ="
          f" {unpinned_pen:.10e}")
    print(f"    (a solid cell's own share rides the staircase flux instead,"
          f" so this is")
    print(f"     only part of the fluid -> body heat; the sum of the two"
          f" is gate [3])")

    ok = d1 <= a.tol
    if a.heat:
        row = np.loadtxt(a.heat, ndmin=2)[-1]
        col = 2 + 3 * (a.index - 1)
        stair, graded, tot = row[col], row[col + 1], row[col + 2]
        dg = abs(graded - unpinned_pen) / max(abs(unpinned_pen), 1e-300)
        d3 = abs(abs(tot) - src_fluid) / src_fluid
        print(f"[3] PHYSICAL GATE -- the solver's runtime heat file"
              f" (last sample, t = {row[1]:.3f}):")
        print(f"    staircase {stair:.10e}   graded {graded:.10e}"
              f"   total {tot:.10e}")
        print(f"    graded vs the snapshot's UNPINNED sum:"
              f" rel dev {dg:.3e}   ({'MATCH' if dg <= a.tol else 'DIFFERS'})")
        print(f"    TOTAL vs source*V_fluid = {src_fluid:.10e}"
              f" (the heat crossing into the body, ONCE):")
        print(f"      rel dev {d3:.3e}"
              f"   ({'MATCH' if d3 <= a.tol else 'DIFFERS'})")
        # Pre-fix, the penalization column was the ALL-cells sum [1], so the
        # reported total was that plus the staircase term -- reconstruct it
        # from the snapshot, not from the (already fixed) total.
        prefix_total = abs(stair) + abs(total_pen)
        print(f"    [before the 2026-08-05 scalar_stats fix this case read"
              f" {prefix_total:.6e} at")
        print(f"     ibm_value = 0 ({100*(prefix_total/src_fluid - 1):.0f} %"
              f" high): the solid cells' heat was counted BOTH")
        print(f"     as penalization and as the staircase flux into them"
              f" -- see README.md]")
        ok = ok and dg <= a.tol and d3 <= a.tol
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


def main() -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("wall")
    p.add_argument("ransgeom")
    p.add_argument("snapshot")
    p.add_argument("--re", type=float, default=180.0)
    p.add_argument("--pr", type=float, default=0.71)
    p.add_argument("--prt", type=float, default=0.85)
    p.add_argument("--tol", type=float, default=1e-12)
    p.set_defaults(func=cmd_wall)

    p = sub.add_parser("budget")
    p.add_argument("snapshot")
    p.add_argument("case")
    p.add_argument("--scalar", default="theta")
    p.add_argument("--index", type=int, default=1)
    p.add_argument("--value", type=float, default=0.0)
    p.add_argument("--source", type=float, default=0.05)
    p.add_argument("--pr", type=float, default=0.71)
    p.add_argument("--lx", type=float, default=np.pi)
    p.add_argument("--ly", type=float, default=2.5)
    p.add_argument("--lz", type=float, default=0.5 * np.pi)
    p.add_argument("--heat", default=None)
    p.add_argument("--tol", type=float, default=1e-8)
    p.set_defaults(func=cmd_budget)

    a = ap.parse_args()
    return a.func(a)


if __name__ == "__main__":
    sys.exit(main())
