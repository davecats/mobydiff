#!/usr/bin/env python3
"""Checker for the passive-scalar STATISTICS gates (increment S4).

The solver accumulates the scalar statistics in scalar_stats.f90 while it
runs; this recomputes THE SAME seven columns from the written snapshots and
compares them row by row. Because the solver samples the end-of-step field
and the snapshot IS the end-of-step field, agreement must be to round-off --
not to statistical tolerance. That makes it a sharp gate on the sampling
kernel, on the row/level bookkeeping and (most of all) on the face fluxes,
which are the transport kernel's own expressions.

Subcommands
  profile  <stats.h5> <field.h5> [field.h5 ...]
           the wall-normal (channel) layout: x-z averaged rows. Several
           snapshots reproduce an ACCUMULATED file (the solver's running sum
           over the same samples).
  plane    <stats.h5> <field.h5> [...]
           the (x,y) (boundary-layer) layout: z averaged rows.
  rows     <stats.h5> <field.h5> [...] --level L
           the per-LEVEL row bookkeeping of a 2:1-refined case, accumulated
           straight from the leaves (no global box), so it gates the
           multi-level statistics files `profile` cannot reach.
  diff     <a.h5> <b.h5>
           two statistics files against each other (rank counts, devices,
           restart continuation). The reduction order differs between them --
           atomics in the sampling kernel, the allreduce tree across ranks --
           so this is a tight tolerance, not an equality.

Assumptions of the recomputation (true for the gate cases, checked here):
single refinement level, and the y walls carry the boundary condition given
by --walls / --wall-types (used to build the ghost row exactly as
apply_scalar_bc_q does: mirror for Dirichlet, copy for Neumann 0). nut ghosts
are zero on a physical face, which is what the solver's halo convention
leaves there, so the wall face diffusivity is D + 0.5 nut_1/Pr_t -- the
kernel's value, NOT the 'nut = 0 at the wall' idealisation.
"""

from __future__ import annotations

import argparse
import sys

import h5py
import numpy as np

from scalar_tools import BlockGeometry
from check_scalar_turb import prt_kays

NSTAT = 7
S, SS, US, CLO, JLO, CHI, JHI = range(7)
NAMES = ["<s>", "<s^2>", "<u_c s>", "<v s>|lo", "J|lo", "<v s>|hi", "J|hi"]


def assemble(f: h5py.File, name: str, geo: BlockGeometry) -> np.ndarray:
    """Single-level block-table dataset -> one global (nz, ny, nx) array."""
    nbx, nby, nbz = geo.nb
    out = np.zeros((int(f.attrs["nz"]), int(f.attrs["ny"]), int(f.attrs["nx"])))
    data = f[name]
    for bid in range(geo.n_blocks):
        ox, oy, oz, lev = (int(t) for t in geo.blocks[bid])
        if lev != 0:
            raise SystemExit("the statistics checker handles single-level files only")
        out[oz:oz + nbz, oy:oy + nby, ox:ox + nbx] = data[bid]
    return out


def widths(line: np.ndarray):
    return 0.5 * (line[:-1] + line[1:]), np.diff(line)


def ghost_row(s: np.ndarray, kind: str, value: float, at_low: bool) -> np.ndarray:
    """The physical ghost apply_scalar_bc_q writes next to the wall."""
    row = s[:, 0, :] if at_low else s[:, -1, :]
    if kind == "dirichlet":
        return 2.0 * value - row
    return row                      # Neumann 0


def face_diffusivity(nutf, pr, prt, model, re):
    d = 1.0 / (re * pr)
    nt = np.maximum(nutf, 0.0)
    if model == "kays":
        return d + nt / prt_kays(nt * re * pr, prt)
    return d + nt / prt


def snapshot_columns(path, a):
    """The seven columns' RAW SUMS and the row weights of one snapshot, in the
    layout `a.layout` -- i.e. exactly what scalar_stats.f90 accumulates."""
    with h5py.File(path, "r") as f:
        geo = BlockGeometry(f)
        s = assemble(f, a.scalar, geo)
        u = assemble(f, "un", geo)
        v = assemble(f, "vn", geo)
        nut = assemble(f, "nut", geo) if "nut" in f else np.zeros_like(s)
        re = float(f.attrs["re"])
        xc, dx = widths(geo.lines[0][0])
        yc, dy = widths(geo.lines[1][0])
        zc, dz = widths(geo.lines[2][0])
    nz, ny, nx = s.shape

    # Ghost rows and the top face values the staggered arrays do not store.
    gs_lo = ghost_row(s, a.wall_types[0], a.walls[0], True)
    gs_hi = ghost_row(s, a.wall_types[1], a.walls[1], False)
    sS = np.concatenate([gs_lo[:, None, :], s[:, :-1, :]], axis=1)
    sN = np.concatenate([s[:, 1:, :], gs_hi[:, None, :]], axis=1)
    vlo = v
    vhi = np.concatenate([v[:, 1:, :], np.zeros((nz, 1, nx))], axis=1)   # wall: v = 0
    # nut halos: the neighbour cell inside the domain, 0 in the physical ghost.
    nS = np.concatenate([np.zeros((nz, 1, nx)), nut[:, :-1, :]], axis=1)
    nN = np.concatenate([nut[:, 1:, :], np.zeros((nz, 1, nx))], axis=1)
    dlo = face_diffusivity(0.5 * (nS + nut), a.pr, a.prt, a.prt_model, re)
    dhi = face_diffusivity(0.5 * (nut + nN), a.pr, a.prt, a.prt_model, re)

    # Inverse centre-to-centre distances; the mirrored ghost centre makes the
    # wall value 1/dy (init.f90 face_at), i.e. an exact wall gradient.
    inv = np.empty(ny + 1)
    inv[1:ny] = 1.0 / np.diff(yc)
    inv[0] = 1.0 / dy[0]
    inv[ny] = 1.0 / dy[-1]
    ilo = inv[:ny][None, :, None]
    ihi = inv[1:][None, :, None]

    uc = 0.5 * (u + np.concatenate([u[:, :, 1:], u[:, :, :1]], axis=2))  # x periodic
    clo = vlo * 0.5 * (sS + s)
    jlo = clo - dlo * (s - sS) * ilo
    chi = vhi * 0.5 * (s + sN)
    jhi = chi - dhi * (sN - s) * ihi
    cols = [s, s * s, uc * s, clo, jlo, chi, jhi]

    if a.layout == "plane":
        w = np.broadcast_to(dz[:, None, None], s.shape)
        weight = w.sum(axis=0).T.reshape(-1)                     # row = (i, j)
        sums = np.stack([(w * c).sum(axis=0).T.reshape(-1) for c in cols], axis=1)
    else:
        w = dz[:, None, None] * dx[None, None, :]
        w = np.broadcast_to(w, s.shape)
        weight = w.sum(axis=(0, 2))
        sums = np.stack([(w * c).sum(axis=(0, 2)) for c in cols], axis=1)
    return sums, weight


def compare(a) -> int:
    with h5py.File(a.stats, "r") as f:
        nstat = int(f.attrs["nstat"])
        prof = f["profile"][...]
        count = f["count"][...]
        step = int(f.attrs["step"])
        coord = f["coord"][...] if "coord" in f else f["ycoord"][...]
    nscal = nstat // NSTAT
    if not 1 <= a.index <= nscal:
        raise SystemExit(f"scalar index {a.index} outside 1..{nscal}")

    total = None
    weight = None
    for path in a.fields:
        sums, w = snapshot_columns(path, a)
        total = sums if total is None else total + sums
        weight = w if weight is None else weight + w
    mine = total / weight[:, None]

    print(f"{a.stats}: step {step}, {nscal} scalar(s), {count.size} rows, "
          f"{len(a.fields)} snapshot(s)")
    dw = float(np.max(np.abs(count / len(a.fields) - weight / len(a.fields))))
    scale = float(np.max(np.abs(weight))) / len(a.fields)
    print(f"row weights (cell areas): max|dev| = {dw:.3e} of {scale:.3e}")
    worst = dw / scale
    for c in range(NSTAT):
        solver = prof[:, NSTAT * (a.index - 1) + c]
        ref = mine[:, c]
        scale = max(float(np.max(np.abs(ref))), 1.0e-300)
        dev = float(np.max(np.abs(solver - ref))) / scale
        worst = max(worst, dev)
        print(f"  {NAMES[c]:>9}  max|solver - recomputed|/max|.| = {dev:.3e}")
    ok = worst <= a.tol
    print(f"worst relative deviation {worst:.3e} (tolerance {a.tol:g}): "
          + ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


def cmd_rows(a) -> int:
    """Per-LEVEL row bookkeeping, for the refined (multi-level) case.

    The `profile` subcommand assembles a global box and therefore needs a
    single level; this one accumulates straight from the leaves, so it works
    on any leaf layout -- which is what the per-level statistics files need
    gating. It checks the row WEIGHTS (cell areas) and the two purely
    cell-centred columns; the face fluxes are gated at single level by
    `profile` (they need halo values no leaf carries in the snapshot).
    """
    with h5py.File(a.stats, "r") as f:
        nstat = int(f.attrs["nstat"])
        prof = f["profile"][...]
        count = f["count"][...]
    acc = np.zeros((count.size, 2))
    weight = np.zeros(count.size)
    for path in a.fields:
        with h5py.File(path, "r") as f:
            geo = BlockGeometry(f)
            data = f[a.scalar][...]
            for bid in range(geo.n_blocks):
                if int(geo.blocks[bid, 3]) != a.level:
                    continue
                (_, dx), _, (_, dz) = geo.block_axes(bid)
                oy = int(geo.blocks[bid, 1])
                w = (dz[:, None, None] * dx[None, None, :])
                w = np.broadcast_to(w, data[bid].shape)
                for j in range(geo.nb[1]):
                    row = oy + j
                    if row >= count.size:
                        raise SystemExit(f"row {row} outside the {count.size}-row table")
                    s = data[bid][:, j, :]
                    ww = w[:, j, :]
                    weight[row] += float(ww.sum())
                    acc[row, 0] += float((ww * s).sum())
                    acc[row, 1] += float((ww * s * s).sum())
    live = weight > 0.0
    print(f"{a.stats}: level {a.level}, {int(live.sum())} of {count.size} rows covered, "
          f"{len(a.fields)} snapshot(s)")
    worst = float(np.max(np.abs(count[live] - weight[live]))) / float(np.max(weight))
    print(f"  row weights (cell areas): max|dev|/max = {worst:.3e}")
    for c, name in ((0, "<s>"), (1, "<s^2>")):
        ref = acc[live, c] / weight[live]
        solver = prof[live, NSTAT * (a.index - 1) + c]
        scale = max(float(np.max(np.abs(ref))), 1.0e-300)
        dev = float(np.max(np.abs(solver - ref))) / scale
        worst = max(worst, dev)
        print(f"  {name:>9}  max|solver - recomputed|/max|.| = {dev:.3e}")
    if not live.any():
        print("FAIL (no rows covered)")
        return 1
    ok = worst <= a.tol
    print(f"worst relative deviation {worst:.3e} (tolerance {a.tol:g}): "
          + ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


def cmd_diff(a) -> int:
    worst = 0.0
    with h5py.File(a.stats, "r") as f, h5py.File(a.other, "r") as g:
        if int(f.attrs["step"]) != int(g.attrs["step"]):
            raise SystemExit(f"different steps: {f.attrs['step']} vs {g.attrs['step']}")
        print(f"{a.stats} vs {a.other}: step {int(f.attrs['step'])}")
        for name in ("count", "raw_sum", "profile"):
            x, y = f[name][...], g[name][...]
            scale = max(float(np.max(np.abs(x))), 1.0e-300)
            dev = float(np.max(np.abs(x - y))) / scale
            worst = max(worst, dev)
            print(f"  {name:>8}  max|a - b|/max|a| = {dev:.3e}")
    ok = worst <= a.tol
    print(f"worst relative deviation {worst:.3e} (tolerance {a.tol:g}): "
          + ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


def cmd_heat(a) -> int:
    """The solver's runtime body-heat samples against check_scalar_ibm.py's
    `surface`, the Python form validated in S3 by the full energy budget
    (d/dt int s dV to 3.9e-4 on a case with no boundary flux). Both are the
    cancellation-free pair (staircase interface flux, graded-cell
    penalization); they share the geometry and the field, so what is measured
    here is the transcription -- it must agree to round-off."""
    from check_scalar_ibm import body_heat_release

    data = np.loadtxt(a.file, ndmin=2)
    col = 2 + 3 * (a.index - 1)
    worst = 0.0
    for row in data:
        step = int(row[0])
        field = a.pattern.replace("STEP", str(step))
        stair, graded, _, _, _ = body_heat_release(field, a.case, a.scalar,
                                                   a.value, a.pr, a.re)
        for label, mine, ref in (("staircase", row[col], stair),
                                 ("graded", row[col + 1], graded),
                                 ("total", row[col + 2], stair + graded)):
            scale = max(abs(ref), 1.0e-300)
            dev = abs(mine - ref) / scale
            worst = max(worst, dev)
            print(f"  step {step:>6} {label:>10}: solver {mine:.12e}  "
                  f"python {ref:.12e}  rel {dev:.3e}")
    ok = worst <= a.tol
    print(f"worst relative deviation {worst:.3e} (tolerance {a.tol:g}): "
          + ("PASS" if ok else "FAIL"))

    if a.balance and data.shape[0] >= 2:
        # The physics statement, with the SOLVER's own Q: on a case with no
        # boundary flux the heat the body releases is exactly the rate at
        # which the domain stores the scalar, to the trapezoid's O(dt^2).
        from scalar_tools import integrate
        (s1, t1), (s2, t2) = ((row[col + 2], row[1]) for row in data[-2:])
        tot = []
        for row in data[-2:]:
            with h5py.File(a.pattern.replace("STEP", str(int(row[0]))), "r") as f:
                tot.append(integrate(f, a.scalar))
        storage = (tot[1] - tot[0]) / (t2 - t1)
        src = 0.5 * (s1 + s2)
        rel = abs(src - storage) / max(abs(storage), 1.0e-300)
        print(f"energy budget: 1/2[Q(t1) + Q(t2)] = {src:.10e}, "
              f"d/dt int {a.scalar} dV = {storage:.10e}, rel {rel:.3e}")
        okb = rel <= a.balance_tol
        print(f"  (tolerance {a.balance_tol:g}): " + ("PASS" if okb else "FAIL"))
        ok = ok and okb
    return 0 if ok else 1


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for layout in ("profile", "plane"):
        p = sub.add_parser(layout)
        p.add_argument("stats")
        p.add_argument("fields", nargs="+")
        p.add_argument("--scalar", default="theta")
        p.add_argument("--index", type=int, default=1)
        p.add_argument("--pr", type=float, default=0.71)
        p.add_argument("--prt", type=float, default=0.85)
        p.add_argument("--prt-model", default="constant", choices=["constant", "kays"])
        p.add_argument("--walls", type=float, nargs=2, default=[1.0, -1.0])
        p.add_argument("--wall-types", nargs=2, default=["dirichlet", "dirichlet"],
                       choices=["dirichlet", "neumann"])
        p.add_argument("--tol", type=float, default=1.0e-12)
        p.set_defaults(layout=layout, fn=compare)

    p = sub.add_parser("rows")
    p.add_argument("stats")
    p.add_argument("fields", nargs="+")
    p.add_argument("--scalar", default="theta")
    p.add_argument("--index", type=int, default=1)
    p.add_argument("--level", type=int, default=0)
    p.add_argument("--tol", type=float, default=1.0e-12)
    p.set_defaults(fn=cmd_rows)

    p = sub.add_parser("heat")
    p.add_argument("file", help="the runtime file [scalar] heat_interval writes")
    p.add_argument("pattern", help="snapshot path with STEP standing for the step number")
    p.add_argument("--case", required=True)
    p.add_argument("--scalar", default="theta")
    p.add_argument("--index", type=int, default=1)
    p.add_argument("--re", type=float, required=True)
    p.add_argument("--pr", type=float, required=True)
    p.add_argument("--value", type=float, default=1.0)
    p.add_argument("--tol", type=float, default=1.0e-12)
    p.add_argument("--balance", action="store_true",
                   help="also close the energy budget with the last two samples")
    p.add_argument("--balance-tol", type=float, default=1.0e-2)
    p.set_defaults(fn=cmd_heat)

    p = sub.add_parser("diff")
    p.add_argument("stats")
    p.add_argument("other")
    p.add_argument("--tol", type=float, default=1.0e-12)
    p.set_defaults(fn=cmd_diff)

    a = ap.parse_args()
    return a.fn(a)


if __name__ == "__main__":
    sys.exit(main())
