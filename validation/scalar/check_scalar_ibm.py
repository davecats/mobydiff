#!/usr/bin/env python3
"""Checkers for the passive-scalar IMMERSED-BODY gates (increment S3).

Subcommands
-----------
solid     <field.h5> --scalar theta --value V
          Dirichlet body mode: every SOLID cell must hold the body value
          EXACTLY. "Solid" is the analytic wavy-wall indicator (the same
          isInBody the coefficient kernel evaluates at the cell centre), or,
          with --case, the cells whose coef_p tile is >= SOLID/Re.

conserve  <before.h5> <after.h5> --scalar phi
          Adiabatic body mode: int phi dV over every stored cell, before and
          after, relative to max|phi| * V. The face mask is symmetric across
          a face, so the flux form telescopes and the drift must be round-off.

coefp     <case.h5>
          The prepared case file's coef_p_blocks tiles against an independent
          Python transcription of the graded sharp-interface formula
          sum((d0-d)/d)/d0^2 / Re over the analytic wavy wall, at EVERY leaf
          level. Reports the worst relative deviation per level.

surface   <field.h5> --case case.h5 --scalar theta --re R --pr P --value V
          The heated body's TOTAL heat release and its Nusselt number,
          Nu = Q / (Lz pi D_th dT) with D_th = 1/(Re Pr) (the diameter
          cancels, so this is Nu for the unit-diameter cylinder). Q is the
          sum of two cancellation-free terms: the flux across every
          staircase face separating a solid cell from a fluid one, and the
          penalization delivered into the GRADED fluid cells.

balance   <before.h5> <after.h5> --case case.h5 --scalar theta --re R --pr P
          The rigorous check on `surface`: with no boundary flux anywhere
          (periodic or Neumann-0 faces) the heat the body releases must be
          exactly the rate at which the domain stores it,
             1/2 [Q(t1) + Q(t2)]  =  [int s dV(t2) - int s dV(t1)] / (t2-t1),
          to the trapezoid's O(dt^2). Independent of any Nusselt convention,
          and it also reports how much of the truth the naive penalization
          integral (`nusselt`) sees.

nusselt   <field.h5> --case case.h5 --scalar theta --re R --pr P --value V
          The A2 penalization integral Q = int coef_p (s_body - s)/Pr dV,
          transposed from the FORCE diagnostic. KEEP FOR THE RECORD ONLY:
          for a DIRICHLET body it is structurally incomplete. A solid cell's
          stored scalar is the body value to the LAST BIT (that is gate (n)),
          so its coef_p (s_body - s) evaluates to 1e28 x 0 = 0, while the
          cell is really re-heated every substage by exactly the flux it
          loses to its fluid neighbours. Measured on the wavy case: it sees
          63 % of the truth. The force version survives the same arithmetic
          only because u_body = 0 and the stored velocity keeps a ~1e-26
          residual whose product with coef_p is O(1). Use `surface` or `cv`.

cv        <field.h5> --scalar theta --re R --pr P --box X0 X1 Y0 Y1
          The INDEPENDENT far-field cross-check: the same heat, as the net
          flux (u s - D_th grad s) . n through a control-volume border in
          the fluid, plus (with --prev) the storage term inside the box. It
          shares nothing with `surface` but the scalar field itself.
"""

from __future__ import annotations

import argparse
import sys

import numpy as np
import h5py

from scalar_tools import BlockGeometry, subdivide

SOLID = 1.0e30
# ibm.f90 set_ibm_geometry_defaults / wavy_wall_height
WAVY_AMP, WAVY_N, WAVY_PHASE, WAVY_OFFSET = 2.5e-2, 1, 0.0, 1.0e-2
# ibm.f90 bisection
BISECT_TOL, BISECT_MAX = 1.0e-10, 200


def wavy_height(x, lx):
    return WAVY_AMP * 0.5 * (1.0 + np.sin(2.0 * np.pi * WAVY_N * x / lx + WAVY_PHASE)) \
        + WAVY_OFFSET


def in_body(x, y, lx):
    return y < wavy_height(x, lx)


# --------------------------------------------------------------------------
# case-file geometry (grid lines are x_nodes/y_nodes/z_nodes there, and the
# coefficient tiles are ghost-inclusive with (block, i, j, k) index order)
# --------------------------------------------------------------------------
class CaseGeometry:
    def __init__(self, h5: h5py.File):
        self.blocks = h5["blocks"][...]
        # A case file stores the CUBIC block size once (block_nb) and the
        # node lines as x_nodes/y_nodes/z_nodes -- a field file's per-
        # direction block_nb_{x,y,z} and x/y/z are the solver-output layout.
        nb = int(h5.attrs["block_nb"])
        self.nb = (nb, nb, nb)
        self.mask = np.asarray(h5.attrs.get("refine_dims", [1, 1, 1]), dtype=int)
        self.leng = (float(h5.attrs["lx"]), float(h5.attrs["ly"]), float(h5.attrs["lz"]))
        lmax = int(self.blocks[:, 3].max())
        base = [h5["x_nodes"][...], h5["y_nodes"][...], h5["z_nodes"][...]]
        self.lines = [[base[d].copy()] for d in range(3)]
        for d in range(3):
            for _ in range(lmax):
                nxt = subdivide(self.lines[d][-1]) if self.mask[d] else self.lines[d][-1]
                self.lines[d].append(nxt)

    def centres(self, bid: int):
        """Interior cell centres (x, y, z) of block bid, one 1-D array each."""
        ox, oy, oz, lev = (int(v) for v in self.blocks[bid])
        out = []
        for d, o in enumerate((ox, oy, oz)):
            line = self.lines[d][lev]
            lo = line[o:o + self.nb[d]]
            hi = line[o + 1:o + self.nb[d] + 1]
            out.append(0.5 * (lo + hi))
        return out


# --------------------------------------------------------------------------
# (n) Dirichlet body mode: solid cells hold the body value exactly
# --------------------------------------------------------------------------
def cmd_solid(a) -> int:
    with h5py.File(a.field, "r") as f:
        geo = BlockGeometry(f)
        data = f[a.scalar]
        worst, nsolid = 0.0, 0
        solid_mask = None
        if a.case:
            solid_mask = case_solid_mask(a.case, a.re)
        for bid in range(geo.n_blocks):
            x, y, _, _ = geo.mesh(bid)
            if solid_mask is not None:
                sel = solid_mask[bid]
            else:
                sel = np.broadcast_to(in_body(x, y, geo.leng[0]), data[bid].shape)
            if not np.any(sel):
                continue
            dev = np.abs(data[bid][sel] - a.value)
            nsolid += int(np.count_nonzero(sel))
            worst = max(worst, float(np.max(dev)))
    print(f"solid cells: {nsolid}")
    print(f"max|{a.scalar} - {a.value}| over solid cells = {worst:.6e}")
    ok = worst == 0.0 if a.exact else worst <= a.tol
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


def case_solid_mask(case_file: str, re: float):
    """Per-block boolean masks (in FIELD dataset order) of coef_p >= SOLID/Re."""
    out = []
    with h5py.File(case_file, "r") as f:
        coefp = f["coef_p_blocks"]
        nbi, nbj, nbk = coefp.shape[1] - 2, coefp.shape[2] - 2, coefp.shape[3] - 2
        for bid in range(coefp.shape[0]):
            tile = coefp[bid, 1:nbi + 1, 1:nbj + 1, 1:nbk + 1]   # (i, j, k)
            out.append(np.transpose(tile, (2, 1, 0)) >= 0.5 * SOLID / re)  # (k, j, i)
    return out


# --------------------------------------------------------------------------
# (o) adiabatic body mode: int s dV conserved
# --------------------------------------------------------------------------
def cmd_conserve(a) -> int:
    tot, peak, vol = [], 0.0, 0.0
    for path in (a.before, a.after):
        with h5py.File(path, "r") as f:
            geo = BlockGeometry(f)
            data = f[a.scalar]
            s, v = 0.0, 0.0
            for bid in range(geo.n_blocks):
                _, _, _, dV = geo.mesh(bid)
                block = data[bid][...]
                s += float(np.sum(block * dV))
                v += float(np.sum(dV))
                peak = max(peak, float(np.max(np.abs(block))))
            tot.append(s)
            vol = v
    drift = (tot[1] - tot[0]) / (peak * vol)
    print(f"int {a.scalar} dV: {tot[0]:.16e} -> {tot[1]:.16e}")
    print(f"relative drift (per max|s| V, V = {vol:g}) = {drift:.3e}")
    ok = abs(drift) <= a.tol

    # Conservation alone would also hold if the mask never fired, so check
    # the POSITIVE half too: a cell whose six staggered faces are all inside
    # the body is sealed on every side and must be frozen EXACTLY.
    frozen, worst = 0, 0.0
    with h5py.File(a.before, "r") as fb, h5py.File(a.after, "r") as fa:
        geo = BlockGeometry(fb)
        for bid in range(geo.n_blocks):
            (xc, dx), (yc, dy), (zc, dz) = geo.block_axes(bid)
            lx = geo.leng[0]
            # faces: x +- dx/2 at (yc, zc), y +- dy/2 at (xc, zc), z at (xc, yc)
            xm = xc[None, None, :] - 0.5 * dx[None, None, :]
            xp = xc[None, None, :] + 0.5 * dx[None, None, :]
            y0 = yc[None, :, None]
            ym = y0 - 0.5 * dy[None, :, None]
            yp = y0 + 0.5 * dy[None, :, None]
            x0 = xc[None, None, :]
            sel = (in_body(xm, y0, lx) & in_body(xp, y0, lx)
                   & in_body(x0, ym, lx) & in_body(x0, yp, lx)
                   & in_body(x0, y0, lx))
            sel = np.broadcast_to(sel, fb[a.scalar][bid].shape)
            if not np.any(sel):
                continue
            d = np.abs(fa[a.scalar][bid][sel] - fb[a.scalar][bid][sel])
            frozen += int(np.count_nonzero(sel))
            worst = max(worst, float(np.max(d)))
    print(f"sealed (all six staggered faces solid) cells: {frozen}, "
          f"max|delta {a.scalar}| = {worst:.3e}")
    ok = ok and frozen > 0 and worst == 0.0
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


# --------------------------------------------------------------------------
# (q) the case file's coef_p tiles vs an independent transcription
# --------------------------------------------------------------------------
def bisect_surface(xa, xb, lx):
    """ibm.f90 bisection, transcribed: returns the surface point between the
    fluid end xa and the solid end xb (3-vectors)."""
    a, b = np.array(xa, float), np.array(xb, float)
    la = in_body(a[0], a[1], lx)
    m = a.copy()
    for _ in range(BISECT_MAX):
        m = 0.5 * (a + b)
        if np.linalg.norm(b - a) < BISECT_TOL:
            break
        if in_body(m[0], m[1], lx) == la:
            a = m
        else:
            b = m
    return m


def graded_coeff(xc, yc, zc, i, j, k, lx, re):
    """coef at the cell-centred point (i,j,k) of the given centre lines."""
    p = np.array([xc[i], yc[j], zc[k]])
    if in_body(p[0], p[1], lx):
        return SOLID / re
    total = 0.0
    for d, (line, idx) in enumerate(((xc, i), (yc, j), (zc, k))):
        for step in (-1, 1):
            n = idx + step
            if n < 0 or n >= line.size:
                return None            # neighbour outside the stored block
            q = p.copy()
            q[d] = line[n]
            if not in_body(q[0], q[1], lx):
                continue
            d0 = np.linalg.norm(q - p)
            s = bisect_surface(p, q, lx)
            dd = np.linalg.norm(s - p)
            total += ((d0 - dd) / dd) / d0**2
    return total / re


def cmd_coefp(a) -> int:
    with h5py.File(a.case, "r") as f:
        geo = CaseGeometry(f)
        coefp = f["coef_p_blocks"][...]
        lx = geo.leng[0]
        levels = {}
        rng = np.random.default_rng(12345)
        for bid in range(geo.blocks.shape[0]):
            lev = int(geo.blocks[bid, 3])
            xc, yc, zc = geo.centres(bid)
            # Only cells whose 6 neighbours are inside this block's INTERIOR
            # are checked: the ghost ring's node lines follow blocks.f90's
            # extension rules, which are not the point of this gate.
            ni, nj, nk = xc.size, yc.size, zc.size
            cand = [(i, j, k)
                    for i in range(1, ni - 1) for j in range(1, nj - 1)
                    for k in range(1, nk - 1)]
            if a.sample and len(cand) > a.sample:
                cand = [cand[t] for t in rng.choice(len(cand), a.sample, replace=False)]
            for (i, j, k) in cand:
                ref = graded_coeff(xc, yc, zc, i, j, k, lx, a.re)
                if ref is None:
                    continue
                got = coefp[bid, i + 1, j + 1, k + 1]
                scale = max(abs(ref), 1.0e-12)
                rel = abs(got - ref) / scale
                worst, n, ngrad, nsol = levels.get(lev, (0.0, 0, 0, 0))
                levels[lev] = (max(worst, rel), n + 1,
                               ngrad + (0.0 < ref < 0.5 * SOLID / a.re),
                               nsol + (ref >= 0.5 * SOLID / a.re))
    ok = True
    for lev in sorted(levels):
        worst, n, ngrad, nsol = levels[lev]
        print(f"level {lev}: {n} cells checked ({ngrad} graded, {nsol} solid), "
              f"worst relative deviation {worst:.3e}")
        ok = ok and worst <= a.tol
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


# --------------------------------------------------------------------------
# (s) Nusselt number: the penalization integral and the CV border flux
# --------------------------------------------------------------------------
def penalization_heat(field: str, case: str, scalar: str, value: float, pr: float):
    """(Q, int s dV, t) of one snapshot: Q = int coef_p (s_body - s)/Pr dV."""
    with h5py.File(case, "r") as c:
        coefp = c["coef_p_blocks"][...]
    with h5py.File(field, "r") as f:
        geo = BlockGeometry(f)
        data = f[scalar]
        if coefp.shape[0] != geo.n_blocks:
            raise SystemExit("case file and field file have different leaf counts")
        nbx, nby, nbz = geo.nb
        q, tot = 0.0, 0.0
        for bid in range(geo.n_blocks):
            _, _, _, dV = geo.mesh(bid)
            # coef_p tile is (i, j, k) ghost-inclusive -> field order (k, j, i)
            cp = np.transpose(coefp[bid, 1:nbx + 1, 1:nby + 1, 1:nbz + 1], (2, 1, 0))
            block = data[bid][...]
            q += float(np.sum(cp * (value - block) * dV))
            tot += float(np.sum(block * dV))
        return q / pr, tot, float(f.attrs["t_current"]), geo


def cmd_balance(a) -> int:
    f1, g1, s1, t1, _ = body_heat_release(a.before, a.case, a.scalar, a.value,
                                          a.pr, a.re)
    f2, g2, s2, t2, _ = body_heat_release(a.after, a.case, a.scalar, a.value,
                                          a.pr, a.re)
    q1, q2 = f1 + g1, f2 + g2
    src = 0.5 * (q1 + q2)
    storage = (s2 - s1) / (t2 - t1)
    rel = abs(src - storage) / max(abs(storage), 1.0e-300)
    print(f"body heat release  1/2 [Q(t1) + Q(t2)] = {src:.10e}"
          f"   (Q: {q1:.6e} -> {q2:.6e})")
    print(f"storage rate  d/dt int {a.scalar} dV     = {storage:.10e}"
          f"   (dt = {t2 - t1:g})")
    print(f"relative difference = {rel:.3e}")
    print(f"   [the penalization integral alone would give {0.5 * (g1 + g2):.6e}"
          f" = {100 * 0.5 * (g1 + g2) / storage:.0f} % of the truth]")
    ok = rel <= a.tol
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


def body_heat_release(field: str, case: str, scalar: str, value: float,
                      pr: float, re: float):
    """The body's heat release, WITHOUT the cancellation cmd_nusselt suffers.

    Total = (the flux across every staircase face separating a solid cell
    from a fluid one) + (the penalization delivered directly into the GRADED
    fluid cells, which is all cmd_nusselt can see). The first term is the one
    the penalization integral loses: a solid cell's stored scalar IS the body
    value to the last bit, so coef_p (s_body - s) reads 1e28 x 0 = 0 there,
    while the cell is in fact re-heated every substage by exactly the flux it
    loses to its fluid neighbours.

    Returns (staircase, graded, int s dV, t, geo)."""
    dth = 1.0 / (re * pr)
    with h5py.File(case, "r") as c:
        coefp = c["coef_p_blocks"][...]
    with h5py.File(field, "r") as f:
        geo = BlockGeometry(f)
        s = assemble(f, scalar, geo)
        u = assemble(f, "un", geo)
        v = assemble(f, "vn", geo)
        w = assemble(f, "wn", geo)
        _, _, dx, dy, dz = uniform_axes(geo)
        solid = assemble_case(coefp, geo, f) >= 0.5 * SOLID / re

    # One term per direction: the face between cell n-1 and cell n carries
    # the staggered velocity of cell n, and the flux counts positive INTO
    # the fluid side.
    staircase = 0.0
    for axis, (vel, d, area) in enumerate(((w, dz, dx * dy),
                                           (v, dy, dx * dz),
                                           (u, dx, dy * dz))):
        lo = [slice(None)] * 3
        hi = [slice(None)] * 3
        lo[axis] = slice(0, -1)
        hi[axis] = slice(1, None)
        lo, hi = tuple(lo), tuple(hi)
        sl, sr = solid[lo], solid[hi]
        face = sl ^ sr                       # exactly one side solid
        flux = vel[hi] * 0.5 * (s[lo] + s[hi]) - dth * (s[hi] - s[lo]) / d
        sign = np.where(sl, 1.0, -1.0)       # solid low => +axis enters the fluid
        staircase += float(np.sum(np.where(face, sign * flux, 0.0))) * area

    graded, total_s, t, _ = penalization_heat(field, case, scalar, value, pr)
    return staircase, graded, total_s, t, geo


def cmd_surface(a) -> int:
    staircase, graded, _, _, geo = body_heat_release(
        a.field, a.case, a.scalar, a.value, a.pr, a.re)
    dth = 1.0 / (a.re * a.pr)
    lz = geo.leng[2]
    total = staircase + graded
    nu = total / (lz * np.pi * dth * a.value)
    print(f"staircase-interface flux/Lz     = {staircase / lz:.6e}")
    print(f"graded fluid-cell penalization  = {graded / lz:.6e}  (all cmd_nusselt sees)")
    print(f"total body heat release/Lz      = {total / lz:.6e}")
    print(f"Nu (surface) = {nu:.4f}   (D_th = {dth:.6e})")
    if a.band:
        lo, hi = a.band
        ok = lo <= nu <= hi
        print(f"literature band [{lo}, {hi}]: " + ("PASS" if ok else "FAIL"))
        return 0 if ok else 1
    return 0


def assemble_case(tiles, geo, field) -> np.ndarray:
    """A per-leaf (block, i, j, k) case-file tile set -> a global (nz,ny,nx)
    array, interior cells only (single level, like assemble())."""
    nbx, nby, nbz = geo.nb
    out = np.zeros((int(field.attrs["nz"]), int(field.attrs["ny"]),
                    int(field.attrs["nx"])))
    for bid in range(geo.n_blocks):
        ox, oy, oz, lev = (int(t) for t in geo.blocks[bid])
        if lev != 0:
            raise SystemExit("surface: single-level fields only")
        out[oz:oz + nbz, oy:oy + nby, ox:ox + nbx] = np.transpose(
            tiles[bid, 1:nbx + 1, 1:nby + 1, 1:nbz + 1], (2, 1, 0))
    return out


def cmd_nusselt(a) -> int:
    """Q = int coef_p (s_body - s)/Pr dV; Nu = Q / (Lz pi D_th dT).

    CAVEAT (measured, see README (s)): for a DIRICHLET body this sees only
    the GRADED FLUID band. A solid cell's stored scalar is the body value
    exactly, so its coef_p (s_body - s) is 1e28 x 0 = 0 while the cell is
    really re-heated every substage. Use `surface` or `cv` for the total."""
    q, _, _, geo = penalization_heat(a.field, a.case, a.scalar, a.value, a.pr)
    dth = 1.0 / (a.re * a.pr)
    lz = geo.leng[2]
    nu = q / (lz * np.pi * dth * a.value)
    print(f"penalization heat release Q/Lz = {q / lz:.6e}")
    print(f"Nu = {nu:.4f}   (D_th = {dth:.6e}, Pr = {a.pr}, Re = {a.re})")
    if a.band:
        lo, hi = a.band
        ok = lo <= nu <= hi
        print(f"literature band [{lo}, {hi}]: " + ("PASS" if ok else "FAIL"))
        return 0 if ok else 1
    return 0


def cmd_cv(a) -> int:
    """Net heat leaving a control-volume border in the fluid: the sum of
    (u s - D_th ds/dn) over the four vertical faces of the box, per unit
    span. Cell-centred data on a uniform mesh, so a face between two cells
    carries the average scalar and the two-point gradient.

    With --prev the STORAGE term d/dt int_box s dV is added, which turns the
    balance into an identity that holds whether or not the field has
    converged: outflow + storage = the heat the body released inside the
    box. Without it the comparison assumes a steady state."""
    dth = 1.0 / (a.re * a.pr)
    x0, x1, y0, y1 = a.box
    with h5py.File(a.field, "r") as f:
        geo = BlockGeometry(f)
        s = assemble(f, a.scalar, geo)
        u = assemble(f, "un", geo)
        v = assemble(f, "vn", geo)
        xc, yc, dx, dy, dz = uniform_axes(geo)
        t_now = float(f.attrs["t_current"])
    # Face-centred quantities: u lives on the x face i (left face of cell i).
    ix0, ix1 = nearest_face(xc, dx, x0), nearest_face(xc, dx, x1)
    iy0, iy1 = nearest_face(yc, dy, y0), nearest_face(yc, dy, y1)
    # The z faces need no term: the box spans the periodic z direction, so
    # whatever leaves one z face enters the other.

    def face_x(i, sign):
        # cells i-1 (left) and i (right) share the x face i
        sface = 0.5 * (s[:, iy0:iy1, i - 1] + s[:, iy0:iy1, i])
        grad = (s[:, iy0:iy1, i] - s[:, iy0:iy1, i - 1]) / dx
        flux = u[:, iy0:iy1, i] * sface - dth * grad
        return sign * float(np.sum(flux)) * dy * dz

    def face_y(j, sign):
        sface = 0.5 * (s[:, j - 1, ix0:ix1] + s[:, j, ix0:ix1])
        grad = (s[:, j, ix0:ix1] - s[:, j - 1, ix0:ix1]) / dy
        flux = v[:, j, ix0:ix1] * sface - dth * grad
        return sign * float(np.sum(flux)) * dx * dz

    total = face_x(ix1, +1.0) + face_x(ix0, -1.0) + face_y(iy1, +1.0) + face_y(iy0, -1.0)
    lz = geo.leng[2]
    storage = 0.0
    if a.prev:
        with h5py.File(a.prev, "r") as g:
            sp = assemble(g, a.scalar, BlockGeometry(g))
            t_prev = float(g.attrs["t_current"])
        box = (slice(None), slice(iy0, iy1), slice(ix0, ix1))
        cell = dx * dy * dz
        storage = float(np.sum(s[box] - sp[box])) * cell / (t_now - t_prev)
    nu = (total + storage) / (lz * np.pi * dth * a.value)
    print(f"CV box x [{xc[ix0] - 0.5 * dx:.4f}, {xc[ix1] - 0.5 * dx:.4f}] "
          f"y [{yc[iy0] - 0.5 * dy:.4f}, {yc[iy1] - 0.5 * dy:.4f}]")
    print(f"net heat outflow/Lz = {total / lz:.6e}"
          + (f", storage/Lz = {storage / lz:.6e} "
             f"({100 * storage / max(abs(total), 1e-300):.1f} % of it)" if a.prev else ""))
    print(f"Nu (CV) = {nu:.4f}")
    if a.nu_pen is not None:
        rel = abs(nu - a.nu_pen) / abs(a.nu_pen)
        print(f"vs the penalization Nu {a.nu_pen:.4f}: {100 * rel:.2f} %")
        ok = rel <= a.tol
        print("PASS" if ok else "FAIL")
        return 0 if ok else 1
    return 0


def assemble(f: h5py.File, name: str, geo: BlockGeometry) -> np.ndarray:
    """Single-level block-table dataset -> one global (nz, ny, nx) array."""
    nbx, nby, nbz = geo.nb
    nx = int(f.attrs["nx"])
    ny = int(f.attrs["ny"])
    nz = int(f.attrs["nz"])
    out = np.zeros((nz, ny, nx))
    data = f[name]
    for bid in range(geo.n_blocks):
        ox, oy, oz, lev = (int(t) for t in geo.blocks[bid])
        if lev != 0:
            raise SystemExit("cv: single-level fields only")
        out[oz:oz + nbz, oy:oy + nby, ox:ox + nbx] = data[bid]
    return out


def uniform_axes(geo: BlockGeometry):
    x = geo.lines[0][0]
    y = geo.lines[1][0]
    z = geo.lines[2][0]
    dx = float(x[1] - x[0])
    dy = float(y[1] - y[0])
    dz = float(z[1] - z[0])
    return 0.5 * (x[:-1] + x[1:]), 0.5 * (y[:-1] + y[1:]), dx, dy, dz


def nearest_face(centres, d, target):
    """Index of the cell whose LOW face is closest to target."""
    faces = centres - 0.5 * d
    return int(np.argmin(np.abs(faces - target)))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("solid")
    p.add_argument("field")
    p.add_argument("--scalar", default="theta")
    p.add_argument("--value", type=float, default=1.0)
    p.add_argument("--case", default=None)
    p.add_argument("--re", type=float, default=100.0)
    p.add_argument("--tol", type=float, default=0.0)
    p.add_argument("--exact", action="store_true", default=True)
    p.set_defaults(fn=cmd_solid)

    p = sub.add_parser("conserve")
    p.add_argument("before")
    p.add_argument("after")
    p.add_argument("--scalar", default="phi")
    p.add_argument("--tol", type=float, default=1.0e-14)
    p.set_defaults(fn=cmd_conserve)

    p = sub.add_parser("coefp")
    p.add_argument("case")
    p.add_argument("--re", type=float, default=100.0)
    p.add_argument("--tol", type=float, default=1.0e-12)
    p.add_argument("--sample", type=int, default=0,
                   help="check at most this many cells per block (0 = all)")
    p.set_defaults(fn=cmd_coefp)

    p = sub.add_parser("balance")
    p.add_argument("before")
    p.add_argument("after")
    p.add_argument("--case", required=True)
    p.add_argument("--scalar", default="theta")
    p.add_argument("--re", type=float, required=True)
    p.add_argument("--pr", type=float, required=True)
    p.add_argument("--value", type=float, default=1.0)
    p.add_argument("--tol", type=float, default=1.0e-3)
    p.set_defaults(fn=cmd_balance)

    p = sub.add_parser("surface")
    p.add_argument("field")
    p.add_argument("--case", required=True)
    p.add_argument("--scalar", default="theta")
    p.add_argument("--re", type=float, required=True)
    p.add_argument("--pr", type=float, required=True)
    p.add_argument("--value", type=float, default=1.0)
    p.add_argument("--band", type=float, nargs=2, default=None)
    p.set_defaults(fn=cmd_surface)

    p = sub.add_parser("nusselt")
    p.add_argument("field")
    p.add_argument("--case", required=True)
    p.add_argument("--scalar", default="theta")
    p.add_argument("--re", type=float, required=True)
    p.add_argument("--pr", type=float, required=True)
    p.add_argument("--value", type=float, default=1.0)
    p.add_argument("--band", type=float, nargs=2, default=None)
    p.set_defaults(fn=cmd_nusselt)

    p = sub.add_parser("cv")
    p.add_argument("field")
    p.add_argument("--scalar", default="theta")
    p.add_argument("--re", type=float, required=True)
    p.add_argument("--pr", type=float, required=True)
    p.add_argument("--value", type=float, default=1.0)
    p.add_argument("--box", type=float, nargs=4, required=True)
    p.add_argument("--prev", default=None,
                   help="earlier snapshot: adds the storage term d/dt int_box s dV")
    p.add_argument("--nu-pen", type=float, default=None)
    p.add_argument("--tol", type=float, default=0.10)
    p.set_defaults(fn=cmd_cv)

    a = ap.parse_args()
    return a.fn(a)


if __name__ == "__main__":
    sys.exit(main())
