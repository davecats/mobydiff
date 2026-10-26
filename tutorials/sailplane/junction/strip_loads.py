#!/usr/bin/env python3
"""Surface pressure and skin friction on an immersed STL body from a
block-table field snapshot, integrated per spanwise strip.

    strip_loads.py FIELD.h5 --case CASE.h5 --stl BODY.stl [--scale S] [--translate X Y Z]
                   --wing-y0 Y0 --wing-y1 Y1 [--strips W] [--wing-xmax X]
                   [--tri-size 0.01] [--cp-start 1.5] [--cp-depth 10]
                   [--cp-station Y ...] [--pinf-box x0 x1 y0 y1 z0 z1] [--out strips]

Step 2a of docs/next_session_sailplane_junction.md. The sailplane case has
no runtime force output ([case] name = generic), so every load of the
junction study comes from the field file through this script.

Method (per SUB-TRIANGLE of the STL, everything in the solver's grid units):

  Geometry. The STL (ASCII or binary) is transformed as moby_prepare does
  (v*scale + translate), the outward normal is the winding's cross product
  (checked against the stored facet normal, which the sailplane STL and
  make_airfoil_stl agree with), and every triangle is bisected along its
  longest edge until no edge exceeds --tri-size (a wing panel 3 m long
  spanwise becomes a row of pieces, a sliver stays a sliver): one Cp sample
  per piece, pieces of 0.01 c = 5 cells at c/512.

  Cp. pn is sampled trilinearly from the leaf block holding the point at
  centroid + s*h along the normal, s = cp_start, cp_start + 1, ... up to
  cp_depth samples, h the cell size of the leaf holding the centroid (the
  first sample sits just outside the penalization band), fitted linearly in
  s and read at the wall, s = 0 (the NACA script's wall extrapolation; one
  sample with --cp-depth 1). C_p = 2 (p_w - p_inf), p_inf the mean level-0
  pressure of the first cell column behind the inlet (or of --pinf-box).

  Cf. The PENALIZATION-BAND estimator of tutorials/naca/rans/postProcess/
  surface_cp_cf.py: the band faces (0 < coef < solid) carry the wall
  gradient the scheme imposes, so a through-origin fit of each staggered
  velocity component against (dn t_comp) over the band points assigned to
  a piece (nearest centroid, 0 < dn < 2 h, in-plane offset under tri_size)
  gives G = du_t/dn and Cf = 2 nu G, signed along the piece's CHORDWISE
  tangent t = x-direction projected onto the facet, TE-ward (attached flow
  positive, the sign is the separation line).

  Strips. The wing is the set of pieces with wing_y0 <= y < wing_y1 and
  x < wing_xmax; they are binned in y by --strips. Per strip the planform
  area is the projected upper surface A = sum(max(n_z, 0) dA), so
  c_l = sum(-C_p n_z dA)/A needs no chord or width convention; c_dp and
  c_df (pressure, friction) likewise; c_m about the quarter chord of the
  strip's own chord (x extent, LE z) nose-up positive, divided by c A; the
  upper-surface separation as the first / last x/c with C_f < 0 beyond
  x/c 0.02 and the separated area fraction. Pieces outside the wing bins
  (the fuselage) are summed as one "rest" row (force / q only).

VALIDATION (the gate of 2a, 2026-10-26, workstation, this script at its
defaults): the NACA 0012 nose case of tutorials/naca/rans (Re 4e5, inflow at
5 deg, the section unpitched, refine_dims = xz, nose at (50, 48)), snapshot
~/horeka_stage/naca/c11_nose_660000.h5 with assets/geometry/
ibm_coeff_c11_nose.h5 and n0012_b11.stl (binary, span-long slivers:
218 576 pieces inside the domain, levels 11 and 12), run as

    strip_loads.py c11_nose_660000.h5 --case ibm_coeff_c11_nose.h5 --stl n0012_b11.stl
        --wing-y0 0 --wing-y1 0.1875 --strips 0.1875 --alpha 5
        --cp-station 0.09375 --station-width 0.03     (15 s, 1.1 GB)

  c_l 0.5171 (pressure 0.5137 + friction) against the runtime control-volume
    lift 0.5205 +- 0.0005 of the same state (README): -0.65 %; the offline CV
    boxes read 0.5224 +- 0.004. Cp depth 1 / 4 / 10 moves it by +-0.0003.
  c_d 0.0141 = c_dp 0.0053 + c_df 0.0088 against CV 0.0131 +- 0.0002 (the
    offline boxes 0.0133 +- 0.001, 0.0148 on the 2.5 c box): the pressure
    drag is the fore/aft cancellation of a Cp known to 0.003.
  Cp(x/c) against surface_cp_cf.py's shipped curve (assets/reval_2026-10-24/
    cpcf_c11_nose_final_cp.dat, 300 arc-length bins): rms 0.0037 upper /
    0.0034 lower, beyond x/c 0.05 rms 0.002, max 0.009; the largest
    differences (0.02-0.03) sit at the stagnation point and the suction peak,
    where the shipped values are bin averages over 0.0068 c of arc and these
    are point values (the field holds Cp 1.023 at the stagnation point out to
    8 cells from the wall; the shipped 1.0024 is its bin's mean).
  Cf(x/c) against the shipped cf curve, x/c 0.05-0.95: median relative
    difference -0.9 % upper, -0.4 % lower (median |rel| 1 %). No Cf < 0 on
    the upper surface (attached, as the reference).
  Cp depth: depth 4 from 1.5 cells is the default (rms 0.0037 / 0.0034);
    one sample 0.0051 / 0.0037, ten samples 0.0031 / 0.0034 but Cp_min 0.016
    deeper; start 2.5 cells no better.

Writes OUT_strips.csv, OUT_pieces.npz (every piece: centroid, normal, area,
Cp, Cf, strip) and OUT_cp_y<Y>.csv per --cp-station (x/c, Cp, Cf, upper).
"""
import argparse
import struct
import sys

import h5py
import numpy as np
from scipy.spatial import cKDTree

SOLID_MIN = 1.0e10          # a face inside the body carries SOLID/Re, not a graded coef


# ----------------------------------------------------------------------------- STL
def read_stl(path):
    """Triangles (n, 3, 3) and stored facet normals (n, 3), ASCII or binary."""
    with open(path, "rb") as f:
        head = f.read(5)
    if head == b"solid":
        v, nrm = [], []
        with open(path) as f:
            for line in f:
                w = line.split()
                if not w:
                    continue
                if w[0] == "vertex":
                    v.append([float(t) for t in w[1:4]])
                elif w[0] == "facet":
                    nrm.append([float(t) for t in w[2:5]])
        return np.array(v).reshape(-1, 3, 3), np.array(nrm)
    with open(path, "rb") as f:
        f.read(80)
        n = struct.unpack("<I", f.read(4))[0]
        d = np.frombuffer(f.read(), dtype=np.dtype([("n", "<3f4"), ("v", "<9f4"), ("a", "<u2")]), count=n)
    return d["v"].reshape(-1, 3, 3).astype(float), d["n"].astype(float)


def subdivide(tri, size):
    """Split every triangle until no edge exceeds size, winding preserved.
    The two edges meeting at the APEX (the vertex opposite the shortest
    edge) are the long ones; both are divided into n = ceil(longest/size)
    segments and the triangle becomes one row of 2n-1 pieces. A sliver
    (a span-long STL facet of a few cells' width) is thereby one row of
    short pieces, not the n^2 of a bisection; a large triangle needs a few
    passes."""
    done, todo = [], tri
    while len(todo):
        e = np.stack([np.linalg.norm(todo[:, 1] - todo[:, 0], axis=1),
                      np.linalg.norm(todo[:, 2] - todo[:, 1], axis=1),
                      np.linalg.norm(todo[:, 0] - todo[:, 2], axis=1)], 1)
        big = e.max(1) > size
        done.append(todo[~big])
        t, e = todo[big], e[big]
        if not len(t):
            break
        apex = (e.argmin(1) + 2) % 3                    # e0 = (v0,v1) faces v2, ...
        idx = (np.arange(3)[None, :] + apex[:, None]) % 3
        t = np.take_along_axis(t, idx[:, :, None], axis=1)   # apex first
        n = np.ceil(np.maximum(np.linalg.norm(t[:, 1] - t[:, 0], axis=1),
                               np.linalg.norm(t[:, 2] - t[:, 0], axis=1)) / size).astype(int)
        pieces = []
        for nk in np.unique(n):
            u = t[n == nk]
            s = np.arange(nk + 1) / nk
            a = u[:, None, 0] + s[None, :, None] * (u[:, None, 1] - u[:, None, 0])   # (m, nk+1, 3)
            b = u[:, None, 0] + s[None, :, None] * (u[:, None, 2] - u[:, None, 0])
            pieces.append(np.stack([u[:, 0], a[:, 1], b[:, 1]], 1))
            if nk > 1:
                pieces.append(np.stack([a[:, 1:-1], a[:, 2:], b[:, 2:]], 2).reshape(-1, 3, 3))
                pieces.append(np.stack([a[:, 1:-1], b[:, 2:], b[:, 1:-1]], 2).reshape(-1, 3, 3))
        todo = np.concatenate(pieces)
    return np.concatenate(done)


# ----------------------------------------------------------------------------- field
class Field:
    """A block-table snapshot: leaf lookup for points, trilinear sampling of
    a cell-centred dataset, staggered band points with their coefficients."""

    def __init__(self, field, case):
        self.f = h5py.File(field, "r")
        self.c = h5py.File(case, "r")
        f = self.f
        self.blocks = f["blocks"][...].astype(np.int64)
        if not np.array_equal(self.c["blocks"][...], self.blocks):
            raise SystemExit(f"{case}: blocks table differs from {field}")
        self.nb = int(f.attrs["block_nb_x"])
        if int(f.attrs["block_nb_y"]) != self.nb or int(f.attrs["block_nb_z"]) != self.nb:
            raise SystemExit("non-cubic blocks are not handled")
        self.n0 = np.array([int(f.attrs[k]) for k in ("nx", "ny", "nz")])
        self.L = np.array([float(f.attrs[k]) for k in ("lx", "ly", "lz")])
        self.m = np.asarray(f.attrs.get("refine_dims", [1, 1, 1]), dtype=np.int64)
        self.lmax = int(self.blocks[:, 3].max())
        self.re = float(f.attrs["re"])
        # cells per direction at every level, and the leaf key table
        self.n_at = np.array([self.n0 * 2 ** (l * self.m) for l in range(self.lmax + 1)])
        self.key = self._encode(self.blocks[:, 3], self.blocks[:, :3])
        self.order = np.argsort(self.key)
        self.sorted_key = self.key[self.order]
        self._cache = {}

    def h(self, level):
        return self.L / self.n_at[level]

    def _encode(self, level, origin):
        n = self.n_at[self.lmax] + 1
        return ((level * n[0] + origin[:, 0]) * n[1] + origin[:, 1]) * n[2] + origin[:, 2]

    def locate(self, pts):
        """Row of the leaf holding each point (-1 outside / in a removed
        block) and its level."""
        pts = np.asarray(pts, dtype=float)
        row = np.full(len(pts), -1, dtype=np.int64)
        lev = np.full(len(pts), -1, dtype=np.int64)
        inside = np.all((pts >= 0) & (pts < self.L), axis=1)
        todo = inside.copy()
        for l in range(self.lmax, -1, -1):
            if not todo.any():
                break
            h = self.h(l)
            cell = np.floor(pts[todo] / h).astype(np.int64)
            cell = np.minimum(cell, self.n_at[l] - 1)
            org = (cell // self.nb) * self.nb
            k = self._encode(np.full(len(org), l), org)
            pos = np.searchsorted(self.sorted_key, k)
            pos = np.minimum(pos, len(self.sorted_key) - 1)
            hit = self.sorted_key[pos] == k
            at = np.nonzero(todo)[0][hit]
            row[at] = self.order[pos[hit]]
            lev[at] = l
            todo[at] = False
        return row, lev

    def rows(self, name, rows):
        """Dataset rows (any order) with a small cache keyed on the dataset."""
        rows = np.asarray(rows, dtype=np.int64)
        uniq, inv = np.unique(rows, return_inverse=True)
        if name not in self._cache or not np.array_equal(self._cache[name][0], uniq):
            src = self.c if name in ("coef_blocks",) else self.f
            data = src[name][uniq] if len(uniq) else src[name][:0]
            self._cache[name] = (uniq, data)
        return self._cache[name][1][inv]

    def sample(self, name, pts, row=None, lev=None):
        """Trilinear sample of a cell-centred dataset (z, y, x tiles) at
        points, clamped to the holding block's own cells. NaN outside."""
        if row is None:
            row, lev = self.locate(pts)
        out = np.full(len(pts), np.nan)
        ok = row >= 0
        if not ok.any():
            return out
        pts, row, lev = pts[ok], row[ok], lev[ok]
        data = self.rows(name, row)                       # (n, nb, nb, nb) as (z, y, x)
        h = self.L[None, :] / self.n_at[lev]
        u = (pts - self.blocks[row, :3] * h) / h - 0.5    # cell-centre coordinates
        i0 = np.floor(u).astype(np.int64)
        fr = u - i0
        i0c = np.clip(i0, 0, self.nb - 1)
        i1c = np.clip(i0 + 1, 0, self.nb - 1)
        r = np.arange(len(row))
        val = np.zeros(len(row))
        for dx, wx in ((i0c[:, 0], 1 - fr[:, 0]), (i1c[:, 0], fr[:, 0])):
            for dy, wy in ((i0c[:, 1], 1 - fr[:, 1]), (i1c[:, 1], fr[:, 1])):
                for dz, wz in ((i0c[:, 2], 1 - fr[:, 2]), (i1c[:, 2], fr[:, 2])):
                    val += wx * wy * wz * data[r, dz, dy, dx]
        out[ok] = val
        return out

    def band_points(self, rows):
        """Every staggered face of the given blocks that lies in the
        penalization band (0 < coef < SOLID_MIN): position, velocity
        component, component index, cell size. Also the solid-cell flags
        (all six faces solid) of those blocks, for the Cp sanity check."""
        rows = np.unique(rows)
        nb = self.nb
        coef = self.c["coef_blocks"][rows][:, 1:nb + 1, 1:nb + 1, 1:nb + 1, :]   # (n, x, y, z, 3)
        coef = coef.transpose(0, 3, 2, 1, 4)                                       # -> (n, z, y, x, 3)
        lev = self.blocks[rows, 3]
        h = self.L[None, :] / self.n_at[lev]                                      # (n, 3)
        org = self.blocks[rows, :3] * h
        kk, jj, ii = np.mgrid[0:nb, 0:nb, 0:nb]
        pos, val, comp, hh, brow = [], [], [], [], []
        for ci, (name, off) in enumerate((("un", (0.0, 0.5, 0.5)),
                                          ("vn", (0.5, 0.0, 0.5)),
                                          ("wn", (0.5, 0.5, 0.0)))):
            q = self.f[name][rows]                                                 # (n, z, y, x)
            c = coef[..., ci]
            m = (c > 0) & (c < SOLID_MIN)
            n_i, z_i, y_i, x_i = np.nonzero(m)
            p = org[n_i] + np.stack([(x_i + off[0]) * h[n_i, 0],
                                     (y_i + off[1]) * h[n_i, 1],
                                     (z_i + off[2]) * h[n_i, 2]], 1)
            pos.append(p); val.append(q[m]); comp.append(np.full(m.sum(), ci))
            hh.append(h[n_i, 0]); brow.append(rows[n_i])
        solid = np.all(coef >= SOLID_MIN, axis=4)
        # a cell is solid when its six faces are: low faces are its own, the
        # high ones belong to the next cell (shift by one, halo faces unknown)
        s6 = solid.copy()
        s6[:, :, :, :-1] &= solid[:, :, :, 1:]
        s6[:, :, :-1, :] &= solid[:, :, 1:, :]
        s6[:, :-1, :, :] &= solid[:, 1:, :, :]
        return (np.concatenate(pos), np.concatenate(val), np.concatenate(comp),
                np.concatenate(hh), np.concatenate(brow), dict(zip(rows.tolist(), s6)))

    def p_inf(self, box=None):
        """Mean level-0 pressure in the first cell column behind the inlet
        (x_min), or in a box."""
        blocks = self.blocks
        h0 = self.h(0)
        tot, n = 0.0, 0
        for b in np.nonzero(blocks[:, 3] == 0)[0]:
            ox, oy, oz = blocks[b, :3]
            if box is None:
                if ox != 0:
                    continue
                tot += float(self.f["pn"][b][:, :, 0].sum()); n += self.nb * self.nb
            else:
                xc = (ox + 0.5 + np.arange(self.nb)) * h0[0]
                yc = (oy + 0.5 + np.arange(self.nb)) * h0[1]
                zc = (oz + 0.5 + np.arange(self.nb)) * h0[2]
                sx = (xc > box[0]) & (xc < box[1]); sy = (yc > box[2]) & (yc < box[3])
                sz = (zc > box[4]) & (zc < box[5])
                if sx.any() and sy.any() and sz.any():
                    tot += float(self.f["pn"][b][np.ix_(sz, sy, sx)].sum())
                    n += int(sx.sum() * sy.sum() * sz.sum())
        if n == 0:
            raise SystemExit("no level-0 cells for p_inf")
        return tot / n


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("field")
    ap.add_argument("--case", required=True, help="the prepared case file (coef tiles)")
    ap.add_argument("--stl", required=True, help="the STL moby_prepare read (already rotated)")
    ap.add_argument("--scale", type=float, default=1.0, help="[ibm] stl_scale")
    ap.add_argument("--translate", type=float, nargs=3, default=[0.0, 0.0, 0.0], help="[ibm] stl_translate")
    ap.add_argument("--wing-y0", type=float, required=True, help="first strip starts here (fuselage side)")
    ap.add_argument("--wing-y1", type=float, required=True, help="last strip ends here (the cut)")
    ap.add_argument("--strips", type=float, default=0.1, help="strip width")
    ap.add_argument("--wing-xmax", type=float, default=np.inf, help="pieces beyond this x are not wing")
    ap.add_argument("--tri-size", type=float, default=0.01, help="longest edge of a piece")
    ap.add_argument("--cp-start", type=float, default=1.5, help="first Cp sample, in local cells")
    ap.add_argument("--cp-depth", type=int, default=4, help="number of Cp samples one cell apart")
    ap.add_argument("--cp-station", type=float, nargs="*", default=[], help="y of the Cp(x/c) stations")
    ap.add_argument("--station-width", type=float, default=None,
                    help="y window of a station (default tri_size; wider when the spanwise cell is coarser)")
    ap.add_argument("--alpha", type=float, default=0.0,
                    help="free-stream angle in the x-z plane, deg (lift/drag are taken in wind axes; "
                         "0 when the STL itself is pitched, as in the sailplane cases)")
    ap.add_argument("--pinf-box", type=float, nargs=6, default=None)
    ap.add_argument("--out", default="strips")
    a = ap.parse_args()

    F = Field(a.field, a.case)
    nu = 1.0 / F.re

    # ---- geometry -----------------------------------------------------------
    tri, nrm_stored = read_stl(a.stl)
    tri = tri * a.scale + np.asarray(a.translate)
    cr = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    live = np.linalg.norm(cr, axis=1) > 0                   # degenerate facets have no winding
    agree = (cr * nrm_stored).sum(1) > 0
    print(f"{a.stl}: {len(tri)} triangles ({(~live).sum()} degenerate), winding agrees with "
          f"the stored normal on {(agree & live).sum()} of {live.sum()}")
    if (~agree & live).sum() > 0:
        print(f"  WARNING: {(~agree & live).sum()} facets disagree -- the winding is taken as the truth")
    # only facets reaching into the domain (the mirrored half of a symmetric STL is dropped here)
    inbox = np.all((tri.max(1) > 0) & (tri.min(1) < F.L[None, :]), axis=1)
    piece = subdivide(tri[live & inbox], a.tri_size)
    cen = piece.mean(1)
    cr = np.cross(piece[:, 1] - piece[:, 0], piece[:, 2] - piece[:, 0])
    area = 0.5 * np.linalg.norm(cr, axis=1)
    keep = area > 0
    piece, cen, cr, area = piece[keep], cen[keep], cr[keep], area[keep]
    nrm = cr / (2 * area)[:, None]
    # chordwise tangent: the x direction projected onto the facet, TE-ward
    t = np.array([1.0, 0.0, 0.0])[None, :] - nrm[:, 0:1] * nrm
    tn = np.linalg.norm(t, axis=1)
    t = np.where(tn[:, None] > 1e-6, t / np.maximum(tn, 1e-300)[:, None], 0.0)
    print(f"  {len(piece)} pieces of edge <= {a.tri_size}, area {area.sum():.5g}")

    row0, lev0 = F.locate(cen)
    inside = row0 >= 0
    print(f"  pieces inside the domain with a leaf: {inside.sum()} of {len(piece)}")
    piece, cen, nrm, area, t = piece[inside], cen[inside], nrm[inside], area[inside], t[inside]
    row0, lev0 = row0[inside], lev0[inside]
    h = F.h(lev0)[:, 0]
    print(f"  surface levels: {dict(zip(*np.unique(lev0, return_counts=True)))}")

    # ---- Cp: linear wall extrapolation of trilinear samples along the normal
    p_inf = F.p_inf(a.pinf_box)
    nsmp = a.cp_depth
    S = np.arange(nsmp) + a.cp_start                             # in local cells
    P = np.full((nsmp, len(cen)), np.nan)
    first_row = None
    for k in range(nsmp):
        pts = cen + (S[k] * h)[:, None] * nrm
        r, l = F.locate(pts)
        if k == 0:
            first_row, first_pts = r, pts
        P[k] = F.sample("pn", pts, r, l)
    good = np.isfinite(P)
    n_ok = good.sum(0)
    # least squares p = p_w + b s over the valid samples of each piece
    s0 = good.sum(0); s1 = (good * S[:, None]).sum(0); s2 = (good * S[:, None] ** 2).sum(0)
    Pz = np.where(good, P, 0.0)
    t0 = Pz.sum(0); t1 = (Pz * S[:, None]).sum(0)
    det = s0 * s2 - s1 * s1
    pw = np.full(len(cen), np.nan)
    two = (s0 >= 2) & (det > 0)
    pw[two] = (s2[two] * t0[two] - s1[two] * t1[two]) / det[two]
    one = (s0 == 1)
    pw[one] = t0[one]
    cp = 2.0 * (pw - p_inf)
    print(f"  p_inf = {p_inf:.6e}; Cp from {np.isfinite(cp).sum()} pieces "
          f"(samples per piece: min {n_ok.min()}, median {int(np.median(n_ok))})")

    # ---- Cf: through-origin fit on the penalization band ---------------------
    half = cen + (0.5 * h)[:, None] * nrm
    rh, _ = F.locate(half)
    rows = np.unique(np.concatenate([row0, first_row[first_row >= 0], rh[rh >= 0]]))
    pos, val, comp, hb, brow, solid6 = F.band_points(rows)
    print(f"  band: {len(rows)} blocks, {len(pos)} staggered band faces")
    tree = cKDTree(cen)
    _, near = tree.query(pos)
    dvec = pos - cen[near]
    dn = (dvec * nrm[near]).sum(1)
    inplane = np.linalg.norm(dvec - dn[:, None] * nrm[near], axis=1)
    # a component whose tangent weight is roundoff (the spanwise one on a
    # 2D section, t_y ~ 1e-17) would fit its noise as a wall gradient
    m = (dn > 0) & (dn < 2 * hb) & (inplane < a.tri_size) & (np.abs(t[near, comp]) > 1e-3)
    g = dn[m] * t[near[m], comp[m]]
    num = np.bincount(near[m], weights=val[m] * g, minlength=len(cen))
    den = np.bincount(near[m], weights=g * g, minlength=len(cen))
    ncf = np.bincount(near[m], minlength=len(cen))
    cf = np.full(len(cen), np.nan)
    okf = den > (0.25 * h) ** 2                 # at least a quarter cell of leverage
    cf[okf] = 2.0 * nu * num[okf] / den[okf]
    print(f"  Cf from {okf.sum()} pieces (band points per piece: median {int(np.median(ncf[okf]))})")
    # pieces that received no band point (a spanwise cell coarser than the
    # pieces, the NACA's xz mode) take the Cf of the nearest fitted piece on
    # the same side, so the friction integral covers the whole surface
    cf_fit = cf.copy()
    if 0 < okf.sum() < len(cen):
        for side in (nrm[:, 2] > 0, nrm[:, 2] <= 0):
            src = okf & side; dst = ~okf & side
            if src.any() and dst.any():
                _, j = cKDTree(cen[src]).query(cen[dst])
                cf[dst] = cf[src][j]

    # sanity: the first Cp sample must sit in a fluid cell (a wrong normal
    # would read the solid interior)
    n_solid = 0
    for r in np.unique(first_row[first_row >= 0]):
        if r not in solid6:
            continue
        sel = first_row == r
        hh = F.h(F.blocks[r, 3])
        ijk = np.floor((first_pts[sel] - F.blocks[r, :3] * hh) / hh).astype(int).clip(0, F.nb - 1)
        n_solid += int(solid6[r][ijk[:, 2], ijk[:, 1], ijk[:, 0]].sum())
    print(f"  first Cp samples landing in a solid cell: {n_solid} (should be 0)")

    # ---- strips ---------------------------------------------------------------
    y = cen[:, 1]
    wing = (y >= a.wing_y0) & (y < a.wing_y1) & (cen[:, 0] < a.wing_xmax)
    edges = np.arange(a.wing_y0, a.wing_y1 - 1e-9, a.strips)
    edges = np.append(edges, a.wing_y1)
    strip = np.full(len(cen), -1)
    strip[wing] = np.clip(np.searchsorted(edges, y[wing], side="right") - 1, 0, len(edges) - 2)
    # pressure and friction force / q per piece, then in WIND axes
    fp = np.where(np.isfinite(cp), -cp, 0.0)[:, None] * nrm * area[:, None]
    ff = np.where(np.isfinite(cf), cf, 0.0)[:, None] * t * area[:, None]
    ca, sa = np.cos(np.radians(a.alpha)), np.sin(np.radians(a.alpha))
    e_lift, e_drag = np.array([-sa, 0.0, ca]), np.array([ca, 0.0, sa])
    fz = fp @ e_lift; fx = fp @ e_drag; ffx = ff @ e_drag; ffz = ff @ e_lift
    upper = nrm[:, 2] > 0
    hdr = ("strip  y0      y1      c_strip  A_plan   c_l      c_dp     c_df     c_m      "
           "x_sep1  x_sep2  f_sep   n_cp   n_cf")
    print(hdr)
    rows_out = []
    for s in range(len(edges) - 1):
        m = strip == s
        if not m.any():
            continue
        A = (np.maximum(nrm[m, 2], 0) * area[m]).sum()
        xle = cen[m, 0].min(); c_s = cen[m, 0].max() - xle
        zle = cen[m, 2][np.argmin(cen[m, 0])]
        xr, zr = xle + 0.25 * c_s, zle
        Fx = fx[m].sum() + ffx[m].sum(); Fz = fz[m].sum() + ffz[m].sum()
        My = ((cen[m, 2] - zr) * (fp[m, 0] + ff[m, 0]) - (cen[m, 0] - xr) * (fp[m, 2] + ff[m, 2])).sum()
        cl = Fz / A; cdp = fx[m].sum() / A; cdf = ffx[m].sum() / A; cm = My / (c_s * A)
        u = m & upper & np.isfinite(cf) & ((cen[:, 0] - xle) / c_s > 0.02)
        neg = u & (cf < 0)
        xs1 = (cen[neg, 0].min() - xle) / c_s if neg.any() else np.nan
        xs2 = (cen[neg, 0].max() - xle) / c_s if neg.any() else np.nan
        fsep = area[neg].sum() / max(area[u].sum(), 1e-300)
        rows_out.append((s, edges[s], edges[s + 1], c_s, A, cl, cdp, cdf, cm, xs1, xs2, fsep,
                         int(np.isfinite(cp[m]).sum()), int(np.isfinite(cf_fit[m]).sum())))
        print("%3d  %7.4f %7.4f  %7.4f  %7.4f  %+.4f  %+.5f  %+.5f  %+.4f  %6.3f  %6.3f  %5.3f  %6d %6d"
              % rows_out[-1])
    rest = ~wing
    if rest.any():
        print("rest (not wing): %d pieces, pressure force/q = (%+.5f, %+.5f, %+.5f), friction (%+.5f, %+.5f, %+.5f)"
              % ((rest.sum(),) + tuple(fp[rest].sum(0)) + tuple(ff[rest].sum(0))))
    with open(a.out + "_strips.csv", "w") as f:
        f.write("# " + hdr + "\n")
        for r in rows_out:
            f.write(",".join("%d" % v if isinstance(v, (int, np.integer)) else "%.6g" % v for v in r) + "\n")
    np.savez(a.out + "_pieces.npz", cen=cen, nrm=nrm, area=area, cp=cp, cf=cf, cf_fit=cf_fit,
             strip=strip, tangent=t, level=lev0, p_inf=p_inf, re=F.re, ncf=ncf, alpha=a.alpha)
    print(f"wrote {a.out}_strips.csv, {a.out}_pieces.npz")

    # ---- Cp(x/c) stations -----------------------------------------------------
    sw = a.station_width or a.tri_size
    for ys in a.cp_station:
        m = (np.abs(y - ys) < 0.5 * sw) & (cen[:, 0] < a.wing_xmax) & np.isfinite(cp)
        if not m.any():
            print(f"station y = {ys}: no pieces"); continue
        xle = cen[m, 0].min(); c_s = cen[m, 0].max() - xle
        xoc = (cen[m, 0] - xle) / c_s
        o = np.lexsort((xoc, ~upper[m]))
        path = f"{a.out}_cp_y{ys:g}.csv"
        with open(path, "w") as f:
            f.write(f"# y = {ys}, x_le = {xle:.5f}, chord = {c_s:.5f}; x/c, Cp, Cf, upper\n")
            for i in o:
                f.write("%.5f,%.6f,%s,%d\n" % (xoc[i], cp[m][i],
                                               "%.6e" % cf[m][i] if np.isfinite(cf[m][i]) else "nan",
                                               int(upper[m][i])))
        print(f"wrote {path}: {m.sum()} pieces, chord {c_s:.4f}, Cp min {cp[m].min():+.4f} "
              f"at x/c {xoc[np.argmin(cp[m])]:.4f}")


if __name__ == "__main__":
    main()
