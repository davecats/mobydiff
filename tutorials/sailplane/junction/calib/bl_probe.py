#!/usr/bin/env python3
"""Wall-normal profiles and the upstream flow angle at spanwise stations of
a block-table snapshot, for the question of 2026-11-04 (why the c/2048
calibration box carries half the XFOIL lift: incidence, Reynolds number,
or a boundary layer that is not turbulent before it separates).

For every station y and chordwise position x/c the upper-surface piece of
the strip-load run nearest (x_le + x/c * chord, y) gives the wall point
and the outward normal; u, v, w, p, k, nut are sampled trilinearly at
distances s along the normal (in local cells, from 0.5 to --smax), and the
tangential velocity is the projection on the chordwise tangent of that
piece. Upstream of the leading edge the flow angle atan(w/u) is sampled
at --upstream chords on the station's leading-edge height, and at the
inlet column. One csv per station (`<out>_y<Y>.csv`, s in cells, distance
in chords) and one for the angles (`<out>_inflow.csv`).

    bl_probe.py FIELD.h5 --case CASE.h5 --pieces PIECES.npz --station 0.722 1.203 2.9
                [--xc 0.05 0.1 ...] [--smax 64] [--upstream 0.5 1 2 4 6] [--out probe]
"""
import argparse, os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from strip_loads import Field

ap = argparse.ArgumentParser()
ap.add_argument("field"); ap.add_argument("--case", required=True); ap.add_argument("--pieces", required=True)
ap.add_argument("--station", type=float, nargs="+", required=True)
ap.add_argument("--xc", type=float, nargs="+", default=[0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.8])
ap.add_argument("--smax", type=float, default=64.0, help="farthest sample, in local cells")
ap.add_argument("--upstream", type=float, nargs="+", default=[0.5, 1.0, 2.0, 4.0, 6.0])
ap.add_argument("--wing-xmax", type=float, default=8.5)
ap.add_argument("--out", default="probe")
a = ap.parse_args()

F = Field(a.field, a.case)
P = np.load(a.pieces)
cen, nrm, tan, lev = P["cen"], P["nrm"], P["tangent"], P["level"]
S = np.unique(np.concatenate([np.geomspace(0.5, a.smax, 40), [1, 2, 3, 4, 6, 8, 12, 16, 24, 32]]))
FIELDS = ["un", "vn", "wn", "pn", "k", "nut"]

inflow = []
for y in a.station:
    slab = (np.abs(cen[:, 1] - y) < 0.02) & (cen[:, 0] < a.wing_xmax)
    if not slab.any():
        print(f"station y = {y}: no pieces"); continue
    x_le = cen[slab, 0].min(); x_te = cen[slab, 0].max(); chord = x_te - x_le
    z_le = cen[slab][np.argmin(cen[slab, 0]), 2]
    upper = slab & (nrm[:, 2] > 0.2)
    rows = []
    for xc in a.xc:
        target = np.array([x_le + xc * chord, y, 0.0])
        cand = np.nonzero(upper)[0]
        d = np.hypot(cen[cand, 0] - target[0], cen[cand, 1] - target[1])
        i = cand[np.argmin(d)]
        h = F.h(int(lev[i]))[0]
        t = tan[i] - (tan[i] @ nrm[i]) * nrm[i]; t /= np.linalg.norm(t)
        if t[0] < 0: t = -t                                   # chordwise, pointing downstream
        pts = cen[i][None, :] + (S * h)[:, None] * nrm[i][None, :]
        row, lv = F.locate(pts)
        q = {n: F.sample(n, pts, row, lv) for n in FIELDS}
        ut = q["un"] * t[0] + q["vn"] * t[1] + q["wn"] * t[2]
        un_ = q["un"] * nrm[i][0] + q["vn"] * nrm[i][1] + q["wn"] * nrm[i][2]
        for k in range(len(S)):
            rows.append((xc, S[k], S[k] * h / chord, ut[k], un_[k], q["pn"][k], q["k"][k], q["nut"][k], q["nut"][k] * F.re, int(lv[k])))
        # the edge velocity (max |u_t| over the profile) and the 99 % thickness
        ue = np.nanmax(np.abs(ut)); k99 = np.argmax(np.abs(ut) >= 0.99 * ue)
        print(f"y {y:5.3f} x/c {xc:4.2f}: h {h:.2e} ({h/chord:.1e} c), u_e {ue:.3f}, delta99 {S[k99]*h/chord:.4f} c, "
              f"k_max {np.nanmax(q['k']):.2e}, nut/nu max {np.nanmax(q['nut'])*F.re:6.1f} at s {S[np.nanargmax(q['nut'])]:.1f}, "
              f"u_t(0.5 cell) {ut[0]:+.3f}")
    hdr = (f"station y = {y}, x_le = {x_le:.5f}, chord = {chord:.5f}, z_le = {z_le:.5f}; "
           "x/c, s(cells), s/c, u_t, u_n, p, k, nut, nut/nu, level")
    np.savetxt(f"{a.out}_y{y}.csv", np.array(rows), delimiter=",", header=hdr,
               fmt="%.3f,%.4f,%.6e,%+.6e,%+.6e,%+.6e,%.6e,%.6e,%.4f,%d")
    # upstream flow angle on the leading-edge height, and at the inlet
    for d in a.upstream + [None]:
        p = np.array([[0.5 if d is None else x_le - d * chord, y, z_le]])
        if p[0, 0] <= 0: continue
        row, lv = F.locate(p)
        u, v, w = (F.sample(n, p, row, lv)[0] for n in ("un", "vn", "wn"))
        ang = np.degrees(np.arctan2(w, u))
        inflow.append((y, -1.0 if d is None else d, p[0, 0], u, v, w, ang))
        print(f"y {y:5.3f} upstream {'inlet' if d is None else f'{d:3.1f} c'} (x {p[0,0]:.3f}, z {z_le:.3f}): u {u:+.4f} w {w:+.4f} angle {ang:+.3f} deg")
np.savetxt(f"{a.out}_inflow.csv", np.array(inflow), delimiter=",",
           header="y, distance upstream (chords; -1 = inlet), x, u, v, w, angle atan(w/u) deg", fmt="%.4f")
print(f"wrote {a.out}_y*.csv, {a.out}_inflow.csv")
