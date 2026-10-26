#!/usr/bin/env python3
"""Nonlinear lifting line of the FRUE V0 wing with AeroSandbox (step 1 of
docs/next_session_sailplane_junction.md, Route A).

The whole wing, NO fuselage (the DNS box supplies it), the measured planform
of section 1 of the handout, the real sections sliced from the STL by
stl_sections.py, viscous section data from NeuralFoil (XFOIL surrogate; the
XFOIL cross-check is section 3.1 of the handout).  Prints C_L(alpha), the
fuselage-axis alpha for the target C_L values, and at the chosen C_L the
strip table of the root region: c, c_l, alpha_i (two ways) and alpha_eff.
alpha_box = alpha_fuselage - alpha_i(root region) is the number the box takes.

    python3 liftingline.py SECTION_DIR [--cl 1.0] [--re 3e5] [--ncrit 9] [--polars DIR] [--csv OUT]

With --polars DIR the alpha_eff inversion uses the XFOIL polars written by
xfoil_polars.py (DIR/sec_y<Y>_ncrit<N>.csv) instead of NeuralFoil; the
lifting line itself (the loading, the downwash) stays on NeuralFoil, and
alpha_box from the polar route = alpha_eff(c_l root) - incidence does not
depend on the lifting line's own alpha at the target C_L.

SECTION_DIR holds sec_y{0.5,1.0,2.0,3.0,5.0,7.0,9.0,9.6}.dat from
`stl_sections.py section`.  Needs a venv with `pip install aerosandbox[full]`
(4.2.10 / neuralfoil 0.3.3 on 2026-10-26).
"""
import argparse
import numpy as np
import aerosandbox as asb

# (y, x_le, chord, z_le, section file key), metres, measured 2026-10-25;
# the root chord is carried to y = 0 (the fuselage is not modelled);
# the chord-line incidence is +2.3 deg at every station
PLANFORM = [
    (0.0, 2.630, 0.831, 0.445, "0.5"), (0.5, 2.630, 0.831, 0.445, "0.5"),
    (0.75, 2.619, 0.784, 0.459, "0.5"), (1.0, 2.607, 0.768, 0.473, "1.0"),
    (1.5, 2.584, 0.763, 0.501, "1.0"), (2.0, 2.560, 0.759, 0.529, "2.0"),
    (3.0, 2.514, 0.749, 0.586, "3.0"), (4.0, 2.470, 0.730, 0.642, "3.0"),
    (5.0, 2.433, 0.690, 0.696, "5.0"), (6.0, 2.395, 0.651, 0.749, "5.0"),
    (7.0, 2.362, 0.598, 0.802, "7.0"), (8.0, 2.337, 0.518, 0.855, "7.0"),
    (9.0, 2.322, 0.413, 0.908, "9.0"), (9.745, 2.321, 0.308, 0.950, "9.6"),
]
INCIDENCE = 2.3
C_ROOT = 0.831


def build(section_dir):
    keys = sorted({k for *_, k in PLANFORM}, key=float)
    afs = {k: asb.Airfoil(name="sec" + k, coordinates=np.loadtxt("%s/sec_y%s.dat" % (section_dir, k))).repanel(100)
           for k in keys}
    wing = asb.Wing(name="wing", symmetric=True,
                    xsecs=[asb.WingXSec(xyz_le=[x, y, z], chord=c, twist=INCIDENCE, airfoil=afs[k])
                           for (y, x, c, z, k) in PLANFORM])
    plane = asb.Airplane(name="FRUE V0 wing", xyz_ref=[2.84, 0, 0.47], wings=[wing],
                         s_ref=wing.area(), c_ref=C_ROOT, b_ref=wing.span())
    return plane, afs


def solve(plane, atm, U, alpha, res=6):
    ll = asb.LiftingLine(airplane=plane, op_point=asb.OperatingPoint(atmosphere=atm, velocity=U, alpha=alpha),
                         spanwise_resolution=res, model_size="xlarge")
    return ll, ll.run()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("section_dir")
    ap.add_argument("--cl", type=float, default=1.0, help="target C_L for the strip table")
    ap.add_argument("--re", type=float, default=3e5, help="Reynolds number on the root chord")
    ap.add_argument("--ncrit", type=float, default=9.0, help="NeuralFoil n_crit for the polar inversion")
    ap.add_argument("--polars", default=None, help="XFOIL polar csv directory for the alpha_eff inversion")
    ap.add_argument("--csv", default=None, help="write the strip table here")
    a = ap.parse_args()

    plane, afs = build(a.section_dir)
    atm = asb.Atmosphere(altitude=0)
    nu = atm.kinematic_viscosity()
    U = a.re * nu / C_ROOT
    wing = plane.wings[0]
    print("S %.2f m2  AR %.1f  Re_root %.2e  U %.2f m/s" % (wing.area(), wing.aspect_ratio(), a.re, U))

    alphas = np.arange(-3.0, 7.0, 1.0)
    CL = np.array([solve(plane, atm, U, al)[1]["CL"] for al in alphas])
    for al, cl in zip(alphas, CL):
        print("alpha_fuselage %+.0f  C_L %.4f" % (al, cl))
    imax = np.argmax(CL)
    for tgt in (0.6, 0.8, 0.9, 1.0, 1.1):
        if tgt <= CL[imax]:
            print("C_L %.1f -> alpha_fuselage %.2f deg" % (tgt, np.interp(tgt, CL[:imax + 1], alphas[:imax + 1])))
    if a.cl > CL[imax]:
        raise SystemExit("C_L %.2f is above the lifting line's maximum %.3f" % (a.cl, CL[imax]))
    a1 = float(np.interp(a.cl, CL[:imax + 1], alphas[:imax + 1]))

    ll, out = solve(plane, atm, U, a1, res=8)
    yc = np.array(ll.vortex_centers)[:, 1]
    ch = np.array(ll.chords).ravel()
    G = np.array(ll.vortex_strengths, dtype=float).ravel()
    cl = 2 * G / (U * ch)
    # alpha_i from the downwash of the whole vortex system at the bound vortex
    vi = np.array(ll.get_induced_velocity_at_points(np.array(ll.vortex_centers)))
    ai_w = np.degrees(np.arctan2(-vi[:, 2], U + vi[:, 0]))
    # alpha_i from the section polar: alpha_eff such that c_l(alpha_eff) = strip c_l
    aa = np.linspace(-6, 5, 111)
    rows = []
    for i in np.argsort(yc):
        if yc[i] < 0:
            continue
        k = min(afs, key=lambda k: abs(float(k) - max(yc[i], 0.5)))
        if a.polars:
            tab = np.loadtxt("%s/sec_y%s_ncrit%g.csv" % (a.polars, k, a.ncrit), delimiter=",")
            imx = np.argmax(tab[:, 1])
            al, pol = tab[:imx + 1, 0], tab[:imx + 1, 1]   # the rising branch
        else:
            al = aa
            pol = np.ravel(afs[k].get_aero_from_neuralfoil(alpha=aa, Re=U * ch[i] / nu, n_crit=a.ncrit,
                                                            model_size="xlarge")["CL"])
        aeff = float(np.interp(float(cl[i]), pol, al))
        rows.append((yc[i], ch[i], cl[i], ai_w[i], a1 + INCIDENCE - aeff, aeff))
    rows = np.array(rows)
    print("\nC_L %.2f: alpha_fuselage %.2f deg (chord line %.2f); right half, root region and every metre"
          % (a.cl, a1, a1 + INCIDENCE))
    print("   y[m]   c[m]    c_l   alpha_i(downwash)  alpha_i(polar)  alpha_eff(chord line)")
    for y, c, clv, aiw, aip, ae in rows:
        if y < 3.2 or abs(y - round(y)) < 0.06:
            print("  %5.2f  %.3f  %.3f   %+.2f              %+.2f           %+.2f" % (y, c, clv, aiw, aip, ae))
    root = rows[(rows[:, 0] > 0.55) & (rows[:, 0] < 1.05)]
    cut = rows[np.argmin(np.abs(rows[:, 0] - 2.95))]
    print("\nroot region y = 0.55-1.05 m: c_l %.3f, alpha_i %.2f (downwash) / %.2f (polar) deg"
          % (root[:, 2].mean(), root[:, 3].mean(), root[:, 4].mean()))
    print("loading drop from y = 0.6 to the cut at 2.95 m: %.1f %%" % (100 * (1 - cut[2] / root[0, 2])))
    print("alpha_box (fuselage axis): downwash route %.2f - %.2f = %.2f deg; polar route (%s, Ncrit %g) alpha_eff - %.1f = %.2f deg"
          % (a1, root[:, 3].mean(), a1 - root[:, 3].mean(), "XFOIL" if a.polars else "NeuralFoil", a.ncrit, INCIDENCE, root[:, 5].mean() - INCIDENCE))
    if a.csv:
        np.savetxt(a.csv, rows, fmt="%.4f", header="y c cl alpha_i_downwash alpha_i_polar alpha_eff (C_L %.2f, alpha_fus %.3f)" % (a.cl, a1))


if __name__ == "__main__":
    main()
