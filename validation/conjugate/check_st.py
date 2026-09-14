#!/usr/bin/env python3
"""STAGE 0 of docs/next_session_tangential.md: does a ONE-SIDED s_t converge?

    ./check_st.py plane    --tag obl30 --theta 30 --kappa 10
    ./check_st.py cylinder --tag cyl   --radius 0.25 --kappa 10

The shipped estimator (check_oblique.debias) differences T ACROSS the
interface and then removes the straddle in closed form. That de-bias is
derived for a PLANE, which is why it does not converge on a curved one:
43 / 54 / 46 % of the signal at h = 1/64, 1/128, 1/256 (README table (b)).

The candidates never straddle -- see check_oblique.onesided. This scores all
of them against the EXACT s_t on both geometries and reports the observed
order, which is KILL GATE 1: a candidate that does not converge on the
cylinder cannot be the basis of anything downstream.

The plane is the unit test, not the gate: the manufactured field is piecewise
LINEAR and phi is exact, so every one-sided difference and the closure itself
are exact -- a one-sided estimator must return machine precision there, and
anything else is a bug rather than a truncation error.
"""

from __future__ import annotations

import argparse
import sys

import numpy as np

from check_oblique import (Plane, load_phi, face_measure, debias,   # noqa: E402
                           onesided, onesided_spread)
from check_cylinder import Dipole, Multipole                        # noqa: E402

CANDIDATES = [("raw projection  ", lambda m: (m["st_num"], None)),
              ("shipped de-bias ", lambda m: (debias(m), None)),
              ("1-sided flux  F ", lambda m: onesided(m, "flux", "fluid")),
              ("1-sided flux  S ", lambda m: onesided(m, "flux", "solid")),
              ("1-sided flux  FS", lambda m: onesided(m, "flux", "both")),
              ("1-sided flux kmx", lambda m: onesided(m, "flux", "kmax")),
              ("1-sided agree FS", lambda m: onesided(m, "agree", "both"))]


def score(case, model, exact_of):
    """rms error of every candidate against the exact s_t, over all cut faces."""
    phi, centres = load_phi(case)
    acc = {lab: [] for lab, _ in CANDIDATES}
    ex, spread, nfall, ntot = [], [], 0, 0
    for d in range(3):
        m = face_measure(phi, centres, model, d)
        if m is None:
            continue
        e = exact_of(model, phi, centres, m, d)
        ex.append(e)
        ntot += m["n"]
        for lab, fn in CANDIDATES:
            v, ok = fn(m)
            acc[lab].append(v - e)
            if ok is not None and lab.startswith("1-sided flux  S"):
                nfall += int((~ok).sum())
        spread.append(onesided_spread(m))
    escale = float(np.sqrt(np.mean(np.concatenate(ex) ** 2)))
    out = {lab: float(np.sqrt(np.mean(np.concatenate(v) ** 2)))
           for lab, v in acc.items()}
    return escale, out, float(np.sqrt(np.mean(np.concatenate(spread) ** 2))), nfall, ntot


def plane_exact(model, phi, centres, m, d):
    return np.full(m["n"], model.s_t(d))


def dipole_exact(model, phi, centres, m, d):
    """The exact s_t at each cut face's own polar angle (the nearest interface
    point of a face centre is its radial projection, so the face centre's
    angle is exact -- no O(h) error enters the reference)."""
    shape = phi.shape
    axis = 2 - d
    gk, gj, gi = np.meshgrid(*[np.arange(s) for s in shape], indexing="ij")
    xx, yy = centres[0][gi], centres[1][gj]
    hi = [slice(1, shape[q] - 1) for q in range(3)]
    lo = [slice(1, shape[q] - 1) for q in range(3)]
    lo[axis] = slice(0, shape[axis] - 2)
    hi, lo = tuple(hi), tuple(lo)
    cut = (phi[lo] < 0.0) != (phi[hi] < 0.0)
    th = model.theta(0.5 * (xx[hi] + xx[lo])[cut], 0.5 * (yy[hi] + yy[lo])[cut])
    return model.s_t_exact(th, d)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("geometry", choices=["plane", "cylinder", "multipole"])
    ap.add_argument("--tag", required=True, help="case-file stem, e.g. obl30 or cyl")
    ap.add_argument("--grids", type=int, nargs="+", default=[64, 128, 256])
    ap.add_argument("--kappa", type=float, required=True)
    ap.add_argument("--theta", type=float, default=30.0)
    ap.add_argument("--radius", type=float, default=0.25)
    ap.add_argument("--mode", type=int, default=2,
                    help="multipole only: m = 1 reproduces the dipole (LINEAR "
                         "interior, which flatters every solid-side score); "
                         "m >= 2 makes neither side linear")
    ap.add_argument("--q-n", type=float, default=1.0, dest="q_n")
    ap.add_argument("--amp", type=float, default=1.0,
                    help="plane only: the tangential amplitude. s_t = amp*t_d, "
                         "so amp = 0 leaves NOTHING to estimate")
    ap.add_argument("--tolerance", type=float, default=None,
                    help="plane only: the one-sided candidates must be below this")
    a = ap.parse_args()

    print(f"== s_t estimators on the {a.geometry}, kappa_s = {a.kappa:g}"
          f"   (rms error as % of the exact rms)")
    hdr = "   n     " + "".join(f"{lab:>18s}" for lab, _ in CANDIDATES)
    print(hdr)
    table = {lab: [] for lab, _ in CANDIDATES}
    hs = []
    for n in a.grids:
        case = f"{a.tag}_{n}.h5"
        if a.geometry == "plane":
            model, exact_of = Plane(a.theta, 0.5, 0.5117, a.kappa, a.q_n, a.amp), plane_exact
        elif a.geometry == "cylinder":
            model, exact_of = Dipole(0.5, 0.5, a.radius, a.kappa), dipole_exact
        else:
            model, exact_of = Multipole(0.5, 0.5, a.radius, a.kappa, a.mode), dipole_exact
        try:
            esc, out, spread, nfall, ntot = score(case, model, exact_of)
            if esc == 0.0:
                raise SystemExit("the exact s_t is identically zero -- "
                                 "raise --amp, or there is nothing to score")
        except (OSError, SystemExit) as exc:
            print(f"   {n:4d}   -- {case}: {exc}")
            continue
        hs.append(1.0 / n)
        row = f"   {n:4d}  "
        for lab, _ in CANDIDATES:
            table[lab].append(out[lab])
            row += f"{out[lab] / esc * 100.0:17.2f}%"
        print(row)
        print(f"          exact rms {esc:.4e}   f/s spread {spread:.3e}"
              f"   solid-side fallbacks {nfall}/{ntot} faces")
    if len(hs) < 2:
        return 1

    print("   observed order (least squares over the grids above):")
    ok = True
    for lab, _ in CANDIDATES:
        v = np.array(table[lab])
        p = np.polyfit(np.log(hs), np.log(np.maximum(v, 1e-300)), 1)[0]
        print(f"     {lab}  order = {p:6.2f}")
        if lab.startswith("1-sided") and a.geometry != "plane" and p <= 0.0:
            ok = False
    if a.tolerance is not None:
        worst = max(table[lab][-1] for lab, _ in CANDIDATES if lab.startswith("1-sided"))
        ok = worst <= a.tolerance
        print(f"   plane unit test: worst one-sided rms = {worst:.3e}"
              f"   {'PASS' if ok else 'FAIL'} (tolerance {a.tolerance:g})")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
