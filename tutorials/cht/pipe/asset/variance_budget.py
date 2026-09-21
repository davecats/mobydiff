#!/usr/bin/env python3
"""The <theta'^2> budget, ours against the reference, term by term.

In a statistically steady pipe, homogeneous in z and (at bccode 0) in phi,

    0 = P - eps + T_t + T_nu,
    P    = -2 <u_r' t'> d<t>/dr                  (production)
    eps  = 2 alpha <|grad t'|^2>                 (dissipation)
    T_t  = -(1/r) d/dr [ r <u_r' t'^2> ]         (turbulent transport)
    T_nu = (alpha/r) d/dr [ r d<t'^2>/dr ]       (molecular transport)

WHY THIS TEST: the core variance deficit is not sampling, not the projection
and not under-production, and our dissipation per unit variance is only 3 %
high AT THE AXIS (it is 32 % high near the wall, where the variance matches).
That points at the term which feeds the core -- turbulent transport from the
wall region -- so this measures it directly. The reference stores the triple
products (t*t*u, t*t*v), which is what makes the comparison possible.

    ./variance_budget.py [snapshots...] --case pipe_fine.h5
"""
from __future__ import annotations
import argparse, sys
import h5py, numpy as np
import os
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))      # pipe_stats.py, beside the inis
from pipe_stats import load_field, load_solid_mask
from compare_neuhauser import ref_budget

NU, PR = 1/5300., 0.71
AL = NU/PR
DELTA_V = 2.759318e-3


def reference(nbins=None, edges=None):
    """The reference budget, from `asset/neuhauser_profiles.npz`.

    The reduction from the published (time, r, phi) moments lives in
    extract_neuhauser.py -- including the SOURCE term, which is easy to miss:
    the Kasagi heating is proportional to u_z, so it fluctuates and feeds the
    variance directly at 2 S <u_z' theta'>. Their S follows from making
    q_w = 1: integrating the source over the cross-section against the flux
    through the perimeter gives S = 2 q_w/(R u_b).
    """
    R = dict(ref_budget())
    R["src"] = 2.0 * (2.0 / (0.5 * float(np.mean(R["uz"][R["r"] < 0.5])))) \
        * R["tuz"]
    return R


def ours(snaps, case, dr):
    solid = load_solid_mask(case)
    with h5py.File(case, "r") as h:
        node = {k: h[f"{k}_nodes"][...] for k in "xyz"}
    cx = 0.5*(node["x"][:-1] + node["x"][1:]); cy = 0.5*(node["y"][:-1] + node["y"][1:])
    dx = float(node["x"][1] - node["x"][0]); dz = float(node["z"][1] - node["z"][0])
    # (y, x) order, to match the FIELD layout (k, j, i). r is symmetric under
    # x<->y and survives a transpose, but cos(phi) is not: built the other way
    # round it silently becomes sin(phi) and the radial flux comes out 5x too
    # small. Same trap the case-file readers document.
    Y, X = np.meshgrid(cy, cx, indexing="ij")
    rp = np.hypot(X - 0.65, Y - 0.65)
    cs, sn = (X - 0.65)/np.maximum(rp, 1e-30), (Y - 0.65)/np.maximum(rp, 1e-30)
    nbin = int(np.ceil(rp.max()/dr))
    ib = np.minimum((rp/dr).astype(np.int32), nbin - 1)
    acc = {k: np.zeros(nbin) for k in
           ("n", "t", "tt", "ur", "urt", "urtt", "g2", "uz", "uzt")}
    for sp in snaps:
        with h5py.File(sp, "r") as h:
            t = load_field(h, "c0"); u = load_field(h, "un"); v = load_field(h, "vn")
            w_ = load_field(h, "wn")
        uc = 0.5*(u + np.roll(u, -1, 2)); vc = 0.5*(v + np.roll(v, -1, 1))
        wc = 0.5*(w_ + np.roll(w_, -1, 0))
        del u, v, w_
        gx = (np.roll(t, -1, 2) - np.roll(t, 1, 2))/(2*dx)
        gy = (np.roll(t, -1, 1) - np.roll(t, 1, 1))/(2*dx)
        gz = (np.roll(t, -1, 0) - np.roll(t, 1, 0))/(2*dz)
        g2 = gx*gx + gy*gy + gz*gz
        del gx, gy, gz
        fl = ~solid
        ur = uc*cs[None] + vc*sn[None]
        del uc, vc
        idx = np.broadcast_to(ib[None], t.shape)[fl]
        w = lambda q: np.bincount(idx, weights=q[fl], minlength=nbin)
        acc["n"] += np.bincount(idx, minlength=nbin)
        acc["t"] += w(t); acc["tt"] += w(t*t); acc["ur"] += w(ur)
        acc["urt"] += w(ur*t); acc["urtt"] += w(ur*t*t); acc["g2"] += w(g2)
        acc["uz"] += w(wc); acc["uzt"] += w(wc*t)
        del wc
        print(f"   {sp}: binned", flush=True)
        del t, ur, g2
    n = np.maximum(acc["n"], 1)
    m = {k: acc[k]/n for k in acc if k != "n"}
    r = (np.arange(nbin) + 0.5)*dr
    var = m["tt"] - m["t"]**2
    flux = m["urt"] - m["ur"]*m["t"]
    trip = m["urtt"] - 2*m["t"]*m["urt"] - m["ur"]*m["tt"] + 2*m["ur"]*m["t"]**2
    good = acc["n"] > 0
    P = np.where(good, -2*flux*np.gradient(m["t"], r), 0.0)
    eps = 2*AL*(m["g2"] - (np.gradient(m["t"], r))**2)
    Tt = -np.gradient(r*trip, r)/np.maximum(r, 1e-12)
    Tn = AL*np.gradient(r*np.gradient(var, r), r)/np.maximum(r, 1e-12)
    # the fluctuating part of the Kasagi source, at our own source = 1.0
    src = 2.0 * (m["uzt"] - m["uz"]*m["t"])
    return dict(r=r, P=P, eps=eps, Tt=Tt, Tn=Tn, var=var, flux=flux, trip=trip,
                src=src, good=good)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("snaps", nargs="+")
    ap.add_argument("--case", default="../pipe_prod.h5")
    ap.add_argument("--dr", type=float, default=0.01)
    ap.add_argument("--qw", type=float, required=True,
                    help="our wall flux from the interface-heat balance")
    a = ap.parse_args()
    O = ours(a.snaps, a.case, a.dr)
    R = reference(None, None)
    # WALL UNITS on each side, from its own theta_tau = q_w/u_tau: the two runs
    # normalise their heating differently (theirs makes q_w = 1), so the raw
    # budget terms differ by (theta_tau ratio)^2 ~ 17 and cannot be compared.
    U_TAU = 0.068379
    fo = NU/(U_TAU**2*(a.qw/U_TAU)**2)
    fr = NU/(U_TAU**2*(1.00330/U_TAU)**2)
    for k in ("P", "eps", "Tt", "Tn", "src"):
        O[k] = O[k]*fo; R[k] = R[k]*fr
    O["var"] = O["var"]/(a.qw/U_TAU)**2
    R["var"] = R["var"]/(1.00330/U_TAU)**2
    print("\n<theta'^2> budget in WALL UNITS (each side by its own theta_tau)")
    print(f"{'y+':>7}{'':3}{'P':>10}{'eps':>10}{'T_turb':>10}{'T_mol':>10}"
          f"{'source':>10}{'sum':>10}{'variance':>11}")
    for target in (0.02, 0.10, 0.20, 0.30, 0.40, 0.45):
        i = np.argmin(np.abs(O["r"] - target)); j = np.argmin(np.abs(R["r"] - target))
        yp = (0.5 - O["r"][i])/DELTA_V
        for lab, D, k in (("ours", O, i), ("ref ", R, j)):
            print(f"{yp:7.1f} {lab}{D['P'][k]:10.3f}{-D['eps'][k]:10.3f}"
                  f"{D['Tt'][k]:10.3f}{D['Tn'][k]:10.3f}{D['src'][k]:10.3f}"
                  f"{D['P'][k]-D['eps'][k]+D['Tt'][k]+D['Tn'][k]+D['src'][k]:10.3f}"
                  f"{D['var'][k]:11.3f}")
        print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
