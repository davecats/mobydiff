#!/usr/bin/env python3
"""A1/A2 cylinder gate metrics (docs/next_session_airfoil.md).

  steady   forces_re40.txt          mean C_D over the converged tail vs 1.5-1.6, C_L -> 0
  strouhal forces_re100.txt         St from C_L zero crossings, mean C_D, symmetry
  empty    forces_empty.txt         C_L = C_D = 0.0 exactly
  cv       <cyl_re40_*.h5> --re 40 --cd-pen <val>
                                    Gauss/CV outer-box flux drag vs the
                                    penalization C_D (steady snapshot)
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "tools"))


def read_forces(path):
    d = np.loadtxt(path, skiprows=1)
    return d[:, 1], d[:, 2], d[:, 3]        # t, cl, cd


def cmd_steady(a):
    t, cl, cd = read_forces(a.forces)
    n = len(t)
    tail = slice(int(0.8*n), n)
    cdm, cds = float(np.mean(cd[tail])), float(np.std(cd[tail]))
    clm = float(np.mean(np.abs(cl[tail])))
    print(f"C_D (last 20%): {cdm:.4f} +- {cds:.2e}   |C_L|: {clm:.2e}")
    ok = 1.4 <= cdm <= 1.7 and cds < 1e-3 and clm < 5e-2
    print("steady gate:", "PASS" if ok else "FAIL", "(band 1.5-1.6 nominal, 1.4-1.7 hard)")
    return 0 if ok else 1


def cmd_strouhal(a):
    t, cl, cd = read_forces(a.forces)
    n = len(t)
    tail = slice(int(0.5*n), n)
    tt, cc = t[tail], cl[tail] - np.mean(cl[tail])
    # Fundamental of C_L: the confined/penalization C_L carries a 3rd
    # harmonic of comparable power (which also defeats zero-crossing
    # counting), so take the LOWEST spectral peak within 35% of the largest
    # and refine it by a local centroid.
    dt = float(np.mean(np.diff(tt)))
    cc_w = cc*np.hanning(len(cc))
    spec = np.abs(np.fft.rfft(cc_w))
    freq = np.fft.rfftfreq(len(cc_w), dt)
    pk = [i for i in range(1, len(spec) - 1)
          if spec[i] >= spec[i - 1] and spec[i] >= spec[i + 1]
          and spec[i] >= 0.35*np.max(spec[1:])]
    kpk = min(pk) if pk else int(np.argmax(spec[1:])) + 1
    sl = slice(max(1, kpk - 2), kpk + 3)
    st = float(np.sum(freq[sl]*spec[sl])/np.sum(spec[sl]))
    cdm = float(np.mean(cd[tail]))
    clm = float(np.mean(cl[tail]))
    amp = float(np.max(np.abs(cc)))
    print(f"St = {st:.4f} (spectral peak over {tt[-1]-tt[0]:.0f} time units)   "
          f"mean C_D = {cdm:.4f}   mean C_L = {clm:.3e}   C_L amplitude = {amp:.3f}")
    ok = 0.15 <= st <= 0.18 and 1.2 <= cdm <= 1.5 and abs(clm) < 0.1*amp
    print("strouhal gate:", "PASS" if ok else "FAIL", "(St 0.16-0.17 nominal)")
    return 0 if ok else 1


def cmd_empty(a):
    t, cl, cd = read_forces(a.forces)
    mcl, mcd = float(np.max(np.abs(cl))), float(np.max(np.abs(cd)))
    print(f"max|C_L| = {mcl:.3e}   max|C_D| = {mcd:.3e}")
    ok = mcl == 0.0 and mcd == 0.0
    print("empty gate:", "PASS (exact)" if ok else "FAIL")
    return 0 if ok else 1


def cmd_cv(a):
    import h5py
    from compare_fields import load_field
    with h5py.File(a.h5, "r") as f:
        u = load_field(f, "un").mean(axis=0)   # (ny, nx), z-averaged (quasi-2D)
        v = load_field(f, "vn").mean(axis=0)
        p = load_field(f, "pn").mean(axis=0)
        x, y, z = f["x"][...], f["y"][...], f["z"][...]
    lz = float(z[-1] - z[0])
    h = float(x[1] - x[0])
    nu = 1.0/a.re
    # collocate to cell centres (u face i sits between centres i-1, i)
    uc = 0.5*(u + np.roll(u, -1, axis=1))
    vc = 0.5*(v + np.roll(v, -1, axis=0))
    dudx = np.gradient(uc, h, axis=1)
    dudy = np.gradient(uc, h, axis=0)
    dvdx = np.gradient(vc, h, axis=1)
    # control box in cell-centre indices (well inside the domain, around the body)
    i0, i1 = int(a.x0/h), int(a.x1/h)
    j0, j1 = int(a.y0/h), int(a.y1/h)
    rows = slice(j0, j1 + 1)
    cols = slice(i0, i1 + 1)
    f_ew = -p - uc*uc + 2.0*nu*dudx                 # x-flux through x-normal faces
    f_ns = nu*(dudy + dvdx) - uc*vc                 # x-flux through y-normal faces
    fx = (np.sum(f_ew[rows, i1]) - np.sum(f_ew[rows, i0]))*h \
       + (np.sum(f_ns[j1, cols]) - np.sum(f_ns[j0, cols]))*h
    cd_cv = 2.0*fx*lz/(1.0**2*1.0*lz)
    print(f"CV box x=[{i0*h:.2f},{i1*h:.2f}] y=[{j0*h:.2f},{j1*h:.2f}]: C_D = {cd_cv:.4f}")
    if a.cd_pen is not None:
        rel = abs(cd_cv - a.cd_pen)/abs(a.cd_pen)
        print(f"penalization C_D = {a.cd_pen:.4f}   relative difference = {rel:.3f}")
        ok = rel < 0.08
        print("cv gate:", "PASS" if ok else "FAIL", "(<= 8% = discretization error)")
        return 0 if ok else 1
    return 0


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name, fn in (("steady", cmd_steady), ("strouhal", cmd_strouhal), ("empty", cmd_empty)):
        p = sub.add_parser(name)
        p.add_argument("forces")
        p.set_defaults(func=fn)
    p = sub.add_parser("cv")
    p.add_argument("h5")
    p.add_argument("--re", type=float, required=True)
    p.add_argument("--cd-pen", type=float, default=None)
    p.add_argument("--x0", type=float, default=3.0)
    p.add_argument("--x1", type=float, default=11.0)
    p.add_argument("--y0", type=float, default=4.0)
    p.add_argument("--y1", type=float, default=12.0)
    p.set_defaults(func=cmd_cv)
    a = ap.parse_args()
    sys.exit(a.func(a))


if __name__ == "__main__":
    main()
