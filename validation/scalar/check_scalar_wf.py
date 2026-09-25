#!/usr/bin/env python3
"""Gates for increment S5a -- the passive scalars' THERMAL WALL FUNCTION.

Under `[rans] wall_treatment = wall_function` the first cell sits in the log
layer, so the scalar's wall flux is set by the Kader/Jayatilleke closure the
solver installs as a wall-cell eddy diffusivity (scalar.f90's
`wall_diffusivity` on `rans_wall_yplus`'s k-based y+), not by the resolved
gradient. These checks read theta_tau and the wall flux STRAIGHT OUT OF THE
S4 STATISTICS (`[scalar] stats_*`, the same numbers `tools/scalar_stats.py
profile` prints) and cross them with the snapshot the same run wrote.

Subcommands
  wall     <stats.h5> <snapshot.h5>
           the wall-cell gates of ONE case: y+_1 from the snapshot's own k,
           the DELIVERED-FLUX identity (the solver's discrete wall flux must
           BE u_tau* (theta_1 - theta_w)/theta+_wf -- an identity, so it is
           gated tightly), theta+_1 against Kader, and the implied centreline
           theta+ = 1/theta_tau.
  compare  <ref_stats.h5> <wf_stats.h5> ...
           the wall-function cases against the RESOLVED reference: theta_tau
           and the centreline theta+, the T3 shape of the validation
           (wf180_y30 gated against turb180 through the DNS anchor).

The wall-function correlations are transcribed INDEPENDENTLY here (they are
five lines) so the comparison is not the solver checking itself.
"""

from __future__ import annotations

import argparse
import sys

import h5py
import numpy as np

sys.path.insert(0, ".")
from check_scalar_turb import LevelGrid, kader          # noqa: E402

NSTAT = 7
S, SS, US, CLO, JLO, CHI, JHI = range(7)

KAPPA = 0.41
ELOG = 9.8
CMU25 = 0.09 ** 0.25


# ------------------------------------------------------- the correlations ---
def jayatilleke_p(prat: float) -> float:
    """Jayatilleke's sublayer resistance P(Pr/Pr_t)."""
    return 9.24 * (prat ** 0.75 - 1.0) * (1.0 + 0.28 * np.exp(-0.007 * prat))


def thermal_yplus(pr: float, prt: float, p: float) -> float:
    """Where Pr y+ meets Pr_t[ln(E y+)/kappa + P], beyond the minimum of the
    difference at y* = Pr_t/(kappa Pr) (the physical crossing)."""
    f = lambda y: pr * y - prt * (np.log(ELOG * y) / KAPPA + p)
    lo = prt / (KAPPA * pr)
    if f(lo) >= 0.0:
        return lo
    hi = lo
    while f(hi) <= 0.0:
        hi *= 2.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) < 0.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def theta_plus_wf(yplus, pr: float, prt: float):
    """theta+(y+): conduction below the thermal sublayer thickness, the
    log branch above it. Elementwise -- the identity below is checked cell
    by cell, because the solver averages the FLUX over the wall plane and
    theta+(y+) is nonlinear in the cell's own k."""
    yplus = np.asarray(yplus, dtype=float)
    p = jayatilleke_p(pr / prt)
    ypt = thermal_yplus(pr, prt, p)
    with np.errstate(divide="ignore", invalid="ignore"):
        log = prt * (np.log(ELOG * np.maximum(yplus, 1e-300)) / KAPPA + p)
    return np.where(yplus < ypt, pr * yplus, log)


# ------------------------------------------------------------ the readers ---
def read_stats(path: str) -> dict:
    with h5py.File(path, "r") as f:
        return dict(nstat=int(f.attrs["nstat"]), step=int(f.attrs["step"]),
                    t=float(f.attrs["t_current"]), re=float(f.attrs["re"]),
                    y=f["coord"][...], profile=f["profile"][...],
                    count=f["count"][...])


def col(st: dict, index: int, which: int) -> np.ndarray:
    return st["profile"][:, NSTAT * (index - 1) + which]


def wall_row(path: str, scalar: str) -> dict:
    """A single-level snapshot: the full (z,y,x) blocks of theta and k, plus
    the y geometry. The wall gates need the CELLS, not the profile: the
    solver's statistics average the flux, and the flux is nonlinear in k."""
    with h5py.File(path, "r") as f:
        g = LevelGrid(f, 0)
        th = g.assemble(f[scalar])
        kk = g.assemble(f["k"])
        nut = g.assemble(f["nut"])
        un = g.assemble(f["un"])
        t = float(f.attrs["t_current"])
        ly = float(f.attrs["ly"])
    return dict(t=t, ly=ly, y=g.y, ynode=g.ynode, theta3=th, k3=kk,
                theta=th.mean(axis=(0, 2)), k=kk.mean(axis=(0, 2)),
                nut=nut.mean(axis=(0, 2)), u=un.mean(axis=(0, 2)))


# ------------------------------------------------------------ subcommands ---
def cmd_wall(a) -> int:
    st = read_stats(a.stats)
    sn = wall_row(a.snapshot, a.scalar)
    re, pr, prt = st["re"], a.pr, a.prt
    nu = 1.0 / re
    flux_lo = col(st, a.index, JLO)
    flux_hi = col(st, a.index, JHI)
    jlo, jhi = float(flux_lo[0]), float(flux_hi[-1])
    theta_tau = 0.5 * abs(jlo + jhi)

    print(f"{a.stats}: step {st['step']}, t = {st['t']:.4f}, {st['y'].size} rows, "
          f"Re = {re:g}, scalar '{a.scalar}' (Pr {pr}, Pr_t {prt})")
    print(f"wall flux from the S4 statistics: y_min {jlo:+.6e}  y_max {jhi:+.6e}"
          f"  -> theta_tau = {theta_tau:.6f}")

    status = 0
    p = jayatilleke_p(pr / prt)
    ypt = thermal_yplus(pr, prt, p)
    print(f"thermal wall function: P = {p:+.4f}, y+_T = {ypt:.3f} "
          f"(momentum y+_lam = 11.530)")

    # Both walls: cell 0 against y_min, cell -1 against y_max.
    worst_flux = 0.0
    worst_kader = 0.0
    for side, j, jwall, thw in ((" y_min", 0, jlo, a.walls[0]),
                                (" y_max", -1, jhi, a.walls[1])):
        y1 = sn["y"][j] if j == 0 else sn["ly"] - sn["y"][j]
        kcell = sn["k3"][:, j, :]
        thcell = sn["theta3"][:, j, :]
        k1 = float(kcell.mean())
        th1 = float(thcell.mean())
        utau_c = CMU25 * np.sqrt(np.maximum(kcell, 0.0))
        ypcell = utau_c * y1 / nu
        tpcell = theta_plus_wf(ypcell, pr, prt)
        utau_s = CMU25 * np.sqrt(max(k1, 0.0))
        yplus = float(utau_s * y1 / nu)
        tp_wf = float(theta_plus_wf(np.array([yplus]), pr, prt)[0])
        # The VISCOUS y+ = y u_tau/nu, from the delivered wall stress
        # (nu + nut_1) U_1/y_1 -- T3's u_tau. Kader is a function of THIS
        # y+, not of the k-based one the closure switches on: on a fine
        # grid the sublayer is resolved, k_1 is small and y+_k << y+.
        utau = np.sqrt((nu + sn["nut"][j]) * abs(sn["u"][j]) / y1)
        ypv = float(utau * y1 / nu)
        # The identity the wall function is built to satisfy: the DISCRETE
        # wall flux the statistics report must be u_tau* (theta_1 - theta_w)
        # / theta+, because the wall-cell diffusivity is nu y+/theta+ and the
        # mirrored Dirichlet ghost turns the face flux into
        # D (theta_1 - theta_w)/y_1.
        # J is measured in the +y sense on both walls, so the low wall's
        # flux points INTO the domain and the high wall's out of it:
        #   J_lo = u_tau* (theta_w - theta_1)/theta+
        #   J_hi = u_tau* (theta_N - theta_w)/theta+
        jpred = float(np.mean(utau_c * (thw - thcell) / tpcell)) if j == 0 \
            else float(np.mean(utau_c * (thcell - thw) / tpcell))
        jmeas = jwall
        rel = abs(jmeas - jpred) / max(abs(jmeas), 1e-300)
        tp_meas = abs(thw - th1) / theta_tau
        kd = float(kader(np.array([ypv]), pr)[0])
        kdev = abs(tp_meas - kd) / kd
        worst_flux = max(worst_flux, rel)
        worst_kader = max(worst_kader, kdev)
        print(f"{side}: y_1 = {y1:.5f}, u_tau = {float(utau):.4f} -> y+_1 = "
              f"{ypv:7.3f};  k_1 = {k1:.5e} -> y+_k = {yplus:7.3f}"
              f" ({'log' if yplus >= ypt else 'conduction'} branch)")
        print(f"        theta_1 = {th1:+.6f}, theta+_1 = {tp_meas:7.4f}"
              f"   Kader(y+_1) {kd:7.4f} (dev {kdev*100:.1f}%)"
              f"   closure theta+(y+_k) {tp_wf:7.4f}")
        print(f"        delivered flux {jmeas:+.6e} vs the wall function's "
              f"{jpred:+.6e}   rel {rel:.3e}")

    # Total-flux constancy: steady state makes J(y) exactly constant.
    dev = float(np.max(np.abs(flux_lo - jlo))) / theta_tau
    print(f"total flux J(y): max|J - J_wall|/theta_tau = {dev:.3e}")
    # The centreline anchor (theta_centre = 0 by antisymmetry).
    ic = int(np.argmin(np.abs(st["y"] - 0.5 * sn["ly"])))
    thc = float(col(st, a.index, S)[ic])
    print(f"centreline: <theta> = {thc:+.6e}, implied theta+_c = "
          f"{abs(a.walls[0] - thc)/theta_tau:.4f}  (1/theta_tau = "
          f"{1.0/theta_tau:.4f})")

    if worst_flux > a.flux_tol:
        print(f"FAIL: delivered flux differs from the wall function by "
              f"{worst_flux:.3e} > {a.flux_tol:g}")
        status = 1
    if a.kader_tol is not None and worst_kader > a.kader_tol:
        print(f"FAIL: theta+_1 off Kader by {worst_kader*100:.1f}% "
              f"> {a.kader_tol*100:g}%")
        status = 1
    print("PASS" if status == 0 else "FAIL")
    return status


def cmd_compare(a) -> int:
    ref = read_stats(a.reference)
    jlo = float(col(ref, a.index, JLO)[0])
    jhi = float(col(ref, a.index, JHI)[-1])
    tt_ref = 0.5 * abs(jlo + jhi)
    print(f"reference {a.reference}: theta_tau = {tt_ref:.6f}, "
          f"theta+_c = {1.0/tt_ref:.4f}   (resolved walls)")
    head = f"{'case':>28} {'theta_tau':>12} {'theta+_c':>10} {'dev':>9}"
    print(head + ("   s2/s1" if a.pair else ""))
    status = 0
    for path in a.cases:
        st = read_stats(path)
        j0 = float(col(st, a.index, JLO)[0])
        j1 = float(col(st, a.index, JHI)[-1])
        tt = 0.5 * abs(j0 + j1)
        dev = (tt - tt_ref) / tt_ref
        # theta_tau and theta+_c are reciprocal here (antisymmetric walls),
        # so ONE tolerance covers both statements.
        flag = "" if abs(dev) <= a.tolerance else "  <-- outside tolerance"
        if abs(dev) > a.tolerance:
            status = 1
        extra = ""
        if a.pair:
            # theta_tau of scalar 2 over scalar 1 IN THE SAME FILE: under
            # wall functions the wall cell uses the constant Pr_t for both,
            # so the Kays-Crawford scalar can differ only through the
            # INTERIOR -- where Pe_t is large and the correlation tends to
            # its Prt_inf asymptote. The resolved reference, whose wall
            # cells are where Kays-Crawford bites hardest, separates the
            # pair much further.
            t1 = 0.5 * abs(float(col(st, 1, JLO)[0]) + float(col(st, 1, JHI)[-1]))
            t2 = 0.5 * abs(float(col(st, 2, JLO)[0]) + float(col(st, 2, JHI)[-1]))
            extra = f"   {t2/t1:6.4f}"
        print(f"{path:>28} {tt:12.6f} {1.0/tt:10.4f} {dev*100:+8.2f}%{extra}{flag}")
    print(f"tolerance {a.tolerance*100:g}% on theta_tau; "
          + ("PASS" if status == 0 else "FAIL"))
    return status


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("wall")
    p.add_argument("stats")
    p.add_argument("snapshot")
    p.add_argument("--scalar", default="theta")
    p.add_argument("--index", type=int, default=1)
    p.add_argument("--pr", type=float, default=0.71)
    p.add_argument("--prt", type=float, default=0.85)
    p.add_argument("--walls", type=float, nargs=2, default=[1.0, -1.0])
    p.add_argument("--flux-tol", type=float, default=1e-10,
                   help="the delivered-flux identity is exact up to the "
                        "plane averaging of a 1D field")
    p.add_argument("--kader-tol", type=float, default=None)
    p.set_defaults(func=cmd_wall)

    p = sub.add_parser("compare")
    p.add_argument("reference")
    p.add_argument("cases", nargs="+")
    p.add_argument("--index", type=int, default=1)
    p.add_argument("--tolerance", type=float, default=0.10)
    p.add_argument("--pair", action="store_true",
                   help="also print theta_tau(scalar 2)/theta_tau(scalar 1)")
    p.set_defaults(func=cmd_compare)

    a = ap.parse_args()
    return a.func(a)


if __name__ == "__main__":
    sys.exit(main())
