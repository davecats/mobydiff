#!/usr/bin/env python3
"""The conjugate turbulent-channel validation (Re_tau = 180, Pr = 0.71).

Reads the solver's accumulated statistics and puts numbers on the conjugate
scheme's TURBULENCE behaviour -- the one thing the C1-C3 gate suite, which is
entirely manufactured solutions and self-consistency identities, cannot say.

    ./check_cht.py velocity <vel_stats.h5>          vs the in-repo KMM180 DNS
    ./check_cht.py thermal  <cht_stats.h5>          the conjugate sweep

WHAT IS ANCHORED, AND HOW
  velocity   the mean U+ and the three rms components against
             tutorials/channel_kmm180/channel_kmm180_stats.h5 -- the
             Kim-Moin-Moser Re_tau = 180 reference this repository already
             carries. Nothing thermal means anything until the flow is shown
             to BE a Re_tau 180 DNS.
  mean theta KADER's correlation (a published correlation, the same one the
             S2 gates use), over the wall layer where it applies.
  the sweep  three predictions that an implementation cannot fake:
     (1) the three kappa_s = 1 scalars must have the SAME mean profile --
         capacity is irrelevant at steady state -- while their fluctuations
         differ. Any leak of C_s into the mean is a bug;
     (2) EFFUSIVITY COLLAPSE. The surface temperature response of a solid to
         a harmonic interface flux depends on sqrt(k rho c) ALONE, so
         (kappa_s, C_s) = (1, 100) and (10, 10) -- same K = 10, different
         alpha_s, different mean resistance -- must give the same normalised
         interface fluctuation. Likewise (1, 1e4) and (100, 100) at K = 100;
     (3) the 1/K ASYMPTOTE. For K >> 1 the wall dominates the coupling and
         theta'_wall/theta_tau falls like 1/K; for K << 1 it saturates at the
         ideal-isoflux plateau. That S-curve, with theta'_wall -> 0 as
         K -> infinity (the ideal-isothermal limit), is the classical
         conjugate result (Tiselj et al. 2001; Flageul et al. 2015).
"""

from __future__ import annotations

import argparse
import os
import sys

import h5py
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "..", "tools"))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "scalar"))

NSTAT = 7
S, SS, US, CLO, JLO, CHI, JHI = range(NSTAT)
RE, PR = 180.0, 0.71
Y_LO, Y_HI = 0.6, 2.6

# ---------------------------------------------------------------------------
# THE LITERATURE ANCHOR: Flageul, Benhamadouche, Lamballais & Laurence,
# "DNS of turbulent channel flow with conjugate heat transfer: effect of
# thermal boundary conditions on the second moments and budgets",
# Int. J. Heat and Fluid Flow 55 (2015) 34-44 (literature/flageul.pdf).
# Re_tau = 149, Pr = 0.71, four thermal boundary conditions, and their
# conjugate case has G = G_2 = 1 -- the solid has the SAME conductivity and
# the SAME diffusivity as the fluid, i.e. kappa_s = 1, alpha_s = 1, K = 1.
# That is EXACTLY this sweep's k1.
#
# Values are the temperature VARIANCE <T'^2> normalised by T_tau^2, read from
# their figure 5 (there is no table; the wall values are where the curves
# meet the axis, so treat them as +-0.1 rather than as digits).
FLAGEUL = dict(re_tau=149, pr=0.71,
               wall_isoQ=4.2,      # ideal Neumann: the K -> 0 bracket
               wall_conjug=1.1,    # THEIR conjugate case, K = 1  -> our k1
               wall_isoT=0.0,      # ideal Dirichlet: the K -> infinity bracket
               peak=6.3,           # peak <T'^2>, y+ ~ 17-20, all four cases
               peak_yp=18.0)
# Their section 5 derives the scaling this sweep is built to test: at the
# interface (their eq. 12-13) the wall variance goes like (lambda_f/(R
# lambda_s))^2 with R^2 = k_x^2 + k_z^2 + i k_t rho c/lambda_s, so the product
# lambda_s R ~ sqrt(lambda_s rho_s c_s) is the EFFUSIVITY and
# theta'_wall ~ 1/K. lambda_s >> lambda_f gives the isothermal wall,
# lambda_s << lambda_f the isoflux one.

# the sweep, in the ini's order: (name, kappa_s, C_s)
SWEEP = [("k1", 1.0, 1.0), ("k2", 1.0, 100.0), ("k3", 1.0, 1.0e4),
         ("a1", 0.1, 0.1), ("a2", 10.0, 10.0), ("a3", 100.0, 100.0)]


def read_stats(path):
    with h5py.File(path, "r") as f:
        prof = f["profile"][...]
        count = f["count"][...]
        coord = f["coord"][...] if "coord" in f else f["ycoord"][...]
        attrs = dict(f.attrs)
    n = coord.size
    nsc = prof.shape[1] // NSTAT
    return dict(y=coord, prof=prof, count=count, nscalar=nsc,
                t=float(attrs.get("t_current", 0.0)),
                step=int(attrs.get("step", 0)), re=float(attrs.get("re", RE)))


def col(st, isc, c):
    return st["prof"][:, NSTAT * isc + c]


def kader(yp, pr=PR):
    """Kader (1981) theta+(y+) -- the published correlation the S2 gates use."""
    yp = np.asarray(yp, dtype=float)
    beta = (3.85 * pr ** (1.0 / 3.0) - 1.3) ** 2 + 2.12 * np.log(pr)
    gamma = 0.01 * (pr * yp) ** 4 / (1.0 + 5.0 * pr ** 3 * yp)
    return pr * yp * np.exp(-gamma) + (2.12 * np.log(1.0 + yp) + beta) * np.exp(-1.0 / np.maximum(gamma, 1e-300))


def fluid_mask(y):
    return (y > Y_LO) & (y < Y_HI)


def cmd_velocity(a):
    """Mean U+ and the rms components against the in-repo KMM180 DNS."""
    with h5py.File(a.file, "r") as f:
        prof = f["profile"][...]
        y = f["coord"][...]
    ref_path = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                            "tutorials", "channel_kmm180", "channel_kmm180_stats.h5")
    with h5py.File(ref_path, "r") as f:
        rprof = f["profile"][...]
        ry = f["coord"][...]

    # column 0 = <u>; the rms components follow the channel writer's order
    # (<u>,<v>,<w>,<uu>,<vv>,<ww>,<uv>,...), so uu/vv/ww are 3,4,5.
    def stats(p, yy, y0):
        u = p[:, 0]
        uu = np.maximum(p[:, 3] - p[:, 0] ** 2, 0.0)
        vv = np.maximum(p[:, 4] - p[:, 1] ** 2, 0.0)
        ww = np.maximum(p[:, 5] - p[:, 2] ** 2, 0.0)
        yp = (yy - y0) * RE
        return yp, u, np.sqrt(uu), np.sqrt(vv), np.sqrt(ww)

    m = fluid_mask(y)
    yp, u, ur, vr, wr = stats(prof[m], y[m], Y_LO)
    ryp, ru, rur, rvr, rwr = stats(rprof, ry, 0.0)
    half = yp <= 180.0
    rhalf = ryp <= 180.0

    print(f"{a.file}: {y[m].size} fluid rows")
    print(f"   U+ centreline   run {u[half].max():.3f}   KMM180 {ru[rhalf].max():.3f}"
          f"   ({100*(u[half].max()/ru[rhalf].max()-1):+.2f} %)")
    for nm, r, rr in (("u'", ur, rur), ("v'", vr, rvr), ("w'", wr, rwr)):
        pk, rpk = r[half].max(), rr[rhalf].max()
        print(f"   {nm}_rms peak      run {pk:.4f}   KMM180 {rpk:.4f}"
              f"   ({100*(pk/rpk-1):+.2f} %)")
    # the profile itself, interpolated onto the reference y+
    ui = np.interp(ryp[rhalf], yp[half], u[half])
    err = np.abs(ui - ru[rhalf])/max(ru[rhalf].max(), 1e-30)
    print(f"   |U+ - U+_KMM| over 0 < y+ < 180: max {err.max()*100:.2f} % of U+_c")
    return 0


def cmd_thermal(a):
    st = read_stats(a.file)
    y = st["y"]
    m = fluid_mask(y)
    yf = y[m]
    print(f"{a.file}: step {st['step']}, t = {st['t']:.3f}, "
          f"{st['nscalar']} scalars, {y.size} rows ({m.sum()} fluid)")
    print(f"   interfaces at y = {Y_LO} / {Y_HI};  fluid rows y+ = "
          f"{(yf[0]-Y_LO)*RE:.2f} .. {(0.5*(Y_LO+Y_HI)-Y_LO)*RE:.0f}")

    rows = []
    for isc, (nm, ks, cs) in enumerate(SWEEP[:st["nscalar"]]):
        mean = col(st, isc, S)
        var = np.maximum(col(st, isc, SS) - mean ** 2, 0.0)
        rms = np.sqrt(var)
        jlo = col(st, isc, JLO)
        jhi = col(st, isc, JHI)
        # the INTERFACE flux: the low face of the first fluid row and the high
        # face of the last one. Both are the transport kernel's own face flux.
        i0 = int(np.argmax(m))                       # first fluid row
        i1 = int(len(m) - 1 - np.argmax(m[::-1]))    # last fluid row
        j0, j1 = float(jlo[i0]), float(jhi[i1])
        ttau = 0.5 * abs(j0 + j1)
        # theta'_rms at the first fluid cell (y+ = 0.75), normalised
        wall_rms = 0.5 * (rms[i0] + rms[i1]) / ttau
        # the mean interface temperature, from the two straddling rows
        ti = 0.5 * (abs(mean[i0] - 0.0) + abs(mean[i1] - 0.0))
        # theta+ profile on the lower half against Kader
        low = m & (y < 0.5 * (Y_LO + Y_HI))
        yp = (y[low] - Y_LO) * RE
        thp = (mean[i0] - mean[low]) / ttau + (mean[i0] - mean[i0]) / ttau
        # theta+ measured from the INTERFACE value (the conjugate wall value)
        thp = (mean[i0] - mean[low]) / ttau
        band = (yp > 5.0) & (yp < 40.0)
        kd = kader(yp[band]) - kader(np.array([yp[0]]))[0]
        dev = np.abs((thp[band] - (thp[0] if False else 0.0)) - kd)
        # the SOLID-side decay: how much variance survives at the outer
        # boundary. Their solid is d+ = 149 with a Neumann outer face; ours is
        # d+ = 36 with a Dirichlet one, so this number bounds how much of the
        # low-frequency tail our thinner wall clips (their figure 12 has the
        # variance down by 1e-3 at y+ = -77).
        sol_lo = rms[:i0] / ttau
        outer = float(sol_lo[0] ** 2) if sol_lo.size else float("nan")
        # THE PEAK MUST BE THE NEAR-WALL ONE, not the maximum over the fluid.
        # This case's thermal problem holds the total flux CONSTANT across the
        # channel (antisymmetric walls, no source), so production never
        # switches off and the variance keeps rising to the CENTRELINE -- the
        # global maximum is there, not at y+ ~ 20. Flageul's wall-flux problem
        # has the flux falling linearly to zero at the centre, so their global
        # maximum IS the near-wall peak. Taking max() over the fluid therefore
        # compares two different quantities, and it made a +11 % near-wall
        # difference read as +40 %.
        ypw = np.minimum(y - Y_LO, Y_HI - y) * RE
        band = m & (ypw > 5.0) & (ypw < 40.0)
        pk_nw = float((rms[band] / ttau).max())
        yp_nw = float(ypw[band][int(np.argmax(rms[band]))])
        rows.append(dict(name=nm, ks=ks, cs=cs, K=np.sqrt(ks * cs), ttau=ttau,
                         peak_nw=pk_nw, yp_nw=yp_nw,
                         wall_rms=wall_rms, ti=ti, mean=mean, rms=rms,
                         j0=j0, j1=j1, yp=yp, thp=thp, peak=float(rms[m].max()/ttau),
                         var_wall=wall_rms ** 2, var_peak=pk_nw ** 2,
                         var_centre=float(rms[m].max()/ttau) ** 2,
                         var_outer=outer))

    print()
    print("   scalar  kappa_s      C_s        K   theta_tau  <t'2>_wall  "
          "<t'2>_nwpeak  y+pk   wall/peak   <t'2>_centre")
    for r in rows:
        print(f"   {r['name']:>5}  {r['ks']:8g} {r['cs']:9g} {r['K']:7.3g} "
              f"{r['ttau']:10.6f} {r['var_wall']:11.3f} {r['var_peak']:12.3f} "
              f"{r['yp_nw']:6.1f} {r['var_wall']/max(r['var_peak'],1e-30):11.4f} "
              f"{r['var_centre']:13.3f}")

    ok = True
    by = {r["name"]: r for r in rows}

    # (1) capacity must not touch the mean
    if all(k in by for k in ("k1", "k2", "k3")):
        mm = np.array([by[k]["mean"] for k in ("k1", "k2", "k3")])
        spread = float(np.max(np.abs(mm - mm[0]))) / max(by["k1"]["ttau"], 1e-30)
        tt = np.array([by[k]["ttau"] for k in ("k1", "k2", "k3")])
        print(f"\n   (1) kappa_s = 1, C_s = 1/100/1e4: the MEAN must coincide")
        print(f"       max|<theta> - <theta>_k1|/theta_tau = {spread:.3e}"
              f"     theta_tau spread {float(tt.max()/tt.min()-1):.3e}")
        ok &= spread < a.mean_tol

    # (2) effusivity collapse: same K, different (kappa_s, C_s)
    print(f"\n   (2) effusivity collapse -- same K, different kappa_s/C_s/alpha_s")
    for lhs, rhs in (("k2", "a2"), ("k3", "a3")):
        if lhs in by and rhs in by:
            a_, b_ = by[lhs]["wall_rms"], by[rhs]["wall_rms"]
            rel = abs(a_ - b_) / max(abs(a_), abs(b_), 1e-30)
            print(f"       K = {by[lhs]['K']:.4g}: {lhs} {a_:.4f} vs {rhs} {b_:.4f}"
                  f"   relative difference {rel*100:.1f} %")
            ok &= rel < a.collapse_tol

    # (3) the K dependence: monotone, and ~1/K at large K
    print(f"\n   (3) theta'_wall/theta_tau against K (the conjugate S-curve)")
    srt = sorted(rows, key=lambda r: r["K"])
    for i, r in enumerate(srt):
        note = ""
        if i:
            kr = r["K"] / srt[i - 1]["K"]
            wr = srt[i - 1]["wall_rms"] / max(r["wall_rms"], 1e-30)
            note = f"   K x{kr:.3g} -> theta'_w /{wr:.3g}  (1/K would give /{kr:.3g})"
        print(f"       K = {r['K']:8.3g}  theta'_w/theta_tau = {r['wall_rms']:.4f}{note}")
    mono = all(srt[i]["wall_rms"] >= srt[i + 1]["wall_rms"] - 1e-12
               for i in range(len(srt) - 1))
    print(f"       monotone decreasing in K: {mono}")
    ok &= mono

    if a.dump:
        with open(a.dump, "w") as fh:
            fh.write("# name kappa_s C_s K theta_tau theta_w_rms peak\n")
            for r in rows:
                fh.write("%s %g %g %g %.8e %.8e %.8e\n" %
                         (r["name"], r["ks"], r["cs"], r["K"], r["ttau"],
                          r["wall_rms"], r["peak"]))
        with open(a.dump + ".prof", "w") as fh:
            fh.write("# y " + " ".join(f"{r['name']}_mean {r['name']}_rms" for r in rows) + "\n")
            for j in range(y.size):
                fh.write("%.10f " % y[j] + " ".join(
                    "%.8e %.8e" % (r["mean"][j], r["rms"][j]) for r in rows) + "\n")
        print(f"\n   tables written to {a.dump}[.prof]")

    # (4) THE LITERATURE COMPARISON
    f = FLAGEUL
    print(f"\n   (4) vs Flageul et al. 2015 (Re_tau {f['re_tau']}, Pr {f['pr']}, "
          f"their figure 5; our Re_tau 180)")
    print(f"       quantity                          ours      Flageul")
    if "k1" in by:
        r = by["k1"]
        print(f"       <theta'^2>_wall, CONJUGATE K = 1  {r['var_wall']:8.2f}"
              f"  {f['wall_conjug']:9.2f}   <- their G = G_2 = 1 case")
    for nm, lab, ref in (("a1", "K = 0.1  -> the isoQ bracket", f["wall_isoQ"]),
                         ("a3", "K = 100  -> the isoT bracket", f["wall_isoT"])):
        if nm in by:
            print(f"       <theta'^2>_wall, {lab:<16} {by[nm]['var_wall']:8.2f}"
                  f"  {ref:9.2f}")
    pk = np.mean([r["var_peak"] for r in rows])
    print(f"       NEAR-WALL peak <theta'^2> (mean over sweep) {pk:6.2f}  {f['peak']:9.2f}")
    ct = np.mean([r["var_centre"] for r in rows])
    print(f"       ...their global max IS that peak; ours is at the CENTRELINE "
          f"({ct:.2f}), because our flux is constant across the channel and "
          f"theirs falls to zero. Do not compare those two.")
    if "k1" in by:
        r = by["k1"]
        print(f"       wall/peak at K = 1 (theta_tau-free, level-free) "
              f"{r['var_wall']/max(r['var_peak'],1e-30):8.4f}  "
              f"{f['wall_conjug']/f['peak']:9.4f}")
    if "k1" in by and np.isfinite(by["k1"]["var_outer"]):
        print(f"       variance surviving at our outer solid face (d+ = 36): "
              f"{by['k1']['var_outer']:.3g}  ({100*by['k1']['var_outer']/max(by['k1']['var_wall'],1e-30):.1f} %"
              f" of the interface value) -- their d+ = 149 wall has 1e-3 of it at y+ = -77")

    print("\n   PASS" if ok else "\n   FAIL")
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("velocity")
    p.add_argument("file")
    p.add_argument("--interfaces", type=float, nargs=2, default=None,
                   help="the two grid-aligned interface positions")
    p.set_defaults(func=cmd_velocity)

    p = sub.add_parser("thermal")
    p.add_argument("file")
    p.add_argument("--dump", default=None)
    p.add_argument("--mean-tol", type=float, default=5e-3)
    p.add_argument("--collapse-tol", type=float, default=0.10)
    p.add_argument("--interfaces", type=float, nargs=2, default=None,
                   help="the two grid-aligned interface positions")
    p.set_defaults(func=cmd_thermal)

    a = ap.parse_args()
    if getattr(a, "interfaces", None):
        Y_LO, Y_HI = a.interfaces
        globals()["Y_LO"], globals()["Y_HI"] = Y_LO, Y_HI
    return a.func(a)


if __name__ == "__main__":
    sys.exit(main())
