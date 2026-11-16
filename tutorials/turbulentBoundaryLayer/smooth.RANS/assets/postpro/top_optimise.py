#!/usr/bin/env python3
"""Find the top-boundary condition that gives a boundary layer of constant
Clauser parameter beta = delta*/tau_w dP_e/dx downstream of a ramp.

    python3 top_optimise.py CASE.json

THE CONTROL IS THE TOP PRESSURE p_top(x); the top velocity V_top(x) is its
result, recorded every iteration. The top stays an outlet (Dirichlet
pressure, predicted normal velocity, the total-pressure datum where fluid
enters) with `[boundary] y_max_p_profile = table`. Prescribing V_top itself
(a generic patch with v Dirichlet, p Neumann) was tried first and fails with
this solver: the only Dirichlet pressure left is the x = 650 outlet, the
iterative projection cannot move the domain-long pressure mode that carries
the top flux to it, and the field keeps a 10 % mass imbalance and
L2(div) ~ 3e-3 whatever niter (6, 30, 100 tested, 2026-11-14). With a
Dirichlet top the same change converges like the ZPG case.

CASE.json (see ../../../apg.RANS/opt/case.json):
    beta        target beta (> 0 adverse, < 0 favourable)
    x_ramp      [x0, x1]: the target beta ramps from 0 at x0 to beta at x1 (half
                cosine); only its value from x_fit[0] on is enforced
    x_fit       [xa, xb]: the window where beta must match (convergence measure)
    p_ramp      [xs, xe]: the top pressure gradient ramps in over this range
    npoly       number of polynomial terms of dp_top/dx downstream of the ramp
    bump        optional [xs, xb]: one more basis function, S over [xs, xb]
                minus S over p_ramp -- a steeper start of the same ramp. The
                layer lags the gradient, and with the ramp ending at x_fit[0]
                the polynomial alone left beta(200) at 0.43 of a 0.5 target
                (the max error stalled at 0.073, 2026-11-15); its coefficient
                comes last in coeffs.txt (files without it are padded with 0)
    tol         max |beta - beta_target| over x_fit to stop
    t_chunk     each iteration marches in chunks of t_chunk time units ...
    steady_tol  ... until max |dc_f|/c_f over one chunk (x 30-640) < steady_tol
    final_tol   the converged iterate is marched on to this (stricter) tolerance
    t_max       ... or t_max time units in all (then flagged NOT steady)
    max_iter    iteration cap
    stall       stop when max |beta - target| fell by less than this fraction
                over the last 3 iterates (default 0.05): the basis, not the
                iteration, is then the limit (the APG run sat at 0.073 for 8
                iterations before the bump was added)
    relax       under-relaxation of the coefficient update
    zpg_field   the converged zero-pressure-gradient field (start)
    template    the ini every iteration copies (@PREFIX@ @RESTART@ @TFINAL@)
    solver      moby_solve binary;  workdir: where it0, it1, ... are created

The top pressure is p_top(x) = -int_x^650 g, g = dp_top/dx = S(x) sum_m a_m
xi^m [+ a_bump (S_bump(x) - S(x))] (S the half-cosine ramp over p_ramp,
xi = (x - xs)/(650 - xs)), so
p_top = 0 at x = 650 where the top meets the x outlet (held at 0). Each
iteration k (directory it<k>/):
  1. evaluate the previous field with pg_bl.pg_integrals (Griffin edge on
     the local free-stream total pressure);
  2. the model of the beta response to a coefficient change da is
        d beta = kappa (delta*/tau_w) R da,
     R_m the linear inviscid response of the WALL pressure gradient to the
     top pressure of basis m (Laplace in the 650 x 100 box: p' = p_top at the
     top, 0 at the outlet, the uniform p_top(0) at the inlet -- u' = 0 there,
     so p' + U u' is that constant -- and dp'/dy = 0 at the wall), with
     delta*/tau_w of the current field. That alone (kappa = 1) is the
     frozen-layer model, and it under-predicts: the layer answers an adverse
     gradient by thickening and losing friction, which raises beta further
     (the first APG step aimed at 0.5 reached 1.2 at steady state). kappa is
     the secant gain measured on the last step with the model of the CURRENT
     state, <d beta_actual, d beta_model> / <d beta_model, d beta_model> over
     [xa, xb], clamped to [0.3, 5] (with the previous state's model the
     thickening is counted twice: kappa 3.0 and a step 6x too small);
  3. da = the least-squares solution of kappa (delta*/tau_w) R da =
     beta_target - beta on [xa, 640]; a low-order g keeps the inversion of
     the smoothing kernel well conditioned;
  4. a += relax da; it<k>/ptop.txt is written and the solver restarts from the
     previous field, chunk by chunk, until the layer is steady (a pressure
     step travels down the plate at the near-wall speed: 500 t.u. left
     x > 400 still adapting, and an update computed from such a field
     chases the transient); the top-face V_top of the result goes to
     it<k>/vtop_result.txt.
The boundary-layer feedback (delta*, tau_w respond to the pressure gradient)
is what the iteration resolves; R only has to be right to first order.
"""
import glob
import json
import os
import subprocess
import sys

import h5py
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.signal import savgol_filter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rans_bl as rb  # noqa: E402
import pg_bl  # noqa: E402

LX, LY = 650.0, 100.0


def top_plane_v(path):
    """Top-face v (the vn_ymax plane of an outlet top, z-averaged) and the x
    cell centres; NaN when the file has no such plane."""
    with h5py.File(path, "r") as f:
        xn = f["x"][...]
        nbx, nbz = int(f.attrs["block_nb_x"]), int(f.attrs["block_nb_z"])
        v = np.full(int(f.attrs["nx"]), np.nan)
        if "vn_ymax" in f:
            P, bl = f["vn_ymax"][...], f["blocks"][...]
            for b, (ox, oy, oz, lev) in enumerate(bl):
                v[ox:ox + nbx] = P[b].reshape(nbz, nbx).mean(axis=0)
    return 0.5 * (xn[1:] + xn[:-1]), v


class PressureResponse:
    """Wall pressure p'(x, 0) of the perturbation potential flow in the box
    for a top pressure p_top(x): Laplace, second order on ng nodes,
    Dirichlet at the top (p_top), the outlet (0) and the inlet (p_top(0)),
    Neumann (ghost node) at the wall."""

    def __init__(self, ng=(261, 41)):
        nx, ny = ng
        self.nx, self.ny = nx, ny
        hx, hy = LX / (nx - 1), LY / (ny - 1)
        self.xg = np.linspace(0.0, LX, nx)
        N = nx * ny
        A = sp.lil_matrix((N, N))
        self.dir_rows = []
        for j in range(ny):
            for i in range(nx):
                k = j * nx + i
                if j == ny - 1 or i == 0 or i == nx - 1:
                    A[k, k] = 1.0                  # Dirichlet node
                    continue
                for di, dj, h in ((1, 0, hx), (-1, 0, hx), (0, 1, hy), (0, -1, hy)):
                    jj = j + dj
                    jj = 1 if jj < 0 else jj        # wall ghost: mirror (dp/dy = 0)
                    A[k, jj * nx + i + di] += 1.0 / h ** 2
                    A[k, k] -= 1.0 / h ** 2
        self.lu = spla.splu(A.tocsc())

    def __call__(self, ptop):
        nx, ny = self.nx, self.ny
        rhs = np.zeros(nx * ny)
        rhs[(ny - 1) * nx:] = ptop                       # top
        rhs[np.arange(ny) * nx] = ptop[0]                # inlet: the uniform p_top(0)
        rhs[np.arange(ny) * nx + nx - 1] = 0.0           # outlet
        rhs[(ny - 1) * nx + nx - 1] = ptop[-1]
        return self.lu.solve(rhs).reshape(ny, nx)[0]


def ramp(x, x0, x1):
    s = np.clip((x - x0) / (x1 - x0), 0.0, 1.0)
    return 0.5 * (1.0 - np.cos(np.pi * s))


def top_pressure(xq, a, p_ramp, bump=None):
    """p_top(xq) = -int_xq^LX g for the coefficients a (fine quadrature); with
    bump, a[-1] is the coefficient of the steeper-start function."""
    xf = np.linspace(0.0, LX, 6501)
    xi = np.clip((xf - p_ramp[0]) / (LX - p_ramp[0]), 0.0, None)
    S = ramp(xf, *p_ramp)
    poly = a[:-1] if bump else a
    g = S * sum(am * xi ** m for m, am in enumerate(poly))
    if bump:
        g = g + a[-1] * (ramp(xf, *bump) - S)
    cum = np.r_[0.0, np.cumsum(0.5 * (g[1:] + g[:-1]) * np.diff(xf))]
    return np.interp(xq, xf, cum - cum[-1])


def last_field(d, prefix):
    fs = [f for f in glob.glob(os.path.join(d, prefix + "_*.h5"))
          if f.rsplit("_", 1)[1][:-3].isdigit()]
    return max(fs, key=lambda s: int(s.rsplit("_", 1)[1][:-3])) if fs else None


def run_steady(cfg, base, d, prefix, prev, t, log):
    """March the case in directory d from field prev (time t) in chunks of
    t_chunk until c_f changes by less than steady_tol over one chunk; keep
    only the final chunk's field. Returns its path."""
    tmpl = open(os.path.join(base, cfg["template"])).read()
    cf_old, keep, c = None, [], 0
    while True:
        name = f"{prefix}c{c}"
        t += cfg["t_chunk"]
        s = tmpl.replace("@PREFIX@", name).replace("@RESTART@", os.path.relpath(prev, d)) \
                .replace("@TFINAL@", f"{t:.1f}")
        open(os.path.join(d, name + ".ini"), "w").write(s)
        with open(os.path.join(d, name + ".log"), "w") as lg:
            r = subprocess.run(["mpirun", "-n", "1", cfg["solver"], name + ".ini"], cwd=d,
                               stdout=lg, stderr=subprocess.STDOUT)
        if r.returncode != 0:
            sys.exit(f"{name}: solver exit {r.returncode}, see {d}/{name}.log")
        prev = last_field(d, name)
        keep.append(prev)
        q = pg_bl.pg_integrals(rb.load(prev))
        sel = (q["x"] > 30) & (q["x"] < 640)
        if cf_old is not None:
            ch = np.abs(q["cf"] / cf_old - 1.0)[sel].max()
            line = f"   {name}: t {t:.0f}, max |dc_f|/c_f over the chunk {100 * ch:.3f} %"
            print(line, flush=True)
            log.write(line + "\n")
            log.flush()
            if ch < cfg["steady_tol"] or t - keep_t0 >= cfg["t_max"]:
                if ch >= cfg["steady_tol"]:
                    log.write("   NOT steady at t_max\n")
                break
        else:
            keep_t0 = t - cfg["t_chunk"]
        cf_old = q["cf"]
        c += 1
    for fn in keep[:-1]:
        os.remove(fn)
    return prev


def load_coeffs(path, n):
    """coeffs.txt, padded with zeros to n (a bump added to a running case)."""
    a = np.atleast_1d(np.loadtxt(path))
    return np.r_[a, np.zeros(n - len(a))]


def evaluate(path):
    f = rb.load(path)
    return f, pg_bl.pg_integrals(f)


def main():
    cfg = json.load(open(sys.argv[1]))
    base = os.path.dirname(os.path.abspath(sys.argv[1]))
    wd = os.path.join(base, cfg.get("workdir", "."))
    beta_t, (x0, x1), (xa, xb) = cfg["beta"], cfg["x_ramp"], cfg["x_fit"]
    zpg = os.path.join(base, cfg["zpg_field"])
    p_ramp, bump = cfg["p_ramp"], cfg.get("bump")
    npoly = cfg.get("npoly", 4) + (1 if bump else 0)          # number of coefficients
    x = rb.load(zpg)["x"]
    pr = PressureResponse()
    R = []
    for m in range(npoly):
        e = np.zeros(npoly)
        e[m] = 1.0
        pw = pr(top_pressure(pr.xg, e, p_ramp, bump))
        R.append(np.interp(x, pr.xg, np.gradient(pw, pr.xg)))
    R = np.array(R).T                                         # (nx, npoly)
    log = open(os.path.join(wd, "optimise.log"), "a")

    k = 0
    while os.path.isdir(os.path.join(wd, f"it{k}")):
        k += 1
    if k == 0:
        prev, a = zpg, np.zeros(npoly)
    else:
        prev = last_field(os.path.join(wd, f"it{k-1}"), f"it{k-1}*")
        a = load_coeffs(os.path.join(wd, f"it{k-1}", "coeffs.txt"), npoly)
    def model(qq):
        """d beta / d a of the frozen-layer model at field qq: (nx, npoly)."""
        return (savgol_filter(qq["dstar"], 41, 2) / savgol_filter(qq["tauw"], 41, 2))[:, None] * R

    # the previous iterate (for the secant gain): it<k-2>'s field, or the ZPG
    hist = []
    if k >= 1:
        hist.append((zpg if k == 1 else last_field(os.path.join(wd, f"it{k-2}"), f"it{k-2}*"),
                     np.zeros(npoly) if k == 1 else load_coeffs(os.path.join(wd, f"it{k-2}", "coeffs.txt"), npoly)))
    kappa = 1.0
    win = (x >= xa) & (x <= xb)
    fit = (x >= xa) & (x <= 640.0)
    bt = beta_t * ramp(x, x0, x1)
    errs = []
    for it in range(k, cfg["max_iter"]):
        f, q = evaluate(prev)
        err = np.abs(q["beta"] - bt)[win].max()
        errs.append(err)
        if hist:
            pf, pa = hist[-1]
            pq = pg_bl.pg_integrals(rb.load(pf))
            db_act = (q["beta"] - pq["beta"])[win]
            db_mod = (model(q) @ (a - pa))[win]     # the model of THIS state: a true secant
            kappa = float(np.clip(np.dot(db_act, db_mod) / np.dot(db_mod, db_mod), 0.3, 5.0))
        msg = (f"it{it}: field {os.path.relpath(prev, wd)} t {f['t']:.0f}  max|beta-target| on "
               f"[{xa:.0f},{xb:.0f}] = {err:.3f}  beta mean {q['beta'][win].mean():+.3f}  "
               f"min cf {q['cf'][(x > 20) & (x < 640)].min():.5f}  kappa {kappa:.2f}  "
               f"a = {np.array2string(a, precision=6)}")
        print(msg, flush=True)
        log.write(msg + "\n")
        log.flush()
        if err < cfg["tol"] and it > 0:
            print("converged", flush=True)
            log.write("converged; final march to final_tol\n")
            cfg_f = dict(cfg, steady_tol=cfg.get("final_tol", cfg["steady_tol"]))
            d = os.path.dirname(prev)
            final = run_steady(cfg_f, base, d, os.path.basename(d) + "f", prev, f["t"], log)
            f, q = evaluate(final)
            log.write(f"final: {os.path.relpath(final, wd)} t {f['t']:.0f}  max|beta-target| "
                      f"{np.abs(q['beta'] - bt)[win].max():.3f}  beta mean {q['beta'][win].mean():+.3f}\n")
            break
        if len(errs) >= 3 and errs[-1] > (1.0 - cfg.get("stall", 0.05)) * errs[-3]:
            msg = (f"stalled: max|beta-target| {errs[-3]:.3f} -> {errs[-1]:.3f} over two iterations; "
                   "the basis cannot reach tol (add terms or a bump)")
            print(msg, flush=True)
            log.write(msg + "\n")
            break
        da = np.linalg.lstsq(kappa * model(q)[fit], (bt - q["beta"])[fit], rcond=None)[0]
        hist.append((prev, a.copy()))
        a = a + cfg["relax"] * da
        d = os.path.join(wd, f"it{it}")
        os.makedirs(d, exist_ok=True)
        np.savetxt(os.path.join(d, "coeffs.txt"), a, header="dp_top/dx = S(x) sum_m a_m xi^m [+ a_bump (S_bump - S), last]")
        xt = np.linspace(0.0, LX, 1301)
        np.savetxt(os.path.join(d, "ptop.txt"), np.c_[xt, top_pressure(xt, a, p_ramp, bump)], fmt="%.10e",
                   header=f"x  p_top  (iteration {it}, beta target {beta_t:+.3f}; p_top(650) = 0)")
        prev = run_steady(cfg, base, d, f"it{it}", prev, f["t"], log)
        xv, vt = top_plane_v(prev)
        np.savetxt(os.path.join(d, "vtop_result.txt"), np.c_[xv, vt], fmt="%.10e",
                   header="x  V_top: the top-face velocity of this iteration's field")
    log.close()


if __name__ == "__main__":
    main()
