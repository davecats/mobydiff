#!/usr/bin/env python3
"""Convergence history of the top-pressure optimisations (top_optimise.py):

    python3 pg_history.py OUT.png LABEL=OPT_DIR:BETA [...]

Every iterate is re-evaluated with the CURRENT pg_bl.py edge, so the curves
are comparable even where the run itself changed definitions on the way
(the log records each restart). State 0 is the ZPG field the optimisation
starts from; state k + 1 is the steady field of it<k>/ (its final-march
field where there is one). Writes OUT.png and OUT.txt (max |beta - target|
and mean beta on x_fit per state).
"""
import glob
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rans_bl as rb  # noqa: E402
import pg_bl  # noqa: E402
import top_optimise as to  # noqa: E402


def states(opt):
    cfg = json.load(open(os.path.join(opt, "case.json")))
    out = [os.path.join(opt, cfg["zpg_field"])]
    k = 0
    while os.path.isdir(os.path.join(opt, f"it{k}")):
        fs = glob.glob(os.path.join(opt, f"it{k}", f"it{k}*_*.h5"))
        out.append(max(fs, key=lambda s: int(s.rsplit("_", 1)[1][:-3])))
        k += 1
    nc = [len(np.atleast_1d(np.loadtxt(os.path.join(opt, f"it{j}", "coeffs.txt")))) for j in range(k)]
    return cfg, out, [0] + nc


def main():
    out = sys.argv[1]
    cases = []
    for arg in sys.argv[2:]:
        lab, spec = arg.split("=", 1)
        d, b = spec.split(":")
        cases.append((lab, d, float(b)))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import cm
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    fig, ax = plt.subplots(len(cases), 2, figsize=(14, 4.6 * len(cases)), squeeze=False,
                           gridspec_kw={"width_ratios": [1.7, 1]})
    rows = []
    for r, (lab, d, beta) in enumerate(cases):
        cfg, fields, ncoef = states(d)
        xa, xb = cfg["x_fit"]
        err, mean = [], []
        n = len(fields)
        cmap = cm.get_cmap("Blues" if beta < 0 else "Oranges")
        for k, fn in enumerate(fields):
            q = pg_bl.pg_integrals(rb.load(fn))
            x = q["x"]
            bt = beta * to.ramp(x, *cfg["x_ramp"])
            w = (x >= xa) & (x <= xb)
            err.append(np.abs(q["beta"] - bt)[w].max())
            mean.append(q["beta"][w].mean())
            rows.append(f"{lab} {k:2d} {err[-1]:.3f} {mean[-1]:+.3f} {ncoef[k]} {os.path.relpath(fn, d)}")
            s = (x > 5) & (x < 645)
            last = k == n - 1
            ax[r, 0].plot(x[s], q["beta"][s], color="k" if last else cmap(0.25 + 0.6 * k / max(n - 1, 1)),
                          lw=2.2 if last else 1.0, label=f"state {k} (final)" if last else None, zorder=3 if last else 2)
        ax[r, 0].plot(x, bt, color="k", ls=":", lw=1, label="target")
        ax[r, 0].axvspan(xa, xb, color="#ddd", alpha=0.4, lw=0)
        ax[r, 0].set(xlabel="x", ylabel="β", xlim=(0, 650),
                     title=f"({'ac'[r]}) {lab}: β(x) of every iterate (light = early, ZPG = state 0)")
        ax[r, 0].legend(frameon=False, fontsize=8)
        k = np.arange(n)
        bump = np.array(ncoef) > cfg.get("npoly", 4)
        ax[r, 1].semilogy(k, err, color="#555", lw=1)
        ax[r, 1].semilogy(k[~bump], np.array(err)[~bump], "o", color="#999", ms=6, label="4-term basis")
        ax[r, 1].semilogy(k[bump], np.array(err)[bump], "s", color="k", ms=6, label="+ steeper ramp start")
        ax[r, 1].axhline(cfg["tol"], color="k", ls=":", lw=1)
        ax[r, 1].text(0.2, cfg["tol"] * 1.1, f"tol {cfg['tol']}", fontsize=8)
        ax[r, 1].set(xlabel="state (iterate)", ylabel=f"max |β − target|, x {xa:.0f}–{xb:.0f}",
                     title=f"({'bd'[r]}) {lab}: convergence")
        ax[r, 1].legend(frameon=False, fontsize=8)
        for a in ax[r]:
            a.grid(alpha=0.25, lw=0.6)
    fig.tight_layout()
    fig.savefig(out, dpi=100)
    with open(os.path.splitext(out)[0] + ".txt", "w") as fh:
        fh.write("# case state max|beta-target| mean_beta n_coeffs field (pg_bl edge of today)\n")
        fh.write("\n".join(rows) + "\n")
    print("wrote", out)


if __name__ == "__main__":
    main()
