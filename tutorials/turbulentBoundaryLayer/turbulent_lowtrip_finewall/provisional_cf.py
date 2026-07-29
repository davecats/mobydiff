#!/usr/bin/env python3
"""PROVISIONAL c_f check for the reduced-tripping case while phase 1 is still
running (no time-averaged stats yet). Span-averages the latest instantaneous
field snapshot(s) and compares c_f(Re_theta) against the converged turbulent_two
run and the SIMSON reference. NOT time-converged -- indicative only."""
import glob
import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

nu = 1 / 450.


def cf_field(fn):
    f = h5py.File(fn, "r")
    y, x = f["y"][...], f["x"][...]
    u = f["un"][0].mean(axis=0).T                      # span (z) average -> (nx, ny)
    yc = 0.5 * (y[:-1] + y[1:]); Ue = u[:, -1]
    yf = np.empty(len(yc) + 1); yf[1:-1] = 0.5 * (yc[:-1] + yc[1:]); yf[0] = 0
    yf[-1] = yc[-1] + (yc[-1] - yf[-2]); dy = np.diff(yf)
    th = np.sum((u / Ue[:, None]) * (1 - u / Ue[:, None]) * dy, axis=1)
    return Ue * th / nu, 2 * (nu * u[:, 0] / yc[0]) / Ue ** 2


snaps = sorted(glob.glob("production_[0-9]*.h5"))[-2:]      # last two snapshots
rc = [cf_field(s) for s in snaps]
reth = np.mean([r for r, _ in rc], axis=0); cf = np.mean([c for _, c in rc], axis=0)

g = h5py.File("../turbulent_two/production_stats.h5", "r")
nx, ny = int(g.attrs["nx"]), int(g.attrs["ny"])
p = g["profile"][...].reshape(nx, ny, -1); yg = g["ycoord"][...]
dyg = np.diff(np.concatenate(([0], 0.5 * (yg[:-1] + yg[1:]), [yg[-1]])))
Um = p[:, :, 0]; Ue2 = Um[:, -1]
th2 = np.sum((Um / Ue2[:, None]) * (1 - Um / Ue2[:, None]) * dyg, axis=1)
reth2 = Ue2 * th2 / nu; cf2 = 2 * (nu * Um[:, 0] / yg[0]) / Ue2 ** 2

s = h5py.File("passivewall.hdf5", "r"); ys = s["mesh"]["y"][...]; us = s["mean"]["u"][...]
Ues = us[:, -1]; ths = np.trapz((us / Ues[:, None]) * (1 - us / Ues[:, None]), ys, axis=1)
rs = Ues * ths / nu
h1, h2 = ys[1], ys[2]
du = (h2 / (h1 * (h2 - h1))) * us[:, 1] - (h1 / (h2 * (h2 - h1))) * us[:, 2]
cfs = 2 * nu * du / Ues ** 2
val = np.zeros(len(rs), bool); val[: np.argmax(rs) + 1] = True; o = np.argsort(rs[val])

fig, ax = plt.subplots(1, 2, figsize=(14, 5.5))
for a in ax:
    m2 = reth2 > 150; m1 = reth > 150
    a.plot(reth2[m2], cf2[m2], "-", c="C0", lw=1.5, label="turbulent_two (trip 0.15, converged)")
    a.plot(rs[val][o], cfs[val][o], "-", c="C3", lw=1.5, label="SIMSON (trip x=10)")
    a.plot(reth[m1], cf[m1], ".", c="C2", ms=4, label="low-trip 0.03 (PROVISIONAL)")
    a.plot(reth2[m2], 0.024 * reth2[m2] ** -0.25, "k--", lw=0.8, label=r"$0.024\,Re_\theta^{-1/4}$")
    a.set_xlabel(r"$Re_\theta$"); a.set_ylabel(r"$c_f$"); a.legend(fontsize=8)
ax[0].set_ylim(0.002, 0.019); ax[0].set_title("full (transition overshoot)")
ax[1].set_ylim(0.004, 0.0058); ax[1].set_xlim(400, 900); ax[1].set_title("developed zoom")
fig.suptitle(f"PROVISIONAL c_f, reduced tripping (span-avg of {', '.join(snaps)}; not time-converged)")
fig.tight_layout(); fig.savefig("provisional_lowtrip_cf.png", dpi=140)
print("wrote provisional_lowtrip_cf.png")
