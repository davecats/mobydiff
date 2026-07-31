#!/usr/bin/env python3
"""Campaign summary: c_f(Re_theta) and H(Re_theta) for the four boundaryLayer
cases (trip x resolution 2x2) vs the SIMSON spectral reference."""
import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

nu = 1 / 450.


def curves(fn):
    f = h5py.File(fn, "r")
    nx, ny = int(f.attrs["nx"]), int(f.attrs["ny"])
    p = f["profile"][...].reshape(nx, ny, -1); y = f["ycoord"][...]
    dy = np.diff(np.concatenate(([0], 0.5 * (y[:-1] + y[1:]), [y[-1]])))
    Um = p[:, :, 0]; Ue = Um[:, -1]
    th = np.sum((Um / Ue[:, None]) * (1 - Um / Ue[:, None]) * dy, axis=1)
    ds = np.sum((1 - Um / Ue[:, None]) * dy, axis=1)
    reth = Ue * th / nu
    cf = 2 * (nu * Um[:, 0] / y[0]) / Ue ** 2
    return reth, cf, ds / th


cases = [
    ("turbulent_two (trip 0.15, std grid)", "turbulent_two/production_stats.h5", "C0"),
    ("case i  (trip 0.03, std grid)", "turbulent_lowtrip/production_stats.h5", "C1"),
    ("case ii (trip 0.15, fine outer)", "turbulent_finewall/production_stats.h5", "C2"),
    ("case iii(trip 0.03, fine outer)", "turbulent_lowtrip_finewall/production_stats.h5", "C4"),
]
s = h5py.File("turbulent_two/passivewall.hdf5", "r")
ys = s["mesh"]["y"][...]; us = s["mean"]["u"][...]; Ues = us[:, -1]
ths = np.trapz((us / Ues[:, None]) * (1 - us / Ues[:, None]), ys, axis=1); rs = Ues * ths / nu
h1, h2 = ys[1], ys[2]
du = (h2 / (h1 * (h2 - h1))) * us[:, 1] - (h1 / (h2 * (h2 - h1))) * us[:, 2]
cfs = 2 * nu * du / Ues; Hs = np.trapz(1 - us / Ues[:, None], ys, axis=1) / ths
val = np.zeros(len(rs), bool); val[: np.argmax(rs) + 1] = True; o = np.argsort(rs[val])

fig, ax = plt.subplots(1, 2, figsize=(14, 5.5))
for lbl, fn, c in cases:
    try:
        r, cf, H = curves(fn); m = (r > 380) & (r < 950)
        ax[0].plot(r[m], cf[m], c=c, lw=1.4, label=lbl)
        ax[1].plot(r[m], H[m], c=c, lw=1.4, label=lbl)
    except Exception as e:
        print(lbl, "->", e)
ax[0].plot(rs[val][o], cfs[val][o], "k-", lw=2, label="SIMSON (spectral)")
ax[1].plot(rs[val][o], Hs[val][o], "k-", lw=2, label="SIMSON (spectral)")
ax[0].set_xlim(400, 900); ax[0].set_ylim(0.0042, 0.0052); ax[0].set_xlabel(r"$Re_\theta$"); ax[0].set_ylabel(r"$c_f$"); ax[0].legend(fontsize=8); ax[0].set_title("skin friction")
ax[1].set_xlim(400, 900); ax[1].set_ylim(1.48, 1.58); ax[1].set_xlabel(r"$Re_\theta$"); ax[1].set_ylabel("H"); ax[1].legend(fontsize=8); ax[1].set_title("shape factor")
fig.suptitle("boundaryLayer campaign: trip x resolution vs SIMSON")
fig.tight_layout(); fig.savefig("campaign_summary.png", dpi=150)
print("wrote campaign_summary.png")
