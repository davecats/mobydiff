#!/usr/bin/env python3
"""Compare the mobydiff ZPG-TBL (finey case) against the reference DNS codes of the
resolution study -- SIMSON (spectral), CaNS and AMPHIBIOUS -- all reduced through
the SAME boundary-layer post-processing so the curves are strictly comparable.

  compare_codes.py [--mobydiff assets/mobydiff/xyz_4096_224_192/data.nc]
                   [--simson assets/postpro/passivewall.hdf5]
                   [--cans <data.nc>] [--amphibious <data.nc>]
                   [--retheta 677] [--out assets/figures/code_comparison.png]

All datasets share the nondimensionalization (length delta*_in, Re_delta*,in=450),
so the mean flow and Reynolds stresses are compared directly at a matched Re_theta.
The reference codes use the strong Schlatter-Orlu trip; mobydiff (finey) uses the
weak trip, so the transition location differs -- the comparison is at matched
Re_theta in the DEVELOPED region, where the trip is forgotten.

Four panels: c_f(Re_theta), H(Re_theta), U+(y+) and the Reynolds stresses
u'/v'/w'_rms+ and -u'v'+ at the matched Re_theta.
"""
import argparse
import os

import numpy as np
import h5py
import netCDF4 as nc

HERE = os.path.dirname(__file__)
REF_DEFAULT = os.path.expanduser("~/sshfsmountpoint/tbl-dns/res_study/data_processed")


def load_nc(path):
    """CaNS/AMPHIBIOUS/mobydiff data.nc -> (nu, x, y, fields[x,y]). Handles the
    (x,y) vs (y,x) dimension order and node- vs cell-centred y transparently."""
    f = nc.Dataset(path, "r")
    nu = float(f.nu)
    x = np.asarray(f["x"][:]); y = np.asarray(f["y"][:])

    def fld(name):
        v = f[name]
        a = np.asarray(v[:])
        return a if v.dimensions == ("x", "y") else a.T   # -> (x, y)

    out = {k: fld(k) for k in ("u_mean", "v_mean", "w_mean",
                               "uu_stress", "vv_stress", "ww_stress", "uv_stress")}
    return nu, x, y, out


def load_simson(path):
    with h5py.File(path, "r") as f:
        nu = 1.0 / float(f["parameters"].attrs["re_d1_0"])
        x = f["mesh"]["x"][...]; y = f["mesh"]["y"][...]
        u = f["mean"]["u"][...]                       # (x, y)
        cov = f["covariance"]
        out = {"u_mean": u, "v_mean": f["mean"]["v"][...], "w_mean": f["mean"]["w"][...],
               "uu_stress": cov["uu"][...], "vv_stress": cov["vv"][...],
               "ww_stress": cov["ww"][...], "uv_stress": cov["uv"][...]}
    return nu, x, y, out


def reduce_bl(nu, x, y, fld):
    """Unified reduction: prepend the wall (y=0, all quantities 0) if the grid does
    not already include it, then compute Re_theta, c_f, H, u_tau on the (x) line.
    Robust to a zeroed top boundary node (some AMPHIBIOUS files) and to all-zero
    inlet/ghost columns (masked out via Ue)."""
    y = np.asarray(y, float)
    fld = {k: np.asarray(v, float) for k, v in fld.items()}
    if y[0] > 1e-9:
        y = np.concatenate(([0.0], y))
        pad = lambda a: np.concatenate((np.zeros((a.shape[0], 1)), a), axis=1)
        fld = {k: pad(v) for k, v in fld.items()}
    # drop non-physical trailing top rows (u collapses to ~0 at a BC node).
    Uref = np.median(np.max(fld["u_mean"], axis=1))
    while y.size > 3 and np.median(fld["u_mean"][:, -1]) < 0.5 * Uref:
        y = y[:-1]; fld = {k: v[:, :-1] for k, v in fld.items()}
    u = fld["u_mean"]
    Ue = np.max(u, axis=1)                      # freestream (robust edge velocity)
    good = Ue > 1e-6
    Ues = np.where(good, Ue, 1.0)[:, None]
    theta = np.trapz((u / Ues) * (1 - u / Ues), y, axis=1)
    dstar = np.trapz(1 - u / Ues, y, axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        reth = np.where(good, Ue * theta / nu, np.nan)
        H = np.where(good & (theta > 0), dstar / theta, np.nan)
    y1, y2 = y[1], y[2]
    dudy = u[:, 1] * y2 / (y1 * (y2 - y1)) - u[:, 2] * y1 / (y2 * (y2 - y1))
    with np.errstate(divide="ignore", invalid="ignore"):
        utau = np.sqrt(np.abs(nu * dudy))
        cf = np.where(good, 2 * nu * dudy / Ue ** 2, np.nan)
    valid = good.copy()
    valid[np.nanargmax(reth) + 1:] = False      # keep the inflow->peak-Re_theta run
    return dict(nu=nu, y=y, fld=fld, Ue=Ue, utau=utau, cf=cf, H=H, reth=reth, valid=valid)


def station(R, retheta):
    idx = np.where(R["valid"] & (R["reth"] > 300))[0]
    i = idx[np.argmin(np.abs(R["reth"][idx] - retheta))]
    ut = R["utau"][i]; yp = R["y"] * ut / R["nu"]
    f = R["fld"]
    prof = dict(reth=R["reth"][i], yp=yp, Up=f["u_mean"][i] / ut,
                urms=np.sqrt(np.maximum(f["uu_stress"][i], 0)) / ut,
                vrms=np.sqrt(np.maximum(f["vv_stress"][i], 0)) / ut,
                wrms=np.sqrt(np.maximum(f["ww_stress"][i], 0)) / ut,
                uv=-f["uv_stress"][i] / ut ** 2)
    return prof


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mobydiff", default=os.path.join(HERE, "..", "mobydiff", "xyz_4096_224_192", "data.nc"))
    ap.add_argument("--simson", default=os.path.join(HERE, "passivewall.hdf5"))
    ap.add_argument("--cans", default=os.path.join(REF_DEFAULT, "cans", "xyz_3200_384_135", "data.nc"))
    ap.add_argument("--amphibious", default=os.path.join(REF_DEFAULT, "amphibious", "xyz_3200_384_135", "data.nc"))
    ap.add_argument("--retheta", type=float, default=677.0)
    ap.add_argument("--out", default=os.path.join(HERE, "..", "figures", "code_comparison.png"))
    a = ap.parse_args()

    datasets = []   # (label, color, style, R)
    def add(label, color, loader, path):
        if path and os.path.exists(path):
            datasets.append((label, color, reduce_bl(*loader(path))))
        else:
            print(f"  (skip {label}: {path} not found)")

    add("SIMSON (spectral)", "k", load_simson, a.simson)
    add("CaNS", "C0", load_nc, a.cans)
    add("AMPHIBIOUS", "C1", load_nc, a.amphibious)
    add("mobydiff (finey)", "C3", load_nc, a.mobydiff)

    print(f"{'code':22s} {'Re_th,max':>9s} {'cf@'+str(int(a.retheta)):>10s} {'H':>7s} "
          f"{'u_tau':>7s} {'u_rms+':>7s} {'-uv+':>7s}")
    for label, _, R in datasets:
        p = station(R, a.retheta)
        i = np.where(R["valid"] & (R["reth"] > 300))[0]
        j = i[np.argmin(np.abs(R["reth"][i] - a.retheta))]
        print(f"{label:22s} {R['reth'][R['valid']].max():9.0f} {1000*R['cf'][j]:9.3f}e-3 "
              f"{R['H'][j]:7.3f} {R['utau'][j]:7.4f} {p['urms'].max():7.3f} {p['uv'].max():7.3f}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 2, figsize=(14, 10))

    for label, c, R in datasets:
        m = R["valid"] & (R["reth"] > 150)
        lw = 1.8 if label.startswith("SIMSON") else 1.3
        ax[0, 0].plot(R["reth"][m], R["cf"][m], c=c, lw=lw, label=label)
        ax[0, 1].plot(R["reth"][m], R["H"][m], c=c, lw=lw, label=label)
    rr = np.linspace(200, 1000, 100)
    ax[0, 0].plot(rr, 0.024 * rr ** -0.25, "k:", lw=1.0, label=r"$0.024\,Re_\theta^{-1/4}$")
    ax[0, 0].axvline(a.retheta, c="0.7", ls="--", lw=0.8)
    ax[0, 0].set_xlabel(r"$Re_\theta$"); ax[0, 0].set_ylabel(r"$c_f$")
    ax[0, 0].set_xlim(150, 1000); ax[0, 0].set_ylim(0.002, 0.006)
    ax[0, 0].legend(fontsize=9); ax[0, 0].set_title("skin friction")

    ax[0, 1].axvline(a.retheta, c="0.7", ls="--", lw=0.8)
    ax[0, 1].set_xlabel(r"$Re_\theta$"); ax[0, 1].set_ylabel("H")
    ax[0, 1].set_xlim(150, 1000); ax[0, 1].set_ylim(1.35, 1.75)
    ax[0, 1].legend(fontsize=9); ax[0, 1].set_title("shape factor")

    for label, c, R in datasets:
        p = station(R, a.retheta)
        lw = 1.8 if label.startswith("SIMSON") else 0
        if label.startswith("SIMSON"):
            ax[1, 0].semilogx(p["yp"], p["Up"], "-", c=c, lw=lw, label=f"{label} ({p['reth']:.0f})")
        else:
            ax[1, 0].semilogx(p["yp"], p["Up"], "o", c=c, ms=3, mfc="none", label=f"{label} ({p['reth']:.0f})")
    ax[1, 0].set_xlabel(r"$y^+$"); ax[1, 0].set_ylabel(r"$U^+$")
    ax[1, 0].set_xlim(1, 500); ax[1, 0].set_ylim(0, 22)
    ax[1, 0].legend(loc="upper left", fontsize=9); ax[1, 0].set_title(f"mean velocity at Re$_\\theta$≈{a.retheta:.0f}")

    for label, c, R in datasets:
        p = station(R, a.retheta)
        for key, sgn in [("urms", 1), ("vrms", 1), ("wrms", 1), ("uv", -1)]:
            if label.startswith("SIMSON"):
                ax[1, 1].plot(p["yp"], sgn * p[key], "-", c=c, lw=1.6)
            else:
                ax[1, 1].plot(p["yp"], sgn * p[key], "o", c=c, ms=2.5, mfc="none")
        ax[1, 1].plot([], [], ("-" if label.startswith("SIMSON") else "o"),
                      c=c, lw=1.6, ms=4, mfc="none", label=label)
    ax[1, 1].annotate(r"$u'_{rms}$", (11, 2.75), fontsize=10)
    ax[1, 1].annotate(r"$w'_{rms}$", (45, 1.35), fontsize=10)
    ax[1, 1].annotate(r"$-\overline{u'v'}$", (70, 0.9), fontsize=10)
    ax[1, 1].annotate(r"$v'_{rms}$", (110, 1.02), fontsize=10)
    ax[1, 1].set_xlabel(r"$y^+$"); ax[1, 1].set_xlim(0, 300); ax[1, 1].set_ylim(0, 3)
    ax[1, 1].legend(fontsize=9, loc="upper right"); ax[1, 1].set_title(f"Reynolds stresses at Re$_\\theta$≈{a.retheta:.0f}")

    fig.suptitle(r"ZPG-TBL code comparison — mobydiff (finey) vs SIMSON / CaNS / AMPHIBIOUS "
                 f"($Re_\\theta$≈{a.retheta:.0f})", fontsize=13)
    fig.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    fig.savefig(a.out, dpi=150)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
