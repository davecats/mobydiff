#!/usr/bin/env python3
"""Compare the conjugate-pipe campaign against Neuhauser's NekRS DNS.

    ./compare_neuhauser.py velocity p_coarse_vel.npz
    ./compare_neuhauser.py thermal  p_coarse_stats.npz --scalar c0

THEIR DEFINITIONS, NOT OURS (docs/next_session_pipe_cht.md Section 2). The
dataset is polar and already z-averaged, with `time` a BATCH axis kept for
uncertainty estimation, so every reference number here is a batch mean and
carries the batch scatter as an error bar. u_tau is taken the way
`dataset/datawrapper.py` takes it -- from d<w>/dr at r -> 0.5 from the fluid
side -- and y+ = (0.5 - r)/delta_v with delta_v = nu/u_tau.

Their (u, v) are CARTESIAN and w is axial; the radial and azimuthal
components are ur = u cos(phi) + v sin(phi), ut = -u sin(phi) + v cos(phi),
which is the rotation pipe_stats.py applies to our own cells.

THEIR RADIAL GRID IS FAR FINER THAN OURS (226 points, clustered at the
interfaces, dr down to 1e-4 = 0.04 wall units), so THEIR profile is
interpolated onto OUR bin centres -- never the reverse.

The isc index is resolved from (type, K, lambda_sf, bccode), i.e. by the case
the row IS, so the mapping to our scalar names cannot silently rot:

    c0 -> CHT K=1    lambda=1        c3 -> CHT K=4    lambda=1
    c1 -> CHT K=1    lambda=2        c4 -> CHT K=0.25 lambda=1
    c2 -> CHT K=1    lambda=0.5      mbc -> MBC       isof -> IF
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np

ZARR = os.environ.get(
    "NEUHAUSER_ZARR",
    "/tmp/claude-215842/neuhauser/joined_datasets.interp.zarr")

# our scalar name -> (type, K, lambda_sf) in their coordinates
CASE_MAP = {
    "c0":   ("CHT", 1.0, 1.0),
    "c1":   ("CHT", 1.0, 2.0),
    "c2":   ("CHT", 1.0, 0.5),
    "c3":   ("CHT", 4.0, 1.0),
    "c4":   ("CHT", 0.25, 1.0),
    "mbc":  ("MBC", None, None),
    "isof": ("IF", None, None),
}


def open_case(name, pr, bccode=0, zarr=ZARR):
    import xarray
    ds = xarray.open_zarr(zarr, consolidated=False)
    typ, K, lam = CASE_MAP[name]
    sel = (np.array([str(x) for x in np.asarray(ds["type"])]) == typ) & \
          (np.asarray(ds["bccode"]) == bccode)
    if K is not None:
        sel &= np.isclose(np.asarray(ds["K"]), K) & \
               np.isclose(np.asarray(ds["lambda_sf"]), lam)
    idx = int(np.flatnonzero(sel)[0])
    return ds.isel(isc=idx).sel(Pr=pr)


CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                     "neuhauser_profiles.npz")


def ref_case(name, pr=0.71, thermal=True, zarr=None):
    """The reference radial profiles for one case: batch mean + standard error.

    PREFERS THE SHIPPED CACHE (`neuhauser_profiles.npz`, written by
    extract_neuhauser.py), so the comparison runs without the 11.8 GB archive.
    Point NEUHAUSER_ZARR at an extracted zarr store -- or pass `zarr` -- to go
    back to the archive; the two paths produce the same numbers, which is what
    makes the cache checkable rather than merely convenient.
    """
    if zarr is None and pr == 0.71 and os.path.exists(CACHE):
        z = np.load(CACHE)
        pref = f"{name}/"
        out = {k[len(pref):]: z[k] for k in z.files if k.startswith(pref)}
        if out:
            return out
    return batch(polar_moments(open_case(name, pr, zarr=zarr or ZARR),
                               thermal=thermal))


def ref_budget():
    """The <theta'^2> budget terms of case c0 (see variance_budget.py)."""
    z = np.load(CACHE)
    pref = "budget/"
    return {k[len(pref):]: z[k] for k in z.files if k.startswith(pref)}


def polar_moments(d, thermal=False):
    """Phi-averaged radial profiles, and the batch scatter of each.

    Second moments are formed at every (r, phi) BEFORE the azimuthal average,
    so the azimuthal variation of the mean does not leak into the rms.
    """
    c, s = np.cos(d["phi"].values), np.sin(d["phi"].values)

    def arr(name):
        return np.asarray(d[name])           # (time, r, phi)

    u, v, w = arr("u"), arr("v"), arr("w")
    ur = u * c + v * s
    ut = -u * s + v * c
    urr = arr("u*u") * c * c + 2 * arr("u*v") * c * s + arr("v*v") * s * s
    utt = arr("u*u") * s * s - 2 * arr("u*v") * c * s + arr("v*v") * c * c
    urw = arr("u*w") * c + arr("v*w") * s

    out = {
        "r": d["r"].values,
        "uz": w.mean(axis=2),
        "ur": ur.mean(axis=2),
        "uz_rms": np.sqrt(np.maximum(arr("w*w") - w * w, 0.0)).mean(axis=2),
        "ur_rms": np.sqrt(np.maximum(urr - ur * ur, 0.0)).mean(axis=2),
        "ut_rms": np.sqrt(np.maximum(utt - ut * ut, 0.0)).mean(axis=2),
        "urz": (urw - ur * w).mean(axis=2),
    }
    if thermal:
        t = arr("t")
        tur = arr("t*u") * c + arr("t*v") * s
        out["t"] = t.mean(axis=2)
        out["t_rms"] = np.sqrt(np.maximum(arr("t*t") - t * t, 0.0)).mean(axis=2)
        out["tur"] = (tur - t * ur).mean(axis=2)
        out["tuz"] = (arr("t*w") - t * w).mean(axis=2)
    return out                                   # each (time, r)


def batch(profile):
    """Batch mean and standard error over the `time` axis."""
    out = {"r": profile["r"]}
    for k, v in profile.items():
        if k == "r":
            continue
        out[k] = v.mean(axis=0)
        out[k + "_err"] = v.std(axis=0, ddof=1) / np.sqrt(v.shape[0])
    return out


def wall_flux(r, t, nu, pr, interface=0.5):
    """q_w = alpha |dT/dr| at the wall, taken on the FLUID side with
    alpha = nu/Pr -- the same definition on both sides of the comparison, so
    the two scales are formed alike even though the two runs normalise their
    heating differently (theirs makes q_w = 1; ours sets source = 1 and lets
    q_w follow the balance, which gives q_w ~ 0.244)."""
    m = (r < interface) & np.isfinite(t)
    rr, tt = r[m][-6:], t[m][-6:]
    slope = np.polyfit(rr, tt, 1)[0]
    return (nu / pr) * abs(slope)


def wall_value(r, t, interface=0.5):
    """theta at the interface, by linear extrapolation of the last fluid rows.
    THE MEAN MUST BE REFERENCED TO IT: their `t` is wall-referenced (their
    datawrapper forms theta = t - t_wall) while ours is absolute, so comparing
    the raw levels compares two different additive constants."""
    m = (r < interface) & np.isfinite(t)
    rr, tt = r[m][-6:], t[m][-6:]
    a, b = np.polyfit(rr, tt, 1)
    return a * interface + b


def wall_scales(p, nu):
    """u_tau and delta_v their way: -d<w>/dr at the last fluid radius."""
    r, w = p["r"], p["uz"]
    m = r < 0.499999
    dwdr = np.gradient(w[m], r[m])[-1]
    u_tau = np.sqrt(nu * -dwdr)
    return u_tau, nu / u_tau


def our_profiles(npz):
    here = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, here)
    sys.path.insert(0, os.path.dirname(here))   # pipe_stats.py, beside the inis
    from pipe_stats import profiles
    acc = dict(np.load(npz, allow_pickle=True))
    return profiles(acc), acc


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("what", choices=("velocity", "thermal"))
    ap.add_argument("npz", help="pipe_stats.py accumulator of our run")
    ap.add_argument("--scalar", default="c0", help="which scalar (thermal)")
    ap.add_argument("--pr", type=float, default=0.71)
    ap.add_argument("--re", type=float, default=5300.0)
    ap.add_argument("--zarr", default=ZARR)
    ap.add_argument("--plot", help="write a comparison figure here")
    ap.add_argument("--wall-flux", type=float, default=None,
                    help="our q_w, from the interface-heat balance"
                         " (|Q|/(2 pi R L)). PREFER IT: the slope estimate is"
                         " differenced across the CUT CELLS and reads 6.6 %%"
                         " low on the coarse grid, against 0.0 %% on the fine"
                         " one -- so normalising by it makes the coarse grid"
                         " look better than it is, by cancelling a too-small"
                         " theta_tau against a too-low profile.")
    a = ap.parse_args()

    nu = 1.0 / a.re
    name = a.scalar if a.what == "thermal" else "c0"
    ref = ref_case(name, a.pr, thermal=(a.what == "thermal"),
                   zarr=(a.zarr if a.zarr != ZARR else None))
    u_tau, delta_v = wall_scales(ref, nu)
    print(f"# reference: {name}  u_tau = {u_tau:.6f}  delta_v = {delta_v:.6e}"
          f"  Re_tau = {0.5/delta_v:.1f}")

    ours, acc = our_profiles(a.npz)
    r = ours["r"]
    fluid = (ours["count"][:, 0] > 0) & (r < 0.5)

    # our own wall scales, from the same definition on our binned profile
    rr, ww = r[fluid], ours["uz"][fluid, 0]
    # The force balance fixes u_tau exactly in a steady pipe; the profile
    # slope is only a MEASUREMENT of it, and a biased one (it is differenced
    # across the cut cells). An accumulator built from the solver's plane
    # statistics carries no velocity at all, so fall back to the exact value.
    u_tau_force = np.sqrt(1.8702750564e-2 * 0.25)
    if np.any(ww != 0.0):
        u_tau_us = np.sqrt(nu * -np.gradient(ww, rr)[-1])
        print(f"# ours:      u_tau = {u_tau_us:.6f} from the profile slope"
              f"  (force balance gives {u_tau_force:.6f})")
        ub_us = np.trapezoid(ww * rr, rr) / np.trapezoid(rr, rr)
        ubr = np.trapezoid(np.interp(rr, ref["r"], ref["uz"]) * rr, rr) \
            / np.trapezoid(rr, rr)
        print(f"# bulk velocity: ours {ub_us:.4f}   reference {ubr:.4f}")
    else:
        u_tau_us = u_tau_force
        print(f"# ours:      no velocity in this accumulator (plane"
              f" statistics) -- using the force-balance u_tau"
              f" {u_tau_us:.6f}")

    def interp(key):
        return np.interp(rr, ref["r"], ref[key])

    if a.what == "velocity":
        keys = [("uz", "uz"), ("uz_rms", "uz_rms"), ("ur_rms", "ur_rms"),
                ("up_rms", "ut_rms"), ("urz_cov", "urz")]
        head = "<uz>", "uz'", "ur'", "ut'", "<ur uz>"
        scale = {k: 1.0 for k in head}
    else:
        keys = [(f"{a.scalar}_mean", "t"), (f"{a.scalar}_rms", "t_rms"),
                (f"{a.scalar}_flux_ur", "tur"), (f"{a.scalar}_flux_uz", "tuz")]
        head = "<theta>+", "theta'+", "<ur t>+", "<uz t>+"
        # Each side is normalised by ITS OWN theta_tau, and the two runs are
        # then on the same scale even though their heating normalisations
        # differ. The fluxes carry u_tau as well.
        # Wall units on each side from ITS OWN wall flux: theta_tau = q_w/u_tau.
        qwr = wall_flux(ref["r"], ref["t"], nu, a.pr)
        qwo_slope = wall_flux(rr, ours[f"{a.scalar}_mean"][fluid, 0], nu, a.pr)
        qwo = a.wall_flux if a.wall_flux else qwo_slope
        if a.wall_flux:
            print(f"# wall flux: balance {qwo:.5f}, near-wall slope"
                  f" {qwo_slope:.5f} ({100*(qwo_slope/qwo-1):+.1f} %)")
        ttr, tto = qwr / u_tau, qwo / u_tau_us
        # ...and the mean referenced to the interface temperature, because
        # theirs is wall-referenced and ours is absolute.
        twr = wall_value(ref["r"], ref["t"])
        two = wall_value(rr, ours[f"{a.scalar}_mean"][fluid, 0])
        print(f"# theta_tau: ours {tto:.5f}   reference {ttr:.5f}"
              f"   (wall flux ours {qwo:.5f}, reference {qwr:.5f};"
              f" interface theta ours {two:.4g}, reference {twr:.4g})")
        ref = dict(ref)
        ours = dict(ours)
        ref["t"] = (ref["t"] - twr) / ttr
        ref["t_rms"] = ref["t_rms"] / ttr
        ref["tur"] = ref["tur"] / (ttr * u_tau)
        ref["tuz"] = ref["tuz"] / (ttr * u_tau)
        ours[f"{a.scalar}_mean"] = (ours[f"{a.scalar}_mean"] - two) / tto
        ours[f"{a.scalar}_rms"] = ours[f"{a.scalar}_rms"] / tto
        ours[f"{a.scalar}_flux_ur"] = ours[f"{a.scalar}_flux_ur"] / (tto * u_tau_us)
        ours[f"{a.scalar}_flux_uz"] = ours[f"{a.scalar}_flux_uz"] / (tto * u_tau_us)

    if a.what == "velocity":
        # THE EXACT LAW, which needs no reference file: integrating the
        # axial momentum balance of a steady, fully-developed pipe once gives
        #     -nu d<w>/dr + <u_r' u_z'> = u_tau^2 (r/R),
        # so the deviation from that line measures how converged the run is
        # (and how well the immersed wall delivers the stress it is driven
        # with). It is the pipe form of the channel law the conjugate-channel
        # campaign gated on. MIND THE SIGN: with r measured from the AXIS the
        # mean shear is negative, so the down-gradient Reynolds stress
        # <u_r' u_z'> is POSITIVE and ADDS to the viscous term -- the opposite
        # of the channel's <u'v'> < 0. Both this code and the reference data
        # carry that sign.
        tau = -nu * np.gradient(ours["uz"][fluid, 0], rr) \
            + ours["urz_cov"][fluid, 0]
        exact = (1.8702750564e-2 * 0.25) * (rr / 0.5)
        m = rr < 0.5 - 0.02
        dev = np.abs(tau[m] - exact[m])
        print(f"# total stress -nu dW/dr + <ur uz> vs u_tau^2 r/R:"
              f"  max dev {dev.max():.3e}  mean {dev.mean():.3e}"
              f"  (u_tau^2 = {1.8702750564e-2*0.25:.5e})")

    print("#      r       y+  " + "".join(f"{h:>12}{'ref':>12}{'dev%':>8}"
                                          for h in head))
    for i in range(len(rr)):
        yp = (0.5 - rr[i]) / delta_v
        row = f"{rr[i]:8.5f} {yp:8.2f}"
        for ours_key, ref_key in keys:
            o = ours[ours_key][fluid, 0][i]
            q = interp(ref_key)[i]
            # A column the accumulator does not carry (the plane statistics
            # have no axial flux) is reported as absent, not as -100 %.
            if not np.any(ours[ours_key][fluid, 0] != 0.0):
                row += f"{'n/a':>12}{q:12.4e}{'':>8}"
                continue
            dev = 100.0 * (o - q) / q if q != 0 else float("nan")
            row += f"{o:12.4e}{q:12.4e}{dev:8.1f}"
        print(row)

    if a.what == "thermal":
        # THE CONJUGATE SIGNATURE (Section 6.3): theta' through the SOLID,
        # as a fraction of its interface value. It does not decay to zero --
        # it flattens over the last 30 %, the signature of the adiabatic
        # (constant-flux) outer surface reflecting the fluctuation.
        solid = (ours["count"][:, 1] > 0) & (r > 0.5) & (r < 0.5 + 0.1)
        rs = r[solid]
        ours_s = ours[f"{a.scalar}_rms"][solid, 1]
        iface_o = ours_s[0] if len(ours_s) else float("nan")
        ref_iface = np.interp(0.5001, ref["r"], ref["t_rms"])
        print("\n# through the solid:  (r-R)/d    ours/iface   ref/iface")
        for i, rv in enumerate(rs):
            q = np.interp(rv, ref["r"], ref["t_rms"]) / ref_iface
            print(f"{(rv - 0.5)/0.1:24.4f} {ours_s[i]/iface_o:12.4f} {q:12.4f}")

    if a.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, len(keys), figsize=(4 * len(keys), 3.2))
        for k, (ours_key, ref_key) in enumerate(keys):
            yp = (0.5 - rr) / delta_v
            ax[k].plot(yp, ours[ours_key][fluid, 0], "o-", ms=3, label="mobydiff")
            ax[k].plot((0.5 - ref["r"]) / delta_v, ref[ref_key], "k-",
                       lw=1, label="Neuhauser")
            ax[k].set_xscale("log")
            ax[k].set_xlabel("$y^+$")
            ax[k].set_title(head[k])
            ax[k].set_xlim(0.5, 200)
        ax[0].legend()
        fig.tight_layout()
        fig.savefig(a.plot, dpi=130)
        print(f"   wrote {a.plot}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
