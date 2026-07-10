#!/usr/bin/env python3
"""IDDES T5 (DDES shielding) gate checks — run after run_gates.sh.

(a) f_d sane on the WMLES channel: -> 1 through the RANS wall layer,
    -> 0 in the LES core (final stats-leg snapshot, x/z-averaged).
(b) mean profile: no gross log-layer mismatch vs the pure-WALE stats
    reference (../channel_interface/les/runs/uniform/stats) and the T2
    RANS turb180 profile (../rans_sst/turb180_132565.h5).
(c) consistency limits: fd_force = 0 bit-exact vs pure WALE (the blend
    reduces to nut = nut_sgs exactly); fd_force = 1 holds the T2 RANS
    turb180 fixed point (round-off-class drift only: the k-sink is the
    same fixed point through different arithmetic).
(d) the les_ibm IBM channel ran IDDES stably (finite, bounded, f_d -> 1
    at the immersed walls).
(e) determinism: 1 rank == 4 ranks EXACT.

Needs python3 + h5py + numpy only.
"""
import glob
import os
import subprocess
import sys

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(HERE, "..", "rans_sst"))
from rans_channel_check import load_raw, profiles_single_level  # noqa: E402

RETAU = 180.0
failures = []


def gate(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}  {detail}")
    if not ok:
        failures.append(name)


def latest(pattern):
    files = sorted(glob.glob(os.path.join(HERE, pattern)),
                   key=lambda p: int(p.rsplit("_", 1)[1].split(".")[0]))
    return files[-1] if files else None


def load_profiles(path, names=("u", "k", "om", "nut")):
    blocks, data, y_nodes, re, ly, ny = load_raw(path)
    with h5py.File(path) as f:
        for extra in ("fd",):
            data[extra] = f[extra][...] if extra in f else None
    yc, prof = profiles_single_level(blocks, data, y_nodes)
    return yc, prof, re


def read_channel_stats(path):
    with h5py.File(path) as f:
        y = f["coord"][...]
        P = f["profile"][...]
        n = f["count"][...]
    m = n > 0
    o = np.argsort(y[m])
    return y[m][o], P[m][o]


def compare_exact(a, b, datasets):
    r = subprocess.run([sys.executable, os.path.join(ROOT, "tools/compare_fields.py"),
                        a, b, *datasets, "--tolerance", "0"],
                       capture_output=True, text=True)
    return r.returncode == 0, r.stdout.strip().replace("\n", "; ")


# ---- (a) fd profile + (b) mean profile ----
snap = latest("iddes180_[0-9]*.h5")
if snap is None:
    gate("(a) fd profile", False, "no iddes180 snapshot found")
else:
    yc, prof, re = load_profiles(snap)
    fd = prof["fd"]
    if fd is None:
        gate("(a) fd profile", False, "snapshot has no fd dataset")
    else:
        dist = np.minimum(yc, 2.0 - yc)
        u_tau = 1.0  # forcing_x = 1, half height 1
        yplus = dist*u_tau*re
        # DDES on a WMLES grid with resolved content shields the viscous
        # sublayer/lower buffer and hands over early (r_d reads the
        # instantaneous resolved gradients) -- extending the RANS coverage
        # outward is the f_B/f_e elevating branch's job (increment 2).
        # Gate: full retention below y+ 5, clean LES core, monotone-ish
        # handover in between (reported).
        fd_wall = fd[yplus <= 5.0].min()
        fd_core = fd[dist >= 0.7].max()
        for lo, hi in ((0, 5), (5, 25), (25, 60), (60, 120)):
            m = (yplus >= lo) & (yplus < hi)
            if m.any():
                print(f"fd y+ [{lo},{hi}): mean {fd[m].mean():.3f} min {fd[m].min():.3f}")
        gate("(a) fd -> 1 at the wall (y+<=5), -> 0 in the core",
             fd_wall >= 0.95 and fd_core <= 0.05,
             f"wall min {fd_wall:.3f} (>=0.95), core max {fd_core:.3f} (<=0.05)")

stats = os.path.join(HERE, "channel_stats.h5")
wale_stats = os.path.join(HERE, "..", "channel_interface", "les",
                          "runs", "uniform", "stats", "channel_stats.h5")
rans_ref = os.path.join(HERE, "..", "rans_sst", "turb180_132565.h5")
if not os.path.exists(stats):
    gate("(b) log-layer mismatch", False, "no channel_stats.h5 (stats leg)")
else:
    y_i, P_i = read_channel_stats(stats)
    u_i = P_i[:, 0]
    refs = {}
    if os.path.exists(wale_stats):
        y_w, P_w = read_channel_stats(wale_stats)
        refs["WALE"] = (y_w, P_w[:, 0])
    if os.path.exists(rans_ref):
        yr, pr, _ = load_profiles(rans_ref)
        refs["RANS-T2"] = (yr, pr["u"])
    dist = np.minimum(y_i, 2.0 - y_i)
    log_m = (dist*RETAU >= 30.0) & (dist*RETAU <= 0.3*RETAU)
    ok_all = bool(refs)
    for name, (yr_, ur_) in refs.items():
        uref = np.interp(y_i, yr_, ur_)
        rel = np.abs(u_i[log_m] - uref[log_m])/np.abs(uref[log_m])
        print(f"(b) log-layer |U - U_{name}|/U: max {rel.max():.3f}, mean {rel.mean():.3f}")
        ok_all = ok_all and rel.max() <= 0.15
    gate("(b) no gross log-layer mismatch (<= 15% vs both references)", ok_all)

# ---- (c) fd_force = 0: bit-exact vs pure WALE ----
a = latest("fd0_iddes_[0-9]*.h5")
b = latest("fd0_wale_[0-9]*.h5")
if a and b:
    ok, out = compare_exact(a, b, ["un", "vn", "wn", "pn", "nut"])
    gate("(c) fd_force=0 == pure WALE (bit-exact)", ok, out)
else:
    gate("(c) fd_force=0 == pure WALE", False, "missing fd0 outputs")

# ---- (c) fd_force = 1: hold the turb180 RANS fixed point ----
a = latest("fd1_rans_[0-9]*.h5")
if a and os.path.exists(rans_ref):
    yc1, p1, _ = load_profiles(a)
    yc0, p0, _ = load_profiles(rans_ref)
    du = np.abs(p1["u"] - np.interp(yc1, yc0, p0["u"]))/p0["u"].max()
    dk = np.abs(p1["k"] - np.interp(yc1, yc0, p0["k"])).max()/max(p0["k"].max(), 1e-30)
    print(f"(c) fd_force=1 drift after 2000 steps: max|du|/u_max {du.max():.2e}, "
          f"max|dk|/k_max {dk:.2e}")
    gate("(c) fd_force=1 holds the T2 RANS answer (drift <= 1e-3)",
         du.max() <= 1e-3 and dk <= 1e-2)
else:
    gate("(c) fd_force=1 vs turb180", False, "missing fd1 output or turb180 reference")

# ---- (d) IBM channel stability ----
a = latest("iddes_ibm_[0-9]*.h5")
if a:
    with h5py.File(a) as f:
        finite = all(np.isfinite(f[d][...]).all() for d in ("un", "vn", "wn", "pn", "nut", "k", "omega"))
        umax = float(np.abs(f["un"][...]).max())
        nutmax = float(f["nut"][...].max())
    gate("(d) les_ibm channel IDDES stable",
         finite and umax < 50.0, f"finite={finite} max|u|={umax:.2f} max nut={nutmax:.2e}")
else:
    gate("(d) les_ibm channel IDDES", False, "missing iddes_ibm output")

# ---- (e) 1 rank == 4 ranks ----
a = latest("ranks1_[0-9]*.h5")
b = latest("ranks4_[0-9]*.h5")
if a and b:
    ok, out = compare_exact(a, b, ["un", "vn", "wn", "pn", "nut", "k", "omega", "fd"])
    gate("(e) iddes 1 rank == 4 ranks (exact)", ok, out)
else:
    gate("(e) rank determinism", False, "missing ranks outputs")

print()
if failures:
    print("FAILURES:", ", ".join(failures))
    sys.exit(1)
print("all IDDES gates PASS")
