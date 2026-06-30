#!/usr/bin/env python3
"""RETIRED: the stability/banding metrics below parse the per-step MOBY_STEPDIV
divergence monitor, which was removed in the code cleanup. Kept as a record of
the method; the divergence peak will read inf until that hook is reinstated.

2:1-INTERFACE benchmark (~250 steps, fast): checks BOTH

  (A) STABILITY -- the energy-conserving constant-1/2 interface (default) keeps the
      refined channel BOUNDED where the old metric/cubic interface BLOWS UP
      (~step 200). Hard PASS/FAIL (+ baseline blow-up as the discriminator).
  (B) INTERFACE BANDING -- the spurious rms EXCESS the coarse cell at each 2:1
      interface develops (u' streak band; v' wall-normal STEP). Reported as TRACKED
      metrics from the energy case's final field -- the quantities a better
      interface treatment must REDUCE (they are not yet zero; the constant-1/2
      default only stabilises the band, it does not remove it).

Fast (~2.5 min/case on one GPU). Run after any interface change to confirm it stays
stable AND to see whether the band metrics moved.

Usage:
  python3 run_benchmark.py [--arch gpu|cpu] [--ranks N] [--energy-only]
PASS = energy bounded (div_max < 0.5, finite) AND baseline blows up (>100 / NaN).
Band metrics are informational (lower = better); current default values printed.
"""
from __future__ import annotations
import argparse
import glob
import math
import os
import re
import subprocess
import sys

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
SRC = os.path.join(ROOT, "tutorials/channel_kmm180/channel_kmm180_restart.h5")
IC = os.path.join(HERE, "IC.h5")
BOUNDED = 0.5      # div_max ceiling for "stable"
BLOWUP = 100.0     # div_max floor for "blew up"
IFACE_Y = (0.643, 1.357)   # the two 2:1 wall-band interfaces (coarse side)


def run_case(binary, ranks, const_half, tag):
    d = os.path.join(HERE, "runs", tag)
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(HERE, "benchmark.ini")) as f:
        ini = f.read()
    ini = re.sub(r"^dims = .*$", f"dims = {ranks} 1 1", ini, flags=re.M)
    ini = re.sub(r"^file = .*$", f"file = {IC}", ini, flags=re.M)
    flag = "true" if const_half else "false"
    ini = re.sub(r"^; interface_constant_half.*$",
                 f"interface_constant_half = {flag}", ini, flags=re.M)
    with open(os.path.join(d, "input.ini"), "w") as f:
        f.write(ini)
    env = dict(os.environ, MOBY_STEPDIV="1")
    log = os.path.join(d, "run.log")
    with open(log, "w") as out:
        subprocess.run(["mpirun", "-n", str(ranks), binary, "input.ini"],
                       cwd=d, env=env, stdout=out, stderr=subprocess.STDOUT)
    dmax, nan = [], False
    for line in open(log):
        m = re.search(r"STEPDIV step\s+(\d+).*div_max=\s*(\S+)", line)
        if m:
            try:
                v = float(m.group(2))
            except ValueError:
                nan = True; continue
            if not math.isfinite(v):
                nan = True
            elif int(m.group(1)) > 30:
                dmax.append(v)
        if re.search(r"NaN|nan|Inf", line):
            nan = True
    return dict(dir=d, peak=max(dmax) if dmax else float("inf"), nan=nan, nsteps=len(dmax))


def _level_line(base, lev):
    line = base.copy()
    for _ in range(lev):
        mid = 0.5 * (line[:-1] + line[1:])
        new = np.empty(2 * len(line) - 1); new[0::2] = line; new[1::2] = mid; line = new
    return line


def rms_rows(field):
    """Per (y-centre) rms of u',v',w' fluctuations (x,z mean removed), all levels.
    h5py block axes = [block, k(z), j(y), i(x)]. Returns sorted [(yc, lev, {u,v,w})]."""
    with h5py.File(field, "r") as f:
        nb = int(f.attrs["block_nb_x"]); yb = f["y"][...]; blocks = f["blocks"][...]
        D = {v: f[{"u": "un", "v": "vn", "w": "wn"}[v]][...] for v in "uvw"}
    rows = {}
    for bid, (ox, oy, oz, lev) in enumerate(blocks):
        for jj in range(nb):
            rows.setdefault((int(lev), int(oy) + jj), []).append(
                {v: D[v][bid, :, jj, :].ravel() for v in "uvw"})
    lines = {L: _level_line(yb, L) for L in {k[0] for k in rows}}
    recs = []
    for (lev, gj), pl in rows.items():
        yc = 0.5 * (lines[lev][gj] + lines[lev][gj + 1])
        c = {v: np.concatenate([p[v] for p in pl]) for v in "uvw"}
        recs.append((yc, lev, {v: float(np.sqrt(((c[v] - c[v].mean()) ** 2).mean())) for v in "uvw"}))
    recs.sort()
    return recs


def band_metrics(field):
    """For each 2:1 interface, the deviation of the rms profile from smooth across it:
      u' EXCESS  = q'(iface) / mean(q' of the two straddling rows)  (>1 = streak band)
      KINK q'    = |q'(i-1) - 2 q'(i) + q'(i+1)|  (2nd diff; 0 for a smooth/linear
                   profile) -- catches the v' wall-normal STEP, where the coarse
                   interface cell sits ABNORMALLY LOW between a suppressed fine side
                   and a raised coarse core (a dip-then-jump the ratio misses).
    Returns (printable summary, max u' excess, max v' kink)."""
    recs = rms_rows(field)
    ys = np.array([r[0] for r in recs])
    lines, u_exc, v_kink = [], [], []
    for yi in IFACE_Y:
        ic = int(np.argmin(np.abs(ys - yi)))                 # coarse interface row
        lo = max(ic - 1, 0); hi = min(ic + 1, len(recs) - 1)
        q = lambda v: (recs[lo][2][v], recs[ic][2][v], recs[hi][2][v])
        ue = recs[ic][2]["u"] / (0.5 * (recs[lo][2]["u"] + recs[hi][2]["u"]))
        vk = abs(q("v")[0] - 2 * q("v")[1] + q("v")[2])
        u_exc.append(ue); v_kink.append(vk)
        lines.append(
            f"   y={recs[ic][0]:.3f}: u'={recs[ic][2]['u']:.3f} (excess x{ue:.2f})  "
            f"v'=[{q('v')[0]:.3f} {q('v')[1]:.3f} {q('v')[2]:.3f}] (kink {vk:.3f})  "
            f"w'={recs[ic][2]['w']:.3f}")
    return "\n".join(lines), max(u_exc), max(v_kink)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arch", default="gpu", choices=["gpu", "cpu"])
    ap.add_argument("--ranks", type=int, default=1)
    ap.add_argument("--energy-only", action="store_true")
    a = ap.parse_args()

    binary = os.path.join(ROOT, f"build_{a.arch}", "main")
    if not os.path.isfile(binary):
        sys.exit(f"binary {binary} not found -- build with ./compile.sh {a.arch}")
    if not os.path.isfile(IC):
        if not os.path.isfile(SRC):
            sys.exit(f"source restart {SRC} not found (see tutorials/channel_kmm180)")
        print("== generating IC.h5")
        subprocess.run([sys.executable, os.path.join(ROOT, "tools/make_channel_restart.py"),
                        "--mode", "refined", "--band-cells", "24",
                        "--source", SRC, "--out", IC], check=True)

    print("== (A) STABILITY: energy case (constant-1/2, default) -- expect BOUNDED")
    e = run_case(binary, a.ranks, True, "energy")
    e_ok = e["nsteps"] > 0 and not e["nan"] and e["peak"] < BOUNDED
    print(f"   div_max peak (step>30) = {e['peak']:.4g}, nan={e['nan']}  "
          f"-> {'STABLE ✓' if e_ok else 'NOT STABLE ✗'}")

    b_ok = True
    if not a.energy_only:
        print("== (A) STABILITY: baseline (interface_constant_half=false) -- expect BLOW-UP")
        b = run_case(binary, a.ranks, False, "baseline")
        b_ok = b["nan"] or b["peak"] > BLOWUP
        print(f"   div_max peak (step>30) = {b['peak']:.4g}, nan={b['nan']}  "
              f"-> {'BLEW UP (as expected) ✓' if b_ok else 'DID NOT BLOW UP ✗'}")

    fields = sorted(glob.glob(os.path.join(e["dir"], "channel_field_*.h5")),
                    key=lambda p: int(re.search(r"\d+", os.path.basename(p)).group()))
    print("\n== (B) INTERFACE BANDING (energy case final field; lower = better, to REDUCE):")
    if fields and e_ok:
        summary, u_exc, v_kink = band_metrics(fields[-1])
        print(summary)
        print(f"   tracked: max u' excess = x{u_exc:.2f}   max v' kink = {v_kink:.3f}   (lower = better)")
    else:
        print("   (no usable field -- energy case did not finish cleanly)")

    ok = e_ok and b_ok
    print(f"\n== STABILITY {'PASS' if ok else 'FAIL'} (banding is tracked, not pass/fail)")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
