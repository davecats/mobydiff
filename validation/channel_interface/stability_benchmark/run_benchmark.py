#!/usr/bin/env python3
"""2:1-interface STABILITY benchmark: the energy-conserving constant-1/2 interface
must keep the refined channel BOUNDED for ~250 steps where the old accuracy-optimal
interface BLOWS UP (~step 200). Fast (~2.5 min/case on one GPU).

Runs two cases of benchmark.ini with MOBY_STEPDIV on, parses div_max(step):
  energy   (interface_constant_half = true,  default) -> EXPECT bounded  (PASS)
  baseline (interface_constant_half = false)          -> EXPECT blow-up  (control)
PASS iff the energy case stays bounded AND the baseline blows up (the benchmark
discriminates). With --energy-only, only the energy case is run/checked.

Usage:
  python3 run_benchmark.py [--arch gpu|cpu] [--ranks N] [--energy-only]
Thresholds: bounded = max div_max (step>30) < 0.5 and finite; blow-up = >100 or NaN.
"""
from __future__ import annotations
import argparse
import math
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
SRC = os.path.join(ROOT, "tutorials/channel_kmm180/channel_kmm180_restart.h5")
IC = os.path.join(HERE, "IC.h5")
BOUNDED = 0.5      # div_max ceiling for "stable"
BLOWUP = 100.0     # div_max floor for "blew up"


def run_case(binary, ranks, const_half, tag):
    d = os.path.join(HERE, "runs", tag)
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(HERE, "benchmark.ini")) as f:
        ini = f.read()
    ini = re.sub(r"^dims = .*$", f"dims = {ranks} 1 1", ini, flags=re.M)
    ini = re.sub(r"^file = .*$", f"file = {IC}", ini, flags=re.M)   # absolute restart
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
    # parse div_max over steps (post-transient)
    dmax = []
    nan = False
    for line in open(log):
        m = re.search(r"STEPDIV step\s+(\d+).*div_max=\s*(\S+)", line)
        if m:
            step = int(m.group(1))
            try:
                v = float(m.group(2))
            except ValueError:
                nan = True
                continue
            if not math.isfinite(v):
                nan = True
            elif step > 30:
                dmax.append(v)
        if "NaN" in line or "nan" in line or "Inf" in line:
            nan = True
    peak = max(dmax) if dmax else float("inf")
    return dict(tag=tag, peak=peak, nan=nan, nsteps=len(dmax), log=log)


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

    print("== energy case (constant-1/2, default) -- expect BOUNDED")
    e = run_case(binary, a.ranks, True, "energy")
    e_ok = e["nsteps"] > 0 and not e["nan"] and e["peak"] < BOUNDED
    print(f"   div_max peak (step>30) = {e['peak']:.4g}, nan={e['nan']}  "
          f"-> {'STABLE ✓' if e_ok else 'NOT STABLE ✗'}")

    b_ok = True
    if not a.energy_only:
        print("== baseline case (interface_constant_half=false) -- expect BLOW-UP")
        b = run_case(binary, a.ranks, False, "baseline")
        b_ok = b["nan"] or b["peak"] > BLOWUP
        print(f"   div_max peak (step>30) = {b['peak']:.4g}, nan={b['nan']}  "
              f"-> {'BLEW UP (as expected) ✓' if b_ok else 'DID NOT BLOW UP ✗'}")

    ok = e_ok and b_ok
    print(f"\n== BENCHMARK {'PASS' if ok else 'FAIL'}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
