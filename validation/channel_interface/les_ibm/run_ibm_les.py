#!/usr/bin/env python3
"""LES<->IBM coupling validation -- run driver (designed to run on ANY machine).

A plane-wall channel whose two flat walls are described by the FILE-BASED immersed
boundary and do NOT coincide with grid nodes (uniform y, walls mid-cell). Re_tau
~180. Tests the ibm_aware solid-cell nut masking in src/modules/les.f90 -- the only
untested piece of the LES path (validated for grid-aligned walls in ../les/).

Cases:
  a_wale  : single-level, IBM wall, WALE LES         (the coupling under test)
  b_none  : single-level, IBM wall, LES OFF          (the control for the mean flow)
  c_refine: refine_body at the IBM wall, WALE LES     (2:1-interface x IBM x LES)

Per case, two legs (mirrors ../les/run_les.py):
  1. transient (t=0..T_TRANSIENT)            stats OFF, no snapshots -- discarded.
  2. statistics (t=T_TRANSIENT..+T_AVERAGE)  channel_stats ON (per-y mean + stresses,
     per level) PLUS field snapshots every --snap-interval steps (snapshots carry
     `nut`, used for the wall-nut measurement and the 2:1-interface nut step).

All prerequisite data files (grid.h5, ibm_coeff.h5, ibm_coeff_blocks.h5, IC.h5,
IC_refine.h5) are committed alongside this script -- they were built with the
geometry venv + mobygrid (see setup.sh / README). This driver needs only the
solver binary, numpy is not required to RUN (only to ANALYSE).

Usage (pick the right mpirun for the host; on the dev box it is the hpcx one):
  MP=/opt/.../hpcx-2.25.1/ompi/bin/mpirun        # or just "mpirun" elsewhere
  python3 run_ibm_les.py --arch gpu --case all --mpirun "$MP"
  python3 measure_nut.py        # gates 1-2 (solid nut==0, no spurious wall spike)
  python3 ibm_les_stats.py      # gate 3 (mean-U log law), gate 4 (interface nut)
"""
from __future__ import annotations
import argparse
import glob
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))

# (ini template, restart IC, coefficient file) per case. Paths relative to HERE.
CASES = {
    "a_wale":   ("channel_ibm.ini",        "IC.h5",        "ibm_coeff.h5"),
    "b_none":   ("channel_ibm.ini",        "IC.h5",        "ibm_coeff.h5"),
    "c_refine": ("channel_ibm_refine.ini", "IC_refine.h5", "ibm_coeff_blocks.h5"),
}


def sh(cmd, cwd=None):
    print("  $", " ".join(cmd), flush=True)
    subprocess.run(cmd, cwd=cwd, check=True)


def make_ini(template, dest, subs):
    t = open(template).read()
    for pat, repl in subs:
        t = re.sub(pat, repl, t, flags=re.MULTILINE)
    open(dest, "w").write(t)


def final_field(d, prefix):
    fs = sorted((p for p in os.listdir(d) if re.match(rf"{prefix}_\d+\.h5", p)),
                key=lambda p: int(re.search(r"\d+", p).group()))
    if not fs:
        sys.exit(f"no {prefix}_*.h5 in {d}")
    return os.path.join(d, fs[-1])


def run_case(name, a, binary):
    template, ic_name, coef_name = CASES[name]
    ini = os.path.join(HERE, template)
    ic = os.path.join(HERE, ic_name)
    coef = os.path.join(HERE, coef_name)
    for p in (ini, ic, coef):
        if not os.path.isfile(p):
            sys.exit(f"missing {p} -- run setup.sh first (needs geometry venv + mobygrid)")
    prefix = "channel_ibm"  # field_prefix in the .ini
    runs = os.path.join(HERE, "runs", name)
    dA, dB = os.path.join(runs, "transient"), os.path.join(runs, "stats")
    os.makedirs(dA, exist_ok=True)
    os.makedirs(dB, exist_ok=True)
    mpirun = a.mpirun.split() + ["-n", str(a.ranks)]
    common = [(r"^dims = .*$", f"dims = {a.ranks} 1 1"),
              (r"^coeff_file = .*$", f"coeff_file = {coef}")]
    model = [(r"^model = .*$", "model = none")] if name == "b_none" else []

    # 1. transient: stats OFF, no snapshots
    print(f"== {name}: transient t=0..{a.t_transient}", flush=True)
    make_ini(ini, os.path.join(dA, "input.ini"),
             [(r"^t_final = .*$", f"t_final = {a.t_transient}"),
              (r"^stats_sample_interval = .*$", "stats_sample_interval = -1"),
              (r"^stats_write_interval = .*$", "stats_write_interval = -1"),
              (r"^field_interval = .*$", "field_interval = 0"),
              (r"^file = .*$", f"file = {ic}")] + common + model)
    sh(mpirun + [binary, "input.ini"], cwd=dA)
    restart = final_field(dA, prefix)

    # 2. statistics: channel_stats ON + nut snapshots
    print(f"== {name}: stats t={a.t_transient}..{a.t_transient + a.t_average}, "
          f"snapshots every {a.snap_interval}", flush=True)
    for stale in glob.glob(os.path.join(dB, f"{prefix}_*.h5")):
        os.remove(stale)
    for stale in glob.glob(os.path.join(dB, "channel_ibm_stats*.h5")):
        os.remove(stale)
    make_ini(ini, os.path.join(dB, "input.ini"),
             [(r"^t_final = .*$", f"t_final = {a.t_transient + a.t_average}"),
              (r"^field_interval = .*$", f"field_interval = {a.snap_interval}"),
              (r"^file = .*$", f"file = {restart}")] + common + model)
    sh(mpirun + [binary, "input.ini"], cwd=dB)
    nsnap = len(glob.glob(os.path.join(dB, f"{prefix}_*.h5")))
    print(f"== {name} done: channel_stats + {nsnap} snapshots in {dB}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arch", default="gpu", choices=["gpu", "cpu"])
    ap.add_argument("--ranks", type=int, default=1)
    ap.add_argument("--case", default="all",
                    choices=["all", "a_wale", "b_none", "c_refine"])
    ap.add_argument("--t-transient", type=float, default=5.0)
    ap.add_argument("--t-average", type=float, default=20.0)
    ap.add_argument("--snap-interval", type=int, default=800,
                    help="dump a field (with nut) every N steps in the stats leg")
    ap.add_argument("--mpirun", default="mpirun")
    a = ap.parse_args()
    binary = os.path.join(ROOT, f"build_{a.arch}", "main")
    if not os.path.isfile(binary):
        sys.exit(f"binary {binary} not found -- ./compile.sh {a.arch}")
    order = ["a_wale", "b_none", "c_refine"]
    for name in (order if a.case == "all" else [a.case]):
        run_case(name, a, binary)
    print("\n== all done. analyse with: python3 measure_nut.py ; python3 ibm_les_stats.py")


if __name__ == "__main__":
    main()
