#!/usr/bin/env python3
"""Gate 5 -- DEVELOPED edge/corner banding validation for the embedded core patch.

For the CORE patch the per-y-row channel_stats are useless (the patch breaks x,z
homogeneity, so a row average smears its localised effect). Instead this driver
dumps 3D field SNAPSHOTS over the averaging window for BOTH the patch run and a
base-128 control, and tools/patch_interface_stats.py turns those into per-cell
time-mean + fluctuation-rms, matched cell-by-cell (patch and base share the base
lattice outside the patch) and aggregated by interface-adjacency class
(interior/face/edge/corner).

Two legs per case (mirrors developed/run_developed.py):
  1. transient (t = 0 .. T_TRANSIENT) -- no dumps, discarded; final field = restart
  2. stats     (t = T_TRANSIENT .. +T_AVERAGE) -- dump every --snap-interval steps

Cases: `patch` (core_patch.ini, the refined run) and `base` (base.ini, the
control). The uniform-256 reference already lives in ../reference.ini /
../run_reference.sh; the uniform-128 isolation control is this `base` case.

Usage (on the run machine):
  python3 run_gate5.py --arch gpu --ranks 2                 # both cases
  python3 run_gate5.py --arch gpu --ranks 2 --case patch    # one case
Then:
  python3 ../../../tools/patch_interface_stats.py \
      --patch 'runs/gate5/patch/stats/channel_field_*.h5' \
      --base  'runs/gate5/base/stats/base_field_*.h5'

Decomposition splits x (dims = N 1 1) so the patch stays whole on each rank and
1==2 rank is comparable; do NOT split y.
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
SRC_RESTART = os.path.join(ROOT, "tutorials/channel_kmm180/channel_kmm180_restart.h5")

# refine box of core_patch.ini (kept in sync; the IC must use the same numbers)
REFINE_BOX = ["4.76147636559703", "7.804894248762142",
              "0.3202674414922831", "1.679732558507717",
              "2.380738182798515", "3.902447124381071"]

CASES = {
    "patch": {"template": "core_patch.ini", "ic": "IC.h5", "prefix": "channel_field",
              "ic_args": ["--mode", "patch", "--refine-box", *REFINE_BOX]},
    "base":  {"template": "base.ini", "ic": "base_IC.h5", "prefix": "base_field",
              "ic_args": ["--mode", "base"]},
}


def sh(cmd, cwd=None):
    print("  $", " ".join(cmd) if isinstance(cmd, list) else cmd)
    subprocess.run(cmd, cwd=cwd, check=True, shell=isinstance(cmd, str))


def make_ini(template, dest, subs):
    with open(template) as f:
        text = f.read()
    for pattern, repl in subs:
        text = re.sub(pattern, repl, text, flags=re.MULTILINE)
    with open(dest, "w") as f:
        f.write(text)


def final_field(d, prefix):
    fields = sorted(
        (p for p in os.listdir(d) if re.match(rf"{prefix}_\d+\.h5", p)),
        key=lambda p: int(re.search(r"\d+", p).group()),
    )
    if not fields:
        sys.exit(f"no {prefix}_*.h5 produced in {d}")
    return os.path.join(d, fields[-1])


def ensure_ic(case, ic_path):
    if os.path.isfile(ic_path):
        return
    if not os.path.isfile(SRC_RESTART):
        sys.exit(f"source restart {SRC_RESTART} not found (see tutorials/channel_kmm180)")
    print(f"== generating {os.path.basename(ic_path)}")
    sh([sys.executable, os.path.join(ROOT, "tools/make_channel_restart.py"),
        *case["ic_args"], "--source", SRC_RESTART, "--out", ic_path])


def run_case(name, a, binary):
    case = CASES[name]
    template = os.path.join(HERE, case["template"])
    ic = os.path.join(HERE, case["ic"])
    ensure_ic(case, ic)

    runs = os.path.join(HERE, "runs", "gate5", name)
    dirA, dirB = os.path.join(runs, "transient"), os.path.join(runs, "stats")
    os.makedirs(dirA, exist_ok=True)
    os.makedirs(dirB, exist_ok=True)

    mpirun = ["mpirun", "-n", str(a.ranks)]
    big_nsteps = [(r"^nsteps = .*$", "nsteps = 100000000")]   # t_final governs
    dims = [(r"^dims = .*$", f"dims = {a.ranks} 1 1")]
    t_end = a.t_transient + a.t_average

    # 1. transient leg -- no dumps
    print(f"== {name}: transient leg t = 0 .. {a.t_transient}")
    make_ini(template, os.path.join(dirA, "input.ini"),
             [(r"^t_final = .*$", f"t_final = {a.t_transient}"),
              (r"^field_interval = .*$", "field_interval = 0"),
              (r"^file = .*$", f"file = {ic}")] + big_nsteps + dims)
    sh(mpirun + [binary, "input.ini"], cwd=dirA)
    restart = final_field(dirA, case["prefix"])
    print(f"   transient final field: {restart}")

    # 2. stats leg -- snapshot dumps over the averaging window
    print(f"== {name}: stats leg t = {a.t_transient} .. {t_end}, "
          f"dump every {a.snap_interval} steps")
    for stale in glob.glob(os.path.join(dirB, f"{case['prefix']}_*.h5")):
        os.remove(stale)
    make_ini(template, os.path.join(dirB, "input.ini"),
             [(r"^t_final = .*$", f"t_final = {t_end}"),
              (r"^field_interval = .*$", f"field_interval = {a.snap_interval}"),
              (r"^file = .*$", f"file = {restart}")] + big_nsteps + dims)
    sh(mpirun + [binary, "input.ini"], cwd=dirB)
    nsnap = len(glob.glob(os.path.join(dirB, f"{case['prefix']}_*.h5")))
    print(f"== {name} done: {nsnap} snapshots in {dirB}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arch", default="gpu", choices=["gpu", "cpu"])
    ap.add_argument("--ranks", type=int, default=2)
    ap.add_argument("--case", default="both", choices=["both", "patch", "base"])
    ap.add_argument("--t-transient", type=float, default=5.0)
    ap.add_argument("--t-average", type=float, default=20.0)
    ap.add_argument("--snap-interval", type=int, default=320,
                    help="dump a field every N steps during the stats leg "
                         "(~200+ snapshots over the window at dtmax)")
    a = ap.parse_args()

    binary = os.path.join(ROOT, f"build_{a.arch}", "main")
    if not os.path.isfile(binary):
        sys.exit(f"binary {binary} not found -- build with ./compile.sh {a.arch}")

    names = ["patch", "base"] if a.case == "both" else [a.case]
    for name in names:
        run_case(name, a, binary)

    print("\n== all done. compute the banding metric with:")
    print(f"   python3 {os.path.join(ROOT, 'tools/patch_interface_stats.py')} \\")
    print(f"       --patch 'runs/gate5/patch/stats/channel_field_*.h5' \\")
    print(f"       --base  'runs/gate5/base/stats/base_field_*.h5'")


if __name__ == "__main__":
    main()
