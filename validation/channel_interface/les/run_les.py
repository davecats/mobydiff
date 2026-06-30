#!/usr/bin/env python3
"""LES validation -- developed-turbulence statistics campaign.

Four cases on Re_tau 180, same code/scheme, WALE + momentum_reflux=false +
interface_constant_half=true (LES cases):
  reference : uniform 128x64x128, NO LES (DNS-adequate) -- the filtered-DNS target
  uniform   : uniform 64x48x64, LES-active (the LES baseline)
  slab      : 64^3 + symmetric wall bands (flat y-interfaces, no edges)
  patch     : 64^3 + embedded core box (edges + corners)

Per case, two legs (mirrors ../developed/run_developed.py):
  1. transient (t=0..T_TRANSIENT)         -- stats OFF, no snapshots, discarded;
                                             its final field restarts leg 2.
  2. statistics (t=T_TRANSIENT..+T_AVERAGE) -- channel_stats ON (per-y mean +
     resolved stresses, time-averaged, per level), PLUS field snapshots every
     --snap-interval steps (snapshots carry `nut`; used for the nut(y) profile and,
     for patch, the edge/corner band metric).

channel_stats gives the homogeneous-case profiles (uniform/slab/reference);
snapshots give nut(y) and the patch band metric (patch_interface_stats, with the
uniform LES run as the matched base control). Analyse with les_stats.py.

IC generation (tools/make_channel_restart.py): reference/uniform = --mode base at
their resolution; slab/patch = --mode patch with the .ini's refine boxes (so the IC
leaf table is bit-identical to the solver's enumeration).

Usage (this host needs the hpcx mpirun, not /usr/bin/mpirun):
  MP=/opt/nvidia/hpc_sdk/Linux_x86_64/26.3/comm_libs/13.1/hpcx/hpcx-2.25.1/ompi/bin/mpirun
  python3 run_les.py --arch gpu --case all --mpirun "$MP"
  # then: python3 les_stats.py
"""
from __future__ import annotations
import argparse, glob, os, re, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
SRC = os.path.join(ROOT, "tutorials/channel_kmm180/channel_kmm180_restart.h5")
MKR = os.path.join(ROOT, "tools/make_channel_restart.py")

# grid per case (nx, ny, nz)
GRID = {"reference": (128, 64, 128), "uniform": (64, 48, 64),
        "slab": (64, 48, 64), "patch": (64, 48, 64)}


def sh(cmd, cwd=None):
    print("  $", " ".join(cmd))
    subprocess.run(cmd, cwd=cwd, check=True)


def refine_boxes(ini):
    out = []
    for line in open(ini):
        m = re.match(r"\s*refine\s*=\s*(.+)", line)
        if m and len(m.group(1).split()) == 6:
            out.append(m.group(1).split())
    return out


def ensure_ic(case, ini, ic):
    if os.path.isfile(ic):
        return
    if not os.path.isfile(SRC):
        sys.exit(f"source restart {SRC} not found (see tutorials/channel_kmm180)")
    nx, ny, nz = GRID[case]
    args = [sys.executable, MKR, "--nx", str(nx), "--ny", str(ny), "--nz", str(nz),
            "--source", SRC, "--out", ic]
    if case in ("reference", "uniform"):
        args += ["--mode", "base"]
    else:
        args += ["--mode", "patch"]
        for b in refine_boxes(ini):
            args += ["--refine-box", *b]
    print(f"== generating {os.path.basename(ic)}")
    sh(args)


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
    ini = os.path.join(HERE, f"{name}.ini")
    ic = os.path.join(HERE, f"{name}_IC.h5")
    ensure_ic(name, ini, ic)
    runs = os.path.join(HERE, "runs", name)
    dA, dB = os.path.join(runs, "transient"), os.path.join(runs, "stats")
    os.makedirs(dA, exist_ok=True); os.makedirs(dB, exist_ok=True)
    mpirun = a.mpirun.split() + ["-n", str(a.ranks)]
    dims = [(r"^dims = .*$", f"dims = {a.ranks} 1 1")]

    # 1. transient: stats OFF, no snapshots
    print(f"== {name}: transient t=0..{a.t_transient}")
    make_ini(ini, os.path.join(dA, "input.ini"),
             [(r"^t_final = .*$", f"t_final = {a.t_transient}"),
              (r"^stats_sample_interval = .*$", "stats_sample_interval = -1"),
              (r"^stats_write_interval = .*$", "stats_write_interval = -1"),
              (r"^field_interval = .*$", "field_interval = 0"),
              (r"^file = .*$", f"file = {ic}")] + dims)
    sh(mpirun + [binary, "input.ini"], cwd=dA)
    restart = final_field(dA, name)

    # 2. statistics: channel_stats ON (kept from the ini) + snapshots for nut/patch
    print(f"== {name}: stats t={a.t_transient}..{a.t_transient + a.t_average}, "
          f"snapshots every {a.snap_interval}")
    for stale in glob.glob(os.path.join(dB, f"{name}_*.h5")):
        os.remove(stale)
    for stale in glob.glob(os.path.join(dB, "channel_stats*.h5")):
        os.remove(stale)
    make_ini(ini, os.path.join(dB, "input.ini"),
             [(r"^t_final = .*$", f"t_final = {a.t_transient + a.t_average}"),
              (r"^field_interval = .*$", f"field_interval = {a.snap_interval}"),
              (r"^file = .*$", f"file = {restart}")] + dims)
    sh(mpirun + [binary, "input.ini"], cwd=dB)
    nsnap = len(glob.glob(os.path.join(dB, f"{name}_*.h5")))
    print(f"== {name} done: channel_stats + {nsnap} snapshots in {dB}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arch", default="gpu", choices=["gpu", "cpu"])
    ap.add_argument("--ranks", type=int, default=1)
    ap.add_argument("--case", default="all",
                    choices=["all", "reference", "uniform", "slab", "patch"])
    ap.add_argument("--t-transient", type=float, default=5.0)
    ap.add_argument("--t-average", type=float, default=20.0)
    ap.add_argument("--snap-interval", type=int, default=800,
                    help="dump a field (with nut) every N steps in the stats leg")
    ap.add_argument("--mpirun", default="mpirun",
                    help="mpirun launcher (use the hpcx full path on this host)")
    a = ap.parse_args()
    binary = os.path.join(ROOT, f"build_{a.arch}", "main")
    if not os.path.isfile(binary):
        sys.exit(f"binary {binary} not found -- ./compile.sh {a.arch}")
    order = ["reference", "uniform", "slab", "patch"]
    for name in (order if a.case == "all" else [a.case]):
        run_case(name, a, binary)
    print("\n== all done. analyse with: python3 les_stats.py")


if __name__ == "__main__":
    main()
