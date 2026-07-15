#!/usr/bin/env python3
"""Run the DEVELOPED 2:1 wall-band channel and collect time-averaged statistics.

Two-leg driver (mirrors validation/channel_interface/run_validation.sh, in Python):
  1. transient leg  (t = 0 .. T_TRANSIENT)  -- statistics OFF, discarded; its
     final field is the restart for leg 2.
  2. statistics leg (t = T_TRANSIENT .. T_TRANSIENT+T_AVERAGE) -- statistics ON,
     accumulated FRESH into channel_stats.h5 (+ _l1) in its own directory.

The energy-conserving constant-1/2 interface is the default. Pass --skew to also
add the skew band correction (interface_skew = true). Decomposition splits x
(dims = N 1 1) so the y wall-bands stay whole on each rank (see README.md).

Usage:
  python3 run_developed.py [--arch gpu|cpu] [--ranks N] [--skew]
                           [--t-transient 5] [--t-average 20] [--name default]
Then:
  python3 ../../../tools/plot_channel_stats.py stats.png \
      runs/<name>/stats/channel_stats.h5:<label>
"""
from __future__ import annotations
import argparse
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
SRC_RESTART = os.path.join(ROOT, "tutorials/channel_kmm180/channel_kmm180_restart.h5")
IC = os.path.join(HERE, "IC.h5")


def sh(cmd, cwd=None):
    print("  $", " ".join(cmd) if isinstance(cmd, list) else cmd)
    subprocess.run(cmd, cwd=cwd, check=True, shell=isinstance(cmd, str))


def make_ini(template, dest, subs):
    """Copy `template` to `dest` applying line-anchored key replacements."""
    with open(template) as f:
        text = f.read()
    for pattern, repl in subs:
        text = re.sub(pattern, repl, text, flags=re.MULTILINE)
    with open(dest, "w") as f:
        f.write(text)


def final_field(d):
    fields = sorted(
        (p for p in os.listdir(d) if re.match(r"channel_field_\d+\.h5", p)),
        key=lambda p: int(re.search(r"\d+", p).group()),
    )
    if not fields:
        sys.exit(f"no channel_field_*.h5 produced in {d}")
    return os.path.join(d, fields[-1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arch", default="gpu", choices=["gpu", "cpu", "gpu_corax"],
                    help="binary build dir suffix (gpu_corax = the cc120 build on istmcorax)")
    ap.add_argument("--ranks", type=int, default=2)
    ap.add_argument("--skew", action="store_true", help="also add interface_skew = true")
    ap.add_argument("--no-reflux", action="store_true",
                    help="set momentum_reflux = false (the reflux is the u'/v' band source -- "
                         "see docs/next_session_edges_les.md)")
    ap.add_argument("--t-transient", type=float, default=5.0)
    ap.add_argument("--t-average", type=float, default=20.0)
    ap.add_argument("--refine-dims", default="xyz", choices=["xyz", "xz"],
                    help="xz: wall bands refined in x,z only (R2D-2 gate; y keeps "
                         "the stretched base line at level 1; own IC_xz.h5)")
    ap.add_argument("--name", default=None, help="run subdir name (default: default / skew / noreflux)")
    a = ap.parse_args()

    binary = os.path.join(ROOT, f"build_{a.arch}", "main")
    if not os.path.isfile(binary):
        sys.exit(f"binary {binary} not found -- build with ./compile.sh {a.arch}")
    xz = a.refine_dims == "xz"
    name = a.name or ("xz" if xz else "noreflux" if a.no_reflux else "skew" if a.skew else "default")
    t_end = a.t_transient + a.t_average
    skew_sub = [(r"^; interface_skew.*$", "interface_skew = true")] if a.skew else []
    reflux_sub = [(r"^momentum_reflux = .*$", "momentum_reflux = false")] if a.no_reflux else []
    # xz quadtree: the wall bands refine in x,z only (the IC leaf table and
    # the solver builder must agree, so the IC is refine_dims-specific).
    dims_key_sub = [(r"^refine_levels = 1$", "refine_levels = 1\nrefine_dims = xz")] if xz else []
    ic = os.path.join(HERE, "IC_xz.h5") if xz else IC
    mpirun = ["mpirun", "-n", str(a.ranks)]
    dims_sub = [(r"^dims = .*$", f"dims = {a.ranks} 1 1")]

    # 0. initial condition
    if not os.path.isfile(ic):
        if not os.path.isfile(SRC_RESTART):
            sys.exit(f"source restart {SRC_RESTART} not found (see tutorials/channel_kmm180)")
        print(f"== generating {os.path.basename(ic)}")
        sh([sys.executable, os.path.join(ROOT, "tools/make_channel_restart.py"),
            "--mode", "refined", "--band-cells", "24",
            "--refine-dims", a.refine_dims,
            "--source", SRC_RESTART, "--out", ic])

    runs = os.path.join(HERE, "runs", name)
    dirA, dirB = os.path.join(runs, "transient"), os.path.join(runs, "stats")
    os.makedirs(dirA, exist_ok=True)
    os.makedirs(dirB, exist_ok=True)

    # 1. transient leg (stats off)
    print(f"== {name}: transient leg t = 0 .. {a.t_transient}")
    make_ini(os.path.join(HERE, "transient.ini"), os.path.join(dirA, "input.ini"),
             [(r"^t_final = .*$", f"t_final = {a.t_transient}"),
              (r"^file = .*$", f"file = {ic}")] + dims_sub + skew_sub + reflux_sub + dims_key_sub)
    sh(mpirun + [binary, "input.ini"], cwd=dirA)
    restart = final_field(dirA)
    print(f"   transient final field: {restart}")

    # 2. statistics leg (stats on, FRESH dir so the accumulator starts at zero)
    for stale in ("channel_stats.h5", "channel_stats_l1.h5"):
        p = os.path.join(dirB, stale)
        if os.path.exists(p):
            os.remove(p)
    print(f"== {name}: statistics leg t = {a.t_transient} .. {t_end}")
    make_ini(os.path.join(HERE, "developed.ini"), os.path.join(dirB, "input.ini"),
             [(r"^t_final = .*$", f"t_final = {t_end}"),
              (r"^file = .*$", f"file = {restart}")] + dims_sub + skew_sub + reflux_sub + dims_key_sub)
    sh(mpirun + [binary, "input.ini"], cwd=dirB)

    stats = os.path.join(dirB, "channel_stats.h5")
    print(f"\n== done. stats: {stats}")
    print("   plot with:")
    print(f"   python3 {os.path.join(ROOT, 'tools/plot_channel_stats.py')} "
          f"{name}_stats.png {stats}:{name}")


if __name__ == "__main__":
    main()
