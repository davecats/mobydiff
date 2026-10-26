#!/usr/bin/env python3
"""Regenerate the reference data and figure of the rough-wall case.

    python3 reproduce.py

Steps:
  1. Export the span+time statistics (production_stats.h5, window from
     production_stats_t0.h5 if present) to the res-study NetCDF format
     -> assets/mobydiff/xyz_3200_384_136/data.nc (same format as ../smooth).
  2. Smooth against rough (rough_compare.py): the momentum balance per x band,
     theta / H / c_f / k+ along the plate, the Hama roughness function
     -> assets/rough_compare.txt and assets/figures/rough_comparison.png.
     Needs the smooth statistics in ../smooth.

Inputs kept LOCALLY (too large for git -- see README): production_stats.h5,
production_stats_t0.h5, grid_nodes.h5 (any field snapshot of the run has the
node lines), and the same three in ../smooth. Steps whose inputs are missing
are skipped with a note.
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SMOOTH = os.path.join(HERE, "..", "smooth")
PP = os.path.join(HERE, "assets", "postpro")
NC = os.path.join(HERE, "assets", "mobydiff", "xyz_3200_384_136", "data.nc")
FIG = os.path.join(HERE, "assets", "figures", "rough_comparison.png")
TXT = os.path.join(HERE, "assets", "rough_compare.txt")
STATS = [os.path.join(HERE, f) for f in ("production_stats.h5", "production_stats_t0.h5")]
SSTATS = [os.path.join(SMOOTH, f) for f in ("production_stats.h5", "production_stats_t0.h5")]
NODES = os.path.join(HERE, "grid_nodes.h5")
os.makedirs(os.path.dirname(NC), exist_ok=True)
os.makedirs(os.path.dirname(FIG), exist_ok=True)


def need(*paths):
    missing = [p for p in paths if not os.path.exists(p)]
    if missing:
        print("   (skip: missing", ", ".join(os.path.relpath(p, HERE) for p in missing), ")")
    return not missing


def window(files):
    """The cumulative file, plus the earlier copy when there is one."""
    return [f for f in files if os.path.exists(f)]


# 1. export in the res-study NetCDF format (the smooth case's exporter)
if need(STATS[0], NODES):
    since = ["--since", STATS[1]] if os.path.exists(STATS[1]) else []
    subprocess.run([sys.executable, os.path.join(SMOOTH, "assets", "postpro", "make_mobydiff_nc.py"),
                    STATS[0], NODES, NC, "--ini", os.path.join(HERE, "production_stats.ini"),
                    "--repo", os.path.join(HERE, "..", "..", ".."), *since], check=True)

# 2. smooth against rough; --half 50 is the balance c_f's running mean
#    (+-50 cells = 10 delta*_in, about two roughness wavelengths)
if need(STATS[0], SSTATS[0]):
    out = subprocess.run([sys.executable, os.path.join(PP, "rough_compare.py"), "--half", "50",
                          "--smooth", *window(SSTATS), "--rough", *window(STATS), "--plot", FIG],
                         check=True, capture_output=True, text=True).stdout
    with open(TXT, "w") as fh:
        fh.write(out)
    print(out)
