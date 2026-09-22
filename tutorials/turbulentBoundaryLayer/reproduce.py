#!/usr/bin/env python3
"""Regenerate the comparison data and figures for the shipped finey case.

    python3 reproduce.py

Steps:
  1. Export the span+time statistics (production_stats.h5) to the res-study NetCDF
     format -> assets/mobydiff/xyz_4096_224_192/data.nc (the committed mobydiff
     comparison data, produced in the SAME format as the CaNS / AMPHIBIOUS files).
  2. Code comparison figure vs SIMSON / CaNS / AMPHIBIOUS -> code_comparison.png.
  3. Single-code figures vs the SIMSON spectral reference (passivewall.hdf5).

Inputs kept LOCALLY (too large for git -- see README):
  production_stats.h5           the converged (x,y) statistics,
  restart_field.h5              a developed field (grid nodes + a snapshot),
  assets/postpro/passivewall.hdf5   the SIMSON spectral reference.
The CaNS / AMPHIBIOUS data (step 2) live on the group LSDF share; step 2 falls
back to SIMSON-only if they are not mounted. Steps that need a missing input are
skipped with a note rather than failing the whole run.
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PP = os.path.join(HERE, "assets", "postpro")
FIG = os.path.join(HERE, "assets", "figures")
NC = os.path.join(HERE, "assets", "mobydiff", "xyz_4096_224_192", "data.nc")
STATS = os.path.join(HERE, "production_stats.h5")
FIELD = os.path.join(HERE, "restart_field.h5")
PW = os.path.join(PP, "passivewall.hdf5")
os.makedirs(FIG, exist_ok=True)
os.makedirs(os.path.dirname(NC), exist_ok=True)


def run(script, *args):
    print("->", script, *map(str, args))
    subprocess.run([sys.executable, os.path.join(PP, script), *map(str, args)], check=True)


def need(*paths):
    missing = [p for p in paths if not os.path.exists(p)]
    if missing:
        print("   (skip: missing", ", ".join(os.path.relpath(p, HERE) for p in missing), ")")
    return not missing


# 1. export mobydiff data in the res-study NetCDF format
if need(STATS, FIELD):
    run("make_mobydiff_nc.py", STATS, FIELD, NC,
        "--ini", os.path.join(HERE, "production_stats.ini"), "--repo", os.path.join(HERE, "..", ".."))

# 2. code comparison (SIMSON / CaNS / AMPHIBIOUS); compare_codes falls back to
#    whatever datasets are present.
if need(NC, PW):
    run("compare_codes.py", "--mobydiff", NC, "--simson", PW,
        "--out", os.path.join(FIG, "code_comparison.png"))

# 3. SIMSON-only figures
if need(STATS, PW):
    run("compare_passivewall.py", STATS, "--ref", PW,
        "--out", os.path.join(FIG, "passivewall_compare.png"))
    run("bl_stats.py", STATS, "--plot", os.path.join(FIG, "first_stats.png"),
        "--retheta", 677, "--ref", PW)

print("done -> assets/figures/ and", os.path.relpath(NC, HERE))
# (viz_flowfield.py renders an instantaneous field but needs tools/compare_fields
#  on PYTHONPATH; run it by hand from the repo root if you want flowfield.png.)
