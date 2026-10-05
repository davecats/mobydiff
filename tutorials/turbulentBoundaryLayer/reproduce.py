#!/usr/bin/env python3
"""Regenerate the comparison data and figure for the shipped CaNS-match case.

    python3 reproduce.py

Steps:
  1. Export the span+time statistics (production_stats.h5) to the res-study NetCDF
     format -> assets/mobydiff/xyz_3200_384_136/data.nc (the committed mobydiff
     comparison data, in the SAME format as the CaNS / AMPHIBIOUS files).
  2. Code comparison figure vs SIMSON / CaNS / AMPHIBIOUS -> code_comparison.png.

Inputs kept LOCALLY (too large for git -- see README):
  production_stats.h5           the converged (x,y) statistics,
  a field snapshot of the run   for the grid node lines (grid_nodes.h5, or any
                                production_p2_*.h5 / coldstart_100000.h5),
  assets/postpro/passivewall.hdf5   the SIMSON spectral reference.
The CaNS / AMPHIBIOUS data live on the group LSDF share; compare_codes.py skips
any reference whose data.nc is not mounted. Steps that need a missing input are
skipped with a note rather than failing the whole run.
"""
import glob
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PP = os.path.join(HERE, "assets", "postpro")
FIG = os.path.join(HERE, "assets", "figures")
NC = os.path.join(HERE, "assets", "mobydiff", "xyz_3200_384_136", "data.nc")
STATS = os.path.join(HERE, "production_stats.h5")
PW = os.path.join(PP, "passivewall.hdf5")
# CaNS / AMPHIBIOUS reference data on the group LSDF share (read-only). AMPHIBIOUS's
# conservative-convection datasets are not published, so the closest available curve
# is the skew-symmetric one (same centered differencing + numerics).
REF = os.path.expanduser("~/sshfsmountpoint/tbl-dns/res_study/data_processed")
CANS = os.path.join(REF, "cans", "xyz_3200_384_135", "data.nc")
AMPHI = os.path.join(REF, "amphibious", "xyz_3200_384_135_TRIP_SKEWSYM_CENTERED", "data.nc")
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


def field_with_nodes():
    """Any field snapshot of the run carries the grid node lines (x/y/z)."""
    for cand in ["grid_nodes.h5", *sorted(glob.glob(os.path.join(HERE, "production_p2_*.h5"))),
                 "coldstart_100000.h5"]:
        p = cand if os.path.isabs(cand) else os.path.join(HERE, cand)
        if os.path.exists(p):
            return p
    return os.path.join(HERE, "grid_nodes.h5")  # reported as missing by need()


# 1. export mobydiff data in the res-study NetCDF format
FIELD = field_with_nodes()
if need(STATS, FIELD):
    run("make_mobydiff_nc.py", STATS, FIELD, NC,
        "--ini", os.path.join(HERE, "production_stats.ini"), "--repo", os.path.join(HERE, "..", ".."))

# 2. code comparison (SIMSON / CaNS / AMPHIBIOUS); compare_codes skips any
#    reference whose file is absent (e.g. the LSDF share not mounted).
if need(NC):
    args = ["--mobydiff", NC]
    if os.path.exists(PW):
        args += ["--simson", PW]
    args += ["--cans", CANS, "--amphibious", AMPHI,
             "--retheta", "677", "--out", os.path.join(FIG, "code_comparison.png")]
    run("compare_codes.py", *args)

print("done -> assets/figures/code_comparison.png and", os.path.relpath(NC, HERE))
