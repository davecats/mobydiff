#!/usr/bin/env python3
"""Regenerate every figure in assets/figures/ from the shipped statistics.

    python3 reproduce.py

Runs the post-processing scripts in assets/postpro/ on:
  - production_stats.h5   the converged span+time-averaged (x,y) statistics,
  - restart_field.h5      a developed instantaneous field (grid + a snapshot),
and the reference data in assets/postpro/. Writes the PNGs into assets/figures/.

Requires the SIMSON reference assets/postpro/passivewall.hdf5 (kept locally; too
large for git -- see README). Re-running the DNS itself (not this script) uses
production.ini then production_stats.ini with the red-black solver."""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PP = os.path.join(HERE, "assets", "postpro")
FIG = os.path.join(HERE, "assets", "figures")
STATS = os.path.join(HERE, "production_stats.h5")
FIELD = os.path.join(HERE, "restart_field.h5")
PW = os.path.join(PP, "passivewall.hdf5")
MAT = os.path.join(PP, "tbl_uncontrolled.mat")
PROF = os.path.join(PP, "ref_schlatter_orlu_Re670.prof")
os.makedirs(FIG, exist_ok=True)


def run(script, *args):
    print("->", script)
    subprocess.run([sys.executable, os.path.join(PP, script), *map(str, args)], check=True)


def fig(name):
    return os.path.join(FIG, name)


run("bl_stats.py", STATS, "--plot", fig("first_stats.png"), "--retheta", 677, "--ref", PW)
run("compare_passivewall.py", STATS, "--ref", PW, "--out", fig("passivewall_compare.png"))
run("compare_spectral.py", STATS, "--mat", MAT, "--out", fig("spectral_compare.png"))
run("compare_schlatter_orlu.py", STATS, "--ref", PROF, "--out", fig("schlatter_orlu_compare.png"))
run("dpdx.py", STATS, "--out", fig("dpdx.png"))
run("vonkarman.py", STATS, "--ref", PW, "--out", fig("vonkarman.png"))
run("resolution.py", STATS, FIELD, "--out", fig("resolution.png"))
run("dyplus_profiles.py", STATS, FIELD, "--out", fig("dyplus_profiles.png"))
run("viz_flowfield.py", FIELD, "--out", fig("flowfield.png"))
print("done -> assets/figures/")
