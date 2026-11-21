#!/usr/bin/env python3
"""Regenerate the shipped statistics and figures of the APG DNS.

    python3 reproduce.py

Steps:
  1. Export the window t 16 500-43 500 of the cumulative statistics to the
     res-study NetCDF format -> assets/mobydiff/xyz_3200_384_272/data.nc (the
     committed data; same format as ../smooth and the CaNS files).
  2. The read-out: assets/APG_DNS_stations.txt (dns_pg.py).
  3. The figures: dns_vs_rans_{x,prof}.png + dns_vs_rans_x.txt,
     beta_dns_rans.png (+ beta_halves.txt), fluct_retheta.png.

Inputs kept LOCALLY (too large for git, see README "Layout"):
  run/stats/stats_apg_{825000,1475000,2175000}.h5   cumulative statistics
  run/apg_2175000.h5                                 any field (node lines)
  ../apg.RANS/production_170182.h5, rans_check/check_169456.h5,
  ../smooth.RANS/production_62365.h5, ../fpg.RANS/production_170012.h5
                                                     the RANS fields
Step 1 is skipped when the statistics are absent (the committed data.nc is
then the input of 2-3); a figure whose RANS field is absent is skipped.
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
A = os.path.join(HERE, "assets")
TBL = os.path.dirname(HERE)
NC = os.path.join(A, "mobydiff", "xyz_3200_384_272", "data.nc")
STATS = os.path.join(HERE, "run", "stats", "stats_apg_%d.h5")
FIELD = os.path.join(HERE, "run", "apg_2175000.h5")
RANS_APG = [os.path.join(TBL, "apg.RANS", "production_170182.h5"),
            os.path.join(HERE, "rans_check", "check_169456.h5")]
RANS_ZPG_FPG = [os.path.join(TBL, "smooth.RANS", "production_62365.h5"),
                os.path.join(TBL, "fpg.RANS", "production_170012.h5")]


def run(script, *args, out=None):
    print("->", os.path.relpath(script, HERE), *(os.path.relpath(a, HERE) if os.path.isabs(a) else a for a in args),
          flush=True)
    r = subprocess.run([sys.executable, script, *args], check=True, capture_output=out is not None, text=True)
    if out:
        with open(out, "w") as f:
            f.write(f"# apg.DNS station table: {os.path.relpath(args[0], HERE)}, assets/dns_pg.py\n" + r.stdout)


def need(*paths):
    missing = [p for p in paths if not os.path.exists(p)]
    if missing:
        print("   (skip: missing", ", ".join(os.path.relpath(p, HERE) for p in missing), ")")
    return not missing


# 1. the committed window
if need(STATS % 2175000, STATS % 825000, FIELD):
    os.makedirs(os.path.dirname(NC), exist_ok=True)
    run(os.path.join(TBL, "smooth", "assets", "postpro", "make_mobydiff_nc.py"), STATS % 2175000, FIELD, NC,
        "--since", STATS % 825000, "--ini", os.path.join(HERE, "apg_dns.ini"), "--repo", os.path.join(TBL, "..", ".."))

# 2. the station read-out
run(os.path.join(A, "dns_pg.py"), NC, out=os.path.join(A, "APG_DNS_stations.txt"))

# 3. the figures
if need(*RANS_APG):
    run(os.path.join(A, "dns_vs_rans.py"))
    run(os.path.join(A, "beta_dns_rans.py"))
if need(*RANS_ZPG_FPG):
    run(os.path.join(A, "fluct_retheta.py"))
print("done -> assets/")
