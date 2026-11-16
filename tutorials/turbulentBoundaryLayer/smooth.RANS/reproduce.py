#!/usr/bin/env python3
"""Regenerate the RANS-vs-DNS comparison of this tutorial.

    python3 reproduce.py

Reads the final fields of the two stages (the highest-step coarse_*.h5 and
production_*.h5 in this directory, kept locally -- see README) and the DNS
reference ../smooth/assets/matlab/tbl_dns_smooth.mat (written by
../smooth/assets/matlab/export_mat.py from the DNS statistics). Writes
  assets/rans_vs_dns.txt               the station table and c_f at matched Re_theta
  assets/figures/rans_vs_dns.png       c_f(x), c_f(Re_theta), Re_theta(Re_x), H(Re_theta)
  assets/figures/rans_profiles.png     U+, -u'v'+, k+, nu_t/nu at Re_theta 400/500/600/700
  assets/figures/rans_maps.png         log10 omega, k, nu_t/nu in (x, y), production grid
"""
import glob
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "assets", "figures")
os.makedirs(FIG, exist_ok=True)


def last(prefix):
    files = [f for f in glob.glob(os.path.join(HERE, prefix + "_*.h5"))
             if re.fullmatch(prefix + r"_\d+\.h5", os.path.basename(f))]
    if not files:
        sys.exit(f"no {prefix}_<step>.h5 here: run the stage first (README)")
    return max(files, key=lambda f: int(f.rsplit("_", 1)[1][:-3]))


prod, coarse = last("production"), last("coarse")
out = subprocess.run([sys.executable, os.path.join(HERE, "assets", "postpro", "rans_bl.py"),
                      prod, coarse, "--labels", "RANS 1600x384", "RANS 800x192",
                      "--plot", os.path.join(FIG, "rans_vs_dns.png"),
                      "--profiles", os.path.join(FIG, "rans_profiles.png"),
                      "--maps", os.path.join(FIG, "rans_maps.png")],
                     check=True, capture_output=True, text=True).stdout
out = out.replace(HERE + os.sep, "")
with open(os.path.join(HERE, "assets", "rans_vs_dns.txt"), "w") as fh:
    fh.write(out)
print(out)
