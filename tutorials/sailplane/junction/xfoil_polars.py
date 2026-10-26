#!/usr/bin/env python3
"""XFOIL polars of the sliced sections (step 1, handout section 3.1), XFOIL
driven directly (AeroSandbox's wrapper refuses the polar files of Drela's
6.99 source build: the CINC column has no header label and the Itr columns
are absent; the numbers are the same).  One csv per section and Ncrit:
OUT_DIR/sec_y<Y>_ncrit<N>.csv with alpha CL CD CDp CM Top_Xtr Bot_Xtr.

    python3 xfoil_polars.py SECTION_DIR OUT_DIR [--xfoil PATH] [--re 3e5] [--ncrit 9 12]

Re per section = re * chord(y) / 0.831 (root chord); alpha -4..10 by 0.5
(ASEQ from 0 upwards and from 0 downwards, so each branch starts where the
boundary layer initialises); PPAR N 240 panels.  The default XFOIL is the
source build of assets/xfoil/build.sh (the Debian binary needs an X display
and traps a floating-point exception in PACC with graphics off).
"""
import argparse, os, shutil, subprocess, sys, tempfile
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CHORD = {"0.5": 0.831, "1.0": 0.768, "2.0": 0.759, "3.0": 0.749, "5.0": 0.690, "7.0": 0.598, "9.0": 0.413, "9.6": 0.323}


def polar(xfoil, dat, re, ncrit, workdir):
    rows = []
    shutil.copy(dat, os.path.join(workdir, "af.dat"))  # XFOIL truncates long file names
    for seq in ("0 10 0.5", "-0.5 -4 -0.5"):
        pol = "pol.txt"
        if os.path.exists(os.path.join(workdir, pol)):
            os.remove(os.path.join(workdir, pol))
        keys = ["plop", "g", "", "load af.dat", "ppar", "n 240", "", "", "", "oper",
                "visc %g" % re, "mach 0", "iter 300", "vpar", "n %g" % ncrit, "", "pacc", pol, "",
                "aseq %s" % seq, "pacc", "", "quit", ""]
        run = subprocess.run([xfoil], input="\n".join(keys), capture_output=True, text=True, timeout=1200, cwd=workdir)
        if not os.path.exists(os.path.join(workdir, pol)):
            sys.exit("no polar written for %s:\n%s" % (dat, run.stdout[-1500:]))
        with open(os.path.join(workdir, pol)) as f:
            lines = f.read().splitlines()
        for l in lines[lines.index(next(x for x in lines if "alpha" in x)) + 2:]:
            w = l.split()
            if len(w) >= 7:
                rows.append([float(v) for v in w[:7]])
    rows = np.array(sorted(rows))
    return rows[np.r_[True, np.diff(rows[:, 0]) > 1e-6]]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("section_dir")
    ap.add_argument("out_dir")
    ap.add_argument("--xfoil", default=os.path.join(HERE, "assets", "xfoil", "Xfoil", "bin", "xfoil"))
    ap.add_argument("--re", type=float, default=3e5)
    ap.add_argument("--ncrit", type=float, nargs="+", default=[9.0, 12.0])
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    with tempfile.TemporaryDirectory() as wd:
        for k, c in CHORD.items():
            dat = os.path.abspath("%s/sec_y%s.dat" % (a.section_dir, k))
            re = a.re * c / 0.831
            for n in a.ncrit:
                r = polar(a.xfoil, dat, re, n, wd)
                out = "%s/sec_y%s_ncrit%g.csv" % (a.out_dir, k, n)
                np.savetxt(out, r, fmt="%.5f", delimiter=",",
                           header="alpha,CL,CD,CDp,CM,Top_Xtr,Bot_Xtr  (y = %s m, chord %.3f m, Re %.3e, Ncrit %g)" % (k, c, re, n))
                al, cl = r[:, 0], r[:, 1]
                print("sec_y%s ncrit %2g: %2d points, CL(0) %.3f CL(3) %.3f  CLmax %.3f at %.1f  Top_Xtr(3) %.2f  CD(3) %.4f"
                      % (k, n, len(al), np.interp(0, al, cl), np.interp(3, al, cl), cl.max(), al[np.argmax(cl)],
                         np.interp(3, al, r[:, 5]), np.interp(3, al, r[:, 2])))


if __name__ == "__main__":
    main()
