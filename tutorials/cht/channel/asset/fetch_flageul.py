#!/usr/bin/env python3
"""Fetch Flageul et al.'s RAW statistics and convert them to plain columns.

    ./fetch_flageul.py [--outdir flageul_data]

Source: https://repo.ijs.si/CFLAG/incompact3d (folder gxay), the data behind
the paper this validation compares against. It replaces the digitised curves
entirely for every quantity it covers -- and it covers far more than figure 5:
mean and rms velocity, both turbulent heat fluxes, the temperature variance in
the FLUID and in the SOLID, and a sweep of nine (G, alpha) cases.

CONVENTIONS, from their README, mapped onto ours:
    their X in gXaY = G   = fluid-to-solid thermal DIFFUSIVITY ratio -> alpha_s = 1/G
    their Y in gXaY = 1/G2, G2 = solid-to-fluid CONDUCTIVITY ratio   -> kappa_s = 1/Y
    so  C_s = kappa_s/alpha_s = G/Y   and   K = sqrt(kappa_s C_s) = sqrt(G)/Y
Their g1a1 is kappa_s = alpha_s = C_s = K = 1, i.e. exactly our k1.

FILE FORMAT: Scilab `write_csv`, which writes COMMA decimal separators and tab
field separators. Read naively with numpy it silently yields nonsense, so the
conversion is explicit here.

COLUMNS, taken from their conjug.sce (lines 494-502) rather than guessed:
    moy1   y+, U+, theta+ (referenced to the wall), dU+/dy+, dtheta+/dy+
    fluct1 y+, u'2, v'2, w'2, -u'v', u'T', -v'T', T'2          [all in wall units]
    moyb   y+ (NEGATIVE, into the bottom solid), theta+ - theta+(outer face)
    fluctb y+ (NEGATIVE, into the bottom solid), T'2
"""

from __future__ import annotations

import argparse
import os
import subprocess

RAW = "https://repo.ijs.si/CFLAG/incompact3d/-/raw/master"
CASES = ["g1a1", "g05a05", "g05a1", "g05a2", "g1a05", "g1a2",
         "g2a05", "g2a1", "g2a2"]
FILES = ["moy1", "fluct1", "moyb", "fluctb"]


def case_params(name):
    """(kappa_s, alpha_s, C_s, K) from their gXaY naming."""
    g, y = name[1:].split("a")
    val = lambda s: float(s) if "." in s else (float(s) / 10 ** (len(s) - 1)
                                               if s.startswith("0") else float(s))
    G, Y = val(g), val(y)
    kappa_s = 1.0 / Y
    alpha_s = 1.0 / G
    return kappa_s, alpha_s, kappa_s / alpha_s, (kappa_s * kappa_s / alpha_s) ** 0.5


def convert(src, dst, header):
    rows = []
    with open(src) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            rows.append([p.replace(",", ".") for p in line.split()])
    with open(dst, "w") as f:
        f.write(f"# {header}\n")
        for r in rows:
            f.write(" ".join(r) + "\n")
    return len(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("--outdir", default=os.path.join(here, "flageul_data"))
    ap.add_argument("--cases", nargs="*", default=CASES)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    print(f"  case      kappa_s   alpha_s      C_s        K")
    for c in a.cases:
        ks, als, cs, K = case_params(c)
        print(f"  {c:<8} {ks:8.3g} {als:9.3g} {cs:8.3g} {K:8.3g}", end="   ")
        got = []
        for fn in FILES:
            # g1a1 sits at the repo ROOT; gxay/g1a1 is a symlink to it, and
            # the raw endpoint will not follow that. Everything else is under
            # gxay/.
            sub = c if c == "g1a1" else f"gxay/{c}"
            url = f"{RAW}/{sub}/csv/{fn}.csv"
            tmp = os.path.join(a.outdir, f".{c}_{fn}.raw")
            r = subprocess.run(["curl", "-sSLf", "--max-time", "60", url, "-o", tmp])
            if r.returncode != 0:
                continue
            n = convert(tmp, os.path.join(a.outdir, f"{c}_{fn}.dat"),
                        f"Flageul et al., {c} (kappa_s={ks:g}, alpha_s={als:g}, "
                        f"C_s={cs:g}, K={K:g}); {fn}.csv, wall units")
            os.remove(tmp)
            got.append(f"{fn}({n})")
        print(" ".join(got))
    print(f"\n  -> {a.outdir}")


if __name__ == "__main__":
    main()
