#!/usr/bin/env python3
"""Build the y node line for the FLAGEUL-MATCHED conjugate channel.

    ./make_ynodes.py [--out ynodes_f149.dat]

The built-in distributions all cluster at the DOMAIN ENDS. This case needs the
opposite: clustering at the two INTERIOR fluid/solid interfaces, which must
also land EXACTLY on cell faces (the conjugate face coefficient is a same-level
arm between two cell centres straddling the interface -- an interface inside a
cell is not representable). Hence a file line, read by `[grid.y] nodes_file`.

Three pieces, C0-matched at the interfaces:
  solid [0, d]      n_s cells, geometric, FINEST at the interface and growing
                    outward -- the interface is where the gradients are.
  fluid [d, d+2]    n_f cells, symmetric tanh: dy+ = 0.49 at both walls to
                    4.8 at the centreline, i.e. Flageul et al.'s spacing.
  solid [d+2, 2d+2] the mirror of the first.

The tanh family fixes N once both spacings are given: the centre/wall ratio is
cosh^2(a), and then dy_centre = L a / (N tanh a). So n_f is not free -- 118
points is what 0.49 and 4.8 ask for, and asking for more would change the
spacings we were told to match.
"""

from __future__ import annotations

import argparse

import numpy as np
from scipy.optimize import brentq

RE = 149.0
DYW_PLUS, DYC_PLUS = 0.49, 4.8       # Flageul et al., table 1
LF, D = 2.0, 1.0                     # fluid gap; solid thickness (d+ = 149)
NS, NF = 53, 118                     # -> ny = 224 = 7 x nb(32)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="ynodes_f149.dat")
    a = ap.parse_args()

    dyw, dyc = DYW_PLUS / RE, DYC_PLUS / RE

    # Solve alpha on the REALISED first cell, not on the analytic derivative:
    # the two differ by ~3 % at this stretching, and the realised value is
    # what sets both the resolution and the diffusive time-step limit.
    j = np.arange(NF + 1)

    def line(al):
        return D + 0.5 * LF * (1.0 + np.tanh(al * (2.0 * j / NF - 1.0)) / np.tanh(al))

    alpha = brentq(lambda al: (line(al)[1] - line(al)[0]) - dyw, 0.2, 4.0)
    yf = line(alpha)
    dw = yf[1] - yf[0]                     # the realised wall cell

    # ---- solid: geometric, first cell = dw at the interface ----------------
    r = brentq(lambda r: dw * (r ** NS - 1.0) / (r - 1.0) - D, 1.0 + 1e-9, 1.5)
    widths = dw * r ** np.arange(NS)       # from the interface outward
    ys = D - np.concatenate([[0.0], np.cumsum(widths)])[::-1]   # 0 .. D

    node = np.concatenate([ys[:-1], yf, (2 * D + LF) - ys[::-1][1:]])
    # Pin the interfaces and the ends EXACTLY: everything downstream compares
    # against these numbers, and the file carries printed decimals.
    node[0], node[-1] = 0.0, 2 * D + LF
    node[NS], node[NS + NF] = D, D + LF

    n = node.size - 1
    dy = np.diff(node)
    assert n == 2 * NS + NF, n
    assert np.all(dy > 0)
    with open(a.out, "w") as f:
        f.write(f"# y node line, {n} cells, Flageul-matched conjugate channel\n")
        f.write(f"# solid {NS} + fluid {NF} + solid {NS};  interfaces at "
                f"y = {D} (index {NS}) and {D+LF} (index {NS+NF})\n")
        f.write(f"# dy+ : wall {dw*RE:.4f}  centre {dy[NS+NF//2]*RE:.4f}  "
                f"solid outer {widths[-1]*RE:.3f}  (growth {r:.5f})\n")
        for v in node:
            f.write(f"{v:.17g}\n")
    print(f"{a.out}: {n} cells, ny = {n}")
    print(f"   fluid  dy+ wall {dw*RE:.4f}   centre {dy[NS+NF//2]*RE:.4f}"
          f"   (asked {DYW_PLUS} / {DYC_PLUS})")
    print(f"   solid  growth {r:.5f}   outer cell dy+ {widths[-1]*RE:.3f}"
          f"   d+ = {D*RE:.0f}")
    print(f"   interfaces at index {NS} -> {node[NS]!r}  and {NS+NF} -> {node[NS+NF]!r}")
    print(f"   max cell-to-cell growth ratio {np.max(dy[1:]/dy[:-1]):.4f}")


if __name__ == "__main__":
    main()
