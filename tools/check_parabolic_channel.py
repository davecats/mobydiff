#!/usr/bin/env python3
"""Check a channel-flow HDF5 field against the Poiseuille profile.

The field must be an x/z-periodic channel with walls in y, forcing_x set,
and IBM disabled. The script reports errors against both the discrete steady
profile implied by the grid and the continuous parabola.
"""

from __future__ import annotations

import argparse
import sys

import h5py
import numpy as np


def second_derivative_matrix(y_nodes: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    n = len(y_nodes) - 1
    y = 0.5 * (y_nodes[:-1] + y_nodes[1:])
    y_ext = np.empty(n + 2)
    y_ext[1:-1] = y
    y_ext[0] = -y[0]
    y_ext[-1] = 2.0 * y_nodes[-1] - y[-1]

    mat = np.zeros((n, n))
    for j in range(n):
        c = j + 1
        hm = y_ext[c] - y_ext[c - 1]
        hp = y_ext[c + 1] - y_ext[c]
        lm = 2.0 / (hm * (hm + hp))
        lp = 2.0 / (hp * (hm + hp))
        diag = -(lm + lp)

        if j == 0:
            diag -= lm
        else:
            mat[j, j - 1] += lm

        if j == n - 1:
            diag -= lp
        else:
            mat[j, j + 1] += lp

        mat[j, j] += diag

    return y, mat


def norm_report(name: str, error: np.ndarray) -> str:
    return (
        f"{name}: max={np.max(np.abs(error)):.8e} "
        f"l2={np.sqrt(np.mean(error * error)):.8e}"
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("field", help="HDF5 field file written by the solver")
    parser.add_argument("--tolerance", type=float, default=None)
    args = parser.parse_args()

    with h5py.File(args.field, "r") as h5:
        u = h5["un"][...]
        y_nodes = h5["y"][...]
        re = float(h5.attrs["re"])
        forcing_x = float(h5.attrs["forcing_x"])

    profile = u.mean(axis=(0, 2))
    y, lap = second_derivative_matrix(y_nodes)

    rhs = np.full_like(y, -re * forcing_x)
    discrete = np.linalg.solve(lap, rhs)
    analytic = 0.5 * re * forcing_x * y * (y_nodes[-1] - y)

    discrete_error = profile - discrete
    analytic_error = profile - analytic

    print(norm_report("discrete", discrete_error))
    print(norm_report("analytic", analytic_error))

    if args.tolerance is not None and np.max(np.abs(discrete_error)) > args.tolerance:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
