#!/usr/bin/env python3
"""Choose a simple 3D Cartesian MPI decomposition for a global grid."""

from __future__ import annotations

import argparse
import math
from typing import Iterable


def factor_triples(n: int) -> Iterable[tuple[int, int, int]]:
    for a in range(1, n + 1):
        if n % a != 0:
            continue
        m = n // a
        for b in range(1, m + 1):
            if m % b != 0:
                continue
            c = m // b
            yield a, b, c


def choose_dims(nx: int, ny: int, nz: int, ranks: int) -> tuple[int, int, int]:
    best: tuple[float, tuple[int, int, int]] | None = None
    for dims in factor_triples(ranks):
        dx, dy, dz = dims
        lx = nx / dx
        ly = ny / dy
        lz = nz / dz
        nondiv = int(nx % dx != 0) + int(ny % dy != 0) + int(nz % dz != 0)
        surface = ly * lz + lx * lz + lx * ly
        aspect = max(lx, ly, lz) / max(min(lx, ly, lz), 1.0e-12)
        # Prefer exact divisibility, then small halo surface, then compact blocks.
        score = nondiv * 1.0e18 + surface + aspect * 1.0e-6
        if best is None or score < best[0]:
            best = (score, dims)
    if best is None:
        raise SystemExit(f"could not factor rank count {ranks}")
    return best[1]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("nx", type=int)
    parser.add_argument("ny", type=int)
    parser.add_argument("nz", type=int)
    parser.add_argument("ranks", type=int)
    args = parser.parse_args()
    dims = choose_dims(args.nx, args.ny, args.nz, args.ranks)
    print(*dims)


if __name__ == "__main__":
    main()
