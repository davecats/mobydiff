#!/usr/bin/env python3
"""The factor that makes the shell sink balance THIS grid's generation.

The sink is right when `C_s |solid_source| V_shell = int source*w dV` over the
fluid, and both integrals are DISCRETE: the generation is what the transport
kernel actually applies, and the shell volume is a cell count. Campaign 1
measured what a 2.8 % error here does -- nothing to the mean profile, but it
dominates the time-averaged variance. Prints the factor to multiply the
current `solid_source` values by.

    ./measure_generation.py <snapshot> <case file> --current 2.21110542
"""
import argparse, sys
import h5py, numpy as np
sys.path.insert(0, __file__.rsplit("/", 1)[0])
from pipe_stats import load_field, load_solid_mask

ap = argparse.ArgumentParser()
ap.add_argument("snapshot"); ap.add_argument("case")
ap.add_argument("--current", type=float, required=True,
                help="the C_s*|solid_source| the run used")
ap.add_argument("--depth", type=float, default=0.1)
ap.add_argument("--centre", type=float, nargs=2, default=(0.65, 0.65))
a = ap.parse_args()

solid = load_solid_mask(a.case)
with h5py.File(a.case, "r") as h:
    node = {d: h[f"{d}_nodes"][...] for d in "xyz"}
cent = {d: 0.5 * (node[d][:-1] + node[d][1:]) for d in "xyz"}
dV = float(np.prod([node[d][1] - node[d][0] for d in "xyz"]))
Z, Y, X = np.meshgrid(cent["z"], cent["y"], cent["x"], indexing="ij")
r = np.hypot(X - a.centre[0], Y - a.centre[1]); del X, Y, Z
Vs = float((solid & (r < 0.5 + a.depth)).sum()) * dV
with h5py.File(a.snapshot, "r") as h:
    gen = float((load_field(h, "wn") * ~solid).sum()) * dV
need = gen / Vs
print(f"{need / a.current:.6f}")
print(f"# int w dV = {gen:.4f}, V_shell = {Vs:.4f} -> C_s|solid_source| ="
      f" {need:.7f} (was {a.current:.7f})", file=sys.stderr)
