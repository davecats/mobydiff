#!/usr/bin/env python3
"""Uniform flow + a 1e-6 random perturbation, written into a solver-minted
template snapshot (a one-step run of the same ini): seed_ic.py <template> <out> U V W"""
import shutil, sys
import h5py, numpy as np
tmpl, out = sys.argv[1], sys.argv[2]
vel = [float(v) for v in sys.argv[3:6]]
shutil.copy(tmpl, out)
rng = np.random.default_rng(1)
with h5py.File(out, "r+") as f:
    for name, u in zip(("un", "vn", "wn"), vel):
        f[name][...] = u + 1.0e-6*(rng.random(f[name].shape) - 0.5)
    f["pn"][...] = 0.0
