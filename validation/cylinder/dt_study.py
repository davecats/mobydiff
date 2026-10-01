#!/usr/bin/env python3
"""Time-accuracy of the cut-cell penalization (numerics review F7).

  dt_study.py <dir with dt_{ref,new}_{0..3}_*.h5, forces_dt_*.txt and ibm_coeff_re100.h5>

Every run integrates the same clean-p Re 100 shedding state over ONE time
unit with a fixed step dt_k = 5e-3/2^k. Both penalization factors converge
to the same semi-discrete limit (du/dt = R - lambda u) as dt -> 0, so ONE
reference serves both: the finest run of the new binary. Errors are taken
over the velocity DOFs of three sets, read from the case file's coefficient
tiles: `cut` (graded coefficient: the penalized fluid-side DOFs), `fluid`
(coefficient zero) and `near` (fluid DOFs of the blocks the surface crosses).
"""
import glob
import os
import sys

import h5py
import numpy as np

d = sys.argv[1] if len(sys.argv) > 1 else "."
NAMES = ("un", "vn", "wn")


def final(side, k):
    f = sorted(glob.glob(os.path.join(d, f"dt_{side}_{k}_*.h5")),
               key=lambda p: int(p.rsplit("_", 1)[1].split(".")[0]))[-1]
    with h5py.File(f) as h:
        return np.stack([h[n][...] for n in NAMES], -1), float(h.attrs["t_current"])


with h5py.File(os.path.join(d, "ibm_coeff_re100.h5")) as h:
    coef = h["coef_blocks"][...]                      # (nBlocks, n+2, n+2, n+2, 3), ghost-inclusive
coef = coef[:, 1:-1, 1:-1, 1:-1, :]
# The tiles are stored (i,j,k) or (k,j,i); the fields are (block, k, j, i).
# The cylinder is extruded in z, so the tile axis along which coef is
# constant is z: put it first.
cutb = [b for b in range(coef.shape[0]) if (coef[b] > 0).any() and (coef[b] == 0).any()]
t = coef[cutb[0]]
if np.ptp(t, axis=0).max() != 0.0:
    assert np.ptp(t, axis=2).max() == 0.0, "no homogeneous tile axis found"
    coef = coef.transpose(0, 3, 2, 1, 4)
cut = (coef > 0) & (coef < 1e20)
fluid = coef == 0
near = np.zeros_like(fluid)
near[cutb] = fluid[cutb]
print(f"DOFs: cut {cut.sum()}, near-body fluid {near.sum()}, fluid {fluid.sum()}  "
      f"(lambda*dt_gamma at dt 5e-3: {coef[cut].min()*5e-3*2/15:.1e} .. {coef[cut].max()*5e-3*8/15:.1e})")

ref, tref = final("new", 3)
print(f"reference: new, dt = {5e-3/8:.3e}, t = {tref:.6f}")
print(f"{'side':4s} {'dt':>9s} | {'cut max':>9s} {'cut rms':>9s} | {'near rms':>9s} | {'fluid rms':>9s} {'fluid max':>9s} | C_D(T)")
rows = {}
for side in ("ref", "new"):
    for k in range(4):
        if side == "new" and k == 3:
            continue
        q, t = final(side, k)
        assert abs(t - tref) < 1e-9, f"{side} {k} ends at t = {t}, reference at {tref}"
        e = np.abs(q - ref)
        fr = np.loadtxt(os.path.join(d, f"forces_dt_{side}_{k}.txt"), skiprows=1)
        rows[side, k] = (e[cut].max(), np.sqrt((e[cut]**2).mean()), np.sqrt((e[near]**2).mean()),
                         np.sqrt((e[fluid]**2).mean()), e[fluid].max(), fr[-1, 3])
        r = rows[side, k]
        print(f"{side:4s} {5e-3/2**k:9.3e} | {r[0]:9.2e} {r[1]:9.2e} | {r[2]:9.2e} | {r[3]:9.2e} {r[4]:9.2e} | {r[5]:.8f}")
print("observed order log2(e(dt)/e(dt/2)):   cut rms   near rms   fluid rms")
for side in ("ref", "new"):
    for k in range(3 if side == "ref" else 2):
        if (side, k + 1) not in rows:
            continue
        a, b = rows[side, k], rows[side, k + 1]
        print(f"  {side:3s} {5e-3/2**k:9.3e} -> {5e-3/2**(k+1):9.3e}    "
              f"{np.log2(a[1]/b[1]):6.2f}   {np.log2(a[2]/b[2]):7.2f}   {np.log2(a[3]/b[3]):8.2f}")
print("error ratio ref/new at equal dt:        cut rms   near rms   fluid rms")
for k in range(3):
    a, b = rows["ref", k], rows["new", k]
    print(f"  dt = {5e-3/2**k:9.3e}                     {a[1]/b[1]:7.2f}   {a[2]/b[2]:8.2f}   {a[3]/b[3]:9.2f}")
