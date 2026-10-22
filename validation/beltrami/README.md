# Beltrami / ABC flow: spatial order and the 2:1 interface artifact

The 3D Beltrami (ABC) flow is an exact decaying solution of the
incompressible Navier-Stokes equations (`[flow] initial = beltrami`,
Re = 100, t = 2, F = exp(-nu k^2 t) = 0.9802): every component varies in
all three directions, so unlike the 2D Taylor-Green vortex it exercises
the whole stencil. `tools/check_beltrami.py` prints the L2 and Linf
velocity error against the exact field at the final time.

```bash
[BIN=/abs/moby_solve] ./run_beltrami.sh [cpu|gpu] [niter]     # niter overrides [pressure] niter
```

| case | grid | what it measures |
|---|---|---|
| `uniform.ini` (nx 32 / 64 / 128, the driver sets the size) | uniform periodic box | spatial convergence: the L2 error must fall 4x per halving (order 2) |
| `slab_x.ini` / `slab_y.ini` / `slab_z.ini` | 64^3 base, one 2:1 slab refined (levels 0/1), normal to x / y / z | the one-direction interface artifact; the three must be permutations of each other |
| `refined.ini` | 64^3 base, a central 3D patch refined (edges and corners) | the 3D interface artifact |

`slab_y.ini` is also the `beltrami_yslab` case of the 9-case bit-exactness
suite (5 steps). The inis run plain damped Jacobi at `niter = 200`
(`sor = 0.8`); that is deliberate for an exactness gate but, as the record
below shows, it is NOT a converged projection on 128^3.

## Re-validated 2026-10-22 (the first recorded result)

RE-VALIDATED 2026-10-22, commit 532ef54 (banner), production flags, workstation
CPU 1 rank (`cpu`, 53 min, the refined patch 52 min of it) and RTX 3060
(`gpu`, 29 min; the uniform 128 and refined cases 0.62 / 0.81 s/step), same
numbers on both. No earlier result was recorded for this directory (the
index row was its only description).

| case | L2(vel) | Linf(u) | note |
|---|---|---|---|
| uniform 32 | 1.076e-4 | 1.286e-4 | 86 steps |
| uniform 64 | 3.198e-5 | 6.574e-5 | 171 steps; 32 -> 64 ratio 3.36 (order 1.75) |
| uniform 128 | 3.805e-5 | 5.593e-5 | 343 steps; NO decrease: the error floor of plain Jacobi 200 |
| uniform 64, Chebyshev 200 | 2.726e-5 | 3.144e-5 | |
| uniform 128, Chebyshev 200 | 6.817e-6 | 7.870e-6 | 64 -> 128 ratio **4.00 (order 2.00)** |
| slab_x / slab_y / slab_z | 3.080e-3 (all three) | 3.40e-2 on the normal component | per-component errors permute exactly (u,v,w = 1.99/1.65/1.68e-3 cycled) |
| slab_y, Chebyshev 200 | 3.077e-3 | 3.40e-2 | the slab error is NOT the projection |
| refined patch | 2.662e-3 | 4.390e-2 | edges and corners |

FINDINGS. (1) **The convergence test is projection-bound at 128^3 with
plain damped Jacobi at 200 sweeps**: the 32 -> 64 step reads order 1.75 and
the 64 -> 128 step none at all, while the same inis with
`accel = chebyshev` at the same 200 iterations converge at order 2.00 (and
read 15 % lower at 64). The gate's statement "order 2" is true of the
spatial scheme with a converged projection; the shipped inis are left as
they are (exactness suite cases), the order claim is made with the
Chebyshev rows. (2) **The 2:1 interface error dominates everything**: a
single slab costs 100x the uniform-64 error (3.1e-3 against 3.2e-5), the
three orientations agree to four digits (the exchange is direction-
symmetric), and the converged projection does not move it -- this is the
const-1/2 transfer's own truncation error (the accepted trade of
`docs/block_refinement_strategy.md`), the same mechanism the band-refined
channels pay for in their interface divergence residual. The 3D patch
(2.7e-3) is no worse than a slab. PASS in the sense of the index row:
order 2 on the uniform grid (with a converged projection) and a bounded,
symmetric interface artifact; no tolerance had been declared, these
numbers are the record. Settings changed for the rerun: none (the driver
gained `BIN`, nx 128 and the three slabs on 2026-10-21).
