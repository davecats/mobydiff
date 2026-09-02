# compute_phi re-profiled after the diagonal hoist, and what is left

`ncu` on istmcetus GPU 0, the same reduced case as the 2026-08-27 profile
(`~/.moby_prof/rect_small.ini`: production ini with nx/8, nz/4, same 64×44×48
block shape), so the comparison is like-for-like. ncu locks the SM clock to
1.40 GHz — read ratios, not seconds.

## Controls first

`jacobi_apply` was not touched by the diagonal hoist and reproduces:

| kernel | 2026-08-27 | now | |
|---|---|---|---|
| apply k1 | 160.19 µs | 158.98 | −0.8 % |
| apply k2 | 585.25 µs | 574.34 | −1.9 % |

(k2's small gain and its 86 → 84 registers are R1a, which moved the
`is_interface` branches out into `interface_correct`.) So the conditions are
comparable and what follows is attributable.

## The hoist did what it was aimed at

| | before | after |
|---|---|---|
| duration | 762.56 µs | **562.94 µs (−26.2 %)** |
| registers | 128 | **100** |
| Compute (SM) | 82.23 % | 74.13 % |
| DRAM | 45.74 % | 44.08 % |
| bytes/cell | 58.8 | **41.8** |
| achieved occupancy | 31.54 % | 31.21 % |

−26.2 % here, −26 % in production: **the ncu ratio transferred exactly.** Bytes
fell 17.0 B/cell against the 16 predicted (−24 `mu`, +8 `pdenom`).

Occupancy did **not** move: 100 registers still gives `Block Limit Registers`
4. It would take ≲ 80 to reach 5 blocks/SM.

`compute_denom`, the new once-per-substage kernel: 387.62 µs, 110 registers,
DRAM 51.3 %, SM 56.4 %, 33.5 B/cell. Per substage the arithmetic is
6 × 762.56 = 4575 µs before against 6 × 562.94 + 387.62 = 3765 µs after
(−18 % at niter 6, −22 % at niter 12).

## What is left is one fp64 division, and it is half the kernel

`compute_phi` is still SM-bound (74 % SM against 44 % DRAM) on a body that is
now only 6 `q` loads, 3 metric loads, 1 `pdenom` load, 1 store and ~10 fp64
operations — **one of which is a division**. On GA102 fp64 runs at a small
fraction of the fp32 rate and a division is a long software sequence.

Measured with a throwaway variant that stores the reciprocal of the diagonal
and multiplies (built, profiled, reverted — **not committed**):

| | duration | SM | DRAM | regs | occupancy |
|---|---|---|---|---|---|
| divide (committed) | 562.94 µs | 74.1 % | 44.1 % | 100 | 31.2 % |
| reciprocal-multiply | **278.85 µs** | 59.6 % | **89.1 %** | 96 | 38.2 % |

**The division alone is 50.5 % of the kernel**, and removing it flips
`compute_phi` from SM-bound to **DRAM-bound at 89 % of peak** — i.e. essentially
optimal, with nothing further to win without moving fewer bytes.

If the ratio transfers as the −26.2 % one did, the production sweep goes
0.0901 → 0.0446 s/step: **≈ 9.6 % of the niter-12 two-rank step.**

## Why it is not committed: this one is NOT bit-exact

`a/b` and `a*(1/b)` differ by up to 1 ulp in IEEE. Every bit-exact gate against
an earlier binary would fail, and all reference snapshots would need
regenerating.

The numerical argument that it is harmless: `pdenom` is the Jacobi
preconditioner diagonal, `phi` is an iterative *increment*, and the projection's
own residual after `niter` iterations is many orders of magnitude larger than a
last-bit change in the diagonal. The operator stays SPD — the diagonal remains
strictly positive; only its final bit moves. Nothing about the SPD pairing of
`face_grad_denom` / `face_grad_corr` is affected, because the correction metric
is untouched.

But it *is* a numerics change, not a refactor, so it needs an explicit decision
rather than being folded into a performance commit. Deciding to take it means
accepting that "bit-exact vs the previous binary" stops being available as a
gate for the projection, and that the replacement gate is the physics ladder
(uniform-flow exactness, 1 == 4 ranks, CPU == GPU, convergence).
