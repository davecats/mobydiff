# Pre-registered readings — the `jacobi_apply` register cut

Written **before** jobs `moby_apply` / `moby_apply_ncu` / `moby_apply16` ran, and
after the register count was already measured on the login node — so the register
number below is **not** a prediction, it is the thing the jobs are asked to
confirm on the run arch and then to price.

## What is already known, off the queue

`cuobjdump -res-usage` on a login-node build pinned to `cc80`, the same source
both sides:

| kernel | before | after |
|---|---|---|
| `jacobi_apply` k2 (face correction) | **REG:88** | **REG:64** |
| `jacobi_apply` k1 (`p += phi·idt`) | 59 | 59 |
| `jacobi_compute_phi` (control) | 94 | 94 |
| `step_momentum` (must not regress) | 128 | 128 |

`STACK:0 LOCAL:0` on every one of them, before and after: nothing spilled.

64 registers at 128 threads/block is ⌊65 536/(128·64)⌋ = **8 blocks/SM = 1 024 of
2 048 threads = 50 % occupancy**, against 5 blocks = 31 % at 88. That is the
best threshold in the handout's table, reached by the cheaper of its two levers
(L2, hoisting `face_grad_corr` out of the `collapse(4)` body) alone. L1 (splitting
the high-face plane into its own kernel) was **not** taken: it costs a launch
(~0.27 ms/step) and the threshold it was meant to reach is already crossed.

## What the jobs must return

| result | conclusion |
|---|---|
| registers 88 → 64 on the A100 build too, `STACK`/`LOCAL` 0 | the lever is real; read the rest |
| `sm__warps_active` ~29.6 % → ~45–50 %, DRAM %peak 44 → ~55–65, k2 7 654 → ~5 300–6 000 us | the occupancy hypothesis of `results_ncu_apply_2026-09-14.md` holds |
| registers fall, occupancy rises, **k2 time does not** | k2 was not occupancy-limited; the handout's line is wrong and must be closed, not extended — the cut is then a free simplification, nothing more |
| `proj_timing: apply` falls but the STEP does not | the saving is being eaten elsewhere (`mpi_wait` absorbs it); report the bucket, not a step win |
| any nonzero `max_abs` | not a refactor — the hoist is bit-exact by construction, so a difference means the hoisted value is not the one the kernel computed. Revert and find out which face. |

**Sizing.** `apply` is 33.1 % of the refined 16-rank step and k2 is 7 654 of its
9 197 us at 1 rank, i.e. ~83 % of the bucket. If k2's time falls with its DRAM
utilisation to k1's 64.8 %, the handout's estimate is 7 654 → ~5 200 us, about
**2.5 ms/step at 16 ranks = 7–8 % of the step**. Anything much larger than that
should be distrusted, not celebrated: the kernel's traffic is unchanged (1.35x
its source-counted minimum, measured), so only the rate can move.

## What these jobs cannot say

- Nothing about `interface_correct` or the exchange: neither is touched.
- Nothing about the CPU path: the hoist removes work there too, but the CPU
  build is gated for bit-exactness only, not timed.
- The 4-rank timeline and the 8/16-rank timings come from **different
  allocations**, so before/after may be compared within each job and never
  across them.
