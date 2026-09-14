# Pre-registered readings — the same hoist on the other two projection kernels

Written **before** job `moby_rdenom` ran. The register counts below come from the
login-node `cc80` loop and are not predictions; what is being predicted is
whether the time follows, which for these two kernels is a **weaker** case than
it was for `jacobi_apply` k2.

## What is already known, off the queue

| kernel | before | after | blocks/SM | theoretical occupancy |
|---|---|---|---|---|
| `compute_rdenom` | **110** | **80** | 4 → 6 | 25 % → 37.5 % |
| `jacobi_compute_phi` | **94** | **80** | 5 → 6 | 31 % → 37.5 % |
| `jacobi_apply` k2 | 64 | 64 | 8 | 50 % |
| `jacobi_apply` k1 | 59 | 59 | 8 | 50 % |

`STACK:0 LOCAL:0` throughout. The change is the same one that worked on k2:
`face_grad_denom` and the divergence metric `d1?(idx,VAR_P,b)` are functions of
the block and the face-normal index alone, so they move into the static tables
`dnLow`/`dnHigh`/`d1P` beside the existing `cfLow`/`cfHigh`. `compute_phi` gains
nothing from `dnLow`/`dnHigh` but loses `d1x`/`d1y`/`d1z` to the single `d1P`.

## Why this is a weaker case than k2, stated up front

k2 was at **42.5 %** of peak DRAM against its sibling k1's 62.4 % at the same
grid and block size — a 1.5x efficiency gap with a 1.5x occupancy gap, which is
what made the mechanism unambiguous. `compute_phi` sits at **51.9 %**, only 1.2x
below k1, and `compute_rdenom` **has never been profiled at all**. Neither has
k2's headroom. The occupancy step is also smaller: 5 → 6 blocks and 4 → 6
blocks, against k2's 5 → 8.

## What the job must return

| result | conclusion |
|---|---|
| `sweep` −5 to −12 %, `setup` −10 to −25 %, step −1 to −2.5 % | the lever generalises at a smaller size; take it and stop |
| `setup` falls but `sweep` does not | `compute_phi` at 51.9 % of peak was already near ITS limit and occupancy was not what held it — say so, keep the change for `compute_rdenom` alone, and do not extend the register line to a fourth kernel |
| neither moves | the k2 result was specific to a kernel that was far from its bandwidth limit. The register line is then closed at one success, and the remaining projection time is exchange and `mpi_wait`, not occupancy |
| `STACK`/`LOCAL` non-zero on the A100 build | spilling; revert regardless of the clock |
| any nonzero `max_abs` | a hoisted value is not the one the kernel computed. Revert and find which metric. |

**Sizing.** At 16 ranks, `sweep` is 12.1 % (`rect`) / 10.3 % (`refined`) of the
step and `setup` 3.3 % / 2.8 %. The upper end of the band above is therefore
worth ~1.8 % / 1.5 % of the step — **an order of magnitude less than the apply
cut**, and that is the honest expectation, not a disappointment to be explained
away afterwards.

## What this job cannot say

- Nothing about 16 ranks: `accelerated` is days out, and job 5145120 (the apply
  16-rank A/B) has priority on that queue.
- Nothing about `redblack_sweep` (96 registers). The same hoist would apply, but
  its `atBnd` predicate is `idx == hi(d)`, the redundant-halo sweep window, not
  `idx == nb(d)` — a different table. It is not the production smoother.
