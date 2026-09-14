# The campaign re-measured after `map(to: c)` — and two published conclusions revised

Job 5142973, 4 nodes (hkn[0620-0621,0627,0632]), 38 min, the unchanged 23-run
matrix run **twice in one allocation**: `ref` = `365af76` (pre-fix), `new` =
`55bee89` (`map(to: c)`). 200 steps, `--map-by numa --bind-to core`. Raw runs in
`horeka/results_job5142973/`.

This exists because the `comm_type` parent-map fix
(`results_kernel_timeline_2026-09-11.md`) removed 6–8 ms from every step of every
blocked configuration, and **every scaling number in
`results_horeka_2026-09-10.md` §9 predated it.** This replaces them.

## 0 — The control: the `ref` column reproduces the published tables

Different nodes, four days later, a rebuilt binary:

| quantity, 16 ranks | published 2026-09-10 | `ref` here |
|---|---|---|
| block tax (1/2/4/8/16) | 1.052 / 1.095 / 1.071 / 1.059 / 1.130 | **1.051 / 1.094 / 1.070 / 1.055 / 1.130** |
| strong scaling, `base` / `rect` / `refined` / `redblack` / `big` | 71 / 66 / 49 / 44 / 81 % | **71 / 66 / 49 / 44 / 80 %** |

So the matrix is sound and the `new` column can be read as the fix and nothing
else.

## 1 — What the fix is worth, across the matrix

| config | n=1 | n=2 | n=4 | n=8 | n=16 |
|---|---|---|---|---|---|
| `base_jacobi` | +0.4 % | +1.7 | +3.4 | +5.2 | **+8.7** |
| `rect_jacobi` | +0.2 | +1.7 | +3.3 | +6.4 | **+11.2** |
| `refined_yp82_rect_jacobi` | +1.1 | +4.5 | +8.9 | +14.1 | **+18.9** |
| `refined_yp82_rect_redblack` | +1.4 | +6.3 | +11.3 | +17.5 | **+23.1** |
| `refined_big_rect_jacobi` | – | – | +2.4 | +4.5 | **+9.2** |

In absolute terms it is a constant, which is the point:

| config | n=1 | n=2 | n=4 | n=8 | n=16 |
|---|---|---|---|---|---|
| `rect_jacobi` | 1.59 ms | 6.51 | 6.65 | 6.95 | 7.59 |
| `refined_yp82_rect_jacobi` | 3.59 | 8.00 | 8.87 | 8.44 | 7.80 |
| `refined_yp82_rect_redblack` | 3.49 | 8.54 | 8.78 | 8.49 | 8.11 |
| `refined_big_rect_jacobi` | – | – | 8.36 | 8.38 | 9.94 |

**The n=1 column is the internal check.** With no peers there is no `pack` or
`unpack`, so only the 39 same-level (+24 cross-level) copy launches ever paid the
descriptor blob — about a third of the 141 that pay it from 2 ranks up. The
saving is correspondingly about a third: 1.6–3.6 ms against 6.5–9.9. Nothing was
fitted to make that come out.

## 2 — Block tax and strong scaling

| ranks | Mcell/GPU | tax ref | **tax new** |
|---|---|---|---|
| 1 | 138.41 | 1.051 | 1.053 |
| 2 | 69.21 | 1.094 | 1.095 |
| 4 | 34.60 | 1.070 | 1.070 |
| 8 | 17.30 | 1.055 | **1.043** |
| 16 | 8.65 | 1.130 | **1.100** |

Strong-scaling efficiency (against each config's smallest rank count):

| config | 1 | 2 | 4 | 8 | 16 |
|---|---|---|---|---|---|
| `base_jacobi` | 100 % | 99 | 94 | 87 | **78** (was 71) |
| `rect_jacobi` | 100 | 95 | 93 | 88 | **74** (was 66) |
| `refined_yp82_rect_jacobi` | 100 | 94 | 88 | 78 | **60** (was 49) |
| `refined_yp82_rect_redblack` | 100 | 95 | 88 | 76 | **56** (was 44) |
| `refined_big_rect_jacobi` | – | – | 100 | 96 | **86** (was 80) |

## 3 — REVISED: the 2:1 machinery is exactly free

Cost of the cells `refined_big_rect_jacobi` adds to `rect_jacobi` — same grid,
same block shape, same solver — measured against the coarse cells beside them:

| ranks | ref | **new** |
|---|---|---|
| 4 | 0.982x | **1.004x** |
| 8 | 0.943x | 0.988x |
| 16 | 0.779x | 0.826x |

The 4-rank figure is the honest one (neither twin is starved there), and it moves
from 0.98 to **1.004**: *the added fine cells cost what the coarse cells beside
them cost, to 0.4 %.* The campaign has quoted "0.98 coarse-cell-equivalents" as
the 2:1 headline since 2026-09-07; the true value was very slightly below 1
because the single-level twin was carrying more of the per-launch bill, and with
that gone the number lands on 1.

The sub-1 values at 8 and 16 ranks still mean the single-level twin is penalised
harder there, not that refinement is free-er than free.

## 4 — REVISED: red-black keeps most of its advantage

`refined_yp82_rect_redblack` / `refined_yp82_rect_jacobi`, s/step:

| binary | 1 | 2 | 4 | 8 | 16 |
|---|---|---|---|---|---|
| ref | 0.761 | 0.763 | 0.782 | 0.813 | 0.852 |
| **new** | 0.759 | **0.749** | **0.762** | **0.781** | **0.808** |

`results_horeka_2026-09-10.md` concluded "the advantage still erodes with rank
count, 0.76 → 0.86". **It erodes much less than that: 0.76 → 0.81.** Red-black's
step is the smallest in the matrix, so the fixed launch bill was the largest
share of it, and removing that recovers most of the apparent erosion. What is
left (0.759 → 0.808) is real and still worth watching, but it is half the effect
that was published.

## 5 — What did NOT move

`mpi_wait` per round, as expected of a device-side fix: `rect` 173.9 → 157.0 us
at 16 ranks, `base` 147.9 → 144.0, `refined_big` 132.8 → 101.7, while
`refined_yp82` goes 82.7 → 96.8 and red-black 81.1 → 85.8. Mixed signs and all
inside the factor this campaign has documented for waits below ~100 us.

But its SHARE grows, because the step shrank: at 16 ranks `mpi_wait` is now
**10.2 % (`base`), 10.2 % (`rect`), 11.3 % (`refined_yp82`), 12.4 %
(red-black)** of the step. It is the largest single item in the exchange again,
and — unlike everything else in this report — neither partitioning nor overlap
will touch it (2026-09-10 §3, and A0).

## 6 — What this does not cover

- **One machine, one node type, one allocation**, 200 steps.
- **`nb16_jacobi` and `blk8_jacobi` are not in this matrix** (they belong to the
  mechanism probe, not the campaign), so the fine-granularity end is unmeasured
  after the fix.
- **The per-phase breakdown is not re-measured.** Section 8 of
  `results_kernel_timeline_2026-09-11.md` estimates `jacobi_apply` at ~34 % of the
  post-fix 16-rank step by subtracting known savings from the pre-fix phase
  table. That is arithmetic, not a measurement, and it is the next thing to run.
