# `jacobi_apply` k2: 88 registers became 64, and the time followed

Jobs 5145099 (registers + 4-rank A/B + every gate), 5145100 (ncu, both binaries,
same node), 5145114 (8-rank A/B, dev partition) and 5145120 (8- and 16-rank A/B,
4 nodes). A100-SXM4-40GB, HoreKa.
`ref` = `e07b8d3`, `new` = `25c13ab`. Raw runs in `horeka/results_job5145099/`,
`results_job5145100/`, `results_job5145114/`, `results_job5145120/`. The readings were written down
first, in `horeka/exchange/PREREGISTERED_apply.md`.

`results_ncu_apply_2026-09-14.md` ended with a hypothesis and no plan: k2 is
latency-bound because 88 registers/thread leave only 5 blocks/SM resident, and
the levers for cutting them were "an open question". **The hypothesis holds, and
the cheaper of the two levers was enough.**

## 1 — The change

`face_grad_corr(face kind, at-boundary, d1f, outlet, refined)` depends on the
block and the index along the face NORMAL alone — not on the two tangential
indices the `collapse(4)` body also spans. Evaluating it inside the kernel
recomputed each value nb² times and, the part that cost something, kept ten array
bases live across the whole body: `physLow`, `physHigh`, `d1x`, `d1y`, `d1z`,
`outLow`, `outHigh`, `refd` on top of `phi`, `q` and `mu`.

So it is precomputed: `cfLow(idx,d,b)` for the low faces (carrying the interface
zeroing) and `cfHigh(d,b)` for the one high face this kernel owns, the outlet.
Both are static — face kinds come from the leaf table, metrics from the node
lines — so unlike `rdenom` they are formed once on the host and mapped once.
The kernel reads one double per face. **Bit-exact by construction: identical
expression, identical inputs, evaluated once instead of nb² times.**

The handout's L1 (splitting the high-face planes into their own kernel) was
**not** taken: it costs a launch (~0.27 ms/step) and the threshold it was meant
to reach is already crossed without it.

## 2 — Registers (job 5145099, `cuobjdump -res-usage` on the A100 build)

| kernel | ref | new |
|---|---|---|
| `jacobi_apply` k2 (face correction) | **88** | **64** |
| `jacobi_apply` k1 (`p += phi·idt`) | 59 | 59 |
| `jacobi_compute_phi` (control) | 94 | 94 |
| `compute_rdenom` | 110 | 110 |
| `interface_correct` ×3 | 54 | 54 |
| `redblack_sweep` | 96 | 96 |
| `step_momentum` (the solver's largest) | 128 | 128 |

`STACK:0 LOCAL:0` on every kernel on both sides: **nothing spilled**, which is
the one way this could have been a loss dressed as a win.

64 registers at 128 threads/block is ⌊65 536/(128·64)⌋ = 8 blocks/SM = 1 024 of
2 048 threads, against 5 blocks = 640 at 88.

## 3 — ncu: occupancy and DRAM utilisation followed, traffic did not move

Job 5145100, hkn0401, `refined_yp82_rect_jacobi` at 1 rank (60.56 Mcell), both
binaries on the same node in the same allocation — the 2026-09-14 numbers were
taken elsewhere, so the BEFORE side is re-measured rather than quoted across
allocations.

| kernel | us ref | **us new** | DRAM %peak | occupancy % | SM %peak | doubles/cell |
|---|---|---|---|---|---|---|
| `jacobi_compute_phi` (control) | 3 181 | 3 181 | 51.9 → 51.9 | 29.1 → 29.1 | 47.1 → 47.1 | 5.30 → 5.30 |
| `jacobi_apply` k1 (control) | 1 536 | 1 536 | 62.4 → 62.5 | 45.5 → 45.4 | 66.2 → 66.2 | 3.08 → 3.08 |
| `jacobi_apply` k2 | 7 565 | **5 309** | **42.5 → 60.5** | **29.6 → 46.0** | 32.1 → 40.6 | **10.31 → 10.31** |

**−29.8 % on k2, with the bytes unchanged to the last digit printed.** That is
the cleanest possible confirmation of the mechanism: the kernel moves exactly what
it moved before, and moves it faster because 1.56x as many warps are resident.
Both control kernels are unmoved to ≤0.1 %, on the same node, in the same run.

Achieved occupancy is 46.0 %, not the theoretical 50 % — tail and launch effects
— and DRAM utilisation lands at 60.5 %, just under k1's 62.5 %. The 2026-09-14
report's estimate ("if k2 reaches 60 %: 7 654 → ~5 600 us") was right, and
slightly conservative.

## 4 — Step times: the saving lands on the step

Two allocations, before/after inside each. 100 steps, `--map-by numa
--bind-to core`, `[output] profile = true`.

| job | case | step ms | apply ms | sweep ms (control) |
|---|---|---|---|---|
| 5145099 | `rect_jacobi` 4×1 | 194.01 → 177.39 **−8.6 %** | 80.38 → 62.89 **−21.8 %** | 27.86 → 28.14 +1.0 % |
| 5145099 | `refined_yp82` 4×1 | 91.33 → 83.80 **−8.2 %** | 37.26 → 29.53 **−20.8 %** | 12.41 → 12.47 +0.5 % |
| 5145114 | `rect_jacobi` 4×1 | 194.21 → 177.51 **−8.6 %** | 80.49 → 62.93 **−21.8 %** | 27.90 → 28.17 +1.0 % |
| 5145114 | `refined_yp82` 4×1 | 91.33 → 83.89 **−8.2 %** | 37.28 → 29.50 **−20.9 %** | 12.41 → 12.47 +0.5 % |
| 5145114 | `rect_jacobi` 8×2 | 102.40 → 93.71 **−8.5 %** | 40.57 → 31.74 **−21.8 %** | 14.14 → 14.21 +0.5 % |
| 5145114 | `refined_yp82` 8×2 | 51.79 → 47.74 **−7.8 %** | 19.91 → 15.90 **−20.1 %** | 6.47 → 6.45 −0.3 % |

Three controls, all passed:

- **`sweep` does not move** (+1.0 % to −0.3 % across six pairs). It shares the
  binary, the node and the allocation with `apply` and is untouched by the change.
- **The two allocations agree.** The 4-rank rows were run on different nodes
  hours apart and reproduce each other to 0.1 %.
- **The committed matrix agrees with the `before` column.** Job 5142973's
  post-`map(to: c)` logs give `apply` = 80.70 / 37.30 ms at 4 ranks and 40.60 /
  19.91 at 8 — our `before` to within 0.4 %.

**The saving on `apply` equals the saving on the step**, to within the noise
(4 ranks: 17.5 vs 16.6 ms `rect`, 7.7 vs 7.5 `refined`; 8 ranks: 8.8 vs 8.7 and
4.0 vs 4.1). Nothing absorbs it.

## 5 — Unlike the exchange fix, this saving is per-cell, not per-launch

| config | 4 ranks | 8 ranks |
|---|---|---|
| `rect_jacobi`, `apply` saved | 17.5 ms/step | 8.8 |
| `refined_yp82`, `apply` saved | 7.7 | 4.0 |

It **halves with the rank count**, because it is a rate improvement on a fixed
amount of work per cell. The `map(to: c)` fix was the opposite — a per-launch
constant of 6.5–9.9 ms/step from 2 ranks up — and the two must not be reasoned
about the same way. In percentage terms this one is nearly rank-independent so
far (−7.8 to −8.6 % of the step at both 4 and 8), because the step shrinks at
about the same rate as the saving.

## 6 — 16 ranks: the prediction, and what it turned into

Written as a prediction before job 5145120 ran, and kept in that order so the
forecast can be graded rather than quietly replaced.

**Predicted**, from the committed pre-change phase table (job 5142973, `new`
column: `apply` = 20.77 ms of a 60.31 ms step for `rect`, 11.06 of 33.44 for
`refined_yp82`) and the fractional reductions measured at 4 and 8 ranks:
`apply` → ~16.2 / ~8.8 ms, i.e. **~4.5 / ~2.2 ms/step saved, −7.5 % / −6.6 % of
the step**, share falling to about 28 % — with the caveat that it might come out
lower, because at 16 ranks each GPU holds only 8.6 / 3.8 Mcell and the kernel has
fewer waves to fill.

**Measured** (job 5145120, 4 nodes, `hkn[0436,0530,0533,0630]`, 100 steps,
before and after in one allocation, `new` built from a worktree pinned at
`6708193` so later work could not contaminate it):

| case | step | `apply` | `apply` share | `sweep` (control) |
|---|---|---|---|---|
| `rect_jacobi` 16×4 | 59.91 → 55.02 ms **−8.2 %** | 20.76 → 16.34 **−21.3 %** | 34.7 → **29.7 %** | +0.6 % |
| `refined_yp82` 16×4 | 32.91 → 30.78 ms **−6.5 %** | 11.17 → 9.08 **−18.7 %** | 33.9 → **29.5 %** | −0.1 % |
| `rect_jacobi` 8×2 | 102.54 → 93.90 **−8.4 %** | 40.60 → 31.81 **−21.7 %** | 39.6 → 33.9 % | +0.7 % |
| `refined_yp82` 8×2 | 51.78 → 47.71 **−7.9 %** | 19.90 → 15.94 **−19.9 %** | 38.4 → 33.4 % | +0.0 % |

The forecast was right to **0.3 ms and 0.7 percentage points** on both cases
(4.42 / 2.09 ms saved against 4.5 / 2.2 predicted; −8.2 / −6.5 % against
−7.5 / −6.6 %), and the caveat fired in the direction stated: the refined case's
fractional `apply` reduction softens from 20.9 % at 4 ranks to **18.7 %** at 16,
as its per-GPU 3.8 Mcell stops filling the machine. Its 8-rank rows reproduce job
5145114's to 0.3 % — a third allocation agreeing with the other two.

The `before` column also settles what `results_horeka_2026-09-14.md` §6 said was
owed: `apply` at 16 ranks is **34.7 % / 33.9 %** of the step, against the ~34 %
inferred there by subtraction.

(The first submission of this job, 5145101, would have rebuilt `new` from the
live working tree days later, and was cancelled and replaced for that reason.)


## 7 — Gates

All at `max_abs 0`, deliberately **without** `-Mnofma`: the two binaries differ
in Fortran source, but the hoisted expression is the identical
`face_grad_corr` call on identical inputs, so a production-flag comparison is the
tighter test and any nonzero difference would be real.

| gate | result |
|---|---|
| Pass G, `rect_jacobi` (138 412 032 points) | un/vn/wn/pn **max_abs 0** |
| Pass G, `refined_yp82_rect_jacobi` (60 555 264 points) | un/vn/wn/pn **max_abs 0** |
| 7-case suite, GPU, 200 steps | all **max_abs 0**, including `nut`, `k`, `omega`, `gamma`, `rethetat` |
| `min_channel` 4 ranks (blocks + 2:1 + Chebyshev) | max_abs 0 |
| `min_channel` 1 rank, `beltrami_slaby`, `les_ibm`, `turb180`, `wf180_y30`, `lam30t` | max_abs 0 |
| CPU build, `min_channel` + `beltrami_slaby`, 20 steps, 1 rank | max_abs 0 |
| `STACK`/`LOCAL` in every kernel, both binaries | 0 |

`les_ibm` + `refine_body` remains the one case of the standard list this gate
does not cover — its `IC_refine.h5` is generated by `setup.sh`, not committed,
and HoreKa has no h5py. Unchanged from previous campaigns.

## 8 — What these numbers do NOT support

- **Not that the kernel got cheaper.** It moves the same bytes — 10.31
  doubles/cell, unchanged — and that is still 1.29x its source-counted minimum of
  8. Nothing was reclaimed; the same traffic is merely issued by more warps.
- **Not a per-kernel 16-rank result.** §6's step and bucket times are measured,
  but the ncu occupancy/DRAM pair behind them comes from a single rank holding a
  whole GPU. Whether k2 still reaches 46 % occupancy on 3.8 Mcell is not
  measured — the softer 18.7 % there suggests not entirely.
- **Not a claim about `interface_correct`, the exchange or `mpi_wait`.** None is
  touched, and `mpi_wait`'s *share* will have grown again simply because the step
  is smaller — as it did after the `map(to: c)` fix.
- **Not that 50 % occupancy was reached.** 46.0 % achieved, 60.5 % of peak DRAM.
  Whether the remaining gap to k1's 62.5 % is worth chasing is unmeasured.
- **Not a CPU result.** The hoist removes work on the CPU path too; the CPU build
  was gated for bit-exactness only, never timed.

## 9 — What is next, sized

`jacobi_apply` is no longer the largest single item. In order:

1. **`compute_rdenom` sits at 110 registers** — the highest in the projection —
   and `proj_timing: setup`, which is essentially that kernel, is 2.3 ms/step at
   16 ranks (~4 % of the step). `face_grad_denom` is static in exactly the same
   way `face_grad_corr` was, so the same hoist applies verbatim and would be
   bit-exact by the same argument. **This is the obvious next lever and it is
   cheap.**
2. `interface_correct`'s three kernels per call (54 launches/step, already at the
   per-launch floor) — a kernel-count question worth ~0.5 ms/step.
3. The 39-round exchange cadence — a numerics change needing a convergence
   argument, not a bit-exactness gate.
4. `mpi_wait`, still ~10–12 % of the step and untouchable by partitioning or
   overlap (2026-09-10 §3, and the A0 probe).

The handout `docs/next_session_apply_registers.md` is **closed**: its question
was "can 88 become 64", the answer is yes by L2 alone, and the time followed as
its §6 predicted.
