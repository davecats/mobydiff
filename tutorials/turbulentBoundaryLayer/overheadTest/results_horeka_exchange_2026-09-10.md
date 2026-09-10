# The exchange at scale is a per-LAUNCH cost, not a per-cell or per-point one

Phases 1 and 2 of `docs/next_session_2to1_performance.md`. Jobs 5139461 and
5139581 (HoreKa `dev_accelerated`, hkn[0401,0403], 2 nodes, 1/4/8 ranks, 100
steps, everything below in ONE allocation unless it says otherwise); sections 1
and 2 are arithmetic on the committed 16-rank logs of job 5139351.

The handout set the target as "the refined case does 82 % of the single-level
twin's exchange work while carrying 44 % of its cells — that factor ~1.9 is the
2:1 interface's cost at scale, it is device-local, and it is the target", and
pre-registered three branches: cross-level volume, copy volume, or a per-point
cost difference between op kinds.

**None of the three is what the data says.** The exchange at scale is bought by
the *number of kernel launches* multiplied by a fixed cost each carries before it
moves anything, and by nothing about the mesh at all. The 2:1 interface's own
share of it is 85 % launch and 15 % data.

## 1 — The 1.9x is not extra points, and this was already in the repository

`exchange sizes` is printed by every profiled run and has been since `be48d44`.
Summing its two numbers (local copy points, which never become messages, and
send points, which do) over the whole job:

| config | cells | exchange points | points per cell |
|---|---|---|---|
| `rect_jacobi` (single level) | 138.41 M | 15.390 M | 1.112e-1 |
| `refined_yp82_rect_jacobi` | 60.56 M | 6.425 M | 1.061e-1 |
| `refined_big_rect_jacobi` | 242.22 M | 25.786 M | 1.065e-1 |

**0.954x and 0.957x.** The refined cases exchange slightly FEWER points per cell
than the single-level twin, and the total is independent of rank count (the same
halo points, some of them promoted from a device copy to a message). So the
handout's first two branches cannot fire: there is no extra volume of either
kind to find, and "1.9x exchange per cell" is **1.9x per exchange POINT**.

## 2 — Where it does come from: a fixed cost per round

Fitting microseconds per exchange ROUND against points per rank, over the rank
counts of one configuration (`results_scaling/new/`, job 5139351; the 1-rank runs
are excluded — they have no peers, so `pack`/`unpack` are identically zero, and
their copy kernel is an order of magnitude larger than at any other rank count):

| config | pack fixed | unpack fixed | local_copy fixed | local_copy ns/pt | worst residual |
|---|---|---|---|---|---|
| `base_jacobi` (34.5k–574k send pts/rank) | **95.9 us** | **86.8 us** | — (no local copies) | — | 2.6 us |
| `rect_jacobi` | 101.5 | 93.6 | **56.1** | 0.098 | 1.3 us |
| `refined_yp82_rect_jacobi` | 98.4 | 85.5 | **114.5** | 0.099 | 0.9 us |
| `refined_big_rect_jacobi` | 103.6 | 97.5 | **109.8** | 0.101 | 1.6 us |
| `refined_yp82_rect_redblack` | 97.3 | 85.1 | **117.9** | 0.099 | 1.7 us |

The `local_copy` column is the PRE-SPLIT aggregate — that binary timed the
same-level and cross-level kernels into one bucket — which is why the three
refined rows read roughly double the single-level one. Section 5 separates them
and shows the doubling is the second kernel's launch.

Residuals of 1–3 us on numbers of 60–120 us, across a 17x span in points: the
affine split is not being forced onto the data. Reading it:

- **The marginal point costs the same everywhere.** 0.098–0.101 ns in every
  configuration — single level or refined, `nb = 16` or `nb = 64 44 48`. The 2:1
  interface carries **no per-point penalty**, which retires the handout's third
  branch as it was written.
- **The fixed part is 250–300 us per round**, and there are 39 rounds in every
  step of every configuration ever measured here. (Section 5 sharpens this: it is
  a cost per kernel LAUNCH, not per round — a round whose copy kernel has no work
  is skipped and costs 0.1 us.) That is 9.8 ms/step for
  `rect_jacobi` and 11.6 ms/step for `refined_yp82` — **14 % and 28 %** of their
  16-rank steps.
- At 16 ranks the two configurations spend **the same absolute device-local
  exchange time**, 13.24 against 13.29 ms/step, while carrying 138.41 M and
  60.56 M cells. *That identity is the 1.9x.* The refined case is not doing more
  work; it is paying an identical, mesh-independent bill on 44 % of the cells.
- **The refined copy intercept is double the single-level one** (114.5 / 109.8 /
  117.9 against 56.1) at an identical slope. The candidate is structural: a full
  exchange launches the cross-level kernel as a SECOND kernel, in the 21 of 39
  rounds per step that carry one.

Everything above is arithmetic on logs that were already committed. What the new
runs add is (a) the same fit inside ONE allocation, (b) the cross-level kernel
timed in its own bucket instead of inferred from a difference, (c) the op split
that says how many points are cross-level at all, and (d) A0 at 8 and 16 ranks.
The predictions written before those runs are in
`horeka/exchange/PREREGISTERED.md`.

## 3 — Phase 2: A0 repeated at 8 ranks. It fails again, harder.

A0 (2026-09-02) inserted a compute-bound target kernel between the `Isend`/
`Irecv` posts and the `MPI_Waitall` at **2** GPU ranks and found that nothing was
hidden. Its own text asked for a repeat on a genuinely many-rank machine before
the conclusion was carried, because the message sizes, peer counts and transports
all change. Repeated here on 2 nodes, 8 ranks, 100 steps, job 5139461 — a
5.7 ms kernel against waits of 49 and 75 us, i.e. **76–116x** the quantity it
would have to hide:

| case | probe us/round | `mpi_wait` us/round, no probe | with probe | | s/step no probe | with probe |
|---|---|---|---|---|---|---|
| `refined_yp82_rect_jacobi` 8x2 | 5701.5 | 49.1 | **70.7** | **1.44x** | 0.060138 | 0.283078 |
| `rect_jacobi` 8x2 | 5700.4 | 75.0 | **79.4** | **1.06x** | 0.109144 | 0.331868 |

**The wait goes UP, never down.** And the step time rises by 0.2229 s = 39 x
5.7 ms exactly: the probe kernel is serialised into the step in full, with no
part of it overlapping anything. There is no asynchrony to exploit here at all —
not a partial one, not a small one.

So the 2-rank conclusion carries: **the transfer happens inside `MPI_Waitall`
and nowhere else, and P2 (overlap) is closed** on this stack. The scope caveat
the 2-rank probe attached to itself is now discharged rather than inherited.
The mild degradation (6 % and 44 %) is the same contention the 2-rank probe saw
(9 %).

## 4 — Phase 1: the op split. Cross-level is 16.4 % of the points, at every rank count

Job 5139581, 2 nodes (hkn[0401,0403]), 100 steps, all runs in ONE allocation.
Points summed over ranks, from the new `exchange by op` line:

| config | ranks | same-level | restrict | prolong | cross-level share |
|---|---|---|---|---|---|
| `rect_jacobi` | 1 / 4 / 8 | 15.39 M / 15.39 M / 15.39 M | 0 | 0 | 0 |
| `nb16_jacobi` | 1 / 4 / 8 | 57.02 M / 57.02 M / 57.02 M | 0 | 0 | 0 |
| `refined_yp82_rect_jacobi` | 1 | 5 369 280 | 211 000 | 844 400 | **16.43 %** |
| | 4 | 5 369 280 | 211 000 | 844 400 | **16.43 %** |
| | 8 | 5 369 280 | 211 000 | 844 400 | **16.43 %** |
| `refined_big_rect_jacobi` | 4 / 8 | 21.56 M | 843 200 | 3 378 400 | **16.38 %** |

(The per-rank split between "local copy" and "send" moves with rank count; the
totals above are invariant, as they must be — the same halo points, some promoted
from a device copy to a message.)

**Prediction 1 confirmed** (≤ 20 %) and **prediction 2 confirmed**: 16.43 % at 1,
4 and 8 ranks, identical to four figures, and 16.38 % on a case with 4x the
leaves. The interface is a geometric property of the leaf table and the Morton
split does not cut it. Prolong outnumbers restrict 4:1, as the 2:1 geometry
requires — four fine destinations per coarse source face.

Entries tell a different story from points: **3 288 cross-level entries against
5 780 same-level ones — 36 % of the entries for 16 % of the points** (321 against
929 points each). A coarse face fed by four fine sub-entries is four small
transfers, and that ratio is what the next section makes expensive.

## 5 — Where the refined case's exchange time actually goes

Per exchange ROUND, rank 0, 8 ranks on 2 nodes, one allocation. `local_copy` is
now the same-level kernel alone; `copy_cross` is the 2:1 restrict/prolong kernel,
which launches **24 times per step** (6 full velocity exchanges + 18 scalar) where
everything else launches 39:

| config | s/step | pack | unpack | local_copy | copy_cross | mpi_wait | device-local ms/step | % of step |
|---|---|---|---|---|---|---|---|---|
| `base_jacobi` | 0.103252 | 165.5 | 125.0 | 0.1 | 0.1 | 133.7 | 11.35 | 11.0 % |
| `rect_jacobi` | 0.109302 | 100.3 | 91.6 | 237.3 | 0.1 | 76.9 | 16.72 | 15.3 % |
| `nb16_jacobi` | 0.141611 | 100.9 | 97.0 | 813.0 | 0.2 | 81.3 | 39.36 | 27.8 % |
| `refined_yp82_rect_jacobi` | 0.060292 | 99.7 | 85.7 | 130.0 | **99.6** | 55.6 | 14.71 | **24.4 %** |
| `refined_big_rect_jacobi` | 0.186315 | 104.1 | 95.5 | 334.9 | **148.0** | 88.9 | 24.42 | 13.1 % |

**`rect_jacobi` spends 16.72 ms/step of device-local exchange on 138.41 M cells;
`refined_yp82` spends 14.71 ms on 60.56 M.** 88 % of the work for 44 % of the
cells, while moving 42 % of the points. That is the handout's 1.9x, measured
inside one allocation with the same binary.

### The split that explains it

Two-point fits (4 and 8 ranks) with the RIGHT x for each kernel — same-level
points for `local_copy`, cross-level points for `copy_cross`:

| config | kernel | fixed us per launch | ns per point |
|---|---|---|---|
| `rect_jacobi` | same-level copy | **56.3** | 0.097 |
| `nb16_jacobi` | same-level copy | **58.5** | 0.107 |
| `refined_big` | same-level copy | **58.5** | 0.106 |
| `refined_yp82` | same-level copy | **66.7** | 0.100 |
| `refined_yp82` | **cross-level copy** | **84.6** | **0.114** |
| `refined_big` | **cross-level copy** | **87.3** | **0.115** |
| `base_jacobi` (widest send span) | pack | **114.0** | 0.090 |
| `base_jacobi` | unpack | **82.5** | 0.074 |

- **A cross-level point costs 0.114–0.115 ns against a same-level point's
  0.097–0.107 — 1.15x, not 2–4x.** Prediction 3, and with it the handout's third
  branch as written ("the per-point cost differs between op kinds"), is
  **REFUTED**. The 2:1 transfer is not expensive per point.
- **The cross-level kernel costs 84.6–87.3 us per LAUNCH before it moves
  anything**, and it launches 24 times per step. At 8 ranks that is 2.03 ms/step
  of pure launch against 0.36 ms of actual transfer: **85 % of the 2:1
  interface's exchange cost is the extra kernel launch, not the data.**
- The fixed cost is a per-LAUNCH cost, not a per-round one, and `base_jacobi`
  proves it: at 8 ranks it has exactly **zero** local copy points, the kernel is
  skipped, and the bucket reads **0.1 us** — while the same kernel with 364 722
  points at 4 ranks costs 82.1 us.

Adding it up at 8 ranks:

Volume is priced with `base_jacobi`'s wide-span slopes (0.090 ns/pt for pack,
0.074 for unpack — the narrow-span configs cannot resolve a slope) and each
config's own copy slope:

| | fixed (launch) | volume | total | check vs measured |
|---|---|---|---|---|
| `rect_jacobi` | 9.27 ms/step (55 %) | 7.47 ms | 16.74 ms | 16.72 |
| `refined_yp82` | **11.61 ms/step (79 %)** | 3.09 ms | 14.69 ms | 14.71 |

**The two configurations pay almost the same launch bill — 9.3 against 11.6 ms —
and it is the whole difference.** `refined_yp82` carries 44 % of the cells, so
its volume term collapses from 7.1 to 3.0 ms while the launch term does not move.
`nb16_jacobi` is the control on the other side: `nb = 16` gives 3.8x `rect`'s
copy volume, its fixed bill stays ~9 ms and its volume term is 29 ms — at fine
granularity, volume still dominates.

**So the fixed per-launch bill is ~10–12 ms/step in EVERY configuration, and it
becomes the dominant term exactly when strong scaling has shrunk the per-rank
volume — which is the regime 2:1 refinement puts you in by design.**

## 6 — The pre-registered scorecard

`horeka/exchange/PREREGISTERED.md` was written and committed before the runs.

| # | prediction | outcome |
|---|---|---|
| 1 | cross-level ≤ 20 % of the refined case's exchange points | **CONFIRMED** — 16.43 % |
| 2 | the op split is rank-count independent | **CONFIRMED** — 16.43 % at 1, 4 and 8 ranks, to four figures |
| 3 | `copy_cross` costs 2–4x per point what `local_copy` does | **REFUTED** — 1.15x (0.114 vs 0.100 ns/pt). Its 32 % share of copy time is real but comes from the launch, not the point |
| 4 | pack/unpack per round vary < 2x while their points vary > 10x | **CONFIRMED** — 1.68x against 16.6x |
| 5 | A0 fails again at scale | **CONFIRMED** — waits at 1.06x, 1.14x, 1.44x of baseline; never below 1 |
| 6 | `copy_cross` reads 105–115 us per call | **PARTIAL** — 114.7 at 4 ranks, 99.6 at 8; the prediction ignored that the quantity has a volume term. Its intercept is 84.6 us |

Prediction 3 is the informative failure: it was the handout's own third branch,
and the instrument built to test it (the `copy_cross` bucket) says the per-point
costs are within 15 % of each other.

Prediction 6 was written to test a mechanism — "5.5 us per map-clause item" —
that a control in the second PREREGISTERED addendum had already refuted before
the run, using `jacobi_compute_phi` (10 map items, 19 us intercept) against
`copy_local_same_level` (11 items, 56 us). Counting only NON-scalar map items
tightens it (6 → 3.2 us each for the sweep; 8.5–15 → 5.6–7.4 us each for the four
exchange kernels), but a 2.3x spread over a 2.5x range in count is a correlate,
not a cause, and revising a hypothesis to fit the datapoint that killed it is how
this campaign produced its previous two wrong mechanisms. **The per-launch cost
is measured; what it consists of is not identified here.**

## 7 — What these numbers do NOT support

- **Not "the 2:1 interface transfer is expensive."** It is 16.4 % of the points
  at 1.15x the per-point cost — 0.36 ms/step of actual transfer at 8 ranks.
- **Not per-level block size (P3) as the exchange lever.** The refined case
  already exchanges 0.95x the points per cell of its single-level twin, its
  volume term is the smaller half (3.0 of 14.7 ms), and bigger coarse blocks
  would not change the number of exchange ROUNDS at all — which is what the
  launch bill is multiplied by. P3 may still be right for *compute*; it is not
  the answer to this.
- **Not partitioning.** The op split is invariant across 1, 4 and 8 ranks, so the
  Morton chain is not cutting the interface, and `mpi_wait` is 55–89 us/round
  (2–4 % of the step) with the corrected GPU mapping.
- **Not overlap.** Section 3.
- **Not a mechanism for the per-launch cost.** Correlated with the number of
  mapped arrays; not proportional to it; not traced.
- **Not 16 ranks.** Everything above is 1/4/8 ranks in one allocation, because
  the 4-node queue is three days out. The committed 16-rank logs of job 5139351
  give the same picture through the older aggregate bucket (13.24 vs 13.29 ms/step
  of device-local exchange on 138.41 M and 60.56 M cells), and the split
  projected to 16 ranks gives 11.7 ms/step of launch cost = **28 % of the refined
  step** — but that projection has not been measured. Job 5139583 is queued for it.
- **One machine, one node type**, as always.

## 8 — The lever, and the experiment that would confirm it

There are ~141 exchange kernel launches per step (39 pack, 39 unpack, 39
same-level copy, 24 cross-level copy) costing **10–12 ms/step in every
configuration measured** — 9.3 ms for the single-level case, 11.6 ms for the
refined one, ~9 ms for `nb16`. Against the projection's `jacobi_compute_phi`,
which launches at **19 us**, 7.0 ms/step (`rect`) to 8.9 ms/step (`refined`) of
that is not obviously necessary.

Ranked by size, and by how well this data supports them:

1. **Find and remove the per-launch cost.** 7.0–8.9 ms/step, 15 % of the refined
   8-rank step and ~19 % projected at 16. Not refinement-specific — every blocked
   configuration pays it. **This is not yet actionable: the first step is a
   per-kernel nsys timeline** (nsys works on this solver for CUDA tracing;
   `results_horeka_2026-09-10.md` §5) to see what those 60–100 us are made of,
   plus a map-clause / hoisted-`target data` A/B on ONE kernel. Do the timeline
   first; do not implement from the correlation.
2. **Fold the local copy into the pack kernel.** They are independent and both
   run before the `Waitall`, so one launch could do both: 39 launches × ~60 us =
   **2.3 ms/step**, bit-exact, no numerics. The cheapest concrete win here.
3. **Fewer rounds.** 39/step, each costing ~240 us of launch alone. Halving the
   18 phi exchanges would save ~4.5 ms/step — but that is a numerics change and
   needs a convergence argument, not a bit-exactness gate.
4. **Fusing the cross-level kernel into the same-level one is NOT recommended**
   on this data. It would save 24 × 84.6 us = 2.0 ms/step of launch, but the two
   were split deliberately: the interface gather needs ~128 registers and runs at
   a third of the light kernel's occupancy (measured with ncu, recorded in
   `copy_local_same_level`'s comment). Paying that on all 39 rounds against a
   same-level volume term of only 2.5 ms/step would very likely lose more than
   the 2.0 ms it saves.
