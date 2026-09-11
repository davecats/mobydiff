# Pre-registered readings — the 2:1 exchange at scale

Written **before** job `moby_exch` ran, from the committed logs of job 5139351
(`results_scaling/new/`) and job 5138715 (`results_balance/`). This campaign has
twice published a mechanism that the next measurement overturned, so what the
numbers are expected to say is written down first and the report is graded
against it.

## What is already derivable from the committed logs

These are not predictions about new data — they are arithmetic on data that is
already in the repository, and they are what motivated the volume sweep added to
Pass 1. They are recorded here so the report cannot be read as having discovered
them in the new runs.

**(A) The refined case's exchange POINTS per cell are the same as the
single-level twin's.** At 16 ranks, `exchange sizes` gives 6.42 M points
(5.735 M local + 0.690 M send) for 60.56 M cells, against 15.39 M for 138.41 M
cells: 1.061e-1 vs 1.112e-1 points per cell, i.e. **0.954x**. The handout's
"1.9x exchange per cell" is therefore **1.9x TIME PER POINT**, not extra points.
Pre-registered branch 1 of the handout ("cross-level points are a small fraction
and the 1.9x is COPY volume") cannot fire as stated: there is no extra volume of
either kind.

**(B) `pack` and `unpack` are almost volume-independent.** Per exchange round,
across five configurations at 8 and 16 ranks, send points per rank span 43 k to
575 k (13x) while `pack` spans 99 to 168 us and `unpack` 86 to 125 us. An affine
fit gives pack = **92 us + 0.13 ns/pt**, unpack = **86 us + 0.07 ns/pt**.

**(C) `local_copy` has a fixed part too**, though smaller: `rect` and `nb16` at
8 ranks (1.86 M and 7.05 M local points per rank) give **31 us + 0.111 ns/pt**.

Summed, that is **~210 us per round of cost that no volume change can touch**,
against 39 rounds per step — 8.2 ms/step, or **20 % of the refined 16-rank step
(41.2 ms)** and 12 % of the single-level twin's (67.8 ms). It is the same
absolute number in both, which is exactly why the refined case looks 1.9x more
expensive per cell: it carries 44 % of the cells and pays the identical fixed
bill.

## Predictions for the new runs

| # | prediction | what refutes it |
|---|---|---|
| 1 | Cross-level (`restrict` + `prolong`) points are **≤ 20 %** of the refined case's exchange points at 1, 4 and 16 ranks. | a share above 20 % — then the interface transfer really is the volume, and the handout's branch 2 fires. |
| 2 | The op split is essentially **rank-count independent** (the interface is a geometric property of the leaf table). | a share that grows with rank count — then the Morton split is cutting the interface. |
| 3 | `copy_cross` costs **2–4x per point** what `local_copy` (same level) costs, so it takes 25–45 % of the copy time while carrying ~15 % of the points. | ≈ 1x — then the interface kernel is not special and the cost is granularity. |
| 4 | The volume sweep confirms (B)/(C) **within one allocation**: pack and unpack per round vary by less than 2x while their points vary by more than 10x. | per-round time tracking points — then the exchange is bandwidth-bound after all and the fixed-cost reading is wrong. |
| 5 | **A0 fails again** at 8 and 16 ranks: inserting a multi-millisecond compute kernel between the posts and the `Waitall` leaves `mpi_wait` at ≥ 0.9x its baseline. | `mpi_wait` collapsing — then MPI does progress during a kernel at these message sizes, P2 (overlap) reopens, and the 2-rank A0 conclusion was scope-limited exactly as it warned. |

## What would follow

- 1 + 3 + 4 ⇒ the lever is the **per-round fixed cost** (kernel launch and map-clause
  marshalling: ~1300 host-to-device memcpys per step were already seen by nsys and
  never chased), not block size, not the entry structure, and not the transport.
  That is the handout's third branch, reached by a different route than it
  anticipated.
- 2 refuted ⇒ partitioning is back on the table for the refined case specifically.
- 5 refuted ⇒ P2 reopens and everything else waits.

---

## Addendum, still before the run: the fits, and a sharper set of predictions

Fitting the committed `results_scaling/new/` logs properly — microseconds per
exchange ROUND against points per rank, **excluding the 1-rank runs** (they have
no peers, so `pack`/`unpack` are identically zero and the copy kernel is an order
of magnitude larger than at any other rank count):

| config | pack fixed | unpack fixed | local_copy fixed | local_copy ns/pt | worst residual |
|---|---|---|---|---|---|
| `base_jacobi` (34.5k–574k send pts) | **95.9** | **86.8** | — | — | 2.6 us |
| `rect_jacobi` | 101.5 | 93.6 | **56.1** | 0.098 | 1.3 us |
| `refined_yp82_rect_jacobi` | 98.4 | 85.5 | **114.5** | 0.099 | 0.9 us |
| `refined_big_rect_jacobi` | 103.6 | 97.5 | **109.8** | 0.101 | 1.6 us |
| `refined_yp82_rect_redblack` | 97.3 | 85.1 | **117.9** | 0.099 | 1.7 us |

Residuals of 1–3 us on numbers of 60–120 us across a 17x span of points: the
affine model is not being forced.

Three things fall out, all of them BEFORE any new measurement:

1. **The marginal point costs the same everywhere** — 0.098–0.101 ns in every
   configuration, refined or not, `nb = 16` or `nb = 64 44 48`. There is no
   per-point penalty for the 2:1 interface.
2. **The fixed per-round cost is 250–300 us** (pack + unpack + copy), against
   39 rounds: 9.8 ms/step for `rect`, 11.6 ms/step for `refined_yp82` — **14 %
   and 28 %** of their 16-rank steps. At 16 ranks the two configurations spend
   the SAME absolute device-local exchange time (13.24 vs 13.29 ms/step) while
   carrying 138.41 M and 60.56 M cells. That identity IS the "1.9x per cell".
3. **The refined case's copy intercept is double the single-level one** (114.5 /
   109.8 / 117.9 against 56.1). The obvious candidate is the cross-level kernel:
   a SECOND launch, in the 21 rounds per step that carry one.

### The map-clause count, and the number this session can falsify

Counting `map(...)` items in each kernel and weighting by the 21 velocity /
18 scalar rounds per step:

| bucket | items (weighted) | fixed us | us per map item |
|---|---|---|---|
| `local_copy` (same level) | 10.1 | 56.1 | **5.55** |
| `unpack` | 16.4 | 93.6 | **5.71** |
| `pack` | 18.6 | 101.5 | **5.46** |

Three kernels, map counts spanning 1.8x, and a per-item cost agreeing to 3 %.
Most of those clauses are REDUNDANT: the arrays they name are already resident
from `init_block_exchange`'s `target enter data`, and the scalars would be
firstprivate without a clause. nsys already saw **~1300 host-to-device memcpys
per step at ~2 us** in every configuration and nobody chased them.

**Prediction 6, the sharp one.** The new `copy_cross` bucket weights 21 items
(velocity) and 19 (scalar) over its 21 launches per step, i.e. 20.1 — so at
5.5 us per item it should read **110 us per call**. Independently, the same
number falls out of the old aggregate: `(114.5 − 56.1) × 39/21 = 108.5 us`.
**If `copy_cross` comes back at 105–115 us per call, four independent kernels
agree on ~5.5 us per map-clause item and the mechanism is named. If it comes
back well outside that, the map-clause reading is wrong** and the fixed cost is
something else (kernel-launch latency, register pressure, an implicit sync) that
a timeline would have to find.

Either way the A/B that settles it is one run: strip the redundant clauses from
one kernel and re-fit its intercept. That is a scheduling change, so it must be
bit-exact — and it is NOT part of this session's deliverable.

---

## Second addendum, still before the run: evidence AGAINST my own mechanism

Before treating "5.5 us per map-clause item" as the mechanism, the same fit was
run on the kernels OUTSIDE the exchange, from the same committed logs —
microseconds per call against Mcell per rank:

| config | kernel | calls/step | Mcell/rank span | fixed us | ps/cell | worst residual |
|---|---|---|---|---|---|---|
| `rect_jacobi` | `sweep` (`jacobi_compute_phi`, **10 map items**) | 18 | 8.65–138.41 | **18.9** | 44.4 | 4.3 us |
| `refined_yp82` | `sweep` | 18 | 3.78–60.56 | 21.2 | 44.3 | 6.2 us |
| `refined_big` | `sweep` | 18 | 15.14–60.56 | 15.8 | 44.4 | 3.3 us |
| `rect_jacobi` | `momentum` (**25 map items**) | 3 | 8.65–138.41 | 11.2 | 378 | 52 us (poor) |

`jacobi_compute_phi` carries **10** map items and a ~19 us intercept;
`copy_local_same_level` carries **11** and a ~56 us intercept. Three times the
fixed cost for the same number of clauses. **The map-item count is not the
explanatory variable**, and the 5 % agreement across the three exchange kernels
is very likely a coincidence of three kernels whose map counts happen to track
their complexity.

A second warning sign, from the same data: the affine model does NOT extrapolate
to 1 rank. `rect_jacobi`'s copy kernel at 1 rank moves 15.39 M points in 772 us
(50 ps/pt) where the 2–16-rank fit predicts 1564 us. And at nearly the same size
`nb16_jacobi` at 4 ranks moves 14.19 M points in 1577 us (111 ps/pt) — the same
kernel, 2.2x the rate difference, driven by entry structure and coalescing. So
the "slope" is a property of a config's entry shape and the "intercept" is the
small-size end of a concave curve, not necessarily a launch cost.

**Prediction 6 stands as a number** (`copy_cross` at 105–115 us per call) because
it also follows from the old aggregate by an independent route. But confirming
it would NOT name the mechanism; it would only confirm that the cross-level
kernel is a second per-round call of the same order as the first. The mechanism
needs a timeline, and this file records that I was one step from publishing a
third wrong one.

---

## Third addendum — pre-registered readings for the nsys timeline (job 5141872)

Written before the trace came back. The quantity being explained is the
per-launch fixed cost measured on 2026-09-10: same-level copy 56–67 us,
cross-level copy 84.6–87.3, pack ~92–114, unpack ~83–100, against
`jacobi_compute_phi` at 24–26.

**The control first.** `jacobi_compute_phi` at 4 ranks on `refined_yp82` should
show a device duration of ~665 us against a host bracket of ~690 — i.e. **the
cheapest kernel in the solver should still cost ~25 us of host-side time per
launch**. If it does not — if its host bracket and device duration agree — then
the 24–26 us intercept is not launch overhead at all and everything built on the
"fixed per launch" reading needs re-examining.

**The number that decides the mechanism: H2D copies per launch.** If the fixed
cost is OpenMP `map`-clause marshalling, each launch should be preceded by one
small host-to-device copy per mapped item that needs one, and the count should
track the kernel's map list — *not* the map count I used in the first addendum
(which the sweep control refuted) but the count of items that actually require a
transfer. Counting derived-type components and module arrays separately from
plain scalars:

| kernel | map items | of which derived-type components / module arrays | predicted H2D per launch |
|---|---|---|---|
| `jacobi_compute_phi` | 10 | 4 components + 2 module arrays | **4–6** |
| `copy_local_same_level` | 11 | 9 components | **9–11** |
| `pack_entries` | 20 | 15 | **15–17** |
| `unpack_entries` | 21 | 15 | **15–17** |
| `copy_local_cross_level` | 21 | 16 | **16–18** |

nsys already reported ~1300 H2D copies per step at ~2 us in every configuration
(`results_horeka_2026-09-10.md` §5, "not chased"); ~200 kernel launches per step
makes that ~6.5 per launch on average, which is the right order for this table.

| outcome | conclusion |
|---|---|
| counts track the component column, sweep low and pack/unpack high | **the fixed cost is map-clause marshalling**, and the fix is to stop re-mapping resident arrays every launch (hoisted `target data`, or device pointers). The first addendum's mechanism was right in kind and wrong in the variable it was counted against. |
| counts are flat across kernels | marshalling is NOT it; look at the API summary for what else scales — `cuLaunchKernel` argument setup, or an implicit synchronise per region. |
| device duration ≈ host bracket for the exchange kernels | the "fixed cost" is device work after all, the kernels are simply inefficient at small sizes, and **sections 5 and 8 of the report are wrong**. |
