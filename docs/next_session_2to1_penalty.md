# Next session — reducing the 2:1 penalty

> **STATUS 2026-09-10 (later) — P2 (OVERLAP) IS CLOSED. THE EXCHANGE COST AT
> SCALE IS KERNEL LAUNCHES.** `overheadTest/results_horeka_exchange_2026-09-10.md`,
> jobs 5139461 / 5139581.
>
> **A0 repeated at 8 ranks / 2 nodes and fails again**, harder than at 2: a 5.7 ms
> compute-bound kernel between the Isend/Irecv posts and the Waitall — 76–116x the
> wait it would have to hide — leaves `mpi_wait` at **1.06x / 1.14x / 1.44x** of
> baseline, never below 1, and the step time rises by exactly 39 x 5.7 ms. The
> probe is serialised in full; there is no asynchrony to exploit. The 2-rank
> probe's own scope caveat ("repeat on a genuinely many-rank GPU machine") is
> **discharged, not inherited**. P2 is closed and nothing in this file is left
> standing.
>
> **What replaces it.** The handout's "1.9x exchange per cell" is not volume of
> any kind: the refined case exchanges **0.954x the points PER CELL** of its
> single-level twin, and **cross-level points are 16.43 % of the total at 1, 4 and
> 8 ranks alike** (16.38 % on a case with 4x the leaves) — the Morton split does
> not cut the interface. Timing the cross-level kernel separately (new
> `exch_timing copy_cross` bucket) **refutes the handout's third branch too**: a
> cross-level point costs 0.114–0.115 ns against a same-level point's
> 0.097–0.107, i.e. **1.15x, not 2–4x**.
>
> The cost is **kernel LAUNCHES**. Each launch carries a fixed cost before it
> moves anything — same-level copy 56–67 us, cross-level copy 84.6–87.3 us, pack
> ~92–114, unpack ~83–100, against the projection's `jacobi_compute_phi` at 19 us
> — and there are ~141 of them per step (39 pack + 39 unpack + 39 same-copy +
> **24** cross-copy). That is **9.7 ms/step for `rect_jacobi` and 11.7 ms for
> `refined_yp82`**, i.e. 80 % of the refined case's device-local exchange at 8
> ranks and the entire difference between the two: they pay the same launch bill,
> and the refined case pays it on 44 % of the cells. `base_jacobi` proves the cost
> is per-LAUNCH: with zero local copy points at 8 ranks the kernel is skipped and
> the bucket reads **0.1 us**.
>
> The 2:1 interface's own exchange cost is therefore **2.03 ms/step of launch
> against 0.36 ms of transfer — 85 % launch, 15 % data.**
>
> **16 RANKS MEASURED (job 5139977, 2026-09-11) and everything holds.** Cross-level
> 16.43 % at a fourth rank count; `rect_jacobi` 13.18 ms/step and `refined_yp82`
> 13.23 ms/step of device-local exchange on 138.41 M and 60.56 M cells — the
> headline is now a direct measurement, not an inference. The 4-and-8-rank fits,
> unrefitted, predict the 16-rank projection kernels to **0.3 %** and the exchange
> kernels to **3.4 %**. Launch-fixed bill: 9.25 ms/step for `rect` (13.9 % of the
> step), 11.66 for `refined` (27.9 %), **13.28 with `interface_correct` = 31.8 %**,
> against 1.7 ms/step of actual halo data movement.
>
> **A THIRD launch cost sits outside the exchange.** The same fit on the
> projection buckets: `jacobi_compute_phi` (one kernel, 6 mapped arrays) launches
> at **24–26 us** whatever the mesh with an identical 43–46 ps/cell — the cheapest
> launch in the solver and the reference the exchange kernels should be read
> against — while `jacobi_apply` is 46 us fixed single-level and **137 us
> refined**, because it launches `interface_correct` as an extra kernel whenever
> interfaces are present, and that kernel does nb² work against the sweep's nb³.
> **91 us x 18 calls = 1.63 ms/step.** With the cross-level copy's 2.03 ms that is
> **3.66 ms/step = 6.1 % of the refined 8-rank step of refinement-specific cost,
> almost all of it launches rather than work. That is the 2:1 tax at scale.**
>
> **THE TIMELINE IS RUN (2026-09-11, job 5141872,
> `overheadTest/results_kernel_timeline_2026-09-11.md`) AND THE MECHANISM IS
> NAMED: every launch of an exchange kernel copies the ENTIRE `comm_type` object
> host-to-device — 7 952 bytes, all 45 of its array descriptors — plus 14–23
> separate descriptor and scalar copies, whatever the kernel actually names.** No
> other kernel in the solver pays one: the projection kernels move 8–148 B per
> launch and cost 5–22 us against the exchange kernels' 62–88. The arithmetic
> identifies the block exactly — 26 rank-1 descriptors at 152 B + 19 rank-2 at
> 200 B = 7 752 B, the rest being the type's scalars and handles — and the three
> copy sizes in the trace (8 / 152 / 200 B) are scalars, rank-1 and rank-2
> descriptors. The arrays are ALREADY resident from `init_block_exchange`'s
> `target enter data`; none of the 10 kB per launch is data the device needs.
> Independently, (untraced bracket − traced device time) reproduces the affine-fit
> intercepts to 1 % on two configs whose device times differ by 3x, so the
> "fixed per launch" reading now rests on two methods.
>
> **Estimated recoverable: 7.8–9.4 ms/step = 19–23 % of the refined 16-rank
> step.** Fix direction: stop referencing `c%component` inside the target regions
> — bind what each kernel needs to local arrays outside, or pass them as dummy
> arguments, so the map list holds plain already-present arrays instead of a
> derived type. **Next step is ONE kernel, not six**: hoist
> `copy_local_same_level`'s `c%` references and check that the 7 952 B block
> disappears and its bracket falls from 62 us toward the ~20 us floor. Scheduling
> change ⇒ bit-exact, 7-case suite + production Pass G, CPU and GPU.
>
> Superseded below: an nsys per-kernel timeline of what
> those 60–100 us consist of. The fixed cost CORRELATES with the number of mapped
> arrays but not proportionally (3.2 us each for the 6-array sweep kernel,
> 5.6–7.4 for the 8.5–15-array exchange kernels), and `PREREGISTERED.md` records
> that fitting a mechanism to that correlation was one step from being this
> campaign's third wrong mechanism. Concrete cheap win meanwhile: **fold the local
> copy into the pack kernel** (independent, both before the Waitall, one launch
> instead of two, 2.3 ms/step, bit-exact). Recommended AGAINST: fusing the
> cross-level kernel into the same-level one — it saves 2.0 ms of launch but puts
> the interface gather's ~128 registers on all 39 rounds, and ncu already measured
> that at a third of the light kernel's occupancy.

> **STATUS 2026-09-10 — FIXED IN THE SOLVER. P1 IS OFF.**
> `overheadTest/results_horeka_2026-09-10.md`. The stall was GPU affinity, not the
> exchange: **the two ends of a cross-node link must sit in the same GPU affinity
> class**. On HoreKa Green all three HCAs sit on NUMA 0 with GPU0/GPU1 while
> GPU2/GPU3 reach a NIC only across the inter-socket link, and the old
> `device = local_rank mod num_devices` paired a node's LAST local rank with the
> next node's FIRST (dev3 against dev0) — the mismatched case. Measured, wait per
> round: mixed 604/604 us, both-far-but-matched 81 us, both-near 75 us. **Matching
> is worth ~8x; being on the NIC-affine side a further ~1.3x.** "Both ends need a
> local HCA" was the first reading and is WRONG — both-bad-but-matched is fast.
>
> `comm.f90 select_target_device` now serves a node's first and last local ranks
> first, identically on every node (node-homogeneity is load-bearing: a per-node
> order cost base_jacobi 14 %). `MOBY_GPU_ORDER` overrides the device order;
> single-node runs keep the identity mapping exactly. Measured ref vs new, back to
> back: rect 8x2 **-20.4 %**, rect 16x4 **-26.3 %**, refined_yp82 16x4 **-26.7 %**,
> base_jacobi 16x4 +1.0 % (noise), single node unchanged. Fields bit-exact
> (max_abs 0, 138 M points). **GATE DISCHARGED 2026-09-10** (job 5139976): the
> 7-case suite runs on HoreKa rather than the unreachable workstation, against
> `build_gpu/moby_solve.ref`, and all seven PASS at max_abs 0 incl. every RANS
> scalar — deliberately WITHOUT nofma, because the arithmetic source is
> byte-identical between the two binaries and a production-flag comparison is
> therefore strictly tighter. Not covered: `les_ibm` + `refine_body`, whose
> `IC_refine.h5` `setup.sh` generates and the repo does not carry.
>
> **Block tax rect/base: 1.505 -> 1.098 at 16 ranks, 1.315 -> 1.058 at 8. There is
> no node-boundary block tax**, and the "19-25 % of the step" prize is collected by
> the device mapping. P1's partitioning rewrite and the transport options in this
> file address a cost that does not exist; **do not start either**. P2 (overlap) is
> untouched and is now the only item here the finding leaves standing.
>
> How it was found, and what the earlier instrument missed: the barrier split
> showed the wait was arrival spread, not transport (rect 8x2: 773 us skew, 1.2 us
> transfer). CUDA traces showed the GPU idle 37 %, so it was not unfinished device
> work, and the per-rank memcpy table named the two staging ranks. Also: every
> `mpi_wait` figure in the 09-07/08/09 reports is **rank 0's**, and rank 0 sits at
> or near the cross-rank minimum; and aggregate min/max/argmax is **blind to
> rotating skew** (rect 8x2 reads a healthy max/min = 1.09 while being 100 %
> arrival spread).
>
> `tools/moby_tune.sh` closes the portability half: it sweeps the C(ndev,2)
> candidate device orders on the real case and **recovers the affinity classes
> with no topology input** (calibrated on HoreKa: the four mixed pairs measure
> 730–757 us, the two matched pairs 77–101). Tune once per machine, put the
> printed line in a submit script, never re-tune inside a production run. Its tie
> threshold is 5 %, set from node-to-node variation — a 4 % winner at 16 ranks
> lost at 8, so the ranking between two matched classes is not a stable property.
>
> CAMPAIGN RE-MEASURED (job 5139351, 46 runs, one allocation, ref and new side by
> side): 1/2/4-rank runs move by ±0.2 % (the mapping is inert without cross-node
> links) and `base_jacobi` by ±0.3 % at every rank count. At 8/16 ranks
> `rect_jacobi` +19.4/+24.5 %, `refined_yp82` +21.4/+26.8 %, red-black
> +25.3/+29.1 %, `refined_big` +13.6/+20.7 %. Strong-scaling efficiency at 16
> ranks: `rect` 50 → 66 %, `refined_yp82` 36 → 49 %, `refined_big` 64 → 81 %.
> **Block tax 1.310 → 1.059 (8 ranks) and 1.499 → 1.130 (16).** The 2:1 headline
> is untouched — its honest 4-rank cost is 0.98 coarse-cell-equivalents either
> way — and `results_horeka_2026-09-07.md` now carries a SUPERSEDED-IN-PART note.
>
> NEXT: run moby_tune on any new machine before trusting the built-in order (on HoreKa it happens to be right, which is exactly what would
> hide a bug elsewhere); and the ~10 % of the step still in `mpi_wait` at 16 ranks
> is now the same for blocked and unblocked — a P2 (overlap) target, not a
> partitioning one.


Handout, written 2026-08-28 from the session that measured all of it. Read the
first section before planning anything: it moves the target.

## The penalty is not where it was assumed to be

"The 2:1 penalty" is three different things, and only one of them is worth work.

**(a) Per-cell machinery cost — 1.4 %, effectively closed.** Controlled
measurement (`overheadTest/results_smoother_refined_2026-08-28.md`): the same
base grid and block shape, Jacobi `niter = 6`, 1 GPU rank, with and without
refinement.

| | cells | s/step | ns/cell |
|---|---|---|---|
| unrefined 2048×176×96 | 34.60 M | 0.3121 | 9.020 |
| refined, 1 level | 60.56 M | 0.5539 | **9.147** |

**+1.4 % per cell.** For comparison the refined case runs 2.20× faster than the
equivalent single-level grid while carrying 0.438 of its cells — i.e.
refinement delivers **96.5 % of its ideal speedup**. At one rank there is
almost nothing left to win here, and effort spent on the interface *kernels*
is effort wasted.

The projection is not degraded either: Jacobi's `L2_div` at `niter = 6` is
2.9336e-05 refined against 2.9395e-05 unrefined — **ratio 0.998, no
convergence penalty from the interface at all**. (Red-black *is* degraded,
ratio 1.286 — see (d).)

**(b) Refinement granularity — a cell-count effect, case-dependent, not
measured.** Refinement is per BLOCK, amplified by 2:1 smoothing and, under
`refine_body`, by a one-block 26-neighbour buffer. It does not make cells more
expensive; it makes you carry more of them. For the boundary-layer case it
currently costs nothing (1 of 4 y-tiles ≈ the 3 of 11 the cubic layout used),
but for body-fitted cases (naca, sd7003) it is where cells are wasted.

**(c) Multi-rank exchange — LARGE, and the whole remaining prize.**
`results_exchange_diag_2026-08-28.md`: the refined config sends **12.8× the
peer points** of the single-level config at 2 ranks while carrying 44 % of the
cells; the exchange goes from ~5 % of the step at 1 rank to **17 % at 2 ranks**,
`mpi_wait` alone 7–8 %, and the projection exchange is **71–75 % of all time
lost to imperfect scaling**.

**(d) Red-black's interface convergence loss** — only if red-black is pursued.

## Why (c) happens — the diagnosis is precise

**In `refine_dims = xz`, the 2:1 interface always sits ON a rank boundary.**
Cross-level peer points are *identical* at 2 and 4 ranks — 923 700 vs 923 692,
not a function of rank count at all. The xz mixed Morton key puts the y tile in
bits 42+, so every fine block precedes every coarse block in the leaf table: the
level change is one contiguous cut in Morton index and a linear split cannot
avoid it. Every interface transfer that could be a device-local copy is an MPI
message instead, and peer traffic costs ~6× a local copy per point here.

This is a **partitioning** defect, and P1 below has now measured that it is
**specific to xz mode** — the xyz octree curve does not have it.
It is upstream of the phi exchange (45 % of GPU bytes), of both velocity shells,
and of the pack/unpack that serve them.

## The plan, ranked

### P1 — DONE (`4f04f6b`, timing 2026-09-01): −3.63 % of the 2-GPU-rank step.

The key change is in: `leaf_key` (blocks.f90) puts the y tile in the low 21
bits for xz mode; xyz untouched. Measured on the wire at 2 GPU ranks, refined
production case:

| | old key | new key |
|---|---|---|
| peer send points | 941 700 | **46 000** (20.5×, prediction 20.3×) |
| per-rank send min/max | 193 700 / 748 000 | **23 000 / 23 000** |
| same-level (copy-only) fraction | 1.9 % | **98.7 %** |
| full nv=4 message | 30.134 MB | **1.472 MB** |
| cross-level peer points | 923 700 | **583** |
| local copy points | 5 482 980 | 6 378 680 (+16 %) |

Gates: 7-case suite bit-exact CPU **and** GPU (xyz path untouched); the xz
renumbering is physics-invariant (old vs new key, `max_abs 0` on all fields
reassembled onto the global lattice); 1 rank == 4 ranks EXACT and CPU == GPU
EXACT with the new key; `block_nb`'s xz uniform-flow gate still EXACT on an
unchanged leaf set.

**Timing — DONE 2026-09-01, clean machine**
(`overheadTest/results_xzkey_2026-08-29.md`). GPU, 2 ranks, 400 steps a side,
zero foreign compute apps before and after each run, drift flat:

| | ref | cand | |
|---|---|---|---|
| **step** | 0.323940 | 0.312171 | **−3.63 %** |
| `mpi_wait` | 0.021820 | 0.003003 | **7.3×**, 6.74 % → 0.96 % of step |
| total exchange | 0.053842 | 0.037001 | 16.6 % → 11.9 % of step |

Runtime diagnostics identical between the binaries at every output step — the
renumbering is physics-invariant on the production case at 2 ranks.

**The estimate was 6 %; the answer is 3.63 %, and the gap is instructive.** The
traffic prediction was exact (20.3× predicted, 20.5× measured) but two things
eat the conversion to time: `local_copy` rises +0.97 % of the step (the traffic
moved onto the device, and a local copy is cheap but not free), and the volume
kernels get 0.37 % slower because the new leaf order changes which blocks are
memory-adjacent. Together they give back 26 % of the −5.20 % exchange saving.
**A byte model predicts traffic, not time.**

Earlier, on a contended machine, a CPU 4-rank A/B gave `mpi_wait` 7.7716 →
0.0853 s/step, a **91×** collapse. The gap between 91× (CPU) and 7.3× (GPU) is
the straggler: the old key gave the CPU decomposition a 9.5× send-volume spread
across 3 peers and `Waitall` pays for the slowest, whereas the GPU case has one
peer per rank and only volume to remove. **Leaf-count balance is not exchange
balance** — compute balance was 1.000 throughout; read `send pts/rank min … max`
from the size report on any new decomposition.

The analysis that produced this, kept because it is what makes the change
defensible:

### The analysis

`tools/partition_analysis.py` reads a leaf table (a `blocks` dataset, or
`leaftable_test` stdout), rebuilds the 26-neighbour block graph in finest-cell
index space and scores candidate partitionings by **shared area**, which is
what an exchange entry carries. Validated against the solver: its proxy is a
uniform 1.30× of the measured peer points at 2 ranks, for both the cross-level
and the total figure.

**Result 1 — the defect is specific to `refine_dims = xz`, and it is not
general.** Cross-level shared area cut by the current linear split:

| case | mode | 4 ranks | 8 ranks | 16 ranks |
|---|---|---|---|---|
| boundary layer, wall band | xz | 87.5 % | 100 % | 100 % |
| compact body-like patch | xz | 0.4 % | 27.9 % | 28.8 % |
| **NACA 0012 `refine_body`, 25 418 leaves, 5 levels** | **xyz** | **1.4 %** | **4.2 %** | **10.1 %** |

In **xyz (octree) mode the standard Morton curve already keeps cross-level
pairs local** — there is nothing to fix for the cylinder/naca/sailplane family,
and a column-wise partitioner is actively *worse* there (total peer area 3.24 M
against 2.57 M at 4 ranks, and load imbalance 2.6 at 16). In **xz mode** the
mixed key puts the y tile in bits 42+, so the curve has no locality along y at
all, and the interface — which in xz mode is necessarily a y-plane — is cut
almost entirely.

**Result 2 — the fix is to move the y tile to the LOW bits of the xz key.**
The leaf order then runs down each (x,z) column before moving on, so **the
existing closed-form `zorder_owner` linear split becomes column-wise for
free**: no partitioner, no ownership table, no balance heuristic, and load
imbalance stays exactly 1.000 (an explicit column-wise partitioner reaches
1.20–2.61 and is strictly worse).

Total peer area, current → y-tile-in-low-bits:

| case | 4 ranks | 8 ranks | 16 ranks |
|---|---|---|---|
| boundary layer | 1 546 032 → 106 728 (**14.5×**) | 2 535 192 → 249 032 (**10.2×**) | 2 657 600 → 533 640 (**5.0×**) |
| compact xz patch | 922 112 → 136 976 (6.7×) | 2 119 502 → 333 056 (6.4×) | 4 283 011 → 717 648 (6.0×) |

with cross-level cut falling to 0.1–0.7 % (boundary layer) and 0.4–3.5 %
(patch).

**Expected step-time gain**, from the 2-rank GPU profile: `mpi_wait` is 7.2 % of
the step and would nearly vanish; the traffic moves into `local_copy` at ~1/6
the per-point cost (+0.4 %). Net **≈ 6 % of the 2-rank step, growing with rank
count**, and doubling again if red-black is adopted, since it exchanges per
colour.

**Cost, and it is the real decision.** Changing the key changes the canonical
leaf id order, so the `blocks` table row order changes. What mirrors it:
`moby_prepare` (writes it), the solver (cross-checks it at read),
`make_channel_restart`, and every existing xz restart/case file. mobygeom is
retired, so the mirror set is small — but existing xz files need either
regeneration or a version flag on the `refine_dims` attribute. **Gate it on
`refine_dims = xz` only**: xyz must keep its current key, which Result 1 shows
is already right.

Verification once implemented: the size report (`[output] profile` → `comm.f90
report_exchange_sizes`) makes it directly measurable — watch the copy-only
prefix fraction rise from 1.9 %, and the cross-level peer count stop being
rank-count-independent.

### P2 — overlap the exchange with compute

Not 2:1-specific, but the exchange it hides is precisely what the interface
inflates, so it pays most on refined multi-rank. Gated on the **A0 progress
probe** (`docs/next_session_multirank_exchange.md`): the entire plan rests on
an in-flight `Isend` progressing while a target kernel runs, which is
unverified on GPU and demonstrably false on CPU (0.46 GB/s). Ten throwaway
lines decide it. If it holds, `jacobi_apply` kernel 1 reads no phi halo at all
and kernel 2 needs only three low halo planes, so ~95 % of the work can run
during the exchange.

P1 and P2 compose: P1 removes bytes, P2 hides what is left.

### P3 — refinement granularity: per-level block size

The structural fix for (b): let a coarse zone use big blocks while fine blocks
stay small, so the wall band can be thin without paying the halo tax
everywhere. `docs/next_session_redblack_interface.md` §7 sketches it — one array
(or offset range) per level, volume kernels launched per level, **not** padded
to `nb_max` (fine blocks dominate; padding wastes 8× memory).

Worth it only when a case shows the granularity actually costing cells. Measure
first: for a candidate case, compare the leaf cell count against the cell count
an ideal (non-block-quantised) refinement of the same region would need. For the
boundary layer today that ratio is ~1, so **do not start here**.

### P4 — level-ordered smoothing — RED-BLACK ONLY

Jacobi's convergence is not degraded by the interface (ratio 0.998), so this is
**not** a general 2:1 fix. It matters only if red-black is pursued, where the
interface costs 1.91× → 1.49× of its per-iteration advantage and pushed the
refined-grid win down to 6.2 % at equal residual. §5's ladder: clamp `omega` at
interface-adjacent cells, or relax coarse → patch → fine so cross-level coupling
becomes multiplicative, at one extra phi exchange per colour.

## Closed — do not reopen

- **R0 / redundant phi halo** — reverted in `e77c75e`; block metrics are not
  bitwise equal across a periodic seam, so it is not bit-exact.
- **"Send less phi"** — 98.1 % of phi bytes at 2 GPU ranks are cross-level,
  which is exactly what an interface must transfer. It is a *consequence* of the
  P1 defect, not an independent lever.
- **Message batching** — already done: one `MPI_Isend` per peer carrying all
  entries × all variables, copy-only as a prefix of the same buffer. `mpi_post`
  is 0.03 % of the step.
- **The interface transfer arithmetic** — uniform oblique flow through a
  3-level patch is EXACT (0.0, `pn` spread 0.0) on both smoothers, 1 rank == 4
  ranks is EXACT, CPU == GPU is EXACT. The operators are consistent; the cost is
  in where the data has to travel, not in what is computed.
- **The interface kernels at one rank** — 1.4 % per cell is the whole budget.

## Measurement discipline for any of this

Same-host same-session baselines only (memory `timing-runs-need-drift-check`);
multi-GPU needs `overheadTest/gpu_rank.sh`; and **re-baseline at `niter = 12`**
— every recorded number in `overheadTest` is `niter = 6` Jacobi, which does not
keep the pressure zero-mode away.
