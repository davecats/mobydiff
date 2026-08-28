# Next session — reducing the 2:1 penalty

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

**The 2:1 interface always sits ON a rank boundary.** Cross-level peer points
are *identical* at 2 and 4 ranks — 923 700 vs 923 692, not a function of rank
count at all. `refine_dims = xz` puts the y tile in the high Morton bits, so
every fine block precedes every coarse block in the leaf table: the level change
is one contiguous cut in Morton index and a linear split cannot avoid it. Every
interface transfer that could be a device-local copy is an MPI message instead,
and peer traffic costs ~6× a local copy per point here.

This is a **partitioning** defect. It is upstream of the phi exchange (45 % of
GPU bytes), of both velocity shells, and of the pack/unpack that serve them.

## The plan, ranked

### P1 — partition so the interface is not cut. The only large 2:1-specific lever.

Ceiling: most of `mpi_wait` becomes `local_copy` at ~1/6 the per-point cost,
≈ 6 % of the 2-rank step, growing with rank count — and doubling again if
red-black is ever adopted, since it exchanges per colour.

**Start with an offline analysis that costs no solver work at all.** The leaf
table is in every case file (`blocks`: origin + level per global id). A script
that reads it and counts, for a candidate partitioning, how many cross-level
neighbour pairs would be cross-rank answers the whole design question before a
line of Fortran is written:

- how much of the 923 700 is avoidable at 2, 4, 8, 16 ranks;
- what a level-aware split costs in load imbalance (blocks per rank spread);
- whether the answer depends on `refine_dims`, on `refine_body` vs box
  refinement, and on rank count.

Only then choose the route, because both are structural:

- **Change the split, keep the curve** — assign blocks so cross-level pairs stay
  co-resident. Costs the closed-form `zorder_owner/start/count` O(1) ownership
  lookup that the exchange build relies on; it would become a table.
- **Change the curve** — interleave the y tile instead of putting it in the high
  bits. Costs the canonical leaf-table id order, which `moby_prepare`, restart
  files and `make_channel_restart` all mirror.

The size report added on 2026-08-28 (`[output] profile` → `comm.f90
report_exchange_sizes`) makes the before/after directly measurable: watch the
copy-only prefix fraction rise from 1.9 %.

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
