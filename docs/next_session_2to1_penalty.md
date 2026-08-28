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

### P1 — RUN 2026-08-28. It is an `xz`-only defect and the fix is one key function.

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
