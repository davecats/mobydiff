# Handout — the 2:1 path at scale: it is a DEVICE-LOCAL exchange problem

> **STATUS 2026-09-10 — PHASES 1 AND 2 DONE. THE ANSWER IS NONE OF THE THREE
> PRE-REGISTERED BRANCHES.** `overheadTest/results_horeka_exchange_2026-09-10.md`,
> jobs 5139461 / 5139581.
>
> Phase 2 first: **A0 fails again at 8 ranks** — 5.7 ms of compute between the
> posts and the `Waitall` leaves `mpi_wait` at 1.06–1.44x of baseline and adds
> exactly 39 x 5.7 ms to the step. **P2 (overlap) is closed**, and the 2-rank
> probe's scope caveat with it.
>
> Phase 1 retires all three readings in the table below. Branch 1 and 2 both
> assume there is extra volume to find and there is not: the refined case
> exchanges **0.954x the points PER CELL** of its single-level twin, and
> cross-level points are **16.43 % of the total at 1, 4 and 8 ranks alike**.
> Branch 3 ("the per-point cost differs between op kinds") was tested with the
> instrument it asked for — the cross-level kernel timed in its own `copy_cross`
> bucket — and is **refuted**: 0.114–0.115 ns per cross-level point against
> 0.097–0.107 per same-level one, **1.15x**.
>
> **What it is instead: kernel LAUNCHES.** ~141 per step (39 pack + 39 unpack +
> 39 same-level copy + 24 cross-level copy), each carrying 56–114 us before it
> moves anything, against the projection's `jacobi_compute_phi` at 19 us. That is
> 9.3 ms/step for `rect_jacobi` and 11.6 for `refined_yp82` — 79 % of the refined
> case's device-local exchange, and the entire difference between twins that pay
> the same launch bill on 138 M and 61 M cells. The 2:1 interface's own share is
> 2.03 ms of launch against 0.36 ms of transfer.
>
> **Phase 3 is therefore NOT P3 (per-level block size)**: bigger coarse blocks cut
> the volume term, which is the smaller half, and do not change the number of
> rounds at all. Next step is an **nsys per-kernel timeline** of what the 60–100 us
> consist of — read `overheadTest/horeka/exchange/PREREGISTERED.md` first, which
> records a plausible mechanism ("~5.5 us per map-clause item") being refuted by a
> control before it could be published.
>
> Outstanding: the 16-rank repeat (the 4-node queue was three days out; 1/4/8
> ranks in one allocation is what this rests on, with the committed 16-rank
> aggregates agreeing).


Written 2026-09-10 from the session that fixed the rank-to-GPU mapping. Read the
first section before planning anything: the target has moved again.

## Start here

```bash
cd $WS/moby-2to1-code
git fetch origin && git checkout optimiseBlockRefinement_parentBoundaryLayer
git pull --ff-only          # you need 82b6845 or later
```

Read, in this order — each closes a door you would otherwise re-open:

1. `tutorials/turbulentBoundaryLayer/overheadTest/results_horeka_2026-09-10.md`
   (sections 6-9: the mapping fix, `moby_tune`, the re-measured campaign)
2. `docs/next_session_2to1_penalty.md` STATUS header
3. `docs/next_session_multirank_exchange.md` §A0 — the overlap probe, and the
   scope caveat that Phase 2 below discharges
4. `CLAUDE.md`

**`select_target_device` is assumed to have passed the workstation gate.** If it
has not, run it before trusting any timing here.

## What is closed, and must not be reopened

| | why |
|---|---|
| The 2:1 per-cell machinery cost | 0.98 coarse-cell-equivalents at 4 ranks, unchanged across the mapping fix. It is free. |
| The rank-split / partitioning rewrite (P1) | The stall it targeted was GPU affinity. Block tax at 16 ranks is now 1.130, not 1.505. |
| `mpi_wait` as a 2:1 problem | 7.7 % of the refined step at 16 ranks, and the same for blocked and unblocked. It is a general exchange cost, not a refinement one. |
| The halo copy KERNEL | At its ~51 %-of-peak coalescing floor. Copy FEWER points, do not copy them faster. |
| Redundant-computation schemes | Blocked by the periodic-seam metric asymmetry (`e77c75e`). |
| Placement sweeps | Settled. `moby_tune` answers it in 3 minutes on a new machine. |

## Where the time actually goes now

`refined_yp82_rect_jacobi`, 16 ranks, corrected mapping, job 5139351:

| phase | s/step | % |
|---|---|---|
| projection | 0.03181 | **77.2** |
| — `apply` | 0.01111 | 27.0 |
| — `phi_exchange` | 0.00841 | 20.4 |
| — `vel_exchange` | 0.00721 | 17.5 |
| — `sweep` | 0.00346 | 8.4 |
| momentum | 0.00446 | 10.8 |

and the exchange, summed over every call site:

| bucket | s/step | % of step |
|---|---|---|
| `local_copy` | 0.00585 | **14.2** |
| `pack` | 0.00387 | **9.4** |
| `unpack` | 0.00334 | **8.1** |
| `mpi_wait` | 0.00318 | 7.7 |
| **total exchange** | **0.01648** | **40.0** |

**The exchange is 40 % of the step and four fifths of it is device-local.**
`pack + unpack + local_copy` = 31.7 %; MPI is 7.7 %. This confirms at 16 ranks
what A0 concluded at 2 (`next_session_multirank_exchange.md`): *the exchange is a
device-local copy problem, not a communication one* — and that conclusion was
explicitly scoped to "repeat on a genuinely many-rank GPU machine" before being
carried. It now has been, on the timing side.

### The 2:1 cost at scale, isolated

Same block shape, same solver, same placement, 16 ranks:

| | cells | exchange s/step | exchange per cell |
|---|---|---|---|
| `rect_jacobi` (single level) | 138.41 M | 0.02009 | 1.00x |
| `refined_yp82` (2:1) | 60.56 M | 0.01648 | **1.9x** |

The refined case does 82 % of the single-level case's exchange work while
carrying 44 % of its cells. **That factor ~1.9 is the 2:1 interface's cost at
scale, it is device-local, and it is the target.** Note it does NOT show up as a
per-cell compute cost (0.98x at 4 ranks) — only in the exchange.

## Phase 1 — attribute the 1.9x (measurement, ~1 h)

`report_exchange_sizes` already prints totals. It does not say how many points
belong to same-level COPY entries against cross-level RESTRICT/PROLONG, which is
exactly the split that decides what to do.

Extend it to break `local copy pts`, `send pts` and the ENTRY COUNT down by op
(`OP_COPY` / `OP_RESTRICT` / `OP_PROLONG`), then run `refined_yp82` and
`rect_jacobi` at 1/4/16 ranks. Diagnostics-only, bit-exact by construction, still
gate it.

**Pre-registered reading:**

| result | conclusion |
|---|---|
| cross-level points are a small fraction and the 1.9x is COPY volume | the cost is block-count granularity, not the interface. Go to P3 (per-level block size) — the refined case carries the same 135 k cells/block as the coarse one and does not need to. |
| cross-level points dominate | the interface transfer itself is the cost. Attack the entry structure: 4 fine sub-entries per coarse face, the tq-aware covering rows, the blended ghosts. |
| neither — the points match but the TIME does not | the per-point cost differs between op kinds (gather maps, strided reads). Profile the three kernels separately. |

## Phase 2 — repeat A0 at 8 and 16 ranks (~1 h)

A0 (2026-09-02, 2 GPU ranks) inserted a compute-bound kernel between the
`Isend`/`Irecv` posts and the `Waitall` and found **nothing was hidden** —
9.95 ms of work bought 0.05 ms of wait. Its own text says to repeat it on a
many-rank GPU machine before carrying the conclusion. Do that now: the message
sizes, peer counts and transports at 16 ranks are all different, and `mpi_wait`
is 7.7 % of the step rather than 0.96 %.

Ten throwaway lines, then delete them. **If transfers still do not progress
during a kernel, close P2 (overlap) for good and say so in
`next_session_2to1_penalty.md`** — it is currently the only item that document
leaves standing, and it should not stay open on a 2-rank measurement.

## Phase 3 — only after 1 and 2

Choose the lever the measurements point at. Do not start any of these blind:

- **P3, per-level block size** (`next_session_redblack_interface.md` §7). The
  structural fix if Phase 1 says COPY volume: let coarse zones use big blocks
  while fine blocks stay small, so the wall band is thin without paying the halo
  tax everywhere. Launch volume kernels per level, **not** padded to `nb_max`.
- **Fewer exchange rounds.** 39 per step, invariant across every configuration
  ever measured here: 3 RK stages x (6 phi + 6 vel) + 3. That is a numerics
  change, not a scheduling one — it needs a convergence argument, not just a
  bit-exactness gate.
- **Fusing pack/copy/unpack** across the per-colour exchanges, if Phase 1's
  third branch fires.

## Rules that still bind

- **Bit-exactness for every scheduling change**: `-Mnofma` / `-gpu=nofma`, the
  7-case suite, CPU AND GPU. `tools/h5maxdiff` does the comparison where h5py is
  unavailable (HoreKa). A refactor that is not bit-exact is a bug.
- **State the placement with every number**, and run `moby_tune` once on any
  machine that is not HoreKa before believing a multi-node timing.
- **Same allocation for anything compared.** Node-to-node variation is ~6 %;
  a 4 % winner at one rank count lost at another in this very campaign.
- **Judge a run from `runtime.txt` drift**, not the final cumulative average.
- **`L2_div` identical** across placements and rank counts, or the run is broken.
- Every `exch_timing` number is **rank 0's**; rank 0 sits at or near the
  cross-rank minimum. Use the `exchange balance:` line for the spread, and
  remember it is blind to skew whose late rank rotates.
- **Never quote a step time from a `BARRIER=1` run.**

## Deliverable

A dated `overheadTest/results_horeka_<date>.md` in the repo convention: the
tables, the mechanism each number supports, and explicitly what the numbers do
**not** support. Update the STATUS header of `docs/next_session_2to1_penalty.md`,
and close A0/P2 there one way or the other. Commit the raw run directories
(`git add -f`) and push.

If Phase 1 and Phase 2 together identify the lever, **stop and write it up.**
Implementing it needs a gate and a fresh session; the attribution is the valuable
half, and this campaign has now twice published a mechanism that the next
measurement overturned.
