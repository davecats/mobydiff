# The blocked exchange is not slow: it is out of step

Job 5138715, 4x HoreKa Green nodes (hkn[0518,0604-0605,0708]), commit `df08a60`
(clean, solver source = `be48d44`), 26 runs of 200 steps in ONE allocation plus
a 6-run barrier pass. Raw logs in `moby-2to1-run/results_balance/` and
`results_barrier/`; driver `horeka/mechanism/submit_session.sh`.

**Headline.** The 708-753 us `mpi_wait` that three reports have called "the
blocked exchange's stall" is **not transport**. An `MPI_Barrier` immediately
before the `Waitall` absorbs **100 %** of it: the ranks are arriving at the
exchange up to ~773 us apart, and once they are synchronised the transfer
completes in **1.2 us**. The `base_jacobi` control -- the configuration that has
been held up all along as "fast" at 120 us -- has a **63.8 us** real transfer,
**53x more actual wire time than the case we were trying to speed up.**

## Deviation from the handout, and why

The handout asked for three added runs resumed on top of the 23 in
`results_mechanism/`. Instead the whole 26-spec matrix was re-run into a fresh
directory with the `be48d44` binary. The 23 came from two node sets, two days
and the pre-diagnostics binary; the re-run cost ~21 min of an otherwise idle 2 h
allocation and removes the cross-job caveat from EVERY comparison rather than
from `rect_jacobi:8:2` alone. It also yields the `exchange balance:` line for all
26 runs instead of the four of Phase 2(a).

Two things this exposed that a resume would have hidden:

- The earlier campaign's `build_gpu/moby_solve` was **25 MB**; the rebuild, with
  byte-identical CMake configuration (`USE_OPENMP_OFFLOAD ON`, `-mp=gpu`, same
  compilers and HDF5), is **4.8 MB**. Consistent with the earlier binary having
  been built where no GPU was visible, so nvfortran emitted every CC target
  instead of cc80 alone. Stated as an observation, not a verified cause.
- `setup_mechanism.sh` still pins `COMMIT=a11e355` and would have detached the
  code tree onto the pre-diagnostics commit, and `submit_mechanism.sh` reuses
  `build_gpu/` whenever the executable exists -- which was exactly that stale
  binary. Neither was run.

## Validation

- **Placement honoured** in all 32 runs: requested node count == `hosts.txt`, no
  `VOID`, no failures.
- **`L2_div` identical** across every placement of every config at every rank
  count.
- **39 exchange rounds per step** in every run, as in all previous matrices.
- **Step time reproduces the earlier campaign to 0.5 %** on nearly every spec
  (worst 1.5 %), across two binaries and three node sets -- so the binary change
  did not move the physics timing.
- **The anomaly is reproducible across three jobs, three node sets, two binaries
  and two mapping policies**: `rect` 8x2 = 708.3 (`--map-by numa`, job 5133554),
  752.6 (`ppr:4:node`, job 5135435), **738.5** (here). Large waits reproduce to
  2-4 %; waits below ~100 us are noisy (up to 2.6x between jobs) and should not
  be quoted to better than a factor.

## 1 — The barrier split: it is arrival spread, not transport

`[output] exchange_barrier = true` puts an `MPI_Barrier` before every `Waitall`.
Per round, 8 ranks:

| config | placement | pass-A `mpi_wait` | `skew_barrier` | residual `mpi_wait` | reading |
|---|---|---|---|---|---|
| `rect_jacobi` | 8x2 | 738.5 us | **773.0** | **1.2** | all spread |
| `nb16_jacobi` | 8x2 | 858.5 | **878.6** | **1.25** | all spread |
| `base_jacobi` | 8x2 | 120.1 | 64.8 | **63.8** | half real transfer |
| `blk8_jacobi` | 8x2 | 123.4 | 83.1 | **49.1** | half real transfer |
| `rect_jacobi` | 8x4 | 149.8 | 160.0 | 7.2 | spread |
| `rect_jacobi` | 4x1 | 60.8 | 63.8 | 8.7 | spread |

The barrier is **non-perturbing in total**: `skew + residual` reproduces the
pass-A wait to 5-10 % in every row, so this is a split of the same quantity, not
a different measurement.

**What this pins down and what it does not.** A barrier costing 773 us means the
ranks genuinely arrive ~750 us apart -- an 8-rank barrier over already-arrived
ranks costs tens of microseconds, not hundreds -- so the skew column is a direct
measurement. The residual is only a LOWER bound on wire time: a transfer that
completed during the barrier leaves ~0 behind, so `rect`'s 1.2 us proves the
transfer **is not the critical path**, not that it would be 1.2 us in isolation.
The informative direction is the large residual: `base` still had 63.8 us of
transfer left *after every rank had arrived*, and `rect` had none.

Either way the arrival spread alone accounts for the anomaly. No transport
effect needs to be invoked, and Occam says none should be.

## 2 — Phase 1: the copy question, closed

`nb16_jacobi:8:2` was the pre-registered direct test -- the anomaly's exact
placement, the same 2-peer chain, 3.4x `rect`'s copy load per rank.

| config @ 8x2 | local copy | `local_copy` us/round | `mpi_wait` us/round |
|---|---|---|---|
| `base_jacobi` | 0 | 0.2 | 120.1 |
| `blk8_jacobi` | 0 | 0.1 | 123.4 |
| `rect_jacobi` | 14.87 M pts | 237.3 | **738.5** |
| `nb16_jacobi` | 56.42 M pts | 812.8 | **858.5** |

Copy PRESENCE separates the fast pair from the slow pair perfectly at this
placement. Copy VOLUME does not: 3.4x the copy buys 1.16x the wait. And the same
`rect` config spread 2-per-node carries an **identical** 237.4 us copy load and
waits 149.8 us. So the slow cell needs all three of {4 ranks/node, off-node
traffic, device-local copy present}, and the copy's role is structural, not
volumetric.

This is coherent with the barrier result, because `copy_local_entries` sits
**inside the arrival path**: the exchange is `pack -> post Isend/Irecv ->
local copy (synchronous device kernels) -> Waitall`. A rank reaches the `Waitall`
only after its copy kernel has completed. That is why only the copy-carrying
configs show a large spread -- but it is NOT sufficient, since `rect` 8x4 has the
same copy and one fifth the spread.

**What this does NOT support.** It does not identify what makes co-resident ranks
drift apart. Copy presence, ranks-per-node and off-node traffic are correlated
across only four cells, and a fourth mechanism built on that correlation would be
the same mistake as the three already refuted. The measurement that separates
them is a timeline, not another aggregate.

## 3 — Two facts about the instrument itself

**Every `mpi_wait` figure in all three previous reports is rank 0's**, and rank 0
sits at or near the *minimum* across ranks in most runs. The true cross-rank mean
is 4-25 % higher. The `rect` 8x2 anomaly, for instance, is 738.5 us on rank 0 and
768.0 us on average.

**Aggregate min/max/argmax is blind to rotating skew.** `rect` 8x2 has
`max/min = 1.09` -- every rank waits equally -- which the pre-registered decision
table reads as "transport, partitioning will not help". The barrier says the same
run is 100 % arrival spread. Both are true: the per-round spread is ~773 us but
the *late rank rotates*, so per-rank totals equalise over 7800 rounds. A
reduction over totals cannot see this, and the decision table's first row is
therefore wrong as stated. Corrected: `max/min ~ 1` distinguishes *persistent*
from *rotating* skew, and says nothing about transport without the barrier.

## 4 — Placement, re-measured in one allocation

Per-round `mpi_wait`, pass A, all at one machine state (rank-0 bucket, for
comparability with the earlier reports; add 4-25 % for the cross-rank mean):

| config | peers | 4 ranks/node (8x2) | 2 ranks/node (8x4) | effect of spreading |
|---|---|---|---|---|
| `base_jacobi` | 7 | 120.1 us | 1621.6 us | **13.5x worse** |
| `blk8_jacobi` | 5 | 123.4 | 2512.0 | **20.4x worse** |
| `nb16_jacobi` | 2 | 858.5 | 175.3 | **4.9x better** |
| `rect_jacobi` | 2 | 738.5 | 149.8 | **4.9x better** |

The inversion reported on 2026-09-09 reproduces on a third node set with a
different binary. End-to-end, `rect` at 8 ranks is 0.13622 s/step packed and
0.11176 spread -- **18 % faster on twice the nodes**.

## 5 — The timeline probe: not device work, and every round

Job 5139026, 2 nodes (hkn[0525-0526]), Nsight Systems 2025.1.1, `--trace=cuda`,
40 steps, one report per rank. Reduced by `mechanism/collect_nsys.py`.

**Two things had to be settled that no timer in this project can settle**, because
both stalls are booked as `mpi_wait` either way: whether the `Waitall` is really
blocked on an unfinished CUDA stream, and whether "39 rounds at 740 us" is a fair
description or a handful of catastrophic rounds hiding in a mean over 7800 calls.

| run | GPU busy | idle | gaps > 50 us | median gap | share of idle | slowest 1 % of rounds |
|---|---|---|---|---|---|---|
| `rect` 8x2 (anomaly) | **62.7 / 62.9 %** | 37 % | 1098 / 1877 | **365 / 349 us** | 78 % | **4-5 %** |
| `base` 8x2 (control) | 79.0 / 79.8 % | 21 % | 1597 / 1673 | 71 / 72 us | 52 % | 16-17 % |
| `rect` 4x1 (no crossing) | 90.1 % | 10 % | 554 | 98 us | 33 % | 21 % |

- **The stall is NOT device work.** The GPU is idle 37 % of the time in the
  anomaly regime; a rank blocked on its own CUDA stream would show the GPU busy.
  The hypothesis that "mpi_wait" was really unfinished offload -- which would have
  invalidated both the partitioning and the transport plans -- is refuted.
- **It is every round.** The slowest 1 % of rounds hold only 4-5 % of `rect` 8x2's
  wait (against 16-21 % in the two faster cases, which are the ones with a spiky
  tail). The idle piles up in ~1100 gaps of median 365 us, which is the same order
  as the 773 us arrival spread the barrier measured directly.

**READ SHAPE, NOT MAGNITUDE.** Tracing inflates the wait, and unevenly: `base`
goes 120 -> ~690 us and `rect` 738 -> ~1255 us. No absolute number in this section
may be compared with an untraced run; only ratios within one traced run and the
distribution shape are used above.

### Two instrument failures worth recording

- **nsys `--trace=mpi` captures nothing from this solver.** OpenMPI's Fortran
  `mpi_f08` bindings reach the C layer as `PMPI_*`, so nsys's interception of the
  `MPI_*` symbols never fires: `does not contain MPI event data`. An LD_PRELOAD
  shim misses it for the same reason. Restoring the MPI view needs solver-side
  NVTX, which is a gated change. The kernel signature `pack -> copy_local ->
  [gap] -> unpack` recovers the per-round wait without it.
- **Under nsys the solver segfaults at teardown**, after the main loop, the
  timing lines and the final field write. The report is written intact; the damage
  is that mpirun kills the siblings mid-report, so only some ranks survive (6 of 8,
  4 of 8, 2 of 4 here). `run_nsys.sh` now wraps each rank so it exits 0 and judges
  the run on `main loop ended` rather than on mpirun's status.
- Also visible, and not chased: **~1300 host-to-device memcpys per step**, ~2 us
  each, ~3-4 % of wall in every configuration including the unblocked one. Sub-MB,
  so almost certainly OpenMP target per-launch argument marshalling rather than
  data staging. Same in `base` and `rect`, so it is not the mechanism -- but it is
  a standing 3-4 % that nobody has looked at.

## What is still open

The mechanism is located, bounded and characterised, but its SOURCE is not
identified: **what desynchronises four co-resident ranks by ~750 us per round
when they carry a device-local copy and the node also drives off-node traffic?**

What the source is NOT, on evidence: not transport (barrier), not device work
(GPU idle 37 %), not a few bad rounds (slowest 1 % hold 4-5 %), not a persistent
straggler (per-rank totals equal to 9 %), not copy volume (3.4x buys 1.16x).

The cheapest decisive follow-up after that is a **second barrier site, before the
`pack`**, which would separate skew that the exchange itself creates (the copy)
from skew the rank brings in from upstream. That is a solver change and needs the
workstation bit-exactness gate, so it is not made here.

## What this changes

`docs/next_session_2to1_penalty.md` P1 (repartition the Morton chain) and the
transport-tuning options are both aimed at a transport cost that **does not
exist**: the blocked exchange moves its bytes in ~1 us once the ranks are
together, against 64 us for the unblocked decomposition it was losing to. The
19-25 %-of-step prize quantified on 2026-09-07 is real, but it is a
**synchronisation** prize, and neither a better partition nor a faster transport
will collect it.
