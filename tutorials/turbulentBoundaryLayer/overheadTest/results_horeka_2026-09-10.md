# The block tax was GPU-to-NIC affinity: 20.6 % of the step, from a launch flag

Job 5138715, 4x HoreKa Green nodes (hkn[0518,0604-0605,0708]), commit `df08a60`
(clean, solver source = `be48d44`), 26 runs of 200 steps in ONE allocation plus
a 6-run barrier pass. Raw logs in `moby-2to1-run/results_balance/` and
`results_barrier/`; driver `horeka/mechanism/submit_session.sh`.

**Headline.** The 708-753 us `mpi_wait` that three reports called "the blocked
exchange's stall" is **not transport, and not the exchange's fault at all**. It is
GPU affinity: **the two ends of a cross-node link must sit in the same GPU
affinity class**, and `comm.f90`'s `device = local_rank mod num_devices` pairs a
node's LAST local rank with the next node's FIRST -- dev3 against dev0 -- which is
the mismatched case.

**Placing them consistently is worth 20-27 % of the step**, and it is now
automatic (`select_target_device`, section 7): 8 ranks -20.4 %, 16 ranks -26.3 %,
16 ranks refined -26.7 %, with the fields bit-identical (`max_abs 0`) and
single-node runs untouched.

The consequence for this campaign's central claim: **"the block tax explodes at
the node boundary" was this artifact.** `rect`/`base` at 16 ranks is 1.505 with the
old mapping and **1.098** with the new one.

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

## 6 — The cause: matched affinity classes, not "a NIC-affine GPU"

The GPU traces pointed here. Per-rank memcpy in `rect` 8x2, steady window:

| rank | H2D MB | D2H events | P2P ms | role |
|---|---|---|---|---|
| 0,1,2,5,6,7 | 25.3 | ~3195 | 18-32 | intra-node peers only |
| **3** | **694.2** | **84 945** | 13.6 | holds the cross-node link |
| **4** | **676.3** | 3 201 | 14.3 | holds the cross-node link |

694 MB over the steady window is ~57 MB/step -- that configuration's ENTIRE
exchange volume (1.47 MB/round x 39). The two cross-node ranks do a
device -> host -> network -> host -> device round trip; the six intra-node ranks
move 25 MB and go GPU-direct over NVLink. `nvidia-smi topo -m` shows why the two
cards differ: all three HCAs sit on NUMA 0 with GPU0/GPU1 (`NODE`), while
GPU2/GPU3 reach every NIC only across the inter-socket link (`SYS`).

### A FIRST, WRONG READING, and the experiment that killed it

The obvious conclusion -- "both cross-node ranks need a NIC-affine GPU" -- is
**wrong**, and it was published in the first version of this file. Two things
refuted it.

**The confound.** The sweep that produced the original 10.6x compared a run with
4 devices visible and no pinning against runs that set `CUDA_VISIBLE_DEVICES` per
rank. That changed the permutation AND the device visibility together. Separating
them (2 binaries x 3 pinnings, one allocation, `rect` 8x2):

| binary | pinning | ndev | s/step | wait/round |
|---|---|---|---|---|
| ref | none | 4 | 0.136639 | 753.5 us |
| ref | identity `0,1,2,3` | 1 | 0.131482 | **604.4** |
| ref | boundary `0,2,3,1` | 1 | 0.108237 | **73.4** |

At fixed visibility the permutation is worth 8.2x on the wait and 21.5 % on the
step, so the effect is real -- but visibility is a genuine second effect
(753.5 -> 604.4) that the original number had folded in.

**The override.** `MOBY_GPU_ORDER=3,2,1,0` puts BOTH cross-node ranks on the
supposedly bad cards and ran FAST (72 us). Under "both ends need a local HCA"
that is the worst case.

### The corrected model, pre-registered and confirmed

Predictions written before the runs, then measured (`rect` 8x2, wait per round):

| ends of the cross-node link | class | predicted | measured |
|---|---|---|---|
| dev3 <-> dev0 | mixed | slow | 604.4 us |
| dev2 <-> dev1 | mixed | slow | **603.9** |
| dev3 <-> dev2 | both far from a NIC | fast | **80.8** |
| dev1 <-> dev0 | both beside a NIC | fast | **74.7** |
| new binary, first preference 0 / 1 / 2 / 3 (always matched) | matched | all fast | 81.9 / 70.9 / 95.4 / 101.9 |

**Matching is worth ~8x; being on the NIC-affine side is a further ~1.3x on the
wait, about 1 % of the step.** Both-far-but-matched is fine. The old
`local_rank mod ndev` pairs dev(ppn-1) with dev0 across every node boundary,
which is exactly the mismatched case.

Correctness is untouched throughout: `L2_div` = 1.07282926E-05 in every variant
and the exchange size report is identical (2 peers, 515 200 send pts, 14 874 752
local copy pts). Only the physical card changes.

### Why it hits the BLOCKED decomposition and not the unblocked one

A 2-peer Morton chain puts only two ranks per node on the network and both can be
placed consistently. A 7-peer Cartesian decomposition puts every rank on the
network, so most links cross classes whatever the map -- `base_jacobi` gains
nothing (and loses nothing) from any permutation. That is also the placement
inversion reported on 2026-09-08/09: `rect` improved when spread because 2
ranks/node lands everyone on GPU0/GPU1, hence matched, by accident.

## 7 — Tier 1: automatic, in the solver

`comm.f90 select_target_device` replaces `device = local_rank mod num_devices`.
The ranks holding cross-node links in a linearly split decomposition are a node's
FIRST and LAST local ranks, so they are served first: **local 0, local ppn-1, then
the rest in order, applied identically on every node**. With the default device
order that is `0->dev0, 1->dev2, 2->dev3, 3->dev1` -- the hand-tuned optimum,
derived from a rule rather than hard-coded. `MOBY_GPU_ORDER` overrides the device
order for machines whose good cards are not the first ones; single-node runs keep
the identity mapping exactly.

**Node-homogeneity is load-bearing.** The first version sorted by each node's
MEASURED off-node degree, which gives the first, middle and last nodes DIFFERENT
maps and hands a many-peer decomposition mixed-class pairs on most of its links:
`base_jacobi` at 16 ranks regressed 14 %. The fixed order removes that.

Reference binary vs new, back to back in one allocation:

| case | ref s/step | new s/step | gain | wait/round |
|---|---|---|---|---|
| `rect_jacobi` 8x2 | 0.139811 | 0.111295 | **-20.4 %** | 755.5 -> 77.5 us |
| `rect_jacobi` 16x4 | 0.094321 | 0.069542 | **-26.3 %** | 768.4 -> 159.3 |
| `refined_yp82` 16x4 | 0.060583 | 0.044433 | **-26.7 %** | 485.5 -> 86.2 |
| `base_jacobi` 16x4 | 0.062679 | 0.063319 | +1.0 % | 140.5 -> 154.7 |
| single node 4x1 | 0.221274 | 0.221739 | none | identity map |

Gates: fields **bit-exact** (`max_abs 0` on un/vn/wn/pn, 138 M points) against the
pre-change binary at 8x2 AND single node, compared with `tools/h5maxdiff` (HoreKa
has neither h5diff nor h5py). `base_jacobi`'s +1.0 % is the same order as the
+0.6 % measured for a hand-applied homogeneous map, i.e. noise-level, not the
14 % the degree-sorted rule cost.

**Still owed: the workstation bit-exactness gate** (nofma, 7-case suite). Nothing
in the arithmetic moves -- only which physical card a rank uses -- and the fields
match exactly here, so it is expected to pass, but it has not been run.

### Block tax, corrected

| | 1 rank | 8 ranks | 16 ranks |
|---|---|---|---|
| old mapping | 1.052 | 1.315 | **1.505** |
| new mapping | 1.052 | 1.058 | **1.098** |

## What is still open

The source is identified, corrected and fixed automatically (sections 6-7). What
remains:

- **The workstation gate** for `select_target_device` (nofma, 7-case suite).
- **Re-measure the campaign at the new mapping.** Every scaling number, block tax
  and strong-scaling efficiency in the 2026-09-07/08/09 reports predates it.
- **Decide the device order from the topology** instead of leaving
  `MOBY_GPU_ORDER` unset-and-hope: reading `/sys/class/infiniband/*/device/numa_node`
  against each GPU's PCI `numa_node` gives it without a vendor API, and only GPU
  enumeration needs one vendor call. On this machine the default order happens to
  be right, which is exactly the situation that hides a bug on the next machine.
- **The residual**: ~80-160 us/round remains, and the arrival spread has shrunk,
  not vanished. A second barrier site before `pack` would say whether what is left
  is the exchange's own copy or upstream drift.
- **Only measured here.** One machine, one node type. The rule assumes contiguous
  rank-to-node blocks and r+-1 peers; both hold for the block Z-order chain under
  `--map-by ppr:N:node`, neither is universal.

## What this changes

`docs/next_session_2to1_penalty.md` P1 (repartition the Morton chain) and the
transport-tuning options are both aimed at a transport cost that **does not
exist**: the blocked exchange moves its bytes in ~1 us once the ranks are
together, against 64 us for the unblocked decomposition it was losing to. The
19-25 %-of-step prize quantified on 2026-09-07 is real, but it is a
**synchronisation** prize, and neither a better partition nor a faster transport
will collect it.
