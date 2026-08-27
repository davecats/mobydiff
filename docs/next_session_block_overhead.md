# Next session — making the block/2:1-interface machinery cheaper

Handout. **Measure first, then optimise.** Every phase below is bit-exact-gatable
except the last, which is deliberately last.

---

## STATUS 2026-08-07 — Phase 0 DONE. The tax is EXCHANGE TRAFFIC, not footprint.

The 41 % is attributed, and the answer changes the ranking below. **Read this
section before the phase list; where they disagree, this section wins.**

Per-phase, GPU (A6000, jacobi, 400 cold steps, 1 rank,
`PROFILE=1 ./run_overhead.sh` + `phase_table.py`), seconds per step:

| phase | base | nb16 | ratio | nb8 | ratio |
|---|---|---|---|---|---|
| `momentum` (fused predictor) | 0.3633 | 0.3676 | **1.012** | 0.3818 | **1.051** |
| `proj/sweep` | 0.3077 | 0.3357 | 1.091 | 0.3670 | 1.193 |
| `proj/apply` | 0.4308 | 0.4533 | 1.052 | 0.5308 | 1.232 |
| `ibm_mu` — the control | 0.0371 | 0.0535 | **1.443** | 0.0706 | **1.904** |
| `proj/vel_exchange` | 0.0075 | 0.3100 | **41.4** | 0.6529 | **87.2** |
| `proj/phi_exchange` | 0.0054 | 0.1005 | **18.6** | 0.2083 | **38.5** |
| `step/vel_exchange` | 0.0012 | 0.0440 | **38.3** | 0.0911 | **79.2** |
| `exch/local_copy` | 0.0120 | 0.4522 | **37.8** | 0.9495 | **79.3** |
| **TOTAL** | 1.2412 | 1.7545 | **1.414** | 2.4052 | **1.938** |
| *predicted `(nb+2)³/nb³`* | | | *1.4238* | | *1.9531* |

**Of the +0.513 s/step that blocking costs at `nb = 16`, `local_copy` is +0.440
— 86 %** (at `nb = 8`: +0.938 of +1.164, 81 %). Everything else together is
14 %. The decision rule in Phase 0 below lands squarely on "exchange
pack/unpack/local_copy ⇒ **traffic-bound** ⇒ prioritise Phases 3, 5".

**The two-point test settles it.** `nb = 8` is the second point on the same grid
with the same binary (the grid admits only `nb ∈ {4, 8, 16}` — it must divide
4096, 176 *and* 192, so the handout's suggested `nb = 32` cross-check is
impossible here, and a second point on the SAME case is stronger anyway). The
control row `ibm_mu` follows the volume law at both sizes (1.443 vs 1.4238;
1.904 vs 1.9531) and so does the TOTAL (1.414, 1.938 — both ~0.8 % under), while
`momentum` stays at 1.012/1.051 throughout. The law predicts the total, but only
because halo cells are what gets *copied*; nothing computes on them.

Caveat on the redblack columns of `phase_table.py`: `base_redblack_prof` was
contended during both attempts (see the drift note below), so redblack RATIOS
from the profiled matrix are unreliable. Its total tax is measured cleanly by
the unprofiled matrix instead — 1.414, matching jacobi's 1.412, which is the
handout's original point that the tax is solver-independent.

Three things this pins down that the handout could not:

1. **The volume kernels do NOT pay the volume tax.** `momentum` is 1.012 — free.
   They loop interiors (`1:nb`) only, so the halo cells enlarge the footprint but
   never enter the work. The one kernel that IS pointwise over the full
   halo-carrying array, `update_ibm_mu`, lands on **1.443** against the predicted
   `(18/16)³ = 1.4238`. That control row both validates the instrument and shows
   why everything else escapes.
2. **The exchange barely exists without blocking.** At `nb` unset the whole
   halo machinery costs 0.012 s/step. Every millisecond of the 0.45 s at
   `nb = 16` is *created* by the block lattice. The cost is per halo CELL, and
   `(nb+2)³/nb³ − 1` is exactly the halo-to-interior cell ratio — so the
   volume law still predicts the total correctly, but through traffic, not
   footprint. Do not let the coincidence mislead the choice of fix.
3. **NEW LEVER, not in the handout: the copy kernel is slow in absolute terms.**
   The halo shell is 58.7 M cells (1736 per block × 33792). A velocity exchange
   moves 3 vars × 8 B (4 on the last iteration, which carries p), read+write
   ≈ 2.8 GB, in 17.2 ms ⇒ **~165 GB/s, about 21 % of the A6000's 768 GB/s
   peak**; phi moves 0.94 GB in 5.6 ms ⇒ 168 GB/s — the same number for a
   quarter of the payload, so it is the kernel, not the size of the message.
   Before moving fewer bytes, it is worth
   asking why the bytes move at a fifth of peak. Up to ~4× on 86 % of the tax
   with no change to the algorithm and trivially bit-exact.

Also measured: **the 2:1 interface is cheap.** The refined case scales at
0.494 of `nb16_jacobi` against a cell ratio of 0.455, and every phase scales
0.43–0.50 — i.e. essentially with cell count. Optimising interface transfer
remains chasing the small term.

### Revised phase order

1. ~~**The copy kernel itself**~~ — **DONE 2026-08-27, block tax 1.430 -> 1.370**
   (see the next section).
2. ~~**Phase 1 — per-direction `nb`**~~ — **DONE 2026-08-27, tax 1.370 -> 1.071**
   (see below).
3. **Phase 3 — fewer exchanges per iteration.** 3a (the `phi` halo) ATTEMPTED
   AND REVERTED — not bit-exact, and worth only 1.8 %. The `vel_exchange` half
   is **DONE 2026-08-27**: the minimal divergence sync, 3.4 % on the production
   layout and 9.0 % at `nb = 16`. Both STATUS sections below.
4. **Phase 5 — fuse same-level interior blocks.** The structural fix for the
   dominant term: fused tiles have no internal halos, so the traffic is not
   reduced but *deleted*. Now clearly worth the disruption. 5a (per-level block
   size) first, as written.
5. **Phase 4 — stop allocating unread halos.** DEMOTED to near-worthless.
   Footprint is not the constraint; `qs` is not read by an exchange, so freeing
   it saves allocation, not traffic. Do it only as tidy-up.
6. **Phase 6 — mixed precision.** Unchanged in position, but note it now has a
   second justification: it halves the exchange payload, which is 86 % of the
   tax.

---

## STATUS 2026-08-27 — copy kernel DONE. Block tax 1.430 -> 1.370.

The first item of the revised order is done: `comm.f90`'s same-rank halo gather
is split into a light same-level copy kernel (`copy_local_same_level`,
`copy_local_scalar_same_level`) and the general cross-level one
(`copy_local_cross_level`). **Pure scheduling change, bit-exact CPU AND GPU.**

**Why it was slow — measured with ncu, not guessed.** The old single kernel was
never bandwidth-starved: `Block Limit Registers = 4` gave 33 % theoretical
occupancy (15 of 48 warps/SM) with DRAM at 37 % and SM at 36 % — latency-bound.
~128 registers/thread went on interface machinery (8 weight variables, the
8-point gather, the ghost blend, the owned-face guard) that same-level copies
never use. The entry list was already ordered copies-first
(`init_block_exchange` round 1 emits exactly the `OP_COPY` entries), so the
split needed no new bookkeeping.

**Why the fast path is bit-exact by construction**, not by luck: for `OP_COPY`,
`entry_gather_map` returns `ga = 1, gs = 0, gc = 1` in every dim (one source at
a constant per-dim offset, weight 1.0), `entry_blend` returns 1.0 (`lWpDst = 0`,
no blend) and `interface_normal_dim` returns 0 (the write guard always passes).
The value is copied, never recomputed.

**Kernel** (A6000, 1.18 M halo points). Three formulations were measured:

| formulation | velocity copy | occupancy | DRAM |
|---|---|---|---|
| original single kernel | 350 us | 31.7 % | 37.0 % |
| split, `collapse(2)` over (var, point) | 317 us | 43.4 % | 44.5 % |
| split, one launch per variable | 3 x 97 us | 57.9 % | 48.1 % |
| **split, point-per-thread, vars inside the thread** | **253 us** | 51.7 % | **53.5 %** |

1.39x on the velocity copy, 1.19x on the scalar copy.

**End to end** (A6000, 400 cold steps, 138.4 M cells, SAME DAY A/B against a
`bcc69eb` binary — see the drift warning below):

| config | reference | candidate | gain |
|---|---|---|---|
| `base_jacobi` (nb unset) | 1.1836 | 1.1764 | 0.6 % |
| `nb16_jacobi` | 1.6923 | **1.6118** | **4.8 %** |
| `refined_yp100` | 0.7809 | **0.7495** | **4.0 %** |
| **block tax** | **1.4297** | **1.3701** | 14 % of the tax removed |

The 0.6 % at `nb` unset is the control: no blocks, no halo copies, no gain.
Runtime lines are IDENTICAL between the two binaries after 400 steps on 138 M
cells (L2_div 1.844810E-05, Linf 1.00341600E+00) — bit-exactness far past the
20-step gate cases.

**HOW MUCH IS LEFT HERE: not much.** The payload now moves at 224 GB/s with
1.84x DRAM amplification, and that amplification is STRUCTURAL: 29.5 % of halo
points lie on x-normal faces, whose halo plane is strided by `nb+2` in the
fastest-varying index and can never fill a 32 B sector (4x waste, and
0.295*4 + 0.705 = 1.885 predicts the measured 1.84). **~51 % of peak is the
floor for this data layout**, so the remaining headroom in this kernel is
occupancy-only, perhaps 1.3-1.5x — not the ~4x that "21 % of peak" suggested
before profiling. Bigger wins now have to come from moving less data
(Phases 1, 3, 5), not from moving it faster.

NOT done here, and worth knowing: `pack_entries`/`unpack_entries` (the off-rank
MPI path) still carry the same register-heavy general form. They are dead at
1 rank, which is why they did not show up, but a multi-GPU run would want the
same treatment.

**MEASUREMENT WARNING, the second time this bit.** The candidate first appeared
5.3 % faster at `nb` unset, where the copy kernel can explain ~1 %. The recorded
baseline was three weeks old and the machine had drifted ~5 % faster in the
meantime. Only a SAME-DAY run of the pre-change binary separated the two.
`run_overhead.sh` now takes `SUFFIX=<tag>` for exactly this, and `summarise.py`
pairs a run with a base carrying the same tag so an A/B never compares across
binaries by accident.

---

## STATUS 2026-08-27 — Phase 1 DONE. Block tax 1.370 -> 1.071.

`[blocks] nb` is per direction: `nb = 16` still broadcasts, `nb = 64 44 48` sets
each. Validation is per direction (even, >= 4, divides the grid) plus all-set-or-
none; the parser counts tokens, so `nb = 64 44` is an error rather than a silent
broadcast of 64.

**Payoff** (A6000, 400 cold steps, 138.4 M cells, same-day A/B, one binary):

| layout | s/step | tax | blocks | halo/interior |
|---|---|---|---|---|
| `nb` unset | 1.1771 | — | 1 | — |
| `nb = 16` | 1.6123 | 1.367 | 33792 | 42.4 % |
| **`nb = 64 44 48`** | **1.2622** | **1.071** | 1024 | 12.3 % |

**21.7 % faster than cubic `nb = 16`.** Together with the copy kernel, the
production layout went 1.6923 -> 1.2622 s/step today, **25.4 %**, with the tax
1.4297 -> 1.0722 — i.e. 83 % of the original block overhead removed, all of it
bit-exact.

**The cost model needs a caveat now.** The allocated-volume law predicted 1.1230;
measured 1.0722. It was calibrated when cost tracked halo CELLS through the old
copy kernel. With 1024 blocks instead of 33792 there are 33x fewer exchange
ENTRIES, and per-entry metadata (index decode, slot lookup, gather map) falls
faster than the cell count. **Treat `(1+2/nb_x)(1+2/nb_y)(1+2/nb_z)` as an upper
bound on the tax, not an estimate.**

**Where it reached**: `dns%block_nb` is a 3-vector through config, blocks, ibm
(three routines had `nb` scalar inside per-direction loops, plus
`block_outside_box`), io, moby_solve, moby_prepare, test_leaftable and the C
case-file readers/writer. Field files already carried `block_nb_x/y/z`, so
restart needed nothing. Case files gain a 3-element `block_nb_xyz` attribute;
the legacy scalar is still written for cubic layouts and still read, so every
committed case file keeps working and a non-cubic run against a legacy file
fails the cross-check rather than silently mismatching.

**Gates** (`validation/block_nb/run_gates.sh`, CPU AND GPU): non-cubic == cubic
== `nb` unset at `max_abs 0`; 1 rank == 4 ranks; uniform oblique flow through a
3-level patch EXACT (0.0, pn spread 0.0) in both `xz` and `xyz` at `nb = 8 4 8`;
7-case suite bit-exact vs 6d7d1e5; prepare round-trip 20/20.

### FINDING: the 2:1 interface is not nb-independent (predates Phase 1)

Two layouts refining the IDENTICAL cell range differ by ~1e-4. This is not a
regression: two CUBIC layouts (`nb = 8` vs `nb = 4`) do the same, reproducing to
9 significant digits on a pre-Phase-1 binary. So "results EXACTLY independent of
nb and rank count" holds for the SINGLE-LEVEL lattice (where the redundant
halo-layer sweep makes it exact) and for rank count always, but NOT across a 2:1
interface, where a block boundary on the interface plane is not equivalent to one
away from it. Consistent with the historical record, which always gated "channel
nb = 4 WITHOUT refinement bit-exact". Details and numbers in
`validation/block_nb/README.md`. **No gate may assume refined layouts are
bit-comparable** -- use uniform-flow preservation instead.

### Interface placement, the price of a large nb_y

A 2:1 interface can only sit on a block boundary, so `nb_y = 44` offers y-index
44 / 88 / 132 (y+ 82 / 318 / 613) where `nb_y = 16` offered every 16 cells. The
refined BL case's y+ ~99 interface moves to y+ ~82. `overheadTest/README.md` has
the placement table.

---

## STATUS 2026-08-27 — Phase 3 (vel_exchange) DONE. 3.4 % on the production layout.

`comm.f90 sync_divergence_halos` replaces the mid-iteration velocity exchange on
single-rank runs.

**The observation.** Between projection iterations the ONLY velocity halo
anything reads is the divergence stencil's. `jacobi_compute_phi` forms
`(q(ip,j,k,U)-q(i,j,k,U))*d1x + ...` over `i,j,k = 1..nb`, so it touches exactly
`q(nb+1)` in each dim, and only the component NORMAL to that face.
`jacobi_apply` writes its own interior faces and reads that same high plane for
owned interface/outlet faces; `apply_bc` works on physical ghosts the block owns.
The low halo planes, the tangential components, the edges, the corners and `p`
are NOT read until the next substage, and the last iteration's full exchange
refreshes them all. So 15 of the 18 per-step calls need three face planes of one
component instead of a 26-direction shell of three: at `nb = 64 44 48`, **8000
values against 49896**.

**Measured** (A6000, 400 cold steps, same-day A/B against `e77c75e`):

| config | reference | candidate | gain |
|---|---|---|---|
| `rect_jacobi` (nb = 64 44 48) | 1.2622 | **1.2192** | **3.4 %** |
| `nb16_jacobi` (cubic) | 1.6116 | **1.4673** | **9.0 %** |

The cubic gain is larger because its halo shell is far bigger — i.e. Phase 1
had already removed most of what this saves.

**Two restrictions, both deliberate:**
- **Same-level only**, exactly like the copy-only exchange it replaces. A block
  whose +axis face is a 2:1 interface has `dsSlot = 0` and keeps the halo it
  had; cross-level transfer happens once per substage, not per iteration.
- **Single rank only** (`c%nPeers == 0`). With peers those planes must come over
  MPI, and the message IS the saving, so it needs the entry list partitioned
  into a per-peer suffix first (the `copyOnly` prefix logic, mirrored). Multi-
  rank keeps the old path: correct, just not faster. **That fallback is also the
  gate**: `1 rank == 4 ranks` now compares the NEW path against the OLD one and
  they are bit-identical.

**Gates**: 7-case suite `max_abs 0` CPU AND GPU vs `e77c75e`; `validation/block_nb`
all pass; 1 rank == 4 ranks exact. The geometry is verified from the gather map
rather than assumed: for a same-level `OP_COPY` with `off = (+1,0,0)`,
`entry_gather_map` gives `gb1 = srcLo1 - dstLo1`, so `b1 = (nb+1) + (1-(nb+1)) = 1`
-- the source really is the neighbour's interior plane, which the kernel hardcodes.

**PROFILER CAVEAT**: the sync is a direct kernel, not an exchange, so it counts
in `proj_timing/vel_exchange` but NOT in `exch_timing/local_copy`. On a
single-rank run `exch_timing`'s total is therefore a SUBSET of the exchange
buckets above it, not equal to them. Noted in `profiling.f90`; do not read the
shortfall as lost measurement.

### WHERE THE TIME ACTUALLY IS NOW — read this before choosing the next target

Per-phase at the production layout (`nb = 64 44 48`, 1.2596 s/step, before this
change):

| phase | s/step | share |
|---|---|---|
| **`proj/apply`** | **0.4072** | **32.3 %** |
| `step/momentum` | 0.3505 | 27.8 % |
| `proj/sweep` | 0.2806 | 22.3 % |
| `proj/vel_exchange` | 0.0664 | 5.3 % |
| `proj/phi_exchange` | 0.0229 | 1.8 % |
| `exch/local_copy` (all exchanges) | 0.0976 | 7.8 % |

**The exchanges are no longer the story.** `jacobi_apply` alone is 32 % --
bigger than momentum, 6x bigger than the exchange this phase just optimised --
and nothing has examined it. It is two kernel launches (pressure update, then
the velocity face corrections) over the same interior, so at minimum the `phi`
read is done twice; whether that is what limits it is UNKNOWN and should be
MEASURED, not assumed. Recipe that worked for the copy kernel:

```
mpirun -n 1 ncu --kernel-name regex:jacobi_apply --launch-skip 12 --launch-count 4 \
    --section SpeedOfLight --section MemoryWorkloadAnalysis --section Occupancy <bin> <ini>
```

NOTE: performance counters are permission-blocked on the local workstation
(`ERR_NVGPUCTRPERM`); `ncu` works on istmcetus. And the copy-kernel precedent is
the reason to measure first: the obvious guess there (bandwidth) was wrong, and
the real cause (register-limited occupancy) called for a completely different fix.

---

## STATUS 2026-08-27 — Phase 3a (R0) ATTEMPTED, NOT BIT-EXACT, REVERTED.

R0 (`docs/next_session_redblack_interface.md` §3) is: compute the Jacobi `phi`
halo redundantly instead of exchanging it. It was implemented in full
(`jacobi_compute_phi` and `cheb_combine` extended over the three low halo face
planes behind a shared `phi_cell_owned` predicate; `crossOnly` on
`exchange_scalar_halos`) and it WORKS -- but it is **not bit-exact**, and the
reason is structural, not a coding slip.

**6 of 7 suite cases fail** (only `lam30t` passes), at ~1e-13. Isolated on one
case, changing nothing but the boundary condition:

| same case | max_abs (pn) |
|---|---|
| periodic x/z | 9.5e-13 |
| **walls all round** | **0.0** |

So the redundant computation reproduces the owner's arithmetic EXACTLY at
interior block seams -- the R0 consistency argument is sound -- and diverges
ONLY at a **periodic domain seam**. There, the two blocks sharing the face build
their metrics through different branches of `init.f90 face_at`: the low side via
`node(idx+n) - length`, the high side via `node(idx)`. For the cell-centred
metric those are bitwise equal (exact negation), but for the FACE-STAGGERED
metric the `- length` enters BEFORE the cell-centre averaging in
`cell_center_at`, so `0.5*fl(a-L) - 0.5*(fl(b-L)+fl(a-L))` is not bitwise
`0.5*(a+L) - 0.5*(b+a)`. The halo cell's metric is therefore not bitwise the
owner's, and no arrangement of the KERNEL fixes that.

**Chebyshev amplifies it but does not cause it** (a wrong turn worth recording:
the first, incomplete gate list looked like a perfect `accel = chebyshev`
correlation; the same case with `accel = none` differs too). With Chebyshev the
seam error enters the recursion through `gamma*delta`: niter = 1 gives 1.8e-13,
niter = 2 gives 1.97e-01.

**Note also**: `1 rank == 4 ranks` does NOT gate this. Both rank counts share
the same seam asymmetry, so it cancels; only the reference-binary comparison
exposes it.

**WHY IT WAS NOT PURSUED FURTHER.** Making the seam metrics canonical would fix
R0, but that change is itself NOT bit-exact (it alters the interior metric
arithmetic, since `0.5*(face(g)-face(g-2))` is not bitwise
`0.5*(face(g-1)+face(g)) - 0.5*(face(g-2)+face(g-1))`), so it needs its own
validation campaign. And the prize has shrunk. Measured at the current
production layout (`nb = 64 44 48`, total 1.2596 s/step):

| phase | s/step | share |
|---|---|---|
| `proj/phi_exchange` — R0's ENTIRE target | 0.0229 | **1.8 %** |
| `proj/vel_exchange` | 0.0664 | 5.3 % |
| `exch/local_copy` (all exchanges) | 0.0976 | 7.8 % |

At `nb = 16` before today R0's target was 5.7 % of the step; the copy-kernel
split and per-direction `nb` have already taken most of it. **A non-bit-exact
metric change plus revalidation to buy at most 1.8 % is a bad trade.**

**If a future session wants R0 anyway**, the order is: (1) make the block
metrics canonical at periodic seams and revalidate as a numerics change, THEN
(2) re-apply R0, which then gates bit-exact. A cheaper alternative that keeps
bit-exactness is to exclude periodic-seam halo cells from ownership and let the
exchange fill only those -- but that needs the entry list partitioned so the
same-level exchange can be skipped for everything else, and it buys less than
1.8 %.

**Where the exchange cost actually is now**: `proj/vel_exchange` (5.3 %), which
R0 never touched, and which runs 18 times per step. That is the target worth
sizing next, not `phi`.

---

### What Phase 0 delivered

- `src/modules/profiling.f90` + `[output] profile` (`docs/configuration.md`):
  three nested profilers on the existing `chron.f90` machinery — `step_timing`
  (disjoint, covers the loop), `proj_timing` (inside the projection bucket),
  `exch_timing` (inside every exchange). Zero-cost disabled, clocks only.
- **Coverage is 1.0000** on every run — the honesty check the trap called for.
  Validated at 1 rank and at 8 CPU ranks, where all five exchange buckets
  populate (`pack`/`mpi_post`/`mpi_wait`/`unpack`/`local_copy`).
- `overheadTest/phase_table.py` produces the table above; `--markdown` for docs.
- MPI buckets are exactly zero at 1 rank by construction, which is the
  production configuration for these runs. Multi-GPU is a separate axis and
  `mpi_wait` was NOT the answer here.
- **The profiler is free**: profiled vs unprofiled totals agree to 0.1 % on
  every jacobi run (1.2403/1.2417, 1.7522/1.7533, 0.8709/0.8713) and
  `nb16_redblack_prof` tracks its unprofiled twin step for step (1.990→2.048 vs
  1.994→2.049).
- **A measurement caveat worth inheriting**: the first profiled redblack pass
  read +27 % and was *contended*, not instrumented-slow — the per-step rate
  drifted WITHIN the run (1.45 → 1.91), and on re-measurement rose then fell
  (1.459 → 1.549 → 1.508), which a fixed instrument cost cannot do. istmcetus's
  second GPU is shared with the production campaign and the interference is
  bursty. Always read `runtime.txt`'s drift, not just the final cumulative
  average, before believing a timing delta.

---

## Working setup (do this first)

This work must NOT run in the tree that is driving the boundary-layer production
campaign. Copy the repository to a fresh directory and branch there:

```bash
cp -a mobydiff.bl mobydiff.blk            # or a fresh clone + carry the untracked dirs
cd mobydiff.blk
git checkout -b blockOverhead             # off boundaryLayer (HEAD 7e1e4b3)
git add tutorials/turbulentBoundaryLayer/overheadTest && git commit -m "overheadTest baseline"
module load toolkits/nvhpc/25.9 && ./compile.sh cpu && ./compile.sh gpu
```

`tutorials/turbulentBoundaryLayer/overheadTest/` is currently UNTRACKED and is the
measurement baseline this whole session hangs on — commit it before touching code.
GPU: corax (5090) is occupied by the campaign; use **istmcetus GPU 1** (A6000, free
2026-08-06), `CUDA_VISIBLE_DEVICES=1`, launch recipe in `overheadTest/run_overhead.sh`.
Absolute times differ from corax by ~3.1x; only ratios transfer.

## The measured starting point (2026-08-06)

Boundary-layer production grid, 138.4 M cells, A6000, 400 cold-start steps
(`overheadTest/README.md` has the full table):

| | s/step | ratio |
|---|---|---|
| `nb` unset, redblack (production layout) | 1.456 | 1.000 |
| `nb = 16`, redblack | 2.052 | **1.409** |
| `nb` unset, jacobi | 1.240 | 1.000 |
| `nb = 16`, jacobi | 1.753 | **1.413** |
| 2:1 xz refined (y+ ~100), `nb = 16`, jacobi | 0.859 | 0.490 (cells 0.455) |

Three facts that set the agenda:

1. **The block lattice costs 41%; the 2:1 interface costs 3-8%.** Optimising the
   interface transfer is chasing the small term.
2. **Cost tracks the ALLOCATED volume**: 1.409 measured vs (18/16)^3 = 1.4238
   predicted. So the overhead of any block shape is predictable as
   `(1+2/nb_x)(1+2/nb_y)(1+2/nb_z)` — no need to re-measure per case.
3. **It is not redundant arithmetic.** redblack (which sweeps the halo layer
   redundantly) and jacobi (which does not) show the SAME overhead, 1.409 vs
   1.413. The tax is bytes — footprint and/or exchange traffic — not flops.
   This already rules out "do less work in halos" as a strategy; Phase 0 decides
   between the two remaining candidates.

Also measured and load-bearing for the gates: **the block decomposition is
result-invariant.** `nb = 16` and `nb` unset produced identical runtime lines
(jacobi L2_div 1.73814896E-05 / Linf 1.00322464E+00; redblack 6.37658051E-06 /
1.00344109E+00). This is the designed Phase-1 property ("results EXACTLY
independent of nb and rank count") and it gives every phase below a free,
razor-sharp gate: **any change to blocking or storage must reproduce the same
fields bit-for-bit.**

---

## Phase 0 — timing instrumentation (DO THIS FIRST)

Nothing below should be started until the 41% is attributed. The whole ranking
changes depending on the answer.

**Do not build a profiler — instantiate the existing one.** `chron.f90` already
has `init_profiler` / `reset_profiler` / `profiler_add` / `write_profiler`, used
by `turbulence.f90` (`turb_timing`, `TURB_PROF_*`). Copy that pattern for the
step. (The old `MOBY_PHASETIME` env hook was deleted in the 2026-06-30 cleanup;
do not resurrect an env-var hook — make it a config key or always-on-and-cheap.)

Categories, chosen so the answer is unambiguous:

- `momentum` (fused predictor kernel), `bodyforce`, `turbulence` (already split)
- projection, split: `sweep` (compute-phi / redblack sweep), `apply`,
  `phi_exchange`, `vel_exchange`, `apply_bc`
- exchange internals, split: `pack`, `mpi_wait`, `unpack`, `local_copy`
  — this split is the crux: it separates device-local traffic from MPI
- `io_stats`

Traps:

- **Device synchronization.** The OpenMP target regions have no `nowait`, so the
  host blocks on each — timing is valid — but the MPI path uses
  `use_device_addr` with `Isend`/`Irecv`. Verify the instrumentation by checking
  that the sum of phases matches `write_chron`'s loop time to within a few
  percent. If it does not, the timers are lying and everything downstream is
  worthless.
- Must be **zero-cost when disabled and bit-exact always** (it only reads clocks).

**The decisive experiment.** Run the SAME case with the profiler on at `nb`
unset and `nb = 16`, and diff the per-phase times. Reuse the existing configs —
`overheadTest/singleLevel/{base,nb16}_{redblack,jacobi}.ini`, cold start, 400
steps, ~10-15 min each. Then repeat on `overheadTest/multiLevel_xz/refined_yp100_jacobi.ini`
to see what the interface's +7.7% is made of.

Decision rule, written down in advance:

| where the 41% lands | conclusion | prioritise |
|---|---|---|
| exchange `pack`/`unpack`/`local_copy` | traffic-bound | Phases 3, 5 |
| spread proportionally over volume kernels | footprint/bandwidth-bound | Phases 1, 4 (exchange tuning is futile) |
| `mpi_wait` | irrelevant at 1 rank — re-run at 4 ranks | overlap (see `docs/nonblocking_overlap_strategy.md`, needs rewriting) |

Secondary check: repeat at `nb = 32` on a grid that allows it (e.g. min_channel)
and confirm the per-phase overhead follows `(nb+2)^3/nb^3` phase by phase.

**Deliverable**: a per-phase table x {nb unset, nb=16, refined}, and a written
decision on the order of the phases below. Commit both.

---

## Phase 1 — per-direction `nb`

Highest certain gain per unit of work. `blk%nb(1:3)` is ALREADY a 3-vector
everywhere in `blocks.f90` and the kernels; only `dns%block_nb` (a scalar) and
the config validation force cubic.

- Config: accept `[blocks] nb = 64 44 48` alongside the current scalar form.
  Validation per direction: `>= 4`, even, divides the level-0 global size.
- Payoff: the boundary-layer case goes 1.42 -> 1.12, i.e. the y+ ~100
  refinement saving 31% -> ~45%. In `refine_dims = xz` the y-blocking is PURE
  overhead (it buys no refinement) and at nb=16 alone accounts for 12.5 of the
  41 points — so `nb_y` should be as large as interface placement allows.
- Trade-off to state in the config docs: larger `nb_y` = fewer y-tiles = fewer
  legal interface heights (the interface can only sit on a tile boundary).

**Audit list** (things that may assume cubic — check each, do not assume):
leaf-table builder and `level_cells`; Morton id encode/decode and `zorder_*`;
the `refine_dims` per-level scaling `2**(l*mask(d))`; exchange entry generation
(already per-dim, likely fine); **io: the block-table datasets are
`(nBlocksGlobal, nb^3)`** — this becomes `nb_x*nb_y*nb_z` and the file's `nb`
attribute changes shape, so the restart cross-check and `moby_prepare` must move
together; `tools/compare_fields.py` and the block-table readers.

**Gates**: cubic `nb` bit-exact (nofma, `max_abs 0`, CPU AND GPU) on the standard
7-case suite; **non-cubic `nb` bit-exact vs cubic `nb` on the same case** (this
is the strong one — nb-independence is a designed property, see above);
1 rank == 4 ranks EXACT; CPU == GPU; prepared case files round-trip.

---

## Phase 3 — fewer exchanges per iteration

Only worth doing if Phase 0 says traffic. Two independent increments:

- **3a — Jacobi's `phi` halo by redundant computation.** Full design already
  written in `docs/next_session_redblack_interface.md` §3 (Increment R0): the
  apply reads only the three LOW halo face planes, which are exactly the cells
  the red-black redundant sweep already computes; extend `jacobi_compute_phi`
  over that layer and delete the same-level part of the scalar exchange. Halves
  the projection's exchanges on single-level grids at ~+9% arithmetic (nb=32) —
  and arithmetic is demonstrably not the constraint. Chebyshev needs
  `cheb_combine` extended over the same layer.
- **3b — cross-level-suffix-only scalar exchange.** A `phase`/`crossOnly`
  argument on `exchange_scalar_halos`, mirroring `copy_local_entries(phase=2)`
  and the `copyOnly` MPI prefix logic (`peerSendCopyOff`/`peerRecvCopyOff`).
  With 3a in place the remaining `phi` exchange exists ONLY to cross a level
  jump, which is the design rule this and the red-black work share.

**Gates**: bit-exact on the 7-case suite + a refined case, CPU AND GPU; report
the s/step delta with the Phase-0 profiler, not with wall-clock guesses.

---

## Phase 4 — stop allocating halos nobody reads

Cheap, small, bit-exact by construction. `oldrhs` is already `1:nb` — the
precedent exists (`blocks.f90:310`).

Halo-carrying today: `blk%q` (NVAR), `blk%qs` (NVEL), `ibm%coef` (3), `ibm%mu`
(3), `turb%nut`, `phi` — 14 components of `(nb+2)^3`.

- `mu`, `q`, `nut`, `phi` are read by stencils in halos — keep.
- `coef` needs halos because `update_ibm_mu` is pointwise over the full array
  including halos (`mu = 1/(1+dt*coef)`) — keep unless `mu`'s halo is exchanged
  instead.
- **`qs` is the candidate**: the predictor writes it on interior faces and the
  projection reads `blk%q`. Confirm nothing reads `qs` in a halo (watch the
  interface momentum path, which predicts the shared face on BOTH sides) — if
  clean, that is 3 of 14 components freed, ~21% off the halo excess.

**Gate**: bit-exact by construction; if it is not bit-exact, something DID read
that storage and the audit was wrong.

---

## Phase 5 — stop paying a per-block halo where there is no refinement

The structural item. Motivation: in the boundary-layer test the refinement
topology was **two slabs**, represented by **15360 blocks**. The equal-size
lattice is bookkeeping for a geometry with essentially no complexity.

**Start with the cheap subset, not the general case:**

- **5a — per-LEVEL block size.** Coarse levels use bigger blocks (their regions
  are large and few). Within a level everything stays uniform and rectangular,
  so kernels keep `collapse(4)`; the change is one array (or offset range) per
  level plus per-level launches. Far less disruptive than general fusion and
  captures much of the same benefit. Composes with Phase 1.
- **5b — fuse same-level, all-interior blocks into storage tiles.** Keep the
  Morton lattice, the leaf table, the file format and the exchange-entry
  semantics; let a rectangular run of same-level same-rank neighbours share one
  array with no internal halo, so halos survive only at rank boundaries and
  level jumps. For the BL y+ ~100 case this would take the tax from 1.41 toward
  1.0 and the saving from 31% to ~50%.
  - **The crux**: kernels rely on a uniform `nb` for `collapse(4)`. Variable
    tile extents mean per-tile launches or a flattened cell list — this touches
    the single most load-bearing performance property of the code. Prototype on
    ONE kernel and one direction (fuse along x only, tiles `(k*nb+2)` x
    `(nb+2)` x `(nb+2)`) and measure before committing.
  - Do not start 5b unless Phase 0 attributes the tax to footprint/traffic AND
    5a has been measured.

**Gate**: still bit-exact — by the nb-independence property, a different tiling
of the same leaves must give identical fields. Note that red-black's redundant
halo sweep pattern changes with the tiling, so this gate genuinely exercises the
nb-independence argument rather than trivially passing.

---

## Phase 6 — mixed precision in the projection

Last, and deliberately so: it is the only phase that breaks the `max_abs 0`
discipline every part of this project has relied on.

- Scope: the smoother in fp32 with fp64 residual/correction accumulation
  (standard), and possibly an fp32 exchange payload for the smoother only.
  The momentum path and the stored state stay fp64.
- Payoff: up to ~2x on the projection if it is purely bandwidth-bound — which
  Phase 0 will have established one way or the other.
- **New risk to gate explicitly**: the interface transfer operators are exactly
  conservative in fp64 (`ifaceRow` restrict = mean, prolong = injection,
  mean-preserving). In fp32 the mean-equality holds only to fp32 round-off, so
  the global mass residual (currently ~1e-20 with a patch) must be re-gated.
- **Replacement gates** (since bit-exactness is gone): converged residual and
  final fields agree with the fp64 solver to a stated tolerance on the standard
  suite; developed-channel statistics unchanged (Re_tau 180 mean U, stresses);
  interface band ratios unchanged; global mass residual bounded; long-run
  stability on a Dirichlet-p-outlet case (the 2-dx mode of the boundary-layer
  finding is exactly the kind of thing reduced precision could excite).

---

## Cross-cutting

- **Standard verification suite** (unchanged): min_channel (blocks + 2:1 +
  chebyshev), les_ibm channel +/- refine_body, Beltrami y-slab, turb180,
  wf180_y30, lam30t — `-Mnofma` / `-gpu=nofma`, `tools/compare_fields.py`,
  `max_abs 0`, CPU AND GPU. Plus 1 rank == 4 ranks and CPU == GPU on a refined
  case.
- **Timing suite**: the `overheadTest/` configs. Cold start, 400 steps, no
  output. Report BOTH the `chron` loop rate and the differenced runtime-line
  rate (`overheadTest/summarise.py` does this).
- **Do not** start Phase 5b or 6 in the same session as anything else. Each of
  0, 1, 3, 4 is a self-contained commit with its own gate; the session can stop
  cleanly after any of them.
- Deliberately NOT in this handout: **multigrid on the existing 2:1 level
  hierarchy**. It is the largest single lever on total cost (it attacks `niter`,
  which multiplies both the halo tax and the interface cost) and the
  restriction/prolongation operators already exist, validated and conservative —
  but it is a project of its own, not an optimisation pass.

## References

- `tutorials/turbulentBoundaryLayer/overheadTest/README.md` — the baseline
  measurements, the placement table, and the cost model.
- memory `block-overhead-measured` — the `(nb+2)^3` law and its consequences.
- `docs/next_session_redblack_interface.md` — Phase 3a is its Increment R0;
  also the reason the refined timing runs had to use `solver = jacobi`.
- `docs/next_session_profiling.md` — the older profiling plan (stale numbers,
  but the same phase-timer need; supersede or merge it in Phase 0).
- `docs/nonblocking_overlap_strategy.md` — predates the Chebyshev solver;
  rewrite before using.
- `docs/block_refinement_strategy.md` §5 (exchange entries), §10 (choosing nb).
