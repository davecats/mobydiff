# Next session — the multi-rank exchange

Handout. The single-rank kernel story is written up in
`docs/next_session_jacobi_apply.md`; this is the parallel one, which has a
different answer.

## STATUS 2026-08-28 — diagnostics run, and they reordered the plan

The ranked plan below was written first and then measured against
(`overheadTest/results_exchange_diag_2026-08-28.md`). **Read this header
before the plan: item 1 was demoted by its own diagnostic.**

- **GPU has no arrival imbalance** (timed barrier: 0.065 % of the step). Its
  `mpi_wait` is pure transfer — 300 MB/step at **12.5 GB/s** aggregate over 39
  rounds, with per-round latency accounting for 2 %. **CPU has real skew**
  (barrier 4.2 % of the step) and moves only **0.46 GB/s**, 25× below intra-node
  shared memory: nothing progresses until the `Waitall`.
- **The byte budget is dominated by the phi exchange** (45 % GPU / 33 % CPU,
  18 rounds/step) and the full velocity shells (53 % / 39 %).
- **The multi-rank fallback is NOT a full exchange.** It is
  `exchange_halos(interp=.false.)`, the copy-only prefix — and that prefix holds
  only 1.9 % of entries on the 2-rank GPU decomposition (18000 of 941700). The
  15 mid-iteration refreshes are **2.2 % of GPU bytes**. My "5–6 % of the step"
  estimate for item 1 assumed a full exchange and was wrong; the poor scaling of
  `proj/vel_exchange` is mostly the 3 *full* shells per step.
  It is decomposition-dependent: at 4 CPU ranks the prefix is 33.6 % of entries
  and 27.9 % of bytes, so the same change matters there.

**Revised order: overlap (item 2) first on GPU** — the "if it is imbalance,
overlap will not help" caveat is resolved in overlap's favour. Then **async MPI
progress on CPU** (try `UCX_*` / progress-thread environment before touching
code — one run, no diff). Then the phi exchange. Then item 1, sized against the
decomposition it will actually run on.

The measured message sizes are now printed at init under `[output] profile`
(`comm.f90 report_exchange_sizes`, bit-exact: 7-case suite `max_abs 0` CPU and
GPU). The barrier was temporary and is reverted; the method note at the end of
the results file says how to re-add it.

## What is established

Measured 2026-08-28 on the refined config
(`overheadTest/multiLevel_xz/refined_yp82_rect_jacobi.ini`, 448 leaves, 60.6 M
cells), 2 GPUs and 4 CPU ranks, each against a same-host 1-rank baseline
(`overheadTest/results_multirank_2026-08-28.md`):

- **85 % parallel efficiency** on both platforms, from absolute rates 115×
  apart.
- **Volume kernels scale essentially perfectly** — momentum 1.02, sweep
  0.99–1.02, apply 0.98–1.00.
- **The projection exchange is 71 % (GPU) / 75 % (CPU) of all time lost**, and
  gets *absolutely* more expensive as ranks are added.
- The whole exchange goes from ~5 % of the step to **17 % on both platforms**;
  `mpi_wait` alone is 7–8 %.

## What is NOT the problem, so do not "fix" it

**Message aggregation is already done.** `start_halo_exchange`
(comm.f90:1106–1118) posts exactly one `MPI_Isend` per peer — not per entry,
not per direction — of length `nSendPts*nv`, i.e. every entry toward that peer
for all `nv` active variables in one contiguous buffer. A copy-only exchange is
a *prefix* of that same buffer, so it needs no message of its own.

The profile confirms it: **`mpi_post` is 0.03 % of the step on GPU and 0.0025 %
on CPU.** Message count and posting overhead are not where the time is. Any
proposal that amounts to "batch the messages" is already implemented; check
this section before re-deriving it.

## MEASURE FIRST — `mpi_wait` is three different problems

7–8 % of the step sitting in `MPI_Waitall` can be transfer volume, per-round
latency × 36 exchange rounds per step, or load imbalance absorbed at the wait.
They need different fixes and nothing measured so far distinguishes them. Two
cheap diagnostics, in this order:

1. **Per-peer buffer sizes** (bytes actually on the wire per round, and per
   step). With the round count this gives an effective GB/s during `mpi_wait`:
   near the interconnect's peak ⇒ bandwidth-bound, far below ⇒ latency or
   imbalance.
2. **A timed `MPI_Barrier` at the START of the exchange** (before the pack, so
   it is not contaminated by progress on this round's own transfers). Its time
   is arrival imbalance carried in from the preceding compute; what remains in
   `mpi_wait` is transfer. If the barrier is large the ranks are unbalanced and
   overlap will not help — which would reorder everything below.

## PLAN A — overlap the exchange with compute (GPU)

Warranted by the diagnostics: on GPU the wait is pure transfer with zero
imbalance, so hiding ~24 ms/step of transfer behind ~280 ms/step of compute is
the right shape of fix. Ceiling: the 7.2 % of the step in `mpi_wait`, plus part
of `pack`/`unpack`.

### A0 — RUN 2026-09-02. IT FAILS. Plan A is CLOSED.

`overheadTest/results_a0probe_2026-09-02.md`. A compute-bound target kernel was
inserted between the `Isend`/`Irecv` posts and the `MPI_Waitall` (in the
pre-key-change binary, so the 30 MB rounds are resolvable):

| | ms/round |
|---|---|
| compute-bound work inserted between post and wait | **9.95** |
| `mpi_wait` without the probe | 0.559 |
| `mpi_wait` with the probe | **0.611** |

**Nothing was hidden.** The transfer happens inside `MPI_Waitall` and nowhere
else, so an interior/shell split buys exactly zero. It also retrospectively
explains why `local_copy` — 0.55 ms/round of device work already sitting
between post and wait — never hid the 0.56 ms transfer: not DMA contention, no
progress at all.

No zero-code fallback on this stack: Open MPI 4.1.9a1 (hpcx) exposes only
`btl_tcp_progress_thread` (TCP, irrelevant intra-node CUDA/UCX) and
`orte_progress_thread_debug`. The remaining route — `nowait` on every kernel to
be overlapped plus host-side `MPI_Test` polling — is invasive and fragile under
nvfortran offload.

**And the prize was already taken by P1.** This plan's ceiling was written as
"the 7.2 % of the step now in `mpi_wait`". That was the pre-key-change number;
after the xz key change `mpi_wait` is **0.96 %** of the 2-rank step, while the
device-local part of the exchange (`local_copy` 7.88 % + `pack` 1.58 % +
`unpack` 1.40 %) is 10.9 % — 11× larger. **The exchange is now a device-local
copy problem, not a communication one.**

Three independent reasons not to implement A1–A3, any one sufficient: the
mechanism does not work, the workaround is invasive with no env fallback, and
the ceiling is under 1 %. Scope caveat: 2 GPU ranks is the only multi-GPU
configuration available, and peer traffic grows with rank count — repeat A0 on
a genuinely many-rank GPU machine before carrying this conclusion there.

The original A1–A3 sketch is kept below for that eventuality.

### A1 — hide the phi exchange behind `jacobi_apply`

The cleanest first target because the dependency is already favourable:

- `jacobi_apply` kernel 1 (`p += phi*idt`, 160 of its 745 µs) reads **no phi
  halo at all** — it can run entirely during the exchange, with no splitting.
- kernel 2 reads only `phi(0,j,k)`, `phi(i,0,k)`, `phi(i,j,0)` (the three low
  halo *face planes*, never edges or corners) and the high plane at outlet and
  2:1 faces. So cells with `i,j,k >= 2` are independent of the exchange:
  ~95 % of the block at `nb = 64 44 48`.

Sequence: `pack` → post → **kernel 1 + kernel 2 interior** → `Waitall` →
`unpack` → kernel 2 boundary shell. `start_halo_exchange` /
`finish_halo_exchange` are already separate subroutines, so the split exists at
the call level.

### A2 — hide the velocity exchange behind the next `jacobi_compute_phi`

`jacobi_compute_phi` reads velocity at `(i,j,k)` and the `+1` neighbour per
dim, so cells with `i <= nb-1` (etc.) do not need the high halo — again ~95 %
of the block. Same shape as A1, one iteration later in the loop.

### A3 — the same for the post-predictor shell and `momentum`

Only if A1/A2 pay off. `momentum` is the largest single consumer of velocity
halos and the most invasive to split.

### Costs and gates

Splitting a kernel into interior + shell is a **scheduling change**: identical
arithmetic per cell, so it must be bit-exact (`max_abs 0`, CPU AND GPU, 7-case
suite + `validation/block_nb` + 1 rank == 4 ranks). Expect to *lose* some kernel
efficiency to two launches and a worse shell shape — the shell is a thin slab
with poor occupancy — so each increment needs its own same-day A/B, and A1 must
be shown to win before A2 is written.

## PLAN B — the phi exchange (45 % of GPU bytes, 18 rounds/step)

### B0 — what can honestly be sent less

`jacobi_apply` reads only the three low halo face planes plus the high plane at
interface/outlet faces (§3 of `docs/next_session_redblack_interface.md`). So in
principle the phi exchange could drop its edge and corner entries and its
same-level high planes.

**But the measurement says that saves almost nothing here**: 98.1 % of the phi
bytes on the 2-rank GPU decomposition are *cross-level* entries, which need both
planes and are exactly what a 2:1 interface must transfer. The reducible part is
the 1.9 % same-level prefix. **Do not spend a session on B0 for the refined
case.** It is worth more on single-level multi-rank runs, where the peer traffic
is 100 % same-level — but there the whole exchange is already small (73 600 pts
against the refined case's 941 700).

**R0 (compute the phi halo redundantly instead of exchanging it) stays closed.**
It was attempted and reverted in `e77c75e`: block metrics are not bitwise equal
across a periodic seam, so the redundantly-computed halo `phi` differs from the
owner's in the last ulp and the change is not bit-exact. That is a numerics
change with its own justification burden, not a refactor.

### B1 — the real lever is upstream: the interface sits ON the rank boundary

See `overheadTest/results_exchange_diag_2026-08-28.md`. The cross-level peer
count is **identical at 2 and 4 ranks** (923 700 / 923 692) because
`refine_dims = xz` puts the y tile in the high Morton bits, so every fine block
precedes every coarse block and the level change is one contiguous cut that a
linear split cannot avoid. The refined case therefore sends 12.8× the peer
points of the single-level case at 2 ranks despite having 44 % of the cells.

This inflates the phi exchange, both velocity shells and their pack/unpack at
once, so it is upstream of A and B0 alike. It is also the most structural: it
costs either the closed-form `zorder_owner/start/count` ownership lookup (change
the split) or the canonical leaf-table id order that `moby_prepare`, restart
files and `make_channel_restart` all mirror (change the curve). Honest ceiling:
most of `mpi_wait` becomes `local_copy` at ~1/6 the per-point cost, ≈ 6 % of the
2-rank step. Measure what a balanced interface split costs in load imbalance
before designing either.

## PLAN C — red-black instead of Jacobi: the accounting

The motivation offered was that red-black converges in `niter = 3` where Jacobi
needs 6, so it should halve the communication. **The round count does not
halve** — but the compute might, and that is the better reason.

The design is already written: `docs/next_session_redblack_interface.md` §4
(R1), ~150–250 lines, no new numerics. Its premise is sound and worth repeating
because it is easy to lose: the *operator* is already SPD and consistent at the
interface (composite `face_grad`, low-side-owns-face, symmetric relaxation,
mean-preserving `ifaceRow` transfer). Red-black is a different **splitting** of
that same `L`, so this is plumbing. The §6c–§6g failures in
`interface_projection_derivation.md` were failures of the old operator and must
not be re-fought.

Rounds per substage, with a 2:1 interface present:

| | per iteration | iterations | phi rounds | velocity rounds | total |
|---|---|---|---|---|---|
| Jacobi `niter = 6` | 1 phi + 1 velocity | 6 | 6 | 6 | 12 |
| red-black `niter = 3` | 2 colours × (1 phi + 1 velocity) | 3 | 6 | 6 | **12** |

**Identical.** Red-black needs a communication *per colour*, and with an
interface each colour needs both a cross-level phi exchange and a velocity
exchange — so halving the iterations exactly cancels against doubling the
exchanges per iteration. On a **single-level** grid red-black needs no phi at
all and the count halves (6 rounds against 12) — but single-level is where
red-black already works today.

Bytes do not save either on the current decomposition: red-black's phi exchange
is cross-level-only, which is 98.1 % of the full phi exchange here. It *would*
save if B1 were fixed.

**Where red-black does win is compute**, and that is the case worth making:
`niter = 3` is 6 half-sweeps against Jacobi's 6 full `compute_phi` + 6 `apply`
passes — roughly half the projection arithmetic. `sweep` + `apply` are 46 % of
the 2-rank step, so the prize is of order **20 % of the step**, several times
anything in Plan A. It is also the only item here that helps the single-rank
case equally.

**Convergence is settled (user measurement, 2026-08-28): Jacobi needs
`niter >= 12` to keep the pressure zero-mode away; red-black stays clean below
`niter = 6` at similar residuals.** Two consequences: every baseline in this
repository runs `niter = 6` Jacobi and is therefore **under-iterated** — the
correct-Jacobi refined 2-rank step is ~0.52 s/step, not the recorded 0.324, and
red-black must be compared against that. And the round count obeys "red-black at
`N` = Jacobi at `2N`", so at the equal-quality settings (12 vs 6) the rounds are
*exactly equal* and a round saving needs red-black below 6.

The remaining unknown is cost, not convergence: at equal `niter` red-black
measured **17 % slower** than Jacobi here (1.4499 vs 1.2422 s/step, single-level,
2026-08-07) because a colour sweep wastes bandwidth on the unused colour. The
prize is `(12 / N_redblack)` against a ~1.17 penalty. Measuring it needs no new
code — `run_overhead.sh` now takes `NITER=` and both smoothers run on the
single-level path.

Still unknown: the `omega` stability limit with an interface present; §5 of the
red-black handout has the escalation ladder.

## Recommended order — REVISED TWICE on 2026-08-28; read both revisions

**First revision** sized Plan C single-level
(`overheadTest/results_smoother_2026-08-28.md`): −26 % at `niter = 6`, −51 % at
3 against the honest Jacobi setting, dwarfing everything else.

**Second revision, after R1 was implemented and measured on a REFINED grid**
(`overheadTest/results_smoother_refined_2026-08-28.md`): **the win mostly does
not survive refinement — 6.2 % at equal residual** (21 % under the zero-mode
criterion, which leaves red-black at 1.23× Jacobi's residual). A control at the
same base resolution with refinement removed shows it is the interface, not the
grid: red-black's per-iteration residual advantage falls 1.91× → 1.49× and its
cost premium rises 1.095× → 1.264×, exactly the additive cross-level coupling
§5 of the red-black handout predicted. The `omega` sweep is done: 1.7 beats the
1.5 that was assumed, but diverges at `niter = 3`.

So on refined grids — what the 2:1 machinery exists for — Plan C no longer
dominates Plan A or B, and the ordering below is a genuine three-way choice
again. The item that would restore it is §5's ladder, in particular
**level-ordered smoothing** (coarse, patch, fine — making cross-level coupling
multiplicative for one extra phi exchange per colour).

1. **C-R1** — red-black + 2:1, `docs/next_session_redblack_interface.md` §4.
   The design exists, convergence and cost are both measured, the operator is
   already SPD at the interface. Start with the uniform-flow-through-a-patch
   gate (§6): if that is not EXACT, stop and fix the transfer.
2. **A0**, the progress probe — half a day, and it decides whether Plan A exists
   at all. Worth doing early anyway because it is cheap and it also tells you
   whether red-black's per-colour exchanges will overlap.
3. **A1** if A0 says transfers progress during a kernel.
4. **B1** (partitioning) once someone is willing to touch the Morton order or
   the ownership lookup; it is the largest *exchange* lever but the most
   structural, and it becomes more attractive after C-R1 because red-black's
   per-colour cross-level phi exchange is exactly the traffic B1 would localise.

Note for whoever runs the A/B on any of these: **re-baseline at `niter = 12`
first.** Every recorded number in `overheadTest` is `niter = 6` Jacobi, which
does not keep the pressure zero-mode away, so all of them understate the
production cost.

B0 and R0 are closed for the refined case; do not reopen without a decomposition
where same-level peer traffic dominates.

## The original ranked plan (superseded above, kept for the reasoning)

**1. `sync_divergence_halos` with peers.** ~5–6 % of a multi-rank step at 2–4
ranks, growing with rank count, on the worst-scaling phase (efficiency
0.21–0.25), with the mechanism already named in `pressure_solver.f90`: the
minimal mid-iteration refresh is single-rank only, so with peers 15 of the 18
per-step velocity refreshes fall back to a full three-component
`exchange_halos(interp=.false.)`. The work is to partition the entry list so
the divergence-only planes (one plane per dim, normal component) form their own
message set. Universal — any multi-rank run, refined or not, IBM or not.
Bounded scope. Start here unless the diagnostics say otherwise.

**2. Overlap MPI with compute.** The structural one, and better placed than it
looks: `start_halo_exchange` / `finish_halo_exchange` are already separate
subroutines, so the split was anticipated. Today the only work between the
`Isend` and the `Waitall` is the local device copy (6.6 % of the step); the
volume kernels (69 %) all run after it. Splitting `jacobi_compute_phi` into an
interior pass and a boundary shell would expose ~88 % of its cells to overlap
(interior fraction at `nb = 64 44 48`), enough to hide the 7.2 % `mpi_wait` in
principle. Costs: two launches and worse kernel shapes, and
`docs/nonblocking_overlap_strategy.md` predates the Chebyshev-Jacobi solver and
needs rewriting first. **Conditional on diagnostic 2** — if the wait is
imbalance rather than transfer, overlap buys nothing.

**3. `jacobi_compute_phi` kernel efficiency.** 18–22 % of the step, and the only
item here that also helps single-rank runs. Compute/latency-bound, not
traffic-bound: 82 % SM throughput at 46 % DRAM, 128 registers, 31.5 %
occupancy. Instruction count and register pressure — the six `face_grad_denom`
calls and the division. See `docs/next_session_jacobi_apply.md`.

**4. The `bodyforce` trip mask — NARROW, needs a decision before starting.**
`fill_trip_kernel` rewrites the entire `bf%f` array (1.45 GB per substage) and
`add_bodyforce_correction` reads it all back, for a force that is identically
zero outside `|x − x₀| ≳ 28` of a 750-long domain. A static per-block "the trip
envelope touches this block" mask would skip ~95 % of blocks in both kernels:
0.0376 → ~0.002 s/step, **7–8 % of the step**, and it removes the load
imbalance behind the 0.68 (GPU) / 0.40 (CPU) scaling of that phase. Bit-exact
in practice (`x + 0.0 == x`; only a `-0.0` sign could differ, which compares
equal at tolerance 0).

But it helps ONLY the boundary-layer trip case — the 2:1-with-IBM cases do not
enable `[force]` at all. That is the same narrowness that got the body-free mu
path reverted on 2026-08-28 (see that entry in
`docs/next_session_jacobi_apply.md`). It differs in being a per-block skip
inside an already case-specific feature rather than a second branch through hot
kernels, but it is the user's call, not a default.

## Gates for anything here

The standard discipline. Every item above is a scheduling change and must be
**bit-exact, `max_abs 0`, CPU AND GPU** on the 7-case suite, plus
`validation/block_nb/run_gates.sh` and 1 rank == 4 ranks. Item 1 additionally
has to hold the multi-rank invariance already demonstrated for this config:
1 vs 2 GPU ranks give identical runtime lines after 400 steps, and CPU 1 vs 4
agree with each other and with the GPU.

Timing: same-host, same-session baselines only (memory:
`timing-runs-need-drift-check`). Multi-GPU runs need `overheadTest/gpu_rank.sh`
— the solver never calls `omp_set_default_device`, so without per-rank
`CUDA_VISIBLE_DEVICES` every rank offloads to the same card and it reads as
catastrophic scaling. CPU runs of this case need `NSTEPS=` (20 s/step at 4
ranks, 64 at 1).
