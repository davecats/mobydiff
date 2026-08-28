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

## The ranked plan

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
