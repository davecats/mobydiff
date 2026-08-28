# What `mpi_wait` actually is, istmcetus, 2026-08-28

Two diagnostics on the refined config
(`multiLevel_xz/refined_yp82_rect_jacobi.ini`), run before designing anything,
because 7–8 % of the step sitting in `MPI_Waitall` can be transfer volume,
per-round latency, or arrival imbalance, and those need different fixes.

1. **Message sizes** — now a permanent init line under `[output] profile`
   (`comm.f90 report_exchange_sizes`).
2. **A timed `MPI_Barrier` before the pack** — temporary, reverted; its time is
   arrival imbalance carried in from the preceding compute, uncontaminated by
   progress on this round's own transfers.

The profiler gives the round count directly: **39 exchange rounds per step**
(15600 `pack` calls / 400 steps), = 3 substages × (1 post-predictor velocity
shell + 6 phi + 6 projection velocity).

## Diagnostic 2: imbalance

| | barrier (arrival skew) | `mpi_wait` after | step cost of the barrier |
|---|---|---|---|
| GPU, 2 ranks | 0.00021 s/step = **0.065 %** | 0.0241 | +0.2 % |
| CPU, 4 ranks | 0.804 s/step = **4.2 %** | 1.313 | +3.1 % |

**GPU: there is no imbalance.** The wait is transfer, full stop.
**CPU: there is real skew** (4.2 % of the step), and it is a *separate* problem
from the transfer — putting the barrier in only dropped `mpi_wait` by 13 %
(1.513 → 1.313) while adding 0.804, i.e. it exposed skew rather than removing
cost.

## Diagnostic 1: the byte budget on the wire

```
GPU 2 ranks: peers/rank 1, send pts/rank min 193700 max 748000, total 941700
             copy-only prefix 18000 pts = 1.9 % of entries
CPU 4 ranks: peers/rank 3, send pts/rank min 101300 max 963600, total 1390900
             copy-only prefix 467208 pts = 33.6 % of entries
```

| per step, both/all ranks | GPU 2 ranks | | CPU 4 ranks | |
|---|---|---|---|---|
| phi (18 rounds, nv=1) | 135.6 MB | **45.2 %** | 200.3 MB | **33.3 %** |
| projection velocity, copy-only (15, nv=3) | 6.5 MB | 2.2 % | 168.2 MB | 27.9 % |
| projection velocity, full (3, nv=4) | 90.4 MB | 30.1 % | 133.5 MB | 22.2 % |
| post-predictor velocity shell (3, nv=3) | 67.8 MB | 22.6 % | 100.1 MB | 16.6 % |
| **total** | **300.3 MB** | | **602.2 MB** | |

Effective rates against the measured `mpi_wait`:

| | rate | reading |
|---|---|---|
| GPU 2 ranks | 300 MB / 24.1 ms = **12.5 GB/s** | ~half of PCIe 4.0 x16 practical (~25 GB/s) — consistent with GPU→GPU staged rather than peer-to-peer. Per-round latency at ~15 µs × 39 is 2 % of `mpi_wait`, so this is **transfer, not latency**. |
| CPU 4 ranks | 602 MB / 1.31 s = **0.46 GB/s** | 25× below intra-node shared-memory rates. **Not bandwidth** — the transfers only progress inside the `Waitall` (no asynchronous progress), plus residual skew. |

## This overturns part of the plan written an hour earlier

The ranked plan in `docs/next_session_multirank_exchange.md` put
"`sync_divergence_halos` with peers" first, at an estimated 5–6 % of the step.
That estimate assumed the multi-rank fallback was a *full* exchange. **It is
not** — it is `exchange_halos(interp=.false.)`, the copy-only prefix, and the
prefix is nearly empty on the GPU decomposition: 18000 of 941700 points, 1.9 %.
The 15 mid-iteration refreshes put **2.2 % of the GPU's bytes** on the wire.

So on the 2-rank GPU layout the divergence sync cannot save meaningful MPI
traffic, and the poor scaling of `proj/vel_exchange` (0.01318 → 0.02603) is
mostly the **3 full shells per step** (90.4 MB, ~7.2 ms at 12.5 GB/s) rather
than the 15 cheap ones.

It is decomposition-dependent, which is the part worth remembering: at 4 CPU
ranks the copy-only prefix is 33.6 % of entries and 27.9 % of bytes, so there
the same change *would* matter. The 2-rank Morton split happens to cut almost
entirely through cross-level (2:1 interface) territory, which the same-level
copy-only prefix excludes by construction.

## Revised ranking

**1. Overlap the exchange with compute.** On GPU this is now unambiguous: the
wait is pure transfer (12.5 GB/s, no imbalance, latency 2 %), so hiding 24
ms/step of transfer behind ~280 ms/step of compute is exactly the right shape
of fix, and the earlier caveat ("if it is imbalance, overlap will not help") is
resolved in overlap's favour. `start_halo_exchange` / `finish_halo_exchange`
are already separate, and splitting `jacobi_compute_phi` into interior + shell
would expose ~88 % of its cells. Ceiling: the 7.2 % of the step now in
`mpi_wait`, plus some of `pack`/`unpack`.

**2. Asynchronous progress on CPU.** 0.46 GB/s says nothing moves until the
`Waitall`. This is a different bug from the GPU one and is cheap to test:
`MPI_Test` polling in the local-copy loop, or an MPI progress thread
(`UCX_*`/`MPICH_ASYNC_PROGRESS` environment, no code change) — try the
environment switch first, it costs one run.

**3. The phi exchange is the single biggest byte consumer** (45 % GPU / 33 %
CPU) and scales at 0.32–0.33. 18 rounds/step of the full entry list at nv=1.
Worth asking whether every Jacobi iteration needs the whole entry list or only
the rows the face corrections read — but that is a numerics question, not a
scheduling one, so it needs care and its own gate.

**4. `sync_divergence_halos` with peers** — demoted. Still right in principle
and still worth its Phase-3 *local-copy* saving, but its MPI saving is 2 % of
bytes on the GPU layout and 28 % on the CPU one. Do it after the above, and
size it against the decomposition it will actually run on.

`jacobi_compute_phi` kernel efficiency and the `bodyforce` trip mask are
unchanged from `docs/next_session_multirank_exchange.md`.

## Method note

The barrier trick is worth repeating: park it in the `mpi_post` bucket (normally
~0.03 % of the step) so no new profiler category is needed, put it *before* the
pack so it is not contaminated by this round's own transfers, and read the
result as "barrier = skew, wait = transfer". Revert it — it costs 0.2 % (GPU) to
3.1 % (CPU) of the step and it serialises what the code deliberately pipelines.
