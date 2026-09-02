# A0: does MPI progress while a target kernel runs? No.

The feasibility probe Plan A (`docs/next_session_multirank_exchange.md`) put
ahead of any overlap work, because the whole plan rests on that one assumption.
istmcetus, 2 GPU ranks, refined production case, clean machine.

## The probe

A **compute-bound** target kernel inserted between the `Isend`/`Irecv` posts and
the `MPI_Waitall`, in the pre-key-change binary (30 MB rounds, so the transfer
is resolvable). Compute-bound on purpose: `local_copy` is memory-bound, so a
memory-bound probe could hide a failure to progress behind DMA contention and
give a false negative.

| | ms/round |
|---|---|
| compute-bound work inserted between post and wait | **9.95** |
| `mpi_wait` without the probe | 0.559 |
| `mpi_wait` with the probe | **0.611** |

**Nothing was hidden.** Ten milliseconds of GPU work sat between the post and
the wait and the transfer did not advance by so much as a percent — it got 9 %
worse, presumably minor contention. The transfer happens inside `MPI_Waitall`
and nowhere else.

This also explains the earlier observation that `local_copy` (0.55 ms/round of
device work, already sitting between post and wait) never hid the 0.56 ms
transfer either. That was not bandwidth contention; there is simply no progress.

## No zero-code workaround on this stack

Open MPI **4.1.9a1** (hpcx). The only progress-thread MCA parameters are
`btl_tcp_progress_thread` (TCP only — irrelevant for an intra-node CUDA/UCX
path) and `orte_progress_thread_debug`. There is no general asynchronous
progress thread to switch on, so the "try the environment first, it costs one
run" option in the plan is not available here.

The remaining route would be `!$omp target ... nowait` on every kernel to be
overlapped, plus host-side `MPI_Test` polling while it runs — invasive, fragile
under nvfortran offload, and it would have to be threaded through the
dependency structure of the projection.

## And the prize has already been taken

Plan A ranked overlap with "ceiling: the 7.2 % of the step now in `mpi_wait`".
That was the **pre-key-change** number. After the xz key change
(`results_xzkey_2026-08-29.md`) the 2-rank exchange looks like this:

| bucket | s/step | % of step |
|---|---|---|
| local_copy | 0.024587 | 7.88 % |
| pack | 0.004938 | 1.58 % |
| unpack | 0.004376 | 1.40 % |
| **mpi_wait** | **0.003003** | **0.96 %** |
| mpi_post | 0.000097 | 0.03 % |
| total | 0.037001 | 11.85 % |

**Overlap targets `mpi_wait` alone: 0.96 % of the step.** The device-local work
that makes up the rest is 10.9 % — 11× larger. P1 removed the bytes P2 existed
to hide.

## Recommendation: do not implement P2

Three independent reasons, any one sufficient:

1. the mechanism does not work on this stack (measured above);
2. the only workaround is invasive and fragile, with no environment fallback;
3. the ceiling is now under 1 % of the step at 2 GPU ranks.

**What the exchange has become is a device-local copy problem**, not a
communication one: `local_copy` + `pack` + `unpack` = 10.9 % of the step. That
is kernel-efficiency work on the halo gather, which
`docs/next_session_block_overhead.md` already measured at its ~51 %-of-peak
coalescing floor — so it is not free either, but it is where the remaining
exchange time actually is.

Caveat on scope: this is 2 GPU ranks, the only multi-GPU configuration
available. Peer traffic grows with rank count (the partition analysis puts
16-rank peer area ~5× the 4-rank figure), so `mpi_wait` would grow too. If a
genuinely many-rank GPU campaign is ever run, A0 is worth repeating there
before this conclusion is carried over.
