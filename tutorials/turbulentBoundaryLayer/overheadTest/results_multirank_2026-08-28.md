# The refined config at multiple ranks, istmcetus, 2026-08-28

`multiLevel_xz/refined_yp82_rect_jacobi.ini` on 4 CPU ranks and 2 GPUs, each
against a **1-rank baseline measured on the same host in the same session**.
Rates do not transfer between hosts (istmcorax > istmcetus > the local
workstation), so every table below is internally paired and no number here
should be compared with one from another machine.

Binary `0f91990` (source identical to `2e74cc4`). The GPU runs are the standard
400 steps; the CPU runs are **20 steps** (`NSTEPS=20`) because 400 would be
2.2 h at 4 ranks and 9 h at 1 — phase shares converge long before that, but the
CPU absolute rates are only comparable to each other.

## Correctness first

1-rank and 2-rank GPU produce **identical runtime lines after 400 steps**
(`L2_div 2.93355685E-05`, `global_div -3.24925358E+00`, `Linf 1.00242065E+00`),
and the CPU runs at 1 and 4 ranks agree with each other and with the GPU on the
first steps (`L2_div 1.09634844E-06`, `1.63250315E-06`). Rank-count and
CPU/GPU invariance hold on the multi-level MPI path for this config.

## Scaling

| | 1 rank | N ranks | speedup | efficiency |
|---|---|---|---|---|
| GPU, N = 2 (A6000 ×2) | 0.55400 s/step | 0.32436 | 1.708× | **85.4 %** |
| CPU, N = 4 | 63.789 s/step | 18.769 | 3.399× | **85.0 %** |

The two platforms land on the same efficiency from very different absolute
rates (the CPU step is 115× the GPU step), and for the same reason.

## The volume kernels scale perfectly; the exchange does not

Efficiency per phase — `(1 rank / N) / (N ranks)`, so 1.00 is perfect:

| phase | GPU 1→2 | CPU 1→4 |
|---|---|---|
| momentum | 1.02 | 1.02 |
| proj/sweep | 1.02 | 0.99 |
| proj/apply | 1.00 | 0.98 |
| ibm_mu | 1.03 | 1.06 |
| **proj/phi_exchange** | **0.32** | **0.33** |
| **proj/vel_exchange** | **0.25** | **0.21** |
| step/vel_exchange | 0.31 | 0.80 |
| bodyforce | 0.68 | 0.40 |

Both projection exchanges get *absolutely* more expensive while the work per
rank halves or quarters:

| | GPU 1→2 | CPU 1→4 |
|---|---|---|
| proj/phi_exchange | 0.01328 → 0.02094 (+58 %) | 2.0755 → 1.5534 (−25 %) |
| proj/vel_exchange | 0.01318 → 0.02603 (+97 %) | 1.1593 → 1.3830 (+19 %) |

**The projection exchange is 71 % (GPU) / 75 % (CPU) of all the time lost to
imperfect scaling** — 0.0337 of 0.0474 s/step, and 2.128 of 2.821 s/step.

## Where the exchange time goes

`exch_timing` is exactly zero in the MPI buckets at 1 rank by construction, so
these columns are the multi-rank cost appearing from nothing:

| bucket | GPU 2 ranks | % of step | CPU 4 ranks | % of step |
|---|---|---|---|---|
| pack | 0.00534 | 1.6 % | 0.2404 | 1.3 % |
| mpi_post | 0.00010 | 0.0 % | 0.0005 | 0.0 % |
| **mpi_wait** | **0.02335** | **7.2 %** | **1.5126** | **8.1 %** |
| unpack | 0.00512 | 1.6 % | 0.3923 | 2.1 % |
| local_copy | 0.02146 | 6.6 % | 1.0477 | 5.6 % |
| **total** | **0.05537** | **17.1 %** | **3.1934** | **17.0 %** |
| (1 rank, local_copy only) | 0.02520 | 4.5 % | 3.8489 | 6.0 % |

The exchange grows from ~5 % of the step to **17 % on both platforms**, and
`mpi_wait` alone is 7–8 %. `mpi_post` is nil, so this is not posting overhead:
the ranks are waiting on data.

## The named cause, now measured

`sync_divergence_halos` is **single-rank only**. Between Jacobi iterations a
1-rank run refreshes just one plane per dimension of the normal velocity
component; with MPI peers `pressure_projection` falls back to a full three-
component `exchange_halos(interp=.false.)`, so 15 of the 18 per-step velocity
refreshes revert to the expensive form. That is exactly the phase whose
efficiency is worst (0.21–0.25) and the only one that grows in absolute terms
on both platforms.

Estimating what it costs, from the measured cost of one full exchange
(`step/vel_exchange` ÷ 3) and of one sync (the 1-rank residual):

| | one full exch | one sync | 18 full (now) | 3 full + 15 sync | saving |
|---|---|---|---|---|---|
| GPU 2 ranks | 0.00301 | ~0.0004 | 0.02603 | ~0.0103 | **~0.016 s/step, ~5 %** |
| CPU 4 ranks | 0.0768 | ~0.0059 | 1.383 | ~0.319 | **~1.06 s/step, ~6 %** |

Rough — the per-sync cost at N ranks is extrapolated from the 1-rank ratio, not
measured — but the two platforms agree, and the saving grows with rank count
because the fallback message is 3 components against the sync's 1 plane.

## What to do

**Partition the exchange entry list so `sync_divergence_halos` works with
peers.** It is the top item for any multi-rank 2:1 run: worth ~5–6 % of the step
at these modest rank counts, on the phase that scales worst, with the mechanism
identified and already documented in `pressure_solver.f90`. This supersedes the
ranking in `results_refined_rect_2026-08-28.md`, which could only see the
1-rank picture where the projection exchange is 4.8 % of the step.

Two secondary items this exposed:

- **`bodyforce` scales at 0.68 (GPU) / 0.40 (CPU)** and is 4.8–8.5 % of the
  multi-rank step. `update_bodyforce` fills a per-block array and should be
  embarrassingly parallel, so something in the trip forcing is not — worth a
  look before it is worth optimising.
- **`proj/phi_exchange` is no better (0.32/0.33)**. It is 18 small scalar
  exchanges per step, so it is latency-bound rather than bandwidth-bound; the
  divergence-sync fix does not touch it, and overlap
  (`docs/nonblocking_overlap_strategy.md`, which needs rewriting for the
  Chebyshev-Jacobi solver) is the lever that would.

## Reproducing

```bash
BIN=./gpu_rank.sh NRANKS=2 GPUS=0,1 PROFILE=1 SUFFIX=gpu2 \
    ./run_overhead.sh multiLevel_xz/refined_yp82_rect_jacobi.ini
BIN=../../../build_cpu/moby_solve NRANKS=4 NSTEPS=20 PROFILE=1 SUFFIX=cpu4 \
    ./run_overhead.sh multiLevel_xz/refined_yp82_rect_jacobi.ini
```

`gpu_rank.sh` is required for any multi-GPU run: the solver never calls
`omp_set_default_device`, so without per-rank `CUDA_VISIBLE_DEVICES` every rank
offloads to the same physical card. Check for "OpenMP target devices available:
1" per rank — a 2 means the pinning did not take and the run is two ranks
sharing one GPU.
