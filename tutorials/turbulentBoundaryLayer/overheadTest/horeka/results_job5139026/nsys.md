# Timeline probe — is the wait device work, and is it every round?

No MPI events: OpenMPI's Fortran mpi_f08 bindings reach the C layer as
PMPI_*, so nsys's MPI_* interception never fires. The exchange's kernel
signature carries the same information: `pack -> copy_local -> [gap =
post + Waitall] -> unpack`.

**READ SHAPE, NOT MAGNITUDE.** Tracing inflates the wait unevenly
(measured: base 120 -> ~690 us, rect 738 -> ~1255 us), so no absolute
value here may be quoted against an untraced run. Ratios WITHIN one
traced run, and the distribution shape, are what this file is for.


## `base_jacobi_r8_N2` — 0 node(s): 

| rank | span ms | **GPU busy** | idle ms | gaps>50us | median | p90 | share of idle | rounds | mean wait | p50 | top 1 % | H2D n | P2P ms | copy ms |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 2 | 2751 | **79.0 %** | 568 | 1597 | 71 us | 122 | 52 % | 468 | 711.2 us | 492.0 | 17 % | 49823 | 124.7 | 0.0 |
| 3 | 2859 | **79.8 %** | 578 | 1673 | 72 us | 131 | 51 % | 493 | 686.7 us | 504.4 | 16 % | 52435 | 130.1 | 0.0 |

- GPU busy 79.0–79.8 % across ranks.
  **The GPU is idle through the wait — the stall is NOT unfinished
  device work.** A rank blocked on a CUDA stream would show the GPU busy.
- The slowest 1 % of rounds hold 16–17 % of the wait.

## `rect_jacobi_r4_N1` — 0 node(s): 

| rank | span ms | **GPU busy** | idle ms | gaps>50us | median | p90 | share of idle | rounds | mean wait | p50 | top 1 % | H2D n | P2P ms | copy ms |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 3794 | **90.1 %** | 375 | 554 | 98 us | 241 | 33 % | 689 | 327.4 us | 223.9 | 21 % | 48388 | 11.7 | 250.7 |

- GPU busy 90.1–90.1 % across ranks.
- The slowest 1 % of rounds hold 21–21 % of the wait.

## `rect_jacobi_r8_N2` — 0 node(s): 

| rank | span ms | **GPU busy** | idle ms | gaps>50us | median | p90 | share of idle | rounds | mean wait | p50 | top 1 % | H2D n | P2P ms | copy ms |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 2588 | **62.7 %** | 966 | 1098 | 365 us | 1586 | 78 % | 620 | 1255.3 us | 1069.0 | 5 % | 43517 | 10.7 | 109.2 |
| 1 | 3870 | **62.9 %** | 1435 | 1877 | 349 us | 1647 | 77 % | 934 | 1265.3 us | 1057.4 | 4 % | 65546 | 25.8 | 163.3 |

- GPU busy 62.7–62.9 % across ranks.
  **The GPU is idle through the wait — the stall is NOT unfinished
  device work.** A rank blocked on a CUDA stream would show the GPU busy.
- The slowest 1 % of rounds hold only 4–5 % of the
  wait: it is EVERY round, so a per-round mean is a fair description.
