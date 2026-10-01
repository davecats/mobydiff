# HoreKa, 2026-10-01: the minimum-surface block order — the node-boundary finding of 09-30 is resolved

Job **5174076**, `accelerated`, 2 nodes (hkn[0526,0535]), 8 x A100-40GB,
nvhpc 25.3, `horeka/exchange/submit_order.sh`, 9 min wall. Two binaries
built in the job:

| column | commit | block order of its case files |
|---|---|---|
| ref | `5bdc5eb` | legacy interleave (x lowest, z on top) |
| new | `d2ac839` | `min_surface_key_order`, recorded in the case file |

Shipped configs on both sides, each run directory prepared by its own
column's `moby_prepare`, `--map-by numa --bind-to core` (the DEFAULT
placement of every matrix), 200 steps. Raw material:
`horeka/results_job5174076/`. Pre-registration:
`horeka/exchange/PREREGISTERED_order.md`.

## 1. The step

| config | ranks | ref | new | change | `mpi_wait` ref → new |
|---|---|---|---|---|---|
| `base_jacobi` | 4 | 141.45 ms | 141.24 ms | −0.15 % | 4.20 → 4.34 ms |
| `base_jacobi` | **8** | **130.67 ms** | **75.93 ms** | **−41.9 %** | **59.07 → 3.63 ms** |
| `rect_jacobi` | 4 | 146.10 ms | 145.00 ms | −0.75 % | 3.52 → 3.28 ms |
| `rect_jacobi` | 8 | 76.06 ms | 75.72 ms | −0.44 % | 2.32 → 2.44 ms |
| `refined_yp82_rect_jacobi` | 4 | 69.07 ms | 69.06 ms | −0.01 % | 2.03 → 2.09 ms |
| `refined_yp82_rect_jacobi` | 8 | 38.43 ms | 38.45 ms | +0.07 % | 1.43 → 1.61 ms |

`base_jacobi` at 8 ranks, phase by phase: projection 95.26 → 51.70 ms,
momentum 19.95 → 20.22 ms; the wait balance over ranks goes from 1.03
(every rank waits equally: transfer) to 1.82 (the usual skew). The
`gpu binding` line is the same in both columns, so nothing but the block
order changed.

The solver's init line, new column:

```
base_jacobi  n=8   partition: face cells shared across ranks 2262016 (8 ranks), across nodes 33792 (2 nodes)
rect_jacobi  n=8   partition: face cells shared across ranks 236544 (8 ranks), across nodes 33792 (2 nodes)
refined_yp82 n=8   partition: face cells shared across ranks 147840 (8 ranks), across nodes 21120 (2 nodes)
base_jacobi  n=4   partition: face cells shared across ranks 820224 (4 ranks), across nodes 0 (1 nodes)
```

`tools/partition_analysis.py lattice ... --order legacy` gives 1,441,792
across nodes for the ref column's `base_jacobi` at 8 ranks: 42.7 times the
cells, 16 times the wait.

## 2. The field gate

`base_jacobi` and `rect_jacobi`, 8 ranks, 20 steps, final snapshot of the
two columns through `tools/h5maxdiff` (built in the job from the new tree),
PRODUCTION flags:

```
  (block row orders differ: rows matched on origin + level)
  un           n=138412032    max_abs=0
  vn           n=138412032    max_abs=0
  wn           n=138412032    max_abs=0
  pn           n=138412032    max_abs=0
OK: worst max_abs = 0 over 4 dataset(s)
```

on both configs. The two columns own their blocks differently (and write
them in different row orders); the fields are the same numbers.

## 3. Against the pre-registration

| # | prediction | measured | |
|---|---|---|---|
| 1 | init line: `base_jacobi` n=8 2262016 / 33792 (2 nodes); `rect_jacobi` n=8 236544 / 33792; n=4 across nodes 0 | exactly these | HIT |
| 2 | `base_jacobi` n=8 ref 132 ms ± 1 %, `mpi_wait` ≈ 59 ms | 130.67 ms (−1.0 %), 59.07 ms | HIT (at the edge of the band) |
| 2 | `base_jacobi` n=8 new **77.7 ms ± 1 %**, `mpi_wait` 3.0–3.6 ms | **75.93 ms (−2.3 %)**, 3.63 ms | **MISSED, on the fast side** |
| 3 | `base_jacobi` n=4 unchanged within 0.5 % | −0.15 % | HIT |
| 4 | `rect_jacobi` n=4, n=8 unchanged within 0.5 % | −0.75 %, −0.44 % | n=8 HIT, n=4 just outside |
| 5 | `refined_yp82_rect_jacobi` unchanged within 0.5 % | −0.01 %, +0.07 % | HIT |
| 6 | field gate `max_abs 0` with the row-order note | as predicted, both configs | HIT |
| 7 | 16 ranks | job 5174077 queued | pending |

The miss of row 2 is in the pre-registered number, not in the change: 77.7 ms
was the 09-26 value of the old rank box, measured BEFORE the predictor
register guard (−1.4 to −2.0 % of the step, `results_horeka_2026-09-30.md`
section 2). 77.7 ms less 1.95 % is 76.2 ms; the same allocation's
`rect_jacobi` reads 75.72 ms. A pre-registered absolute time must be taken
from the head it is compared with. `rect_jacobi` at 4 ranks (−0.75 %) has
the same cross-rank cell count in both orders (101,376); the difference is
inside what one pair in one allocation resolves, and it is the only run
whose `mpi_wait` fell (3.52 → 3.28 ms).

## 4. What this restores

The block tax can be formed at 8 ranks again from a default-placed
`base_jacobi`: `rect / base` = 75.72 / 75.93 = **0.997** (published 09-26:
0.996), and 145.00 / 141.24 = 1.027 at 4 ranks (1.025). The warning of
`results_horeka_2026-09-30.md` section 3 ("no block tax at 8 ranks and
above from a post-step-7 `base_jacobi`") applies to case files prepared
before `d2ac839` only: re-prepare, and read the `partition:` line.

16 ranks (4 nodes) is job 5174077, same script with
`RANKS_SMALL=16 RANKS_BIG=16 GATE=0 CONFIGS="base_jacobi rect_jacobi"`, run
directory `order16_run`; `base_jacobi` is 4 x 2 x 2 blocks there, order
x x y z, 101,376 cells across nodes against 1,475,584 in the legacy order.
