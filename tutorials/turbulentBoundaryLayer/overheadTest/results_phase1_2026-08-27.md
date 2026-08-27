# Phase 1 — per-direction nb, A6000 (istmcetus GPU 1), 2026-08-27

Same-day A/B, one binary, 400 cold steps, 4096 x 176 x 192 = 138.4 M cells.

| layout | s/step | tax | blocks | halo/interior |
|---|---|---|---|---|
| `nb` unset (one block per rank box) | 1.1771 | — | 1 | — |
| `nb = 16` (cubic) | 1.6123 | 1.367 | 33792 | 42.4 % |
| **`nb = 64 44 48`** | **1.2622** | **1.071** | 1024 | 12.3 % |

**21.7 % faster than cubic `nb = 16`**; the block tax falls from 37 points
to 7, i.e. 81 % of what remained after the copy-kernel split.

The allocated-volume law predicted 1.1230 and the measured tax is 1.0722 --
the law now OVER-predicts at large nb. It was calibrated against the old
copy kernel, where cost tracked halo CELLS; with 1024 blocks instead of
33792 there are 33x fewer exchange ENTRIES, and the per-entry metadata
(decode, slot lookup, gather map) falls faster than the cell count does.
Treat the law as an upper bound on the tax, not an estimate of it.

nb-independence holds end to end: `rect_jacobi` and `base_jacobi` produce
identical runtime lines (L2_div 1.844810E-05, Linf 1.00341600E+00) after
400 steps on 138 M cells.

```
run                       chron s/step   marginal   ratio    cells   leaves        L2_div     Linf_vel
--------------------------------------------------------------------------------------------------------
base_jacobi_p1                  1.1771     1.1788                         -  1.844810e-05 1.003416e+00
nb16_jacobi_p1                  1.6123     1.6116   1.367    1.000        -  1.844810e-05 1.003416e+00
rect_jacobi_p1                  1.2622     1.2620   1.071    1.000        -  1.844810e-05 1.003416e+00
nb-independence nb16_jacobi_p1 vs base_jacobi_p1: OK (runtime lines identical)
nb-independence rect_jacobi_p1 vs base_jacobi_p1: OK (runtime lines identical)
```
