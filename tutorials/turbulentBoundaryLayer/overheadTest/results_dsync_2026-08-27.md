# Phase 3 (vel_exchange) — divergence sync, A6000, 2026-08-27

Same-day A/B of `sync_divergence_halos` against the pre-change binary
`e77c75e`, both default flags, 400 cold steps, 138.4 M cells.

| layout | reference | candidate | gain | tax before | tax after |
|---|---|---|---|---|---|
| `nb` unset (control) | 1.1787 | 1.1770 | 0.1 % | — | — |
| `nb = 16` (cubic) | 1.6116 | 1.4673 | 9.0 % | 1.367 | 1.246 |
| **`nb = 64 44 48`** (production) | 1.2622 | **1.2192** | **3.4 %** | 1.070 | **1.037** |

The nb-unset row is the control: one block per rank means almost no halo
to sync, so there is almost nothing to gain — 0.1 % is noise.

The cubic gain is 2.6x the production one because its halo shell is far
bigger. Phase 1 had already removed most of what this saves; successive
optimisations on the same cost keep shrinking each other's returns.

nb-independence holds on both binaries (identical runtime lines after 400
steps on 138 M cells), and the 7-case suite is max_abs 0 on CPU and GPU.

```
run                       chron s/step   marginal   ratio    cells   leaves        L2_div     Linf_vel
--------------------------------------------------------------------------------------------------------
base_jacobi_dsync               1.1770     1.1772                         -  1.844810e-05 1.003416e+00
base_jacobi_dsyncref            1.1787     1.1786                         -  1.844810e-05 1.003416e+00
nb16_jacobi_dsync               1.4673     1.4666   1.246    1.000        -  1.844810e-05 1.003416e+00
nb16_jacobi_dsyncref            1.6116     1.6105   1.367    1.000        -  1.844810e-05 1.003416e+00
rect_jacobi_dsync               1.2192     1.2208   1.037    1.000        -  1.844810e-05 1.003416e+00
rect_jacobi_dsyncref            1.2622     1.2615   1.070    1.000        -  1.844810e-05 1.003416e+00
nb-independence nb16_jacobi_dsync vs base_jacobi_dsync: OK (runtime lines identical)
nb-independence nb16_jacobi_dsyncref vs base_jacobi_dsyncref: OK (runtime lines identical)
nb-independence rect_jacobi_dsync vs base_jacobi_dsync: OK (runtime lines identical)
nb-independence rect_jacobi_dsyncref vs base_jacobi_dsyncref: OK (runtime lines identical)
```
