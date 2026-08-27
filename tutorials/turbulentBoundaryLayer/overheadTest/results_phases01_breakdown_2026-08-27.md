# Per-phase breakdown after Phases 0-1 — A6000, 2026-08-27

Post copy-kernel-split and per-direction nb. `rect_jacobi_p3` (nb = 64 44 48)
is the current production candidate; `nb16_jacobi_p3` is the cubic control,
both measured the same day with the same binary.

Where the step goes at rect_jacobi_p3 (total 1.2596 s/step):

| phase | s/step | share |
|---|---|---|
| proj/sweep | 0.2806 | 22.3 % |
| proj/apply | 0.4072 | 32.3 % |
| proj/vel_exchange | 0.0664 | 5.3 % |
| proj/phi_exchange | 0.0229 | 1.8 % |
| exch/local_copy (inside all exchanges) | 0.0976 | 7.8 % |

This is why Phase 3a (R0) was not pursued: its entire target is the 1.8 %.
See docs/next_session_block_overhead.md STATUS 2026-08-27 (Phase 3a).

```
| phase                  |  base_redblack |  nb16_redblack |    base_jacobi |     nb8_jacobi |    nb16_jacobi | nb16_jacobi_p3 | rect_jacobi_p3 | refined_jacobi |     nb16_redblack/base |        nb8_jacobi/base |       nb16_jacobi/base |    refined_jacobi/base |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| step/momentum          |        0.37113 |        0.36614 |        0.36331 |        0.38180 |        0.36760 |        0.36007 |        0.35048 |        0.16990 |                  0.987 |                  1.051 |                  1.012 |                  0.462 |
| step/ibm_mu            |        0.03714 |        0.05109 |        0.03706 |        0.07058 |        0.05348 |        0.05001 |        0.04003 |        0.02322 |                  1.376 |                  1.904 |                  1.443 |                  0.434 |
| step/bodyforce         |        0.06114 |        0.06535 |        0.06442 |        0.06653 |        0.06124 |        0.05963 |        0.05955 |        0.03954 |                  1.069 |                  1.033 |                  0.951 |                  0.646 |
| step/turbulence        |        0.00000 |        0.00000 |        0.00000 |        0.00000 |        0.00000 |        0.00000 |        0.00000 |        0.00000 |                        |                        |                        |                        |
| step/apply_bc          |        0.00154 |        0.00168 |        0.00147 |        0.00219 |        0.00163 |        0.00161 |        0.00139 |        0.00105 |                  1.093 |                  1.488 |                  1.108 |                  0.642 |
| step/vel_exchange      |        0.00166 |        0.04608 |        0.00115 |        0.09111 |        0.04404 |        0.03647 |        0.01045 |        0.02032 |                 27.688 |                 79.209 |                 38.289 |                  0.461 |
| step/projection        |        1.02394 |        1.50398 |        0.76135 |        1.77432 |        1.21213 |        1.09162 |        0.78624 |        0.60560 |                  1.469 |                  2.330 |                  1.592 |                  0.500 |
| step/io_stats          |        0.01171 |        0.01392 |        0.01239 |        0.01870 |        0.01432 |        0.01386 |        0.01141 |        0.00695 |                  1.189 |                  1.510 |                  1.156 |                  0.485 |
| proj/sweep             |        0.97393 |        0.86762 |        0.30773 |        0.36700 |        0.33574 |        0.32948 |        0.28060 |        0.16355 |                  0.891 |                  1.193 |                  1.091 |                  0.487 |
| proj/apply             |        0.00000 |        0.00000 |        0.43079 |        0.53084 |        0.45329 |        0.44564 |        0.40723 |        0.22072 |                        |                  1.232 |                  1.052 |                  0.487 |
| proj/phi_exchange      |        0.00000 |        0.00000 |        0.00541 |        0.20825 |        0.10053 |        0.07396 |        0.02285 |        0.05879 |                        |                 38.489 |                 18.580 |                  0.585 |
| proj/vel_exchange      |        0.03060 |        0.61666 |        0.00748 |        0.65289 |        0.31000 |        0.23187 |        0.06641 |        0.13371 |                 20.155 |                 87.229 |                 41.417 |                  0.431 |
| proj/apply_bc          |        0.01938 |        0.01968 |        0.00755 |        0.01076 |        0.00925 |        0.00866 |        0.00750 |        0.02731 |                  1.016 |                  1.424 |                  1.225 |                  2.951 |
| proj/setup             |        0.00000 |        0.00000 |        0.00236 |        0.00457 |        0.00329 |        0.00199 |        0.00163 |        0.00151 |                        |                  1.932 |                  1.393 |                  0.458 |
| exch/pack              |        0.00000 |        0.00000 |        0.00000 |        0.00000 |        0.00000 |        0.00000 |        0.00000 |        0.00000 |                        |                        |                        |                        |
| exch/mpi_post          |        0.00000 |        0.00000 |        0.00000 |        0.00000 |        0.00000 |        0.00000 |        0.00000 |        0.00000 |                        |                        |                        |                        |
| exch/mpi_wait          |        0.00000 |        0.00000 |        0.00000 |        0.00000 |        0.00000 |        0.00000 |        0.00000 |        0.00000 |                        |                        |                        |                        |
| exch/unpack            |        0.00000 |        0.00000 |        0.00000 |        0.00000 |        0.00000 |        0.00000 |        0.00000 |        0.00000 |                        |                        |                        |                        |
| exch/local_copy        |        0.03214 |        0.66268 |        0.01197 |        0.94953 |        0.45224 |        0.34006 |        0.09764 |        0.21164 |                 20.617 |                 79.312 |                 37.774 |                  0.468 |
| TOTAL/chron            |        1.50827 |        2.04825 |        1.24115 |        2.40523 |        1.75445 |        1.61327 |        1.25957 |        0.86658 |                  1.358 |                  1.938 |                  1.414 |                  0.494 |

allocated-volume prediction, nb16: 1.4238
allocated-volume prediction, nb8: 1.9531
base_redblack_prof: coverage 1.0000
nb16_redblack_prof: coverage 1.0000
base_jacobi_prof: coverage 1.0000
nb8_jacobi_prof: coverage 1.0000
nb16_jacobi_prof: coverage 1.0000
nb16_jacobi_prof_p3: coverage 1.0000
rect_jacobi_prof_p3: coverage 1.0000
refined_yp100_jacobi_prof: coverage 1.0000
```
