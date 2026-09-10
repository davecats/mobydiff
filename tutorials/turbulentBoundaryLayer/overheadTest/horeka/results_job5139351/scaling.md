# The campaign, re-measured at the corrected rank-to-GPU mapping

Same 23-run matrix as `results_horeka_2026-09-07.md` (`run_matrix.sh`,
`--map-by numa --bind-to core`, 400 steps), run twice in ONE allocation:
`ref` is the pre-change binary (`device = local_rank mod ndev`), `new`
the topology-aware mapping. The two columns differ only in the mapping.

## s/step, and what the mapping is worth

| config | ranks | ref | new | gain |
|---|---|---|---|---|
| `base_jacobi` | 1 | 0.68429 | 0.68301 | **+0.2 %** |
| `base_jacobi` | 2 | 0.34889 | 0.34885 | **+0.0 %** |
| `base_jacobi` | 4 | 0.18628 | 0.18628 | **-0.0 %** |
| `base_jacobi` | 8 | 0.10277 | 0.10241 | **+0.3 %** |
| `base_jacobi` | 16 | 0.05994 | 0.06003 | **-0.2 %** |
| `rect_jacobi` | 1 | 0.71874 | 0.71847 | **+0.0 %** |
| `rect_jacobi` | 2 | 0.38212 | 0.38215 | **-0.0 %** |
| `rect_jacobi` | 4 | 0.19961 | 0.19952 | **+0.0 %** |
| `rect_jacobi` | 8 | 0.13460 | 0.10848 | **+19.4 %** |
| `rect_jacobi` | 16 | 0.08985 | 0.06784 | **+24.5 %** |
| `refined_yp82_rect_jacobi` | 1 | 0.32383 | 0.32402 | **-0.1 %** |
| `refined_yp82_rect_jacobi` | 2 | 0.17884 | 0.17882 | **+0.0 %** |
| `refined_yp82_rect_jacobi` | 4 | 0.09947 | 0.09951 | **-0.0 %** |
| `refined_yp82_rect_jacobi` | 8 | 0.07606 | 0.05977 | **+21.4 %** |
| `refined_yp82_rect_jacobi` | 16 | 0.05625 | 0.04119 | **+26.8 %** |
| `refined_yp82_rect_redblack` | 1 | 0.24666 | 0.24670 | **-0.0 %** |
| `refined_yp82_rect_redblack` | 2 | 0.13659 | 0.13654 | **+0.0 %** |
| `refined_yp82_rect_redblack` | 4 | 0.07787 | 0.07800 | **-0.2 %** |
| `refined_yp82_rect_redblack` | 8 | 0.06520 | 0.04868 | **+25.3 %** |
| `refined_yp82_rect_redblack` | 16 | 0.04984 | 0.03535 | **+29.1 %** |
| `refined_big_rect_jacobi` | 4 | 0.34643 | 0.34662 | **-0.1 %** |
| `refined_big_rect_jacobi` | 8 | 0.21437 | 0.18528 | **+13.6 %** |
| `refined_big_rect_jacobi` | 16 | 0.13515 | 0.10722 | **+20.7 %** |

## Per-round `mpi_wait`

| config | n=1 | n=2 | n=4 | n=8 | n=16 |
|---|---|---|---|---|---|
| `base_jacobi` (ref) | - | 29.0 us | 116.9 us | 130.8 us | 147.4 us |
| `base_jacobi` (new) | - | 29.9 us | 119.0 us | 122.7 us | 149.6 us |
| `rect_jacobi` (ref) | - | 36.0 us | 65.2 us | 700.5 us | 707.8 us |
| `rect_jacobi` (new) | - | 36.4 us | 64.7 us | 76.8 us | 175.6 us |
| `refined_yp82_rect_jacobi` (ref) | - | 34.9 us | 52.9 us | 440.2 us | 432.4 us |
| `refined_yp82_rect_jacobi` (new) | - | 38.2 us | 49.0 us | 52.4 us | 81.6 us |
| `refined_yp82_rect_redblack` (ref) | - | 23.0 us | 46.6 us | 449.7 us | 426.6 us |
| `refined_yp82_rect_redblack` (new) | - | 34.6 us | 50.6 us | 44.4 us | 80.5 us |
| `refined_big_rect_jacobi` (ref) | - | - | 79.1 us | 782.2 us | 810.7 us |
| `refined_big_rect_jacobi` (new) | - | - | 80.4 us | 80.6 us | 126.7 us |

## Block tax — inverse per-GPU throughput `t_step * n / cells` (ns)

| ranks | Mcell/GPU | base ref | base new | rect ref | rect new | **tax ref** | **tax new** |
|---|---|---|---|---|---|---|---|
| 1 | 138.41 | 4.944 | 4.935 | 5.193 | 5.191 | **1.050** | **1.052** |
| 2 | 69.21 | 5.041 | 5.041 | 5.522 | 5.522 | **1.095** | **1.095** |
| 4 | 34.60 | 5.383 | 5.383 | 5.769 | 5.766 | **1.072** | **1.071** |
| 8 | 17.30 | 5.940 | 5.919 | 7.780 | 6.270 | **1.310** | **1.059** |
| 16 | 8.65 | 6.929 | 6.939 | 10.386 | 7.842 | **1.499** | **1.130** |

## Strong-scaling efficiency (vs each config's smallest rank count)

| config | n=1 | n=2 | n=4 | n=8 | n=16 |
|---|---|---|---|---|---|
| `base_jacobi` (ref) | 100 % | 98 % | 92 % | 83 % | 71 % |
| `base_jacobi` (new) | 100 % | 98 % | 92 % | 83 % | 71 % |
| `rect_jacobi` (ref) | 100 % | 94 % | 90 % | 67 % | 50 % |
| `rect_jacobi` (new) | 100 % | 94 % | 90 % | 83 % | 66 % |
| `refined_yp82_rect_jacobi` (ref) | 100 % | 91 % | 81 % | 53 % | 36 % |
| `refined_yp82_rect_jacobi` (new) | 100 % | 91 % | 81 % | 68 % | 49 % |
| `refined_yp82_rect_redblack` (ref) | 100 % | 90 % | 79 % | 47 % | 31 % |
| `refined_yp82_rect_redblack` (new) | 100 % | 90 % | 79 % | 63 % | 44 % |
| `refined_big_rect_jacobi` (ref) | - | - | 100 % | 81 % | 64 % |
| `refined_big_rect_jacobi` (new) | - | - | 100 % | 94 % | 81 % |

## Headline 1 — the 2:1 machinery, like for like

`refined_big_rect_jacobi` shares grid, block shape and solver with
`rect_jacobi`; it just refines one y-tile. Cost of the ADDED cells
against the coarse cells beside them:

| ranks | binary | rect s/step | big s/step | ns per EXTRA cell | vs coarse-cell avg |
|---|---|---|---|---|---|
| 4 | ref | 0.19961 | 0.34643 | 1.4143 | **0.981x** |
| 4 | new | 0.19952 | 0.34662 | 1.4171 | **0.983x** |
| 8 | ref | 0.13460 | 0.21437 | 0.7684 | **0.790x** |
| 8 | new | 0.10848 | 0.18528 | 0.7397 | **0.944x** |
| 16 | ref | 0.08985 | 0.13515 | 0.4364 | **0.672x** |
| 16 | new | 0.06784 | 0.10722 | 0.3793 | **0.774x** |

## Headline 3 — red-black against Jacobi

| binary | n=1 | n=2 | n=4 | n=8 | n=16 |
|---|---|---|---|---|---|
| ref | 0.762 | 0.764 | 0.783 | 0.857 | 0.886 |
| new | 0.761 | 0.764 | 0.784 | 0.814 | 0.858 |
