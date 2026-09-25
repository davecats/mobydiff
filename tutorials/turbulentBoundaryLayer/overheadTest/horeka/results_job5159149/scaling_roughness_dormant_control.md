# The campaign, re-measured

The 23-run matrix (`run_matrix.sh`, `--map-by numa --bind-to core`), run
TWICE IN ONE ALLOCATION -- `ref` is the pre-change binary, `new` the
post-change one, and the two columns differ only in that change. WHICH
change is recorded in the run directory's provenance.txt, not here: this
collector has now served two of them (the rank-to-GPU mapping, 2026-09-10,
and the comm_type parent map, 2026-09-14), and hard-coding either into the
header is how a table gets quoted against the wrong question.

```
job    : 5159149
date   : 2026-09-25 11:36:27 CEST
nodes  : 4  (hkn[0521,0530-0531,0623])
new    : e5b5f8df8575bf6555f8ed911c761ae340e4a376  (head_e5b5f8d)
mid    : 8aedbb5d48e5d1e3bacc8d01cba9b3fada51de14  (rdenom_disabled_809759e)
ref    : b9414bd6a1d2b9e29178e5605a9ddd80a7db17b6  (published_point_b9414bd)
gpu    : NVIDIA A100-SXM4-40GB
nsteps : 200
layout : --map-by numa --bind-to core (as every earlier matrix)
```

## s/step, and what the mapping is worth

| config | ranks | ref | new | gain |
|---|---|---|---|---|
| `base_jacobi` | 1 | 0.56940 | 0.57245 | **-0.5 %** |
| `base_jacobi` | 2 | 0.28403 | 0.28570 | **-0.6 %** |
| `base_jacobi` | 4 | 0.14713 | 0.14796 | **-0.6 %** |
| `base_jacobi` | 8 | 0.07858 | 0.07902 | **-0.6 %** |
| `base_jacobi` | 16 | 0.04326 | 0.04351 | **-0.6 %** |
| `rect_jacobi` | 1 | 0.57641 | 0.57913 | **-0.5 %** |
| `rect_jacobi` | 2 | 0.29409 | 0.29535 | **-0.4 %** |
| `rect_jacobi` | 4 | 0.15117 | 0.15184 | **-0.4 %** |
| `rect_jacobi` | 8 | 0.07896 | 0.07937 | **-0.5 %** |
| `rect_jacobi` | 16 | 0.04465 | 0.04490 | **-0.5 %** |
| `refined_yp82_rect_jacobi` | 1 | 0.25781 | 0.25888 | **-0.4 %** |
| `refined_yp82_rect_jacobi` | 2 | 0.13428 | 0.13490 | **-0.5 %** |
| `refined_yp82_rect_jacobi` | 4 | 0.07166 | 0.07189 | **-0.3 %** |
| `refined_yp82_rect_jacobi` | 8 | 0.03997 | 0.04020 | **-0.6 %** |
| `refined_yp82_rect_jacobi` | 16 | 0.02538 | 0.02531 | **+0.3 %** |
| `refined_yp82_rect_redblack` | 1 | 0.21953 | 0.21960 | **-0.0 %** |
| `refined_yp82_rect_redblack` | 2 | 0.11567 | 0.11562 | **+0.0 %** |
| `refined_yp82_rect_redblack` | 4 | 0.06255 | 0.06247 | **+0.1 %** |
| `refined_yp82_rect_redblack` | 8 | 0.03639 | 0.03639 | **+0.0 %** |
| `refined_yp82_rect_redblack` | 16 | 0.02544 | 0.02511 | **+1.3 %** |
| `refined_big_rect_jacobi` | 4 | 0.26400 | 0.26521 | **-0.5 %** |
| `refined_big_rect_jacobi` | 8 | 0.13672 | 0.13736 | **-0.5 %** |
| `refined_big_rect_jacobi` | 16 | 0.07489 | 0.07506 | **-0.2 %** |

## Per-round `mpi_wait`

| config | n=1 | n=2 | n=4 | n=8 | n=16 |
|---|---|---|---|---|---|
| `base_jacobi` (ref) | - | 38.1 us | 100.0 us | 78.4 us | 95.2 us |
| `base_jacobi` (new) | - | 41.3 us | 98.4 us | 77.1 us | 96.1 us |
| `rect_jacobi` (ref) | - | 51.3 us | 63.6 us | 50.6 us | 82.9 us |
| `rect_jacobi` (new) | - | 49.6 us | 65.6 us | 51.3 us | 85.0 us |
| `refined_yp82_rect_jacobi` (ref) | - | 28.2 us | 47.2 us | 33.8 us | 54.6 us |
| `refined_yp82_rect_jacobi` (new) | - | 28.8 us | 42.0 us | 38.5 us | 53.2 us |
| `refined_yp82_rect_redblack` (ref) | - | 28.7 us | 42.3 us | 51.6 us | 106.5 us |
| `refined_yp82_rect_redblack` (new) | - | 29.4 us | 38.0 us | 51.0 us | 101.4 us |
| `refined_big_rect_jacobi` (ref) | - | - | 86.7 us | 57.4 us | 91.6 us |
| `refined_big_rect_jacobi` (new) | - | - | 85.4 us | 57.3 us | 89.0 us |

## Block tax — inverse per-GPU throughput `t_step * n / cells` (ns)

| ranks | Mcell/GPU | base ref | base new | rect ref | rect new | **tax ref** | **tax new** |
|---|---|---|---|---|---|---|---|
| 1 | 138.41 | 4.114 | 4.136 | 4.164 | 4.184 | **1.012** | **1.012** |
| 2 | 69.21 | 4.104 | 4.128 | 4.250 | 4.268 | **1.035** | **1.034** |
| 4 | 34.60 | 4.252 | 4.276 | 4.369 | 4.388 | **1.028** | **1.026** |
| 8 | 17.30 | 4.542 | 4.567 | 4.564 | 4.587 | **1.005** | **1.004** |
| 16 | 8.65 | 5.000 | 5.030 | 5.162 | 5.190 | **1.032** | **1.032** |

## Strong-scaling efficiency (vs each config's smallest rank count)

| config | n=1 | n=2 | n=4 | n=8 | n=16 |
|---|---|---|---|---|---|
| `base_jacobi` (ref) | 100 % | 100 % | 97 % | 91 % | 82 % |
| `base_jacobi` (new) | 100 % | 100 % | 97 % | 91 % | 82 % |
| `rect_jacobi` (ref) | 100 % | 98 % | 95 % | 91 % | 81 % |
| `rect_jacobi` (new) | 100 % | 98 % | 95 % | 91 % | 81 % |
| `refined_yp82_rect_jacobi` (ref) | 100 % | 96 % | 90 % | 81 % | 63 % |
| `refined_yp82_rect_jacobi` (new) | 100 % | 96 % | 90 % | 80 % | 64 % |
| `refined_yp82_rect_redblack` (ref) | 100 % | 95 % | 88 % | 75 % | 54 % |
| `refined_yp82_rect_redblack` (new) | 100 % | 95 % | 88 % | 75 % | 55 % |
| `refined_big_rect_jacobi` (ref) | - | - | 100 % | 97 % | 88 % |
| `refined_big_rect_jacobi` (new) | - | - | 100 % | 97 % | 88 % |

## Headline 1 -- the 2:1 machinery, like for like

`refined_big_rect_jacobi` shares grid, block shape and solver with
`rect_jacobi`; it just refines one y-tile. Cost of the ADDED cells
against the coarse cells beside them:

| ranks | binary | rect s/step | big s/step | ns per EXTRA cell | vs coarse-cell avg |
|---|---|---|---|---|---|
| 4 | ref | 0.15117 | 0.26400 | 1.0869 | **0.995x** |
| 4 | new | 0.15184 | 0.26521 | 1.0921 | **0.996x** |
| 8 | ref | 0.07896 | 0.13672 | 0.5564 | **0.975x** |
| 8 | new | 0.07937 | 0.13736 | 0.5586 | **0.974x** |
| 16 | ref | 0.04465 | 0.07489 | 0.2912 | **0.903x** |
| 16 | new | 0.04490 | 0.07506 | 0.2905 | **0.896x** |

## Headline 3 — red-black against Jacobi

| binary | n=1 | n=2 | n=4 | n=8 | n=16 |
|---|---|---|---|---|---|
| ref | 0.852 | 0.861 | 0.873 | 0.910 | 1.003 |
| new | 0.848 | 0.857 | 0.869 | 0.905 | 0.992 |
