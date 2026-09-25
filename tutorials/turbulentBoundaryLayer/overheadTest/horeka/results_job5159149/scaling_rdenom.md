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
| `base_jacobi` | 1 | 0.57245 | 0.54925 | **+4.1 %** |
| `base_jacobi` | 2 | 0.28570 | 0.27427 | **+4.0 %** |
| `base_jacobi` | 4 | 0.14796 | 0.14230 | **+3.8 %** |
| `base_jacobi` | 8 | 0.07902 | 0.07617 | **+3.6 %** |
| `base_jacobi` | 16 | 0.04351 | 0.04198 | **+3.5 %** |
| `rect_jacobi` | 1 | 0.57913 | 0.55716 | **+3.8 %** |
| `rect_jacobi` | 2 | 0.29535 | 0.28430 | **+3.7 %** |
| `rect_jacobi` | 4 | 0.15184 | 0.14628 | **+3.7 %** |
| `rect_jacobi` | 8 | 0.07937 | 0.07653 | **+3.6 %** |
| `rect_jacobi` | 16 | 0.04490 | 0.04325 | **+3.7 %** |
| `refined_yp82_rect_jacobi` | 1 | 0.25888 | 0.24941 | **+3.7 %** |
| `refined_yp82_rect_jacobi` | 2 | 0.13490 | 0.13017 | **+3.5 %** |
| `refined_yp82_rect_jacobi` | 4 | 0.07189 | 0.06951 | **+3.3 %** |
| `refined_yp82_rect_jacobi` | 8 | 0.04020 | 0.03890 | **+3.2 %** |
| `refined_yp82_rect_jacobi` | 16 | 0.02531 | 0.02474 | **+2.3 %** |
| `refined_yp82_rect_redblack` | 1 | 0.21960 | 0.21947 | **+0.1 %** |
| `refined_yp82_rect_redblack` | 2 | 0.11562 | 0.11570 | **-0.1 %** |
| `refined_yp82_rect_redblack` | 4 | 0.06247 | 0.06252 | **-0.1 %** |
| `refined_yp82_rect_redblack` | 8 | 0.03639 | 0.03640 | **-0.0 %** |
| `refined_yp82_rect_redblack` | 16 | 0.02511 | 0.02529 | **-0.7 %** |
| `refined_big_rect_jacobi` | 4 | 0.26521 | 0.25564 | **+3.6 %** |
| `refined_big_rect_jacobi` | 8 | 0.13736 | 0.13256 | **+3.5 %** |
| `refined_big_rect_jacobi` | 16 | 0.07506 | 0.07269 | **+3.2 %** |

## Per-round `mpi_wait`

| config | n=1 | n=2 | n=4 | n=8 | n=16 |
|---|---|---|---|---|---|
| `base_jacobi` (ref) | - | 41.3 us | 98.4 us | 77.1 us | 96.1 us |
| `base_jacobi` (new) | - | 38.7 us | 97.8 us | 82.1 us | 92.2 us |
| `rect_jacobi` (ref) | - | 49.6 us | 65.6 us | 51.3 us | 85.0 us |
| `rect_jacobi` (new) | - | 46.8 us | 66.7 us | 53.4 us | 75.8 us |
| `refined_yp82_rect_jacobi` (ref) | - | 28.8 us | 42.0 us | 38.5 us | 53.2 us |
| `refined_yp82_rect_jacobi` (new) | - | 27.1 us | 39.5 us | 35.5 us | 55.4 us |
| `refined_yp82_rect_redblack` (ref) | - | 29.4 us | 38.0 us | 51.0 us | 101.4 us |
| `refined_yp82_rect_redblack` (new) | - | 29.5 us | 42.7 us | 52.5 us | 98.7 us |
| `refined_big_rect_jacobi` (ref) | - | - | 85.4 us | 57.3 us | 89.0 us |
| `refined_big_rect_jacobi` (new) | - | - | 77.5 us | 57.7 us | 86.5 us |

## Block tax — inverse per-GPU throughput `t_step * n / cells` (ns)

| ranks | Mcell/GPU | base ref | base new | rect ref | rect new | **tax ref** | **tax new** |
|---|---|---|---|---|---|---|---|
| 1 | 138.41 | 4.136 | 3.968 | 4.184 | 4.025 | **1.012** | **1.014** |
| 2 | 69.21 | 4.128 | 3.963 | 4.268 | 4.108 | **1.034** | **1.037** |
| 4 | 34.60 | 4.276 | 4.112 | 4.388 | 4.227 | **1.026** | **1.028** |
| 8 | 17.30 | 4.567 | 4.403 | 4.587 | 4.423 | **1.004** | **1.005** |
| 16 | 8.65 | 5.030 | 4.852 | 5.190 | 5.000 | **1.032** | **1.030** |

## Strong-scaling efficiency (vs each config's smallest rank count)

| config | n=1 | n=2 | n=4 | n=8 | n=16 |
|---|---|---|---|---|---|
| `base_jacobi` (ref) | 100 % | 100 % | 97 % | 91 % | 82 % |
| `base_jacobi` (new) | 100 % | 100 % | 96 % | 90 % | 82 % |
| `rect_jacobi` (ref) | 100 % | 98 % | 95 % | 91 % | 81 % |
| `rect_jacobi` (new) | 100 % | 98 % | 95 % | 91 % | 81 % |
| `refined_yp82_rect_jacobi` (ref) | 100 % | 96 % | 90 % | 80 % | 64 % |
| `refined_yp82_rect_jacobi` (new) | 100 % | 96 % | 90 % | 80 % | 63 % |
| `refined_yp82_rect_redblack` (ref) | 100 % | 95 % | 88 % | 75 % | 55 % |
| `refined_yp82_rect_redblack` (new) | 100 % | 95 % | 88 % | 75 % | 54 % |
| `refined_big_rect_jacobi` (ref) | - | - | 100 % | 97 % | 88 % |
| `refined_big_rect_jacobi` (new) | - | - | 100 % | 96 % | 88 % |

## Headline 1 -- the 2:1 machinery, like for like

`refined_big_rect_jacobi` shares grid, block shape and solver with
`rect_jacobi`; it just refines one y-tile. Cost of the ADDED cells
against the coarse cells beside them:

| ranks | binary | rect s/step | big s/step | ns per EXTRA cell | vs coarse-cell avg |
|---|---|---|---|---|---|
| 4 | ref | 0.15184 | 0.26521 | 1.0921 | **0.996x** |
| 4 | new | 0.14628 | 0.25564 | 1.0535 | **0.997x** |
| 8 | ref | 0.07937 | 0.13736 | 0.5586 | **0.974x** |
| 8 | new | 0.07653 | 0.13256 | 0.5397 | **0.976x** |
| 16 | ref | 0.04490 | 0.07506 | 0.2905 | **0.896x** |
| 16 | new | 0.04325 | 0.07269 | 0.2836 | **0.908x** |

## Headline 3 — red-black against Jacobi

| binary | n=1 | n=2 | n=4 | n=8 | n=16 |
|---|---|---|---|---|---|
| ref | 0.848 | 0.857 | 0.869 | 0.905 | 0.992 |
| new | 0.880 | 0.889 | 0.899 | 0.936 | 1.022 |
