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
| `base_jacobi` | 1 | 0.56940 | 0.54925 | **+3.5 %** |
| `base_jacobi` | 2 | 0.28403 | 0.27427 | **+3.4 %** |
| `base_jacobi` | 4 | 0.14713 | 0.14230 | **+3.3 %** |
| `base_jacobi` | 8 | 0.07858 | 0.07617 | **+3.1 %** |
| `base_jacobi` | 16 | 0.04326 | 0.04198 | **+3.0 %** |
| `rect_jacobi` | 1 | 0.57641 | 0.55716 | **+3.3 %** |
| `rect_jacobi` | 2 | 0.29409 | 0.28430 | **+3.3 %** |
| `rect_jacobi` | 4 | 0.15117 | 0.14628 | **+3.2 %** |
| `rect_jacobi` | 8 | 0.07896 | 0.07653 | **+3.1 %** |
| `rect_jacobi` | 16 | 0.04465 | 0.04325 | **+3.1 %** |
| `refined_yp82_rect_jacobi` | 1 | 0.25781 | 0.24941 | **+3.3 %** |
| `refined_yp82_rect_jacobi` | 2 | 0.13428 | 0.13017 | **+3.1 %** |
| `refined_yp82_rect_jacobi` | 4 | 0.07166 | 0.06951 | **+3.0 %** |
| `refined_yp82_rect_jacobi` | 8 | 0.03997 | 0.03890 | **+2.7 %** |
| `refined_yp82_rect_jacobi` | 16 | 0.02538 | 0.02474 | **+2.5 %** |
| `refined_yp82_rect_redblack` | 1 | 0.21953 | 0.21947 | **+0.0 %** |
| `refined_yp82_rect_redblack` | 2 | 0.11567 | 0.11570 | **-0.0 %** |
| `refined_yp82_rect_redblack` | 4 | 0.06255 | 0.06252 | **+0.1 %** |
| `refined_yp82_rect_redblack` | 8 | 0.03639 | 0.03640 | **-0.0 %** |
| `refined_yp82_rect_redblack` | 16 | 0.02544 | 0.02529 | **+0.6 %** |
| `refined_big_rect_jacobi` | 4 | 0.26400 | 0.25564 | **+3.2 %** |
| `refined_big_rect_jacobi` | 8 | 0.13672 | 0.13256 | **+3.0 %** |
| `refined_big_rect_jacobi` | 16 | 0.07489 | 0.07269 | **+2.9 %** |

## Per-round `mpi_wait`

| config | n=1 | n=2 | n=4 | n=8 | n=16 |
|---|---|---|---|---|---|
| `base_jacobi` (ref) | - | 38.1 us | 100.0 us | 78.4 us | 95.2 us |
| `base_jacobi` (new) | - | 38.7 us | 97.8 us | 82.1 us | 92.2 us |
| `rect_jacobi` (ref) | - | 51.3 us | 63.6 us | 50.6 us | 82.9 us |
| `rect_jacobi` (new) | - | 46.8 us | 66.7 us | 53.4 us | 75.8 us |
| `refined_yp82_rect_jacobi` (ref) | - | 28.2 us | 47.2 us | 33.8 us | 54.6 us |
| `refined_yp82_rect_jacobi` (new) | - | 27.1 us | 39.5 us | 35.5 us | 55.4 us |
| `refined_yp82_rect_redblack` (ref) | - | 28.7 us | 42.3 us | 51.6 us | 106.5 us |
| `refined_yp82_rect_redblack` (new) | - | 29.5 us | 42.7 us | 52.5 us | 98.7 us |
| `refined_big_rect_jacobi` (ref) | - | - | 86.7 us | 57.4 us | 91.6 us |
| `refined_big_rect_jacobi` (new) | - | - | 77.5 us | 57.7 us | 86.5 us |

## Block tax — inverse per-GPU throughput `t_step * n / cells` (ns)

| ranks | Mcell/GPU | base ref | base new | rect ref | rect new | **tax ref** | **tax new** |
|---|---|---|---|---|---|---|---|
| 1 | 138.41 | 4.114 | 3.968 | 4.164 | 4.025 | **1.012** | **1.014** |
| 2 | 69.21 | 4.104 | 3.963 | 4.250 | 4.108 | **1.035** | **1.037** |
| 4 | 34.60 | 4.252 | 4.112 | 4.369 | 4.227 | **1.028** | **1.028** |
| 8 | 17.30 | 4.542 | 4.403 | 4.564 | 4.423 | **1.005** | **1.005** |
| 16 | 8.65 | 5.000 | 4.852 | 5.162 | 5.000 | **1.032** | **1.030** |

## Strong-scaling efficiency (vs each config's smallest rank count)

| config | n=1 | n=2 | n=4 | n=8 | n=16 |
|---|---|---|---|---|---|
| `base_jacobi` (ref) | 100 % | 100 % | 97 % | 91 % | 82 % |
| `base_jacobi` (new) | 100 % | 100 % | 96 % | 90 % | 82 % |
| `rect_jacobi` (ref) | 100 % | 98 % | 95 % | 91 % | 81 % |
| `rect_jacobi` (new) | 100 % | 98 % | 95 % | 91 % | 81 % |
| `refined_yp82_rect_jacobi` (ref) | 100 % | 96 % | 90 % | 81 % | 63 % |
| `refined_yp82_rect_jacobi` (new) | 100 % | 96 % | 90 % | 80 % | 63 % |
| `refined_yp82_rect_redblack` (ref) | 100 % | 95 % | 88 % | 75 % | 54 % |
| `refined_yp82_rect_redblack` (new) | 100 % | 95 % | 88 % | 75 % | 54 % |
| `refined_big_rect_jacobi` (ref) | - | - | 100 % | 97 % | 88 % |
| `refined_big_rect_jacobi` (new) | - | - | 100 % | 96 % | 88 % |

## Headline 1 -- the 2:1 machinery, like for like

`refined_big_rect_jacobi` shares grid, block shape and solver with
`rect_jacobi`; it just refines one y-tile. Cost of the ADDED cells
against the coarse cells beside them:

| ranks | binary | rect s/step | big s/step | ns per EXTRA cell | vs coarse-cell avg |
|---|---|---|---|---|---|
| 4 | ref | 0.15117 | 0.26400 | 1.0869 | **0.995x** |
| 4 | new | 0.14628 | 0.25564 | 1.0535 | **0.997x** |
| 8 | ref | 0.07896 | 0.13672 | 0.5564 | **0.975x** |
| 8 | new | 0.07653 | 0.13256 | 0.5397 | **0.976x** |
| 16 | ref | 0.04465 | 0.07489 | 0.2912 | **0.903x** |
| 16 | new | 0.04325 | 0.07269 | 0.2836 | **0.908x** |

## Headline 3 — red-black against Jacobi

| binary | n=1 | n=2 | n=4 | n=8 | n=16 |
|---|---|---|---|---|---|
| ref | 0.852 | 0.861 | 0.873 | 0.910 | 1.003 |
| new | 0.880 | 0.889 | 0.899 | 0.936 | 1.022 |
