# The campaign, re-measured

The 23-run matrix (`run_matrix.sh`, `--map-by numa --bind-to core`), run
TWICE IN ONE ALLOCATION -- `ref` is the pre-change binary, `new` the
post-change one, and the two columns differ only in that change. WHICH
change is recorded in the run directory's provenance.txt, not here: this
collector has now served two of them (the rank-to-GPU mapping, 2026-09-10,
and the comm_type parent map, 2026-09-14), and hard-coding either into the
header is how a table gets quoted against the wrong question.

```
job    : 5150030
date   : 2026-09-23 01:02:33 CEST
nodes  : 4  (hkn[0701,0709,0711,0733])
new    : b9414bd6a1d2b9e29178e5605a9ddd80a7db17b6  (index_list_b9414bd)
mid    : 95312d75f867312298955840e446d54546936312  (prefix_form_95312d7)
ref    : 3c2903a62c7e027fb7eb250925f33077ca73d88a  (before_divergence_work_3c2903a)
gpu    : NVIDIA A100-SXM4-40GB
nsteps : 200
layout : --map-by numa --bind-to core (as every earlier matrix)
```

## s/step, and what the mapping is worth

| config | ranks | ref | new | gain |
|---|---|---|---|---|
| `base_jacobi` | 1 | 0.57917 | 0.57723 | **+0.3 %** |
| `base_jacobi` | 2 | 0.28999 | 0.28777 | **+0.8 %** |
| `base_jacobi` | 4 | 0.15282 | 0.14895 | **+2.5 %** |
| `base_jacobi` | 8 | 0.08367 | 0.07993 | **+4.5 %** |
| `base_jacobi` | 16 | 0.04838 | 0.04524 | **+6.5 %** |
| `rect_jacobi` | 1 | 0.58606 | 0.58590 | **+0.0 %** |
| `rect_jacobi` | 2 | 0.30979 | 0.29769 | **+3.9 %** |
| `rect_jacobi` | 4 | 0.15926 | 0.15254 | **+4.2 %** |
| `rect_jacobi` | 8 | 0.08428 | 0.07990 | **+5.2 %** |
| `rect_jacobi` | 16 | 0.05268 | 0.04556 | **+13.5 %** |
| `refined_yp82_rect_jacobi` | 1 | 0.26269 | 0.26172 | **+0.4 %** |
| `refined_yp82_rect_jacobi` | 2 | 0.14084 | 0.13578 | **+3.6 %** |
| `refined_yp82_rect_jacobi` | 4 | 0.07513 | 0.07224 | **+3.8 %** |
| `refined_yp82_rect_jacobi` | 8 | 0.04418 | 0.04110 | **+7.0 %** |
| `refined_yp82_rect_jacobi` | 16 | 0.02997 | 0.02751 | **+8.2 %** |
| `refined_yp82_rect_redblack` | 1 | 0.22161 | 0.22172 | **-0.0 %** |
| `refined_yp82_rect_redblack` | 2 | 0.11640 | 0.11647 | **-0.1 %** |
| `refined_yp82_rect_redblack` | 4 | 0.06284 | 0.06284 | **-0.0 %** |
| `refined_yp82_rect_redblack` | 8 | 0.03795 | 0.03808 | **-0.4 %** |
| `refined_yp82_rect_redblack` | 16 | 0.02671 | 0.02657 | **+0.5 %** |
| `refined_big_rect_jacobi` | 4 | 0.27732 | 0.26708 | **+3.7 %** |
| `refined_big_rect_jacobi` | 8 | 0.14469 | 0.13884 | **+4.0 %** |
| `refined_big_rect_jacobi` | 16 | 0.08341 | 0.07561 | **+9.3 %** |

## Per-round `mpi_wait`

| config | n=1 | n=2 | n=4 | n=8 | n=16 |
|---|---|---|---|---|---|
| `base_jacobi` (ref) | - | 21.9 us | 102.9 us | 117.0 us | 180.7 us |
| `base_jacobi` (new) | - | 19.2 us | 79.2 us | 86.0 us | 130.7 us |
| `rect_jacobi` (ref) | - | 20.3 us | 55.2 us | 69.7 us | 213.3 us |
| `rect_jacobi` (new) | - | 21.0 us | 38.9 us | 47.1 us | 88.7 us |
| `refined_yp82_rect_jacobi` (ref) | - | 21.9 us | 39.7 us | 89.5 us | 134.7 us |
| `refined_yp82_rect_jacobi` (new) | - | 16.3 us | 31.7 us | 49.8 us | 103.3 us |
| `refined_yp82_rect_redblack` (ref) | - | 20.2 us | 39.3 us | 84.3 us | 138.1 us |
| `refined_yp82_rect_redblack` (new) | - | 21.6 us | 39.9 us | 87.7 us | 136.6 us |
| `refined_big_rect_jacobi` (ref) | - | - | 81.4 us | 70.7 us | 205.5 us |
| `refined_big_rect_jacobi` (new) | - | - | 54.8 us | 48.9 us | 84.7 us |

## Block tax — inverse per-GPU throughput `t_step * n / cells` (ns)

| ranks | Mcell/GPU | base ref | base new | rect ref | rect new | **tax ref** | **tax new** |
|---|---|---|---|---|---|---|---|
| 1 | 138.41 | 4.184 | 4.170 | 4.234 | 4.233 | **1.012** | **1.015** |
| 2 | 69.21 | 4.190 | 4.158 | 4.476 | 4.301 | **1.068** | **1.034** |
| 4 | 34.60 | 4.416 | 4.305 | 4.602 | 4.408 | **1.042** | **1.024** |
| 8 | 17.30 | 4.836 | 4.620 | 4.872 | 4.618 | **1.007** | **1.000** |
| 16 | 8.65 | 5.593 | 5.230 | 6.090 | 5.267 | **1.089** | **1.007** |

## Strong-scaling efficiency (vs each config's smallest rank count)

| config | n=1 | n=2 | n=4 | n=8 | n=16 |
|---|---|---|---|---|---|
| `base_jacobi` (ref) | 100 % | 100 % | 95 % | 87 % | 75 % |
| `base_jacobi` (new) | 100 % | 100 % | 97 % | 90 % | 80 % |
| `rect_jacobi` (ref) | 100 % | 95 % | 92 % | 87 % | 70 % |
| `rect_jacobi` (new) | 100 % | 98 % | 96 % | 92 % | 80 % |
| `refined_yp82_rect_jacobi` (ref) | 100 % | 93 % | 87 % | 74 % | 55 % |
| `refined_yp82_rect_jacobi` (new) | 100 % | 96 % | 91 % | 80 % | 59 % |
| `refined_yp82_rect_redblack` (ref) | 100 % | 95 % | 88 % | 73 % | 52 % |
| `refined_yp82_rect_redblack` (new) | 100 % | 95 % | 88 % | 73 % | 52 % |
| `refined_big_rect_jacobi` (ref) | - | - | 100 % | 96 % | 83 % |
| `refined_big_rect_jacobi` (new) | - | - | 100 % | 96 % | 88 % |

## Headline 1 -- the 2:1 machinery, like for like

`refined_big_rect_jacobi` shares grid, block shape and solver with
`rect_jacobi`; it just refines one y-tile. Cost of the ADDED cells
against the coarse cells beside them:

| ranks | binary | rect s/step | big s/step | ns per EXTRA cell | vs coarse-cell avg |
|---|---|---|---|---|---|
| 4 | ref | 0.15926 | 0.27732 | 1.1373 | **0.988x** |
| 4 | new | 0.15254 | 0.26708 | 1.1034 | **1.001x** |
| 8 | ref | 0.08428 | 0.14469 | 0.5819 | **0.956x** |
| 8 | new | 0.07990 | 0.13884 | 0.5678 | **0.984x** |
| 16 | ref | 0.05268 | 0.08341 | 0.2960 | **0.778x** |
| 16 | new | 0.04556 | 0.07561 | 0.2895 | **0.879x** |

## Headline 3 — red-black against Jacobi

| binary | n=1 | n=2 | n=4 | n=8 | n=16 |
|---|---|---|---|---|---|
| ref | 0.844 | 0.826 | 0.836 | 0.859 | 0.891 |
| new | 0.847 | 0.858 | 0.870 | 0.926 | 0.966 |
