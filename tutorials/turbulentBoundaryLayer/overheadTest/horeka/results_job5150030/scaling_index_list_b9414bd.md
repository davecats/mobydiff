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
| `base_jacobi` | 1 | 0.57886 | 0.57723 | **+0.3 %** |
| `base_jacobi` | 2 | 0.28794 | 0.28777 | **+0.1 %** |
| `base_jacobi` | 4 | 0.14908 | 0.14895 | **+0.1 %** |
| `base_jacobi` | 8 | 0.08013 | 0.07993 | **+0.3 %** |
| `base_jacobi` | 16 | 0.04438 | 0.04524 | **-1.9 %** |
| `rect_jacobi` | 1 | 0.59131 | 0.58590 | **+0.9 %** |
| `rect_jacobi` | 2 | 0.30036 | 0.29769 | **+0.9 %** |
| `rect_jacobi` | 4 | 0.15397 | 0.15254 | **+0.9 %** |
| `rect_jacobi` | 8 | 0.08057 | 0.07990 | **+0.8 %** |
| `rect_jacobi` | 16 | 0.04559 | 0.04556 | **+0.1 %** |
| `refined_yp82_rect_jacobi` | 1 | 0.26529 | 0.26172 | **+1.3 %** |
| `refined_yp82_rect_jacobi` | 2 | 0.13758 | 0.13578 | **+1.3 %** |
| `refined_yp82_rect_jacobi` | 4 | 0.07295 | 0.07224 | **+1.0 %** |
| `refined_yp82_rect_jacobi` | 8 | 0.04154 | 0.04110 | **+1.1 %** |
| `refined_yp82_rect_jacobi` | 16 | 0.02718 | 0.02751 | **-1.2 %** |
| `refined_yp82_rect_redblack` | 1 | 0.22958 | 0.22172 | **+3.4 %** |
| `refined_yp82_rect_redblack` | 2 | 0.12024 | 0.11647 | **+3.1 %** |
| `refined_yp82_rect_redblack` | 4 | 0.06476 | 0.06284 | **+3.0 %** |
| `refined_yp82_rect_redblack` | 8 | 0.03874 | 0.03808 | **+1.7 %** |
| `refined_yp82_rect_redblack` | 16 | 0.02611 | 0.02657 | **-1.8 %** |
| `refined_big_rect_jacobi` | 4 | 0.27012 | 0.26708 | **+1.1 %** |
| `refined_big_rect_jacobi` | 8 | 0.14024 | 0.13884 | **+1.0 %** |
| `refined_big_rect_jacobi` | 16 | 0.07606 | 0.07561 | **+0.6 %** |

## Per-round `mpi_wait`

| config | n=1 | n=2 | n=4 | n=8 | n=16 |
|---|---|---|---|---|---|
| `base_jacobi` (ref) | - | 19.4 us | 80.3 us | 84.4 us | 111.0 us |
| `base_jacobi` (new) | - | 19.2 us | 79.2 us | 86.0 us | 130.7 us |
| `rect_jacobi` (ref) | - | 20.5 us | 42.0 us | 45.4 us | 91.9 us |
| `rect_jacobi` (new) | - | 21.0 us | 38.9 us | 47.1 us | 88.7 us |
| `refined_yp82_rect_jacobi` (ref) | - | 16.0 us | 37.2 us | 57.3 us | 95.2 us |
| `refined_yp82_rect_jacobi` (new) | - | 16.3 us | 31.7 us | 49.8 us | 103.3 us |
| `refined_yp82_rect_redblack` (ref) | - | 21.6 us | 40.2 us | 88.7 us | 122.8 us |
| `refined_yp82_rect_redblack` (new) | - | 21.6 us | 39.9 us | 87.7 us | 136.6 us |
| `refined_big_rect_jacobi` (ref) | - | - | 63.1 us | 49.9 us | 78.1 us |
| `refined_big_rect_jacobi` (new) | - | - | 54.8 us | 48.9 us | 84.7 us |

## Block tax — inverse per-GPU throughput `t_step * n / cells` (ns)

| ranks | Mcell/GPU | base ref | base new | rect ref | rect new | **tax ref** | **tax new** |
|---|---|---|---|---|---|---|---|
| 1 | 138.41 | 4.182 | 4.170 | 4.272 | 4.233 | **1.021** | **1.015** |
| 2 | 69.21 | 4.161 | 4.158 | 4.340 | 4.301 | **1.043** | **1.034** |
| 4 | 34.60 | 4.308 | 4.305 | 4.450 | 4.408 | **1.033** | **1.024** |
| 8 | 17.30 | 4.631 | 4.620 | 4.657 | 4.618 | **1.006** | **1.000** |
| 16 | 8.65 | 5.130 | 5.230 | 5.270 | 5.267 | **1.027** | **1.007** |

## Strong-scaling efficiency (vs each config's smallest rank count)

| config | n=1 | n=2 | n=4 | n=8 | n=16 |
|---|---|---|---|---|---|
| `base_jacobi` (ref) | 100 % | 101 % | 97 % | 90 % | 82 % |
| `base_jacobi` (new) | 100 % | 100 % | 97 % | 90 % | 80 % |
| `rect_jacobi` (ref) | 100 % | 98 % | 96 % | 92 % | 81 % |
| `rect_jacobi` (new) | 100 % | 98 % | 96 % | 92 % | 80 % |
| `refined_yp82_rect_jacobi` (ref) | 100 % | 96 % | 91 % | 80 % | 61 % |
| `refined_yp82_rect_jacobi` (new) | 100 % | 96 % | 91 % | 80 % | 59 % |
| `refined_yp82_rect_redblack` (ref) | 100 % | 95 % | 89 % | 74 % | 55 % |
| `refined_yp82_rect_redblack` (new) | 100 % | 95 % | 88 % | 73 % | 52 % |
| `refined_big_rect_jacobi` (ref) | - | - | 100 % | 96 % | 89 % |
| `refined_big_rect_jacobi` (new) | - | - | 100 % | 96 % | 88 % |

## Headline 1 -- the 2:1 machinery, like for like

`refined_big_rect_jacobi` shares grid, block shape and solver with
`rect_jacobi`; it just refines one y-tile. Cost of the ADDED cells
against the coarse cells beside them:

| ranks | binary | rect s/step | big s/step | ns per EXTRA cell | vs coarse-cell avg |
|---|---|---|---|---|---|
| 4 | ref | 0.15397 | 0.27012 | 1.1189 | **1.006x** |
| 4 | new | 0.15254 | 0.26708 | 1.1034 | **1.001x** |
| 8 | ref | 0.08057 | 0.14024 | 0.5748 | **0.987x** |
| 8 | new | 0.07990 | 0.13884 | 0.5678 | **0.984x** |
| 16 | ref | 0.04559 | 0.07606 | 0.2935 | **0.891x** |
| 16 | new | 0.04556 | 0.07561 | 0.2895 | **0.879x** |

## Headline 3 — red-black against Jacobi

| binary | n=1 | n=2 | n=4 | n=8 | n=16 |
|---|---|---|---|---|---|
| ref | 0.865 | 0.874 | 0.888 | 0.932 | 0.961 |
| new | 0.847 | 0.858 | 0.870 | 0.926 | 0.966 |
