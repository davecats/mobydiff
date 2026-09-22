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
| `base_jacobi` | 1 | 0.57917 | 0.57886 | **+0.1 %** |
| `base_jacobi` | 2 | 0.28999 | 0.28794 | **+0.7 %** |
| `base_jacobi` | 4 | 0.15282 | 0.14908 | **+2.4 %** |
| `base_jacobi` | 8 | 0.08367 | 0.08013 | **+4.2 %** |
| `base_jacobi` | 16 | 0.04838 | 0.04438 | **+8.3 %** |
| `rect_jacobi` | 1 | 0.58606 | 0.59131 | **-0.9 %** |
| `rect_jacobi` | 2 | 0.30979 | 0.30036 | **+3.0 %** |
| `rect_jacobi` | 4 | 0.15926 | 0.15397 | **+3.3 %** |
| `rect_jacobi` | 8 | 0.08428 | 0.08057 | **+4.4 %** |
| `rect_jacobi` | 16 | 0.05268 | 0.04559 | **+13.5 %** |
| `refined_yp82_rect_jacobi` | 1 | 0.26269 | 0.26529 | **-1.0 %** |
| `refined_yp82_rect_jacobi` | 2 | 0.14084 | 0.13758 | **+2.3 %** |
| `refined_yp82_rect_jacobi` | 4 | 0.07513 | 0.07295 | **+2.9 %** |
| `refined_yp82_rect_jacobi` | 8 | 0.04418 | 0.04154 | **+6.0 %** |
| `refined_yp82_rect_jacobi` | 16 | 0.02997 | 0.02718 | **+9.3 %** |
| `refined_yp82_rect_redblack` | 1 | 0.22161 | 0.22958 | **-3.6 %** |
| `refined_yp82_rect_redblack` | 2 | 0.11640 | 0.12024 | **-3.3 %** |
| `refined_yp82_rect_redblack` | 4 | 0.06284 | 0.06476 | **-3.1 %** |
| `refined_yp82_rect_redblack` | 8 | 0.03795 | 0.03874 | **-2.1 %** |
| `refined_yp82_rect_redblack` | 16 | 0.02671 | 0.02611 | **+2.3 %** |
| `refined_big_rect_jacobi` | 4 | 0.27732 | 0.27012 | **+2.6 %** |
| `refined_big_rect_jacobi` | 8 | 0.14469 | 0.14024 | **+3.1 %** |
| `refined_big_rect_jacobi` | 16 | 0.08341 | 0.07606 | **+8.8 %** |

## Per-round `mpi_wait`

| config | n=1 | n=2 | n=4 | n=8 | n=16 |
|---|---|---|---|---|---|
| `base_jacobi` (ref) | - | 21.9 us | 102.9 us | 117.0 us | 180.7 us |
| `base_jacobi` (new) | - | 19.4 us | 80.3 us | 84.4 us | 111.0 us |
| `rect_jacobi` (ref) | - | 20.3 us | 55.2 us | 69.7 us | 213.3 us |
| `rect_jacobi` (new) | - | 20.5 us | 42.0 us | 45.4 us | 91.9 us |
| `refined_yp82_rect_jacobi` (ref) | - | 21.9 us | 39.7 us | 89.5 us | 134.7 us |
| `refined_yp82_rect_jacobi` (new) | - | 16.0 us | 37.2 us | 57.3 us | 95.2 us |
| `refined_yp82_rect_redblack` (ref) | - | 20.2 us | 39.3 us | 84.3 us | 138.1 us |
| `refined_yp82_rect_redblack` (new) | - | 21.6 us | 40.2 us | 88.7 us | 122.8 us |
| `refined_big_rect_jacobi` (ref) | - | - | 81.4 us | 70.7 us | 205.5 us |
| `refined_big_rect_jacobi` (new) | - | - | 63.1 us | 49.9 us | 78.1 us |

## Block tax — inverse per-GPU throughput `t_step * n / cells` (ns)

| ranks | Mcell/GPU | base ref | base new | rect ref | rect new | **tax ref** | **tax new** |
|---|---|---|---|---|---|---|---|
| 1 | 138.41 | 4.184 | 4.182 | 4.234 | 4.272 | **1.012** | **1.021** |
| 2 | 69.21 | 4.190 | 4.161 | 4.476 | 4.340 | **1.068** | **1.043** |
| 4 | 34.60 | 4.416 | 4.308 | 4.602 | 4.450 | **1.042** | **1.033** |
| 8 | 17.30 | 4.836 | 4.631 | 4.872 | 4.657 | **1.007** | **1.006** |
| 16 | 8.65 | 5.593 | 5.130 | 6.090 | 5.270 | **1.089** | **1.027** |

## Strong-scaling efficiency (vs each config's smallest rank count)

| config | n=1 | n=2 | n=4 | n=8 | n=16 |
|---|---|---|---|---|---|
| `base_jacobi` (ref) | 100 % | 100 % | 95 % | 87 % | 75 % |
| `base_jacobi` (new) | 100 % | 101 % | 97 % | 90 % | 82 % |
| `rect_jacobi` (ref) | 100 % | 95 % | 92 % | 87 % | 70 % |
| `rect_jacobi` (new) | 100 % | 98 % | 96 % | 92 % | 81 % |
| `refined_yp82_rect_jacobi` (ref) | 100 % | 93 % | 87 % | 74 % | 55 % |
| `refined_yp82_rect_jacobi` (new) | 100 % | 96 % | 91 % | 80 % | 61 % |
| `refined_yp82_rect_redblack` (ref) | 100 % | 95 % | 88 % | 73 % | 52 % |
| `refined_yp82_rect_redblack` (new) | 100 % | 95 % | 89 % | 74 % | 55 % |
| `refined_big_rect_jacobi` (ref) | - | - | 100 % | 96 % | 83 % |
| `refined_big_rect_jacobi` (new) | - | - | 100 % | 96 % | 89 % |

## Headline 1 -- the 2:1 machinery, like for like

`refined_big_rect_jacobi` shares grid, block shape and solver with
`rect_jacobi`; it just refines one y-tile. Cost of the ADDED cells
against the coarse cells beside them:

| ranks | binary | rect s/step | big s/step | ns per EXTRA cell | vs coarse-cell avg |
|---|---|---|---|---|---|
| 4 | ref | 0.15926 | 0.27732 | 1.1373 | **0.988x** |
| 4 | new | 0.15397 | 0.27012 | 1.1189 | **1.006x** |
| 8 | ref | 0.08428 | 0.14469 | 0.5819 | **0.956x** |
| 8 | new | 0.08057 | 0.14024 | 0.5748 | **0.987x** |
| 16 | ref | 0.05268 | 0.08341 | 0.2960 | **0.778x** |
| 16 | new | 0.04559 | 0.07606 | 0.2935 | **0.891x** |

## Headline 3 — red-black against Jacobi

| binary | n=1 | n=2 | n=4 | n=8 | n=16 |
|---|---|---|---|---|---|
| ref | 0.844 | 0.826 | 0.836 | 0.859 | 0.891 |
| new | 0.865 | 0.874 | 0.888 | 0.932 | 0.961 |
