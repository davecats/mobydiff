# The campaign, re-measured

The 23-run matrix (`run_matrix.sh`, `--map-by numa --bind-to core`), run
TWICE IN ONE ALLOCATION -- `ref` is the pre-change binary, `new` the
post-change one, and the two columns differ only in that change. WHICH
change is recorded in the run directory's provenance.txt, not here: this
collector has now served two of them (the rank-to-GPU mapping, 2026-09-10,
and the comm_type parent map, 2026-09-14), and hard-coding either into the
header is how a table gets quoted against the wrong question.

```
job    : 5147466
date   : 2026-09-17 11:08:38 CEST
nodes  : 4  (hkn[0532,0627,0720,0733])
new    : 95312d75f867312298955840e446d54546936312  (+ the divergence-halo exchange)
mid    : 3c2903a62c7e027fb7eb250925f33077ca73d88a  (+ the register cuts and step work)
ref    : 55bee89aec5d8a29321a11a59b360915b97ed1d7  (map(to: c); = the 2026-09-14 'new' column)
gpu    : NVIDIA A100-SXM4-40GB
nsteps : 200
layout : --map-by numa --bind-to core (as every earlier matrix)
```

## s/step, and what the mapping is worth

| config | ranks | ref | new | gain |
|---|---|---|---|---|
| `base_jacobi` | 1 | 0.69433 | 0.56942 | **+18.0 %** |
| `base_jacobi` | 2 | 0.34785 | 0.28420 | **+18.3 %** |
| `base_jacobi` | 4 | 0.18119 | 0.14706 | **+18.8 %** |
| `base_jacobi` | 8 | 0.09785 | 0.07886 | **+19.4 %** |
| `base_jacobi` | 16 | 0.05496 | 0.04327 | **+21.3 %** |
| `rect_jacobi` | 1 | 0.72837 | 0.58972 | **+19.0 %** |
| `rect_jacobi` | 2 | 0.37947 | 0.29917 | **+21.2 %** |
| `rect_jacobi` | 4 | 0.19378 | 0.15263 | **+21.2 %** |
| `rect_jacobi` | 8 | 0.10192 | 0.07966 | **+21.8 %** |
| `rect_jacobi` | 16 | 0.05979 | 0.04470 | **+25.2 %** |
| `refined_yp82_rect_jacobi` | 1 | 0.32073 | 0.26217 | **+18.3 %** |
| `refined_yp82_rect_jacobi` | 2 | 0.17116 | 0.13612 | **+20.5 %** |
| `refined_yp82_rect_jacobi` | 4 | 0.09092 | 0.07250 | **+20.3 %** |
| `refined_yp82_rect_jacobi` | 8 | 0.05132 | 0.04021 | **+21.7 %** |
| `refined_yp82_rect_jacobi` | 16 | 0.03315 | 0.02539 | **+23.4 %** |
| `refined_yp82_rect_redblack` | 1 | 0.24341 | 0.22810 | **+6.3 %** |
| `refined_yp82_rect_redblack` | 2 | 0.12802 | 0.11957 | **+6.6 %** |
| `refined_yp82_rect_redblack` | 4 | 0.06923 | 0.06453 | **+6.8 %** |
| `refined_yp82_rect_redblack` | 8 | 0.04013 | 0.03700 | **+7.8 %** |
| `refined_yp82_rect_redblack` | 16 | 0.02665 | 0.02522 | **+5.4 %** |
| `refined_big_rect_jacobi` | 4 | 0.33910 | 0.26738 | **+21.2 %** |
| `refined_big_rect_jacobi` | 8 | 0.17733 | 0.13854 | **+21.9 %** |
| `refined_big_rect_jacobi` | 16 | 0.09847 | 0.07487 | **+24.0 %** |

## Per-round `mpi_wait`

| config | n=1 | n=2 | n=4 | n=8 | n=16 |
|---|---|---|---|---|---|
| `base_jacobi` (ref) | - | 37.1 us | 132.0 us | 122.4 us | 148.7 us |
| `base_jacobi` (new) | - | 30.0 us | 96.9 us | 88.3 us | 97.4 us |
| `rect_jacobi` (ref) | - | 37.5 us | 66.5 us | 72.8 us | 145.8 us |
| `rect_jacobi` (new) | - | 28.0 us | 55.3 us | 52.5 us | 80.0 us |
| `refined_yp82_rect_jacobi` (ref) | - | 34.0 us | 45.3 us | 44.0 us | 85.9 us |
| `refined_yp82_rect_jacobi` (new) | - | 25.1 us | 42.7 us | 36.9 us | 54.9 us |
| `refined_yp82_rect_redblack` (ref) | - | 22.5 us | 36.8 us | 42.7 us | 77.2 us |
| `refined_yp82_rect_redblack` (new) | - | 26.0 us | 41.1 us | 50.2 us | 97.4 us |
| `refined_big_rect_jacobi` (ref) | - | - | 92.9 us | 82.2 us | 115.5 us |
| `refined_big_rect_jacobi` (new) | - | - | 87.5 us | 65.0 us | 71.4 us |

## Block tax — inverse per-GPU throughput `t_step * n / cells` (ns)

| ranks | Mcell/GPU | base ref | base new | rect ref | rect new | **tax ref** | **tax new** |
|---|---|---|---|---|---|---|---|
| 1 | 138.41 | 5.016 | 4.114 | 5.262 | 4.261 | **1.049** | **1.036** |
| 2 | 69.21 | 5.026 | 4.107 | 5.483 | 4.323 | **1.091** | **1.053** |
| 4 | 34.60 | 5.236 | 4.250 | 5.600 | 4.411 | **1.069** | **1.038** |
| 8 | 17.30 | 5.656 | 4.558 | 5.891 | 4.604 | **1.042** | **1.010** |
| 16 | 8.65 | 6.354 | 5.001 | 6.911 | 5.167 | **1.088** | **1.033** |

## Strong-scaling efficiency (vs each config's smallest rank count)

| config | n=1 | n=2 | n=4 | n=8 | n=16 |
|---|---|---|---|---|---|
| `base_jacobi` (ref) | 100 % | 100 % | 96 % | 89 % | 79 % |
| `base_jacobi` (new) | 100 % | 100 % | 97 % | 90 % | 82 % |
| `rect_jacobi` (ref) | 100 % | 96 % | 94 % | 89 % | 76 % |
| `rect_jacobi` (new) | 100 % | 99 % | 97 % | 93 % | 82 % |
| `refined_yp82_rect_jacobi` (ref) | 100 % | 94 % | 88 % | 78 % | 60 % |
| `refined_yp82_rect_jacobi` (new) | 100 % | 96 % | 90 % | 82 % | 65 % |
| `refined_yp82_rect_redblack` (ref) | 100 % | 95 % | 88 % | 76 % | 57 % |
| `refined_yp82_rect_redblack` (new) | 100 % | 95 % | 88 % | 77 % | 57 % |
| `refined_big_rect_jacobi` (ref) | - | - | 100 % | 96 % | 86 % |
| `refined_big_rect_jacobi` (new) | - | - | 100 % | 96 % | 89 % |

## Headline 1 -- the 2:1 machinery, like for like

`refined_big_rect_jacobi` shares grid, block shape and solver with
`rect_jacobi`; it just refines one y-tile. Cost of the ADDED cells
against the coarse cells beside them:

| ranks | binary | rect s/step | big s/step | ns per EXTRA cell | vs coarse-cell avg |
|---|---|---|---|---|---|
| 4 | ref | 0.19378 | 0.33910 | 1.3999 | **1.000x** |
| 4 | new | 0.15263 | 0.26738 | 1.1054 | **1.002x** |
| 8 | ref | 0.10192 | 0.17733 | 0.7264 | **0.986x** |
| 8 | new | 0.07966 | 0.13854 | 0.5672 | **0.986x** |
| 16 | ref | 0.05979 | 0.09847 | 0.3726 | **0.863x** |
| 16 | new | 0.04470 | 0.07487 | 0.2907 | **0.900x** |

## Headline 3 — red-black against Jacobi

| binary | n=1 | n=2 | n=4 | n=8 | n=16 |
|---|---|---|---|---|---|
| ref | 0.759 | 0.748 | 0.761 | 0.782 | 0.804 |
| new | 0.870 | 0.878 | 0.890 | 0.920 | 0.993 |
