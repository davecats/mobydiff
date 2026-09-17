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
| `base_jacobi` | 1 | 0.69433 | 0.58017 | **+16.4 %** |
| `base_jacobi` | 2 | 0.34785 | 0.28835 | **+17.1 %** |
| `base_jacobi` | 4 | 0.18119 | 0.15108 | **+16.6 %** |
| `base_jacobi` | 8 | 0.09785 | 0.08254 | **+15.6 %** |
| `base_jacobi` | 16 | 0.05496 | 0.04736 | **+13.8 %** |
| `rect_jacobi` | 1 | 0.72837 | 0.58410 | **+19.8 %** |
| `rect_jacobi` | 2 | 0.37947 | 0.30698 | **+19.1 %** |
| `rect_jacobi` | 4 | 0.19378 | 0.15785 | **+18.5 %** |
| `rect_jacobi` | 8 | 0.10192 | 0.08342 | **+18.2 %** |
| `rect_jacobi` | 16 | 0.05979 | 0.05065 | **+15.3 %** |
| `refined_yp82_rect_jacobi` | 1 | 0.32073 | 0.25868 | **+19.3 %** |
| `refined_yp82_rect_jacobi` | 2 | 0.17116 | 0.13973 | **+18.4 %** |
| `refined_yp82_rect_jacobi` | 4 | 0.09092 | 0.07453 | **+18.0 %** |
| `refined_yp82_rect_jacobi` | 8 | 0.05132 | 0.04238 | **+17.4 %** |
| `refined_yp82_rect_jacobi` | 16 | 0.03315 | 0.02887 | **+12.9 %** |
| `refined_yp82_rect_redblack` | 1 | 0.24341 | 0.21970 | **+9.7 %** |
| `refined_yp82_rect_redblack` | 2 | 0.12802 | 0.11558 | **+9.7 %** |
| `refined_yp82_rect_redblack` | 4 | 0.06923 | 0.06252 | **+9.7 %** |
| `refined_yp82_rect_redblack` | 8 | 0.04013 | 0.03642 | **+9.2 %** |
| `refined_yp82_rect_redblack` | 16 | 0.02665 | 0.02503 | **+6.1 %** |
| `refined_big_rect_jacobi` | 4 | 0.33910 | 0.27459 | **+19.0 %** |
| `refined_big_rect_jacobi` | 8 | 0.17733 | 0.14316 | **+19.3 %** |
| `refined_big_rect_jacobi` | 16 | 0.09847 | 0.08248 | **+16.2 %** |

## Per-round `mpi_wait`

| config | n=1 | n=2 | n=4 | n=8 | n=16 |
|---|---|---|---|---|---|
| `base_jacobi` (ref) | - | 37.1 us | 132.0 us | 122.4 us | 148.7 us |
| `base_jacobi` (new) | - | 27.2 us | 121.8 us | 116.2 us | 168.6 us |
| `rect_jacobi` (ref) | - | 37.5 us | 66.5 us | 72.8 us | 145.8 us |
| `rect_jacobi` (new) | - | 44.5 us | 76.0 us | 78.2 us | 170.6 us |
| `refined_yp82_rect_jacobi` (ref) | - | 34.0 us | 45.3 us | 44.0 us | 85.9 us |
| `refined_yp82_rect_jacobi` (new) | - | 39.8 us | 42.9 us | 49.9 us | 105.2 us |
| `refined_yp82_rect_redblack` (ref) | - | 22.5 us | 36.8 us | 42.7 us | 77.2 us |
| `refined_yp82_rect_redblack` (new) | - | 24.7 us | 39.5 us | 46.8 us | 96.2 us |
| `refined_big_rect_jacobi` (ref) | - | - | 92.9 us | 82.2 us | 115.5 us |
| `refined_big_rect_jacobi` (new) | - | - | 118.5 us | 92.0 us | 199.2 us |

## Block tax — inverse per-GPU throughput `t_step * n / cells` (ns)

| ranks | Mcell/GPU | base ref | base new | rect ref | rect new | **tax ref** | **tax new** |
|---|---|---|---|---|---|---|---|
| 1 | 138.41 | 5.016 | 4.192 | 5.262 | 4.220 | **1.049** | **1.007** |
| 2 | 69.21 | 5.026 | 4.167 | 5.483 | 4.436 | **1.091** | **1.065** |
| 4 | 34.60 | 5.236 | 4.366 | 5.600 | 4.562 | **1.069** | **1.045** |
| 8 | 17.30 | 5.656 | 4.771 | 5.891 | 4.821 | **1.042** | **1.011** |
| 16 | 8.65 | 6.354 | 5.475 | 6.911 | 5.855 | **1.088** | **1.069** |

## Strong-scaling efficiency (vs each config's smallest rank count)

| config | n=1 | n=2 | n=4 | n=8 | n=16 |
|---|---|---|---|---|---|
| `base_jacobi` (ref) | 100 % | 100 % | 96 % | 89 % | 79 % |
| `base_jacobi` (new) | 100 % | 101 % | 96 % | 88 % | 77 % |
| `rect_jacobi` (ref) | 100 % | 96 % | 94 % | 89 % | 76 % |
| `rect_jacobi` (new) | 100 % | 95 % | 93 % | 88 % | 72 % |
| `refined_yp82_rect_jacobi` (ref) | 100 % | 94 % | 88 % | 78 % | 60 % |
| `refined_yp82_rect_jacobi` (new) | 100 % | 93 % | 87 % | 76 % | 56 % |
| `refined_yp82_rect_redblack` (ref) | 100 % | 95 % | 88 % | 76 % | 57 % |
| `refined_yp82_rect_redblack` (new) | 100 % | 95 % | 88 % | 75 % | 55 % |
| `refined_big_rect_jacobi` (ref) | - | - | 100 % | 96 % | 86 % |
| `refined_big_rect_jacobi` (new) | - | - | 100 % | 96 % | 83 % |

## Headline 1 -- the 2:1 machinery, like for like

`refined_big_rect_jacobi` shares grid, block shape and solver with
`rect_jacobi`; it just refines one y-tile. Cost of the ADDED cells
against the coarse cells beside them:

| ranks | binary | rect s/step | big s/step | ns per EXTRA cell | vs coarse-cell avg |
|---|---|---|---|---|---|
| 4 | ref | 0.19378 | 0.33910 | 1.3999 | **1.000x** |
| 4 | new | 0.15785 | 0.27459 | 1.1246 | **0.986x** |
| 8 | ref | 0.10192 | 0.17733 | 0.7264 | **0.986x** |
| 8 | new | 0.08342 | 0.14316 | 0.5755 | **0.955x** |
| 16 | ref | 0.05979 | 0.09847 | 0.3726 | **0.863x** |
| 16 | new | 0.05065 | 0.08248 | 0.3066 | **0.838x** |

## Headline 3 — red-black against Jacobi

| binary | n=1 | n=2 | n=4 | n=8 | n=16 |
|---|---|---|---|---|---|
| ref | 0.759 | 0.748 | 0.761 | 0.782 | 0.804 |
| new | 0.849 | 0.827 | 0.839 | 0.859 | 0.867 |
