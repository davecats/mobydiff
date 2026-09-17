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
| `base_jacobi` | 1 | 0.58017 | 0.56942 | **+1.9 %** |
| `base_jacobi` | 2 | 0.28835 | 0.28420 | **+1.4 %** |
| `base_jacobi` | 4 | 0.15108 | 0.14706 | **+2.7 %** |
| `base_jacobi` | 8 | 0.08254 | 0.07886 | **+4.5 %** |
| `base_jacobi` | 16 | 0.04736 | 0.04327 | **+8.6 %** |
| `rect_jacobi` | 1 | 0.58410 | 0.58972 | **-1.0 %** |
| `rect_jacobi` | 2 | 0.30698 | 0.29917 | **+2.5 %** |
| `rect_jacobi` | 4 | 0.15785 | 0.15263 | **+3.3 %** |
| `rect_jacobi` | 8 | 0.08342 | 0.07966 | **+4.5 %** |
| `rect_jacobi` | 16 | 0.05065 | 0.04470 | **+11.8 %** |
| `refined_yp82_rect_jacobi` | 1 | 0.25868 | 0.26217 | **-1.3 %** |
| `refined_yp82_rect_jacobi` | 2 | 0.13973 | 0.13612 | **+2.6 %** |
| `refined_yp82_rect_jacobi` | 4 | 0.07453 | 0.07250 | **+2.7 %** |
| `refined_yp82_rect_jacobi` | 8 | 0.04238 | 0.04021 | **+5.1 %** |
| `refined_yp82_rect_jacobi` | 16 | 0.02887 | 0.02539 | **+12.0 %** |
| `refined_yp82_rect_redblack` | 1 | 0.21970 | 0.22810 | **-3.8 %** |
| `refined_yp82_rect_redblack` | 2 | 0.11558 | 0.11957 | **-3.5 %** |
| `refined_yp82_rect_redblack` | 4 | 0.06252 | 0.06453 | **-3.2 %** |
| `refined_yp82_rect_redblack` | 8 | 0.03642 | 0.03700 | **-1.6 %** |
| `refined_yp82_rect_redblack` | 16 | 0.02503 | 0.02522 | **-0.8 %** |
| `refined_big_rect_jacobi` | 4 | 0.27459 | 0.26738 | **+2.6 %** |
| `refined_big_rect_jacobi` | 8 | 0.14316 | 0.13854 | **+3.2 %** |
| `refined_big_rect_jacobi` | 16 | 0.08248 | 0.07487 | **+9.2 %** |

## Per-round `mpi_wait`

| config | n=1 | n=2 | n=4 | n=8 | n=16 |
|---|---|---|---|---|---|
| `base_jacobi` (ref) | - | 27.2 us | 121.8 us | 116.2 us | 168.6 us |
| `base_jacobi` (new) | - | 30.0 us | 96.9 us | 88.3 us | 97.4 us |
| `rect_jacobi` (ref) | - | 44.5 us | 76.0 us | 78.2 us | 170.6 us |
| `rect_jacobi` (new) | - | 28.0 us | 55.3 us | 52.5 us | 80.0 us |
| `refined_yp82_rect_jacobi` (ref) | - | 39.8 us | 42.9 us | 49.9 us | 105.2 us |
| `refined_yp82_rect_jacobi` (new) | - | 25.1 us | 42.7 us | 36.9 us | 54.9 us |
| `refined_yp82_rect_redblack` (ref) | - | 24.7 us | 39.5 us | 46.8 us | 96.2 us |
| `refined_yp82_rect_redblack` (new) | - | 26.0 us | 41.1 us | 50.2 us | 97.4 us |
| `refined_big_rect_jacobi` (ref) | - | - | 118.5 us | 92.0 us | 199.2 us |
| `refined_big_rect_jacobi` (new) | - | - | 87.5 us | 65.0 us | 71.4 us |

## Block tax — inverse per-GPU throughput `t_step * n / cells` (ns)

| ranks | Mcell/GPU | base ref | base new | rect ref | rect new | **tax ref** | **tax new** |
|---|---|---|---|---|---|---|---|
| 1 | 138.41 | 4.192 | 4.114 | 4.220 | 4.261 | **1.007** | **1.036** |
| 2 | 69.21 | 4.167 | 4.107 | 4.436 | 4.323 | **1.065** | **1.053** |
| 4 | 34.60 | 4.366 | 4.250 | 4.562 | 4.411 | **1.045** | **1.038** |
| 8 | 17.30 | 4.771 | 4.558 | 4.821 | 4.604 | **1.011** | **1.010** |
| 16 | 8.65 | 5.475 | 5.001 | 5.855 | 5.167 | **1.069** | **1.033** |

## Strong-scaling efficiency (vs each config's smallest rank count)

| config | n=1 | n=2 | n=4 | n=8 | n=16 |
|---|---|---|---|---|---|
| `base_jacobi` (ref) | 100 % | 101 % | 96 % | 88 % | 77 % |
| `base_jacobi` (new) | 100 % | 100 % | 97 % | 90 % | 82 % |
| `rect_jacobi` (ref) | 100 % | 95 % | 93 % | 88 % | 72 % |
| `rect_jacobi` (new) | 100 % | 99 % | 97 % | 93 % | 82 % |
| `refined_yp82_rect_jacobi` (ref) | 100 % | 93 % | 87 % | 76 % | 56 % |
| `refined_yp82_rect_jacobi` (new) | 100 % | 96 % | 90 % | 82 % | 65 % |
| `refined_yp82_rect_redblack` (ref) | 100 % | 95 % | 88 % | 75 % | 55 % |
| `refined_yp82_rect_redblack` (new) | 100 % | 95 % | 88 % | 77 % | 57 % |
| `refined_big_rect_jacobi` (ref) | - | - | 100 % | 96 % | 83 % |
| `refined_big_rect_jacobi` (new) | - | - | 100 % | 96 % | 89 % |

## Headline 1 -- the 2:1 machinery, like for like

`refined_big_rect_jacobi` shares grid, block shape and solver with
`rect_jacobi`; it just refines one y-tile. Cost of the ADDED cells
against the coarse cells beside them:

| ranks | binary | rect s/step | big s/step | ns per EXTRA cell | vs coarse-cell avg |
|---|---|---|---|---|---|
| 4 | ref | 0.15785 | 0.27459 | 1.1246 | **0.986x** |
| 4 | new | 0.15263 | 0.26738 | 1.1054 | **1.002x** |
| 8 | ref | 0.08342 | 0.14316 | 0.5755 | **0.955x** |
| 8 | new | 0.07966 | 0.13854 | 0.5672 | **0.986x** |
| 16 | ref | 0.05065 | 0.08248 | 0.3066 | **0.838x** |
| 16 | new | 0.04470 | 0.07487 | 0.2907 | **0.900x** |

## Headline 3 — red-black against Jacobi

| binary | n=1 | n=2 | n=4 | n=8 | n=16 |
|---|---|---|---|---|---|
| ref | 0.849 | 0.827 | 0.839 | 0.859 | 0.867 |
| new | 0.870 | 0.878 | 0.890 | 0.920 | 0.993 |
