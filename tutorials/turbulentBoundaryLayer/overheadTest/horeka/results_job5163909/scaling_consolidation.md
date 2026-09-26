# The campaign, re-measured

The 23-run matrix (`run_matrix.sh`, `--map-by numa --bind-to core`), run
TWICE IN ONE ALLOCATION -- `ref` is the pre-change binary, `new` the
post-change one, and the two columns differ only in that change. WHICH
change is recorded in the run directory's provenance.txt, not here: this
collector has now served two of them (the rank-to-GPU mapping, 2026-09-10,
and the comm_type parent map, 2026-09-14), and hard-coding either into the
header is how a table gets quoted against the wrong question.

```
job    : 5163909
date   : 2026-09-26 07:04:09 CEST
nodes  : 2  (hkn[0401,0403])
new    : 2a0ed47ff993a0bcd132b10715f9262cb201de88  (new, configs as shipped)
ref    : 8fa0fc2c38982e13c0a64a7bf773bd5d043e03b0  (ref, configs + convection = skew)
gpu    : NVIDIA A100-SXM4-40GB
nsteps : 200
layout : --map-by numa --bind-to core (as every earlier matrix)
```

## s/step, and what the mapping is worth

| config | ranks | ref | new | gain |
|---|---|---|---|---|
| `base_jacobi` | 1 | 0.55167 | 0.56280 | **-2.0 %** |
| `base_jacobi` | 2 | 0.27447 | 0.27975 | **-1.9 %** |
| `base_jacobi` | 4 | 0.14200 | 0.14447 | **-1.7 %** |
| `base_jacobi` | 8 | 0.07642 | 0.07775 | **-1.7 %** |
| `rect_jacobi` | 1 | 0.55819 | 0.56844 | **-1.8 %** |
| `rect_jacobi` | 2 | 0.28428 | 0.28866 | **-1.5 %** |
| `rect_jacobi` | 4 | 0.14589 | 0.14812 | **-1.5 %** |
| `rect_jacobi` | 8 | 0.07637 | 0.07746 | **-1.4 %** |
| `refined_yp82_rect_jacobi` | 1 | 0.24969 | 0.25369 | **-1.6 %** |
| `refined_yp82_rect_jacobi` | 2 | 0.13007 | 0.13180 | **-1.3 %** |
| `refined_yp82_rect_jacobi` | 4 | 0.06933 | 0.07034 | **-1.4 %** |
| `refined_yp82_rect_jacobi` | 8 | 0.03891 | 0.03951 | **-1.5 %** |
| `refined_yp82_rect_redblack` | 1 | 0.22009 | 0.22403 | **-1.8 %** |
| `refined_yp82_rect_redblack` | 2 | 0.11564 | 0.11779 | **-1.9 %** |
| `refined_yp82_rect_redblack` | 4 | 0.06257 | 0.06366 | **-1.7 %** |
| `refined_yp82_rect_redblack` | 8 | 0.03641 | 0.03701 | **-1.6 %** |
| `refined_big_rect_jacobi` | 4 | 0.25526 | 0.25891 | **-1.4 %** |
| `refined_big_rect_jacobi` | 8 | 0.13254 | 0.13435 | **-1.4 %** |

## Per-round `mpi_wait`

| config | n=1 | n=2 | n=4 | n=8 | n=16 |
|---|---|---|---|---|---|
| `base_jacobi` (ref) | - | 19.3 us | 81.5 us | 79.9 us | - |
| `base_jacobi` (new) | - | 18.7 us | 81.1 us | 80.9 us | - |
| `rect_jacobi` (ref) | - | 23.5 us | 47.4 us | 41.3 us | - |
| `rect_jacobi` (new) | - | 20.6 us | 44.8 us | 43.7 us | - |
| `refined_yp82_rect_jacobi` (ref) | - | 17.9 us | 31.8 us | 32.7 us | - |
| `refined_yp82_rect_jacobi` (new) | - | 16.9 us | 29.5 us | 36.9 us | - |
| `refined_yp82_rect_redblack` (ref) | - | 22.5 us | 41.1 us | 50.4 us | - |
| `refined_yp82_rect_redblack` (new) | - | 23.9 us | 39.4 us | 50.2 us | - |
| `refined_big_rect_jacobi` (ref) | - | - | 55.0 us | 50.9 us | - |
| `refined_big_rect_jacobi` (new) | - | - | 55.0 us | 48.8 us | - |

## Block tax — inverse per-GPU throughput `t_step * n / cells` (ns)

| ranks | Mcell/GPU | base ref | base new | rect ref | rect new | **tax ref** | **tax new** |
|---|---|---|---|---|---|---|---|
| 1 | 138.41 | 3.986 | 4.066 | 4.033 | 4.107 | **1.012** | **1.010** |
| 2 | 69.21 | 3.966 | 4.042 | 4.108 | 4.171 | **1.036** | **1.032** |
| 4 | 34.60 | 4.104 | 4.175 | 4.216 | 4.281 | **1.027** | **1.025** |
| 8 | 17.30 | 4.417 | 4.494 | 4.414 | 4.477 | **0.999** | **0.996** |
| 16 | - | - | - | - | - | - | - |

## Strong-scaling efficiency (vs each config's smallest rank count)

| config | n=1 | n=2 | n=4 | n=8 | n=16 |
|---|---|---|---|---|---|
| `base_jacobi` (ref) | 100 % | 100 % | 97 % | 90 % | - |
| `base_jacobi` (new) | 100 % | 101 % | 97 % | 90 % | - |
| `rect_jacobi` (ref) | 100 % | 98 % | 96 % | 91 % | - |
| `rect_jacobi` (new) | 100 % | 98 % | 96 % | 92 % | - |
| `refined_yp82_rect_jacobi` (ref) | 100 % | 96 % | 90 % | 80 % | - |
| `refined_yp82_rect_jacobi` (new) | 100 % | 96 % | 90 % | 80 % | - |
| `refined_yp82_rect_redblack` (ref) | 100 % | 95 % | 88 % | 76 % | - |
| `refined_yp82_rect_redblack` (new) | 100 % | 95 % | 88 % | 76 % | - |
| `refined_big_rect_jacobi` (ref) | - | - | 100 % | 96 % | - |
| `refined_big_rect_jacobi` (new) | - | - | 100 % | 96 % | - |

## Headline 1 -- the 2:1 machinery, like for like

`refined_big_rect_jacobi` shares grid, block shape and solver with
`rect_jacobi`; it just refines one y-tile. Cost of the ADDED cells
against the coarse cells beside them:

| ranks | binary | rect s/step | big s/step | ns per EXTRA cell | vs coarse-cell avg |
|---|---|---|---|---|---|
| 4 | ref | 0.14589 | 0.25526 | 1.0536 | **1.000x** |
| 4 | new | 0.14812 | 0.25891 | 1.0672 | **0.997x** |
| 8 | ref | 0.07637 | 0.13254 | 0.5411 | **0.981x** |
| 8 | new | 0.07746 | 0.13435 | 0.5481 | **0.979x** |

## Headline 3 — red-black against Jacobi

| binary | n=1 | n=2 | n=4 | n=8 | n=16 |
|---|---|---|---|---|---|
| ref | 0.881 | 0.889 | 0.902 | 0.936 | - |
| new | 0.883 | 0.894 | 0.905 | 0.937 | - |
