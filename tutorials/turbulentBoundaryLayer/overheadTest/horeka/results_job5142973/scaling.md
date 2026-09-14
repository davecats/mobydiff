# The campaign, re-measured

The 23-run matrix (`run_matrix.sh`, `--map-by numa --bind-to core`), run
TWICE IN ONE ALLOCATION -- `ref` is the pre-change binary, `new` the
post-change one, and the two columns differ only in that change. WHICH
change is recorded in the run directory's provenance.txt, not here: this
collector has now served two of them (the rank-to-GPU mapping, 2026-09-10,
and the comm_type parent map, 2026-09-14), and hard-coding either into the
header is how a table gets quoted against the wrong question.

```
job    : 5142973
date   : 2026-09-14 03:26:27 CEST
nodes  : 4  (hkn[0620-0621,0627,0632])
new    : 55bee89aec5d8a29321a11a59b360915b97ed1d7  (map(to: c))
ref    : 365af767dd290bcdaf951a4409fb9acd9edd2154  (pre-fix)
gpu    : NVIDIA A100-SXM4-40GB
nsteps : 200
layout : --map-by numa --bind-to core (as every earlier matrix)
```

## s/step, and what the mapping is worth

| config | ranks | ref | new | gain |
|---|---|---|---|---|
| `base_jacobi` | 1 | 0.68328 | 0.68020 | **+0.4 %** |
| `base_jacobi` | 2 | 0.34926 | 0.34317 | **+1.7 %** |
| `base_jacobi` | 4 | 0.18634 | 0.18003 | **+3.4 %** |
| `base_jacobi` | 8 | 0.10277 | 0.09737 | **+5.2 %** |
| `base_jacobi` | 16 | 0.06008 | 0.05485 | **+8.7 %** |
| `rect_jacobi` | 1 | 0.71790 | 0.71631 | **+0.2 %** |
| `rect_jacobi` | 2 | 0.38221 | 0.37570 | **+1.7 %** |
| `rect_jacobi` | 4 | 0.19936 | 0.19271 | **+3.3 %** |
| `rect_jacobi` | 8 | 0.10847 | 0.10152 | **+6.4 %** |
| `rect_jacobi` | 16 | 0.06790 | 0.06031 | **+11.2 %** |
| `refined_yp82_rect_jacobi` | 1 | 0.32386 | 0.32027 | **+1.1 %** |
| `refined_yp82_rect_jacobi` | 2 | 0.17872 | 0.17072 | **+4.5 %** |
| `refined_yp82_rect_jacobi` | 4 | 0.09962 | 0.09075 | **+8.9 %** |
| `refined_yp82_rect_jacobi` | 8 | 0.05981 | 0.05137 | **+14.1 %** |
| `refined_yp82_rect_jacobi` | 16 | 0.04124 | 0.03344 | **+18.9 %** |
| `refined_yp82_rect_redblack` | 1 | 0.24658 | 0.24309 | **+1.4 %** |
| `refined_yp82_rect_redblack` | 2 | 0.13640 | 0.12786 | **+6.3 %** |
| `refined_yp82_rect_redblack` | 4 | 0.07790 | 0.06912 | **+11.3 %** |
| `refined_yp82_rect_redblack` | 8 | 0.04862 | 0.04013 | **+17.5 %** |
| `refined_yp82_rect_redblack` | 16 | 0.03514 | 0.02703 | **+23.1 %** |
| `refined_big_rect_jacobi` | 4 | 0.34612 | 0.33776 | **+2.4 %** |
| `refined_big_rect_jacobi` | 8 | 0.18515 | 0.17677 | **+4.5 %** |
| `refined_big_rect_jacobi` | 16 | 0.10760 | 0.09766 | **+9.2 %** |

## Per-round `mpi_wait`

| config | n=1 | n=2 | n=4 | n=8 | n=16 |
|---|---|---|---|---|---|
| `base_jacobi` (ref) | - | 42.7 us | 116.5 us | 128.9 us | 147.9 us |
| `base_jacobi` (new) | - | 32.1 us | 118.7 us | 115.4 us | 144.0 us |
| `rect_jacobi` (ref) | - | 45.6 us | 64.4 us | 74.5 us | 173.9 us |
| `rect_jacobi` (new) | - | 35.7 us | 54.8 us | 63.8 us | 157.0 us |
| `refined_yp82_rect_jacobi` (ref) | - | 38.9 us | 59.3 us | 51.6 us | 82.7 us |
| `refined_yp82_rect_jacobi` (new) | - | 27.3 us | 43.7 us | 46.1 us | 96.8 us |
| `refined_yp82_rect_redblack` (ref) | - | 33.4 us | 54.4 us | 53.4 us | 81.1 us |
| `refined_yp82_rect_redblack` (new) | - | 25.1 us | 36.6 us | 45.2 us | 85.8 us |
| `refined_big_rect_jacobi` (ref) | - | - | 67.8 us | 85.4 us | 132.8 us |
| `refined_big_rect_jacobi` (new) | - | - | 58.3 us | 66.5 us | 101.7 us |

## Block tax — inverse per-GPU throughput `t_step * n / cells` (ns)

| ranks | Mcell/GPU | base ref | base new | rect ref | rect new | **tax ref** | **tax new** |
|---|---|---|---|---|---|---|---|
| 1 | 138.41 | 4.937 | 4.914 | 5.187 | 5.175 | **1.051** | **1.053** |
| 2 | 69.21 | 5.047 | 4.959 | 5.523 | 5.429 | **1.094** | **1.095** |
| 4 | 34.60 | 5.385 | 5.203 | 5.761 | 5.569 | **1.070** | **1.070** |
| 8 | 17.30 | 5.940 | 5.628 | 6.269 | 5.868 | **1.055** | **1.043** |
| 16 | 8.65 | 6.946 | 6.341 | 7.849 | 6.972 | **1.130** | **1.100** |

## Strong-scaling efficiency (vs each config's smallest rank count)

| config | n=1 | n=2 | n=4 | n=8 | n=16 |
|---|---|---|---|---|---|
| `base_jacobi` (ref) | 100 % | 98 % | 92 % | 83 % | 71 % |
| `base_jacobi` (new) | 100 % | 99 % | 94 % | 87 % | 78 % |
| `rect_jacobi` (ref) | 100 % | 94 % | 90 % | 83 % | 66 % |
| `rect_jacobi` (new) | 100 % | 95 % | 93 % | 88 % | 74 % |
| `refined_yp82_rect_jacobi` (ref) | 100 % | 91 % | 81 % | 68 % | 49 % |
| `refined_yp82_rect_jacobi` (new) | 100 % | 94 % | 88 % | 78 % | 60 % |
| `refined_yp82_rect_redblack` (ref) | 100 % | 90 % | 79 % | 63 % | 44 % |
| `refined_yp82_rect_redblack` (new) | 100 % | 95 % | 88 % | 76 % | 56 % |
| `refined_big_rect_jacobi` (ref) | - | - | 100 % | 93 % | 80 % |
| `refined_big_rect_jacobi` (new) | - | - | 100 % | 96 % | 86 % |

## Headline 1 -- the 2:1 machinery, like for like

`refined_big_rect_jacobi` shares grid, block shape and solver with
`rect_jacobi`; it just refines one y-tile. Cost of the ADDED cells
against the coarse cells beside them:

| ranks | binary | rect s/step | big s/step | ns per EXTRA cell | vs coarse-cell avg |
|---|---|---|---|---|---|
| 4 | ref | 0.19936 | 0.34612 | 1.4138 | **0.982x** |
| 4 | new | 0.19271 | 0.33776 | 1.3973 | **1.004x** |
| 8 | ref | 0.10847 | 0.18515 | 0.7387 | **0.943x** |
| 8 | new | 0.10152 | 0.17677 | 0.7248 | **0.988x** |
| 16 | ref | 0.06790 | 0.10760 | 0.3824 | **0.779x** |
| 16 | new | 0.06031 | 0.09766 | 0.3598 | **0.826x** |

## Headline 3 — red-black against Jacobi

| binary | n=1 | n=2 | n=4 | n=8 | n=16 |
|---|---|---|---|---|---|
| ref | 0.761 | 0.763 | 0.782 | 0.813 | 0.852 |
| new | 0.759 | 0.749 | 0.762 | 0.781 | 0.808 |
