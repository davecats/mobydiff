# The campaign, re-measured

The 23-run matrix (`run_matrix.sh`, `--map-by numa --bind-to core`), run
TWICE IN ONE ALLOCATION -- `ref` is the pre-change binary, `new` the
post-change one, and the two columns differ only in that change. WHICH
change is recorded in the run directory's provenance.txt, not here: this
collector has now served two of them (the rank-to-GPU mapping, 2026-09-10,
and the comm_type parent map, 2026-09-14), and hard-coding either into the
header is how a table gets quoted against the wrong question.

```
job    : 5174076
date   : 2026-10-01 17:06:01 CEST
nodes  : 2  (hkn[0526,0535])
new    : d2ac8397684f4ac786c3343d954475bb595e58bd  (minimum-surface order)
ref    : 5bdc5eb4ceb6cce423789cf360563169bd49bb3a  (legacy order)
gpu    : NVIDIA A100-SXM4-40GB
nsteps : 200
layout : --map-by numa --bind-to core (the default placement of every matrix)
```

## s/step, and what the mapping is worth

| config | ranks | ref | new | gain |
|---|---|---|---|---|
| `base_jacobi` | 4 | 0.14145 | 0.14124 | **+0.2 %** |
| `base_jacobi` | 8 | 0.13067 | 0.07593 | **+41.9 %** |
| `rect_jacobi` | 4 | 0.14610 | 0.14500 | **+0.8 %** |
| `rect_jacobi` | 8 | 0.07606 | 0.07572 | **+0.4 %** |
| `refined_yp82_rect_jacobi` | 4 | 0.06907 | 0.06906 | **+0.0 %** |
| `refined_yp82_rect_jacobi` | 8 | 0.03843 | 0.03845 | **-0.1 %** |

## Per-round `mpi_wait`

| config | n=1 | n=2 | n=4 | n=8 | n=16 |
|---|---|---|---|---|---|
| `base_jacobi` (ref) | - | - | 107.8 us | 1514.6 us | - |
| `base_jacobi` (new) | - | - | 111.3 us | 93.2 us | - |
| `rect_jacobi` (ref) | - | - | 90.2 us | 59.5 us | - |
| `rect_jacobi` (new) | - | - | 84.0 us | 62.7 us | - |
| `refined_yp82_rect_jacobi` (ref) | - | - | 51.9 us | 36.7 us | - |
| `refined_yp82_rect_jacobi` (new) | - | - | 53.6 us | 41.2 us | - |

## Block tax — inverse per-GPU throughput `t_step * n / cells` (ns)

| ranks | Mcell/GPU | base ref | base new | rect ref | rect new | **tax ref** | **tax new** |
|---|---|---|---|---|---|---|---|
| 1 | - | - | - | - | - | - | - |
| 2 | - | - | - | - | - | - | - |
| 4 | 34.60 | 4.088 | 4.082 | 4.222 | 4.190 | **1.033** | **1.027** |
| 8 | 17.30 | 7.553 | 4.389 | 4.396 | 4.377 | **0.582** | **0.997** |
| 16 | - | - | - | - | - | - | - |

## Strong-scaling efficiency (vs each config's smallest rank count)

| config | n=1 | n=2 | n=4 | n=8 | n=16 |
|---|---|---|---|---|---|
| `base_jacobi` (ref) | - | - | 100 % | 54 % | - |
| `base_jacobi` (new) | - | - | 100 % | 93 % | - |
| `rect_jacobi` (ref) | - | - | 100 % | 96 % | - |
| `rect_jacobi` (new) | - | - | 100 % | 96 % | - |
| `refined_yp82_rect_jacobi` (ref) | - | - | 100 % | 90 % | - |
| `refined_yp82_rect_jacobi` (new) | - | - | 100 % | 90 % | - |

## Headline 1 -- the 2:1 machinery, like for like

`refined_big_rect_jacobi` shares grid, block shape and solver with
`rect_jacobi`; it just refines one y-tile. Cost of the ADDED cells
against the coarse cells beside them:

| ranks | binary | rect s/step | big s/step | ns per EXTRA cell | vs coarse-cell avg |
|---|---|---|---|---|---|

## Headline 3 — red-black against Jacobi

| binary | n=1 | n=2 | n=4 | n=8 | n=16 |
|---|---|---|---|---|---|
| ref | - | - | - | - | - |
| new | - | - | - | - | - |
