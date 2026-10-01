# The campaign, re-measured

The 23-run matrix (`run_matrix.sh`, `--map-by numa --bind-to core`), run
TWICE IN ONE ALLOCATION -- `ref` is the pre-change binary, `new` the
post-change one, and the two columns differ only in that change. WHICH
change is recorded in the run directory's provenance.txt, not here: this
collector has now served two of them (the rank-to-GPU mapping, 2026-09-10,
and the comm_type parent map, 2026-09-14), and hard-coding either into the
header is how a table gets quoted against the wrong question.

```
job    : 5174077
date   : 2026-10-01 17:47:42 CEST
nodes  : 4  (hkn[0502,0526,0535,0617])
new    : d2ac8397684f4ac786c3343d954475bb595e58bd  (minimum-surface order)
ref    : 5bdc5eb4ceb6cce423789cf360563169bd49bb3a  (legacy order)
gpu    : NVIDIA A100-SXM4-40GB
nsteps : 200
layout : --map-by numa --bind-to core (the default placement of every matrix)
```

## s/step, and what the mapping is worth

| config | ranks | ref | new | gain |
|---|---|---|---|---|
| `base_jacobi` | 16 | 0.06740 | 0.04152 | **+38.4 %** |
| `rect_jacobi` | 16 | 0.04270 | 0.04250 | **+0.5 %** |

## Per-round `mpi_wait`

| config | n=1 | n=2 | n=4 | n=8 | n=16 |
|---|---|---|---|---|---|
| `base_jacobi` (ref) | - | - | - | - | 715.3 us |
| `base_jacobi` (new) | - | - | - | - | 92.4 us |
| `rect_jacobi` (ref) | - | - | - | - | 81.0 us |
| `rect_jacobi` (new) | - | - | - | - | 82.1 us |

## Block tax — inverse per-GPU throughput `t_step * n / cells` (ns)

| ranks | Mcell/GPU | base ref | base new | rect ref | rect new | **tax ref** | **tax new** |
|---|---|---|---|---|---|---|---|
| 1 | - | - | - | - | - | - | - |
| 2 | - | - | - | - | - | - | - |
| 4 | - | - | - | - | - | - | - |
| 8 | - | - | - | - | - | - | - |
| 16 | 8.65 | 7.791 | 4.799 | 4.936 | 4.913 | **0.634** | **1.024** |

## Strong-scaling efficiency (vs each config's smallest rank count)

| config | n=1 | n=2 | n=4 | n=8 | n=16 |
|---|---|---|---|---|---|
| `base_jacobi` (ref) | - | - | - | - | 100 % |
| `base_jacobi` (new) | - | - | - | - | 100 % |
| `rect_jacobi` (ref) | - | - | - | - | 100 % |
| `rect_jacobi` (new) | - | - | - | - | 100 % |

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
