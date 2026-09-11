# Exchange op split, per-round cost, and the A0 probe at scale

```
job    : 5139977
date   : 2026-09-11 15:52:49 CEST
nodes  : 4  (hkn[0420,0428,0435,0516])
commit : e0542a53a25db0e453dd0f2cf9afba360821ed02
dirty  : 0
ref    : 703e62d4595dfc25c593d7a8db99035a5f9af6aa
gpu    : NVIDIA A100-SXM4-40GB
```

## Gate -- fields against the pre-change binary

**gate_rect_jacobi**

```
comparing /hkfs/work/workspace/scratch/xt8786-roughBoundaryLayer/optimiseBlockRefinement/moby-2to1-run/results_exchange16/gate_rect_jacobi_ref/overhead_20.h5 vs /hkfs/work/workspace/scratch/xt8786-roughBoundaryLayer/optimiseBlockRefinement/moby-2to1-run/results_exchange16/gate_rect_jacobi_new/overhead_20.h5
  un           n=138412032    max_abs=0
  vn           n=138412032    max_abs=0
  wn           n=138412032    max_abs=0
  pn           n=138412032    max_abs=0
OK: worst max_abs = 0 over 4 dataset(s)
```

**gate_refined_yp82_rect_jacobi**

```
comparing /hkfs/work/workspace/scratch/xt8786-roughBoundaryLayer/optimiseBlockRefinement/moby-2to1-run/results_exchange16/gate_refined_yp82_rect_jacobi_ref/overhead_20.h5 vs /hkfs/work/workspace/scratch/xt8786-roughBoundaryLayer/optimiseBlockRefinement/moby-2to1-run/results_exchange16/gate_refined_yp82_rect_jacobi_new/overhead_20.h5
  un           n=60555264     max_abs=0
  vn           n=60555264     max_abs=0
  wn           n=60555264     max_abs=0
  pn           n=60555264     max_abs=0
OK: worst max_abs = 0 over 4 dataset(s)
```


## Pass 1a -- exchange volume by op (points summed over ranks)

| run | peers/rank | local copy | local restrict | local prolong | send copy | send restrict | send prolong | cross-level % |
|---|---|---|---|---|---|---|---|---|
| base_jacobi_n16 | 11 | 0 | 0 | 0 | 4,735,872 | 0 | 0 | 0.00 % |
| nb16_jacobi_n16 | 2 | 55,738,656 | 0 | 0 | 1,283,040 | 0 | 0 | 0.00 % |
| rect_jacobi_n16 | 2 | 14,285,952 | 0 | 0 | 1,104,000 | 0 | 0 | 0.00 % |
| refined_big_rect_jacobi_n16 | 2 | 20,200,752 | 838,400 | 3,366,400 | 1,362,000 | 6,000 | 12,000 | 16.38 % |
| refined_yp82_rect_jacobi_n16 | 2 | 4,688,280 | 208,000 | 838,400 | 681,000 | 3,000 | 6,000 | 16.43 % |

## Pass 1b -- entries, and points per entry

| run | copy entries | cross entries | pts/entry copy | pts/entry cross |
|---|---|---|---|---|
| base_jacobi_n16 | 224 | 0 | - | - |
| nb16_jacobi_n16 | 821,064 | 0 | - | - |
| rect_jacobi_n16 | 21,776 | 0 | - | - |
| refined_big_rect_jacobi_n16 | 23,336 | 13,232 | 924 | 319 |
| refined_yp82_rect_jacobi_n16 | 5,780 | 3,288 | 929 | 321 |

## Pass 1c -- microseconds per exchange ROUND, against points per round

Points are PER RANK (the printed totals divided by the rank count);
times are rank 0's. `local_copy` is same-level only, `copy_cross` is
the 2:1 restrict/prolong kernel (exactly zero on a single-level grid).

| run | s/step | local pts/rank | send pts/rank | pack us | unpack us | local_copy us | copy_cross us | mpi_wait us | rounds/step |
|---|---|---|---|---|---|---|---|---|---|
| base_jacobi_n16 | 0.061008 | 0 | 295,992 | 130.5 | 105.6 | 0.1 | 0.1 | 165.9 | 39 |
| nb16_jacobi_n16 | 0.083823 | 3,483,666 | 80,190 | 101.7 | 97.6 | 431.7 | 0.2 | 160.6 | 39 |
| rect_jacobi_n16 | 0.066560 | 892,872 | 69,000 | 100.6 | 91.6 | 145.1 | 0.2 | 138.7 | 39 |
| refined_big_rect_jacobi_n16 | 0.106302 | 1,525,347 | 86,250 | 105.7 | 95.0 | 192.0 | 116.1 | 105.8 | 39 |
| refined_yp82_rect_jacobi_n16 | 0.041729 | 358,417 | 43,125 | 100.3 | 86.9 | 93.1 | 95.3 | 86.0 | 39 |

### The same, as a share of the step

| run | exchange total s/step | % of step | device-local % | mpi_wait % |
|---|---|---|---|---|
| base_jacobi_n16 | 0.016844 | 27.6 % | 15.1 % | 10.6 % |
| nb16_jacobi_n16 | 0.031107 | 37.1 % | 29.4 % | 7.5 % |
| rect_jacobi_n16 | 0.018805 | 28.3 % | 19.8 % | 8.1 % |
| refined_big_rect_jacobi_n16 | 0.022541 | 21.2 % | 17.0 % | 3.9 % |
| refined_yp82_rect_jacobi_n16 | 0.016819 | 40.3 % | 31.7 % | 8.0 % |

## Pass 1d -- fixed vs per-point cost, per config

Least squares of microseconds-per-round on points-per-rank over the
rank counts of one config. The INTERCEPT is what one exchange round
costs with no points to move; the SLOPE is the marginal point. Only
configs measured at 3+ rank counts are fitted.

| config | bucket | rank counts | pts/rank span | fixed us/round | ns/pt | worst residual |
|---|---|---|---|---|---|---|

## Pass 2 -- the A0 overlap probe

If an in-flight transfer progressed while the probe kernel ran, `mpi_wait`
would collapse toward zero as `a0_probe` grows.

| case | probe us/round | mpi_wait us/round | vs its baseline | s/step |
|---|---|---|---|---|

## L2_div -- the correctness check

| run | last L2_div |
|---|---|
| gate_rect_jacobi_new | 1.06618625E-05 |
| gate_rect_jacobi_ref | 1.06618625E-05 |
| gate_refined_yp82_rect_jacobi_new | 1.43477585E-05 |
| gate_refined_yp82_rect_jacobi_ref | 1.43477585E-05 |
| op_base_jacobi_n16 | 1.11230086E-05 |
| op_nb16_jacobi_n16 | 1.11230086E-05 |
| op_rect_jacobi_n16 | 1.11230086E-05 |
| op_refined_big_rect_jacobi_n16 | 9.03333418E-06 |
| op_refined_yp82_rect_jacobi_n16 | 2.13035041E-05 |
