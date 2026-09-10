# Exchange op split, per-round cost, and the A0 probe at scale

```
job    : 5139581
date   : 2026-09-10 16:22:56 CEST
nodes  : 2  (hkn[0401,0403])
commit : 5b6b813940e951bc20bd5bdf07412be758c190a6
dirty  : 0
ref    : 703e62d4595dfc25c593d7a8db99035a5f9af6aa
gpu    : NVIDIA A100-SXM4-40GB
```

## Gate -- fields against the pre-change binary

**gate_rect_jacobi**

```
comparing /hkfs/work/workspace/scratch/xt8786-roughBoundaryLayer/optimiseBlockRefinement/moby-2to1-run/results_exchange/gate_rect_jacobi_ref/overhead_20.h5 vs /hkfs/work/workspace/scratch/xt8786-roughBoundaryLayer/optimiseBlockRefinement/moby-2to1-run/results_exchange/gate_rect_jacobi_new/overhead_20.h5
  un           n=138412032    max_abs=0
  vn           n=138412032    max_abs=0
  wn           n=138412032    max_abs=0
  pn           n=138412032    max_abs=0
OK: worst max_abs = 0 over 4 dataset(s)
```

**gate_refined_yp82_rect_jacobi**

```
comparing /hkfs/work/workspace/scratch/xt8786-roughBoundaryLayer/optimiseBlockRefinement/moby-2to1-run/results_exchange/gate_refined_yp82_rect_jacobi_ref/overhead_20.h5 vs /hkfs/work/workspace/scratch/xt8786-roughBoundaryLayer/optimiseBlockRefinement/moby-2to1-run/results_exchange/gate_refined_yp82_rect_jacobi_new/overhead_20.h5
  un           n=60555264     max_abs=0
  vn           n=60555264     max_abs=0
  wn           n=60555264     max_abs=0
  pn           n=60555264     max_abs=0
OK: worst max_abs = 0 over 4 dataset(s)
```


## Pass 1a -- exchange volume by op (points summed over ranks)

| run | peers/rank | local copy | local restrict | local prolong | send copy | send restrict | send prolong | cross-level % |
|---|---|---|---|---|---|---|---|---|
| base_jacobi_n1 | 0 | 1,458,888 | 0 | 0 | 0 | 0 | 0 | 0.00 % |
| base_jacobi_n4 | 3 | 1,458,888 | 0 | 0 | 1,659,864 | 0 | 0 | 0.00 % |
| base_jacobi_n8 | 7 | 0 | 0 | 0 | 4,594,752 | 0 | 0 | 0.00 % |
| nb16_jacobi_n1 | 0 | 57,021,696 | 0 | 0 | 0 | 0 | 0 | 0.00 % |
| nb16_jacobi_n4 | 2 | 56,765,088 | 0 | 0 | 256,608 | 0 | 0 | 0.00 % |
| nb16_jacobi_n8 | 2 | 56,422,944 | 0 | 0 | 598,752 | 0 | 0 | 0.00 % |
| rect_jacobi_n1 | 0 | 15,389,952 | 0 | 0 | 0 | 0 | 0 | 0.00 % |
| rect_jacobi_n4 | 2 | 15,169,152 | 0 | 0 | 220,800 | 0 | 0 | 0.00 % |
| rect_jacobi_n8 | 2 | 14,874,752 | 0 | 0 | 515,200 | 0 | 0 | 0.00 % |
| refined_big_rect_jacobi_n4 | 2 | 21,290,352 | 843,200 | 3,376,000 | 272,400 | 1,200 | 2,400 | 16.38 % |
| refined_big_rect_jacobi_n8 | 2 | 20,927,152 | 841,600 | 3,372,800 | 635,600 | 2,800 | 5,600 | 16.38 % |
| refined_yp82_rect_jacobi_n1 | 0 | 5,369,280 | 211,000 | 844,400 | 0 | 0 | 0 | 16.43 % |
| refined_yp82_rect_jacobi_n4 | 2 | 5,233,080 | 210,400 | 843,200 | 136,200 | 600 | 1,200 | 16.43 % |
| refined_yp82_rect_jacobi_n8 | 2 | 5,051,480 | 209,600 | 841,600 | 317,800 | 1,400 | 2,800 | 16.43 % |

## Pass 1b -- entries, and points per entry

| run | copy entries | cross entries | pts/entry copy | pts/entry cross |
|---|---|---|---|---|
| base_jacobi_n1 | 2 | 0 | - | - |
| base_jacobi_n4 | 44 | 0 | - | - |
| base_jacobi_n8 | 88 | 0 | - | - |
| nb16_jacobi_n1 | 821,064 | 0 | - | - |
| nb16_jacobi_n4 | 821,064 | 0 | - | - |
| nb16_jacobi_n8 | 821,064 | 0 | - | - |
| rect_jacobi_n1 | 21,776 | 0 | - | - |
| rect_jacobi_n4 | 21,776 | 0 | - | - |
| rect_jacobi_n8 | 21,776 | 0 | - | - |
| refined_big_rect_jacobi_n4 | 23,336 | 13,232 | 924 | 319 |
| refined_big_rect_jacobi_n8 | 23,336 | 13,232 | 924 | 319 |
| refined_yp82_rect_jacobi_n1 | 5,780 | 3,288 | 929 | 321 |
| refined_yp82_rect_jacobi_n4 | 5,780 | 3,288 | 929 | 321 |
| refined_yp82_rect_jacobi_n8 | 5,780 | 3,288 | 929 | 321 |

## Pass 1c -- microseconds per exchange ROUND, against points per round

Points are PER RANK (the printed totals divided by the rank count);
times are rank 0's. `local_copy` is same-level only, `copy_cross` is
the 2:1 restrict/prolong kernel (exactly zero on a single-level grid).

| run | s/step | local pts/rank | send pts/rank | pack us | unpack us | local_copy us | copy_cross us | mpi_wait us | rounds/step |
|---|---|---|---|---|---|---|---|---|---|
| base_jacobi_n1 | 0.684979 | 1,458,888 | 0 | 0.0 | 0.0 | 111.0 | 0.1 | 0.0 | 0 |
| base_jacobi_n4 | 0.187168 | 364,722 | 414,966 | 151.2 | 113.2 | 82.1 | 0.2 | 121.0 | 39 |
| base_jacobi_n8 | 0.103252 | 0 | 574,344 | 165.5 | 125.0 | 0.1 | 0.1 | 133.7 | 39 |
| nb16_jacobi_n1 | 0.901933 | 57,021,696 | 0 | 0.0 | 0.0 | 4987.6 | 0.2 | 0.0 | 0 |
| nb16_jacobi_n4 | 0.264418 | 14,191,272 | 64,152 | 99.7 | 96.2 | 1576.7 | 0.1 | 73.5 | 39 |
| nb16_jacobi_n8 | 0.141611 | 7,052,868 | 74,844 | 100.9 | 97.0 | 813.0 | 0.2 | 81.3 | 39 |
| rect_jacobi_n1 | 0.720748 | 15,389,952 | 0 | 0.0 | 0.0 | 1251.4 | 0.1 | 0.0 | 0 |
| rect_jacobi_n4 | 0.200528 | 3,792,288 | 55,200 | 101.3 | 92.4 | 425.5 | 0.2 | 69.1 | 39 |
| rect_jacobi_n8 | 0.109302 | 1,859,344 | 64,400 | 100.3 | 91.6 | 237.3 | 0.1 | 76.9 | 39 |
| refined_big_rect_jacobi_n4 | 0.348129 | 6,377,388 | 69,000 | 103.9 | 96.2 | 620.9 | 208.8 | 86.2 | 39 |
| refined_big_rect_jacobi_n8 | 0.186315 | 3,142,694 | 80,500 | 104.1 | 95.5 | 334.9 | 148.0 | 88.9 | 39 |
| refined_yp82_rect_jacobi_n1 | 0.325160 | 6,424,680 | 0 | 0.0 | 0.0 | 492.2 | 209.8 | 0.0 | 0 |
| refined_yp82_rect_jacobi_n4 | 0.100293 | 1,571,670 | 34,500 | 98.5 | 85.6 | 197.8 | 114.7 | 61.6 | 39 |
| refined_yp82_rect_jacobi_n8 | 0.060292 | 762,835 | 40,250 | 99.7 | 85.7 | 130.0 | 99.6 | 55.6 | 39 |

### The same, as a share of the step

| run | exchange total s/step | % of step | device-local % | mpi_wait % |
|---|---|---|---|---|
| base_jacobi_n1 | 0.002665 | 0.4 % | 0.4 % | 0.0 % |
| base_jacobi_n4 | 0.018809 | 10.0 % | 7.2 % | 2.5 % |
| base_jacobi_n8 | 0.017727 | 17.2 % | 11.0 % | 5.1 % |
| nb16_jacobi_n1 | 0.119704 | 13.3 % | 13.3 % | 0.0 % |
| nb16_jacobi_n4 | 0.072231 | 27.3 % | 26.1 % | 1.1 % |
| nb16_jacobi_n8 | 0.042872 | 30.3 % | 27.8 % | 2.2 % |
| rect_jacobi_n1 | 0.030035 | 4.2 % | 4.2 % | 0.0 % |
| rect_jacobi_n4 | 0.027104 | 13.5 % | 12.0 % | 1.3 % |
| rect_jacobi_n8 | 0.019998 | 18.3 % | 15.3 % | 2.7 % |
| refined_big_rect_jacobi_n4 | 0.040664 | 11.7 % | 10.6 % | 1.0 % |
| refined_big_rect_jacobi_n8 | 0.028149 | 15.1 % | 13.1 % | 1.9 % |
| refined_yp82_rect_jacobi_n1 | 0.016846 | 5.2 % | 5.2 % | 0.0 % |
| refined_yp82_rect_jacobi_n4 | 0.020299 | 20.2 % | 17.6 % | 2.4 % |
| refined_yp82_rect_jacobi_n8 | 0.017109 | 28.4 % | 24.4 % | 3.6 % |

## Pass 1d -- fixed vs per-point cost, per config

Least squares of microseconds-per-round on points-per-rank over the
rank counts of one config. The INTERCEPT is what one exchange round
costs with no points to move; the SLOPE is the marginal point. Only
configs measured at 3+ rank counts are fitted.

| config | bucket | rank counts | pts/rank span | fixed us/round | ns/pt | worst residual |
|---|---|---|---|---|---|---|
| base_jacobi | pack | 1,4,8 | 0-574,344 | 5.5 | 0.303 | 19.8 us |
| base_jacobi | unpack | 1,4,8 | 0-574,344 | 4.0 | 0.229 | 14.3 us |
| base_jacobi | local_copy | 1,4,8 | 0-1,458,888 | 25.1 | 0.065 | 33.4 us |
| base_jacobi | copy_cross | 1,4,8 | 0-1,458,888 | 0.1 | 0.000 | 0.0 us |
| nb16_jacobi | pack | 1,4,8 | 0-74,844 | 1.1 | 1.420 | 7.5 us |
| nb16_jacobi | unpack | 1,4,8 | 0-74,844 | 1.1 | 1.367 | 7.4 us |
| nb16_jacobi | local_copy | 1,4,8 | 7,052,868-57,021,696 | 315.1 | 0.082 | 95.3 us |
| nb16_jacobi | copy_cross | 1,4,8 | 7,052,868-57,021,696 | 0.2 | -0.000 | 0.0 us |
| rect_jacobi | pack | 1,4,8 | 0-64,400 | 1.3 | 1.655 | 8.8 us |
| rect_jacobi | unpack | 1,4,8 | 0-64,400 | 1.1 | 1.510 | 7.9 us |
| rect_jacobi | local_copy | 1,4,8 | 1,859,344-15,389,952 | 121.5 | 0.074 | 24.7 us |
| rect_jacobi | copy_cross | 1,4,8 | 1,859,344-15,389,952 | 0.2 | -0.000 | 0.0 us |
| refined_yp82_rect_jacobi | pack | 1,4,8 | 0-40,250 | 1.1 | 2.610 | 7.4 us |
| refined_yp82_rect_jacobi | unpack | 1,4,8 | 0-40,250 | 1.0 | 2.253 | 6.9 us |
| refined_yp82_rect_jacobi | local_copy | 1,4,8 | 762,835-6,424,680 | 90.0 | 0.063 | 9.1 us |
| refined_yp82_rect_jacobi | copy_cross | 1,4,8 | 762,835-6,424,680 | 84.4 | 0.020 | 0.4 us |

## Pass 2 -- the A0 overlap probe

If an in-flight transfer progressed while the probe kernel ran, `mpi_wait`
would collapse toward zero as `a0_probe` grows.

| case | probe us/round | mpi_wait us/round | vs its baseline | s/step |
|---|---|---|---|---|
| rect_jacobi_r8N2_base | 0.0 | 75.0 | - | 0.109144 |
| rect_jacobi_r8N2_probe | 5700.4 | 79.4 | 1.06x | 0.331868 |
| refined_yp82_rect_jacobi_r4N1_base | 0.0 | 57.3 | - | 0.100154 |
| refined_yp82_rect_jacobi_r4N1_probe | 5695.1 | 65.5 | 1.14x | 0.322845 |
| refined_yp82_rect_jacobi_r8N2_base | 0.0 | 49.1 | - | 0.060138 |
| refined_yp82_rect_jacobi_r8N2_probe | 5701.5 | 70.7 | 1.44x | 0.283078 |

## L2_div -- the correctness check

| run | last L2_div |
|---|---|
| a0_rect_jacobi_r8N2_base | 1.11230086E-05 |
| a0_rect_jacobi_r8N2_probe | 1.11230086E-05 |
| a0_refined_yp82_rect_jacobi_r4N1_base | 2.13035041E-05 |
| a0_refined_yp82_rect_jacobi_r4N1_probe | 2.13035041E-05 |
| a0_refined_yp82_rect_jacobi_r8N2_base | 2.13035041E-05 |
| a0_refined_yp82_rect_jacobi_r8N2_probe | 2.13035041E-05 |
| gate_rect_jacobi_new | 1.06618625E-05 |
| gate_rect_jacobi_ref | 1.06618625E-05 |
| gate_refined_yp82_rect_jacobi_new | 1.43477585E-05 |
| gate_refined_yp82_rect_jacobi_ref | 1.43477585E-05 |
| op_base_jacobi_n1 | 1.11230086E-05 |
| op_base_jacobi_n4 | 1.11230086E-05 |
| op_base_jacobi_n8 | 1.11230086E-05 |
| op_nb16_jacobi_n1 | 1.11230086E-05 |
| op_nb16_jacobi_n4 | 1.11230086E-05 |
| op_nb16_jacobi_n8 | 1.11230086E-05 |
| op_rect_jacobi_n1 | 1.11230086E-05 |
| op_rect_jacobi_n4 | 1.11230086E-05 |
| op_rect_jacobi_n8 | 1.11230086E-05 |
| op_refined_big_rect_jacobi_n4 | 9.03333418E-06 |
| op_refined_big_rect_jacobi_n8 | 9.03333418E-06 |
| op_refined_yp82_rect_jacobi_n1 | 2.13035041E-05 |
| op_refined_yp82_rect_jacobi_n4 | 2.13035041E-05 |
| op_refined_yp82_rect_jacobi_n8 | 2.13035041E-05 |
