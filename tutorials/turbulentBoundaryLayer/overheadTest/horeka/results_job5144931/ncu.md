# What limits `jacobi_apply`

```
job    : 5144931
date   : 2026-09-14 10:03:33 CEST
node   : hkn0403
commit : ca0ea3d60f6a5cec16a3f68888b17aac580644a4
gpu    : NVIDIA A100-SXM4-40GB
```

## refined_yp82_rect_jacobi  (60.56 Mcell)

| kernel | launches | us | doubles/cell | min | vs min | ld sect/req | st sect/req | L2 hit % | DRAM %peak | SM %peak | occupancy % | registers |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `jacobi_compute_phi__F1L418_10` | 4 | 3,204 | 5.48 | 5 | **1.10x** | 2.49 | 9.00 | 51.45 | 53.27 | 47.01 | 29.07 | 94.00 |
| `jacobi_apply__F1L470_14` | 4 | 1,543 | 3.21 | 3 | **1.07x** | 1.90 | 9.00 | 49.86 | 64.80 | 66.12 | 45.34 | 59.00 |
| `jacobi_apply__F1L502_16` | 4 | 7,654 | 10.82 | 8 | **1.35x** | 2.13 | 8.97 | 51.90 | 44.03 | 32.11 | 29.55 | 88.00 |

