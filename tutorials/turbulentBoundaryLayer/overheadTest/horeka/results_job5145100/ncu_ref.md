# What limits `jacobi_apply`

## refined_yp82_rect_jacobi  (60.56 Mcell)

| kernel | launches | us | doubles/cell | min | vs min | ld sect/req | st sect/req | L2 hit % | DRAM %peak | SM %peak | occupancy % | registers |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `jacobi_compute_phi__F1L418_10` | 4 | 3,181 | 5.30 | 5 | **1.06x** | 2.49 | 9.00 | 50.88 | 51.93 | 47.07 | 29.12 | 94.00 |
| `jacobi_apply__F1L470_14` | 4 | 1,536 | 3.08 | 3 | **1.03x** | 1.90 | 9.00 | 49.13 | 62.40 | 66.21 | 45.45 | 59.00 |
| `jacobi_apply__F1L502_16` | 4 | 7,565 | 10.31 | 8 | **1.29x** | 2.13 | 8.97 | 50.86 | 42.45 | 32.14 | 29.56 | 88.00 |

