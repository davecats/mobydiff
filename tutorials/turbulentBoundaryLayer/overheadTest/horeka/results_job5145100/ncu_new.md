# What limits `jacobi_apply`

## refined_yp82_rect_jacobi  (60.56 Mcell)

| kernel | launches | us | doubles/cell | min | vs min | ld sect/req | st sect/req | L2 hit % | DRAM %peak | SM %peak | occupancy % | registers |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `jacobi_compute_phi__F1L499_10` | 4 | 3,181 | 5.30 | 5 | **1.06x** | 2.49 | 9.00 | 50.91 | 51.92 | 47.07 | 29.12 | 94.00 |
| `jacobi_apply__F1L550_14` | 4 | 1,536 | 3.08 | - | - | 1.90 | 9.00 | 49.06 | 62.46 | 66.23 | 45.44 | 59.00 |
| `jacobi_apply__F1L585_16` | 4 | 5,309 | 10.31 | - | - | 2.36 | 8.97 | 50.86 | 60.48 | 40.55 | 46.04 | 64.00 |

