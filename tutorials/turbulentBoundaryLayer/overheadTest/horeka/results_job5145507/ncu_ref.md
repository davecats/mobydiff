# What limits `jacobi_apply`

## refined_yp82_rect_jacobi  (60.56 Mcell)

| kernel | launches | us | doubles/cell | min | vs min | ld sect/req | st sect/req | L2 hit % | DRAM %peak | SM %peak | occupancy % | registers |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `jacobi_compute_phi__F1L499_10` | 6 | 3,207 | 5.54 | 5 | **1.11x** | 2.49 | 9.00 | 51.62 | 53.80 | 47.08 | 29.14 | 94.00 |
| `jacobi_apply__F1L550_14` | 6 | 1,541 | 3.21 | 3 | **1.07x** | 1.90 | 9.00 | 49.90 | 64.86 | 66.22 | 45.47 | 59.00 |
| `jacobi_apply__F1L585_16` | 6 | 5,352 | 10.68 | 8 | **1.33x** | 2.36 | 8.97 | 51.59 | 62.14 | 40.59 | 46.05 | 64.00 |
| `compute_rdenom__F1L447_6` | 1 | 4,608 | 4.27 | 4 | **1.07x** | 2.36 | 9.00 | 54.98 | 28.85 | 41.72 | 23.73 | 110.00 |

