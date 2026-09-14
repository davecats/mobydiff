# What limits `jacobi_apply`

## refined_yp82_rect_jacobi  (60.56 Mcell)

| kernel | launches | us | doubles/cell | min | vs min | ld sect/req | st sect/req | L2 hit % | DRAM %peak | SM %peak | occupancy % | registers |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `jacobi_compute_phi__F1L543_10` | 6 | 2,754 | 5.51 | 5 | **1.10x** | 2.89 | 9.00 | 51.39 | 62.32 | 52.44 | 35.01 | 80.00 |
| `jacobi_apply__F1L594_14` | 6 | 1,542 | 3.21 | 3 | **1.07x** | 1.90 | 9.00 | 49.91 | 64.84 | 66.22 | 45.47 | 59.00 |
| `jacobi_apply__F1L629_16` | 6 | 5,355 | 10.68 | 8 | **1.33x** | 2.36 | 8.97 | 51.59 | 62.10 | 40.59 | 46.05 | 64.00 |
| `compute_rdenom__F1L493_6` | 1 | 3,342 | 4.27 | 4 | **1.07x** | 2.77 | 9.00 | 55.74 | 39.78 | 49.58 | 35.01 | 80.00 |

