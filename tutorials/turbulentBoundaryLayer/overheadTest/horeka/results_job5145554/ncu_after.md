# What limits `jacobi_apply`

## refined_yp82_rect_jacobi  (60.56 Mcell)

| kernel | launches | us | doubles/cell | min | vs min | ld sect/req | st sect/req | L2 hit % | DRAM %peak | SM %peak | occupancy % | registers |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `nvkernel_step_momentum__F1L154_6` | 3 | 15,221 | 16.87 | 16 | **1.05x** | 3.62 | 8.67 | 68.16 | 34.52 | 40.33 | 23.77 | 128.00 |
| `nvkernel_step_momentum__F1L366_8` | 3 | 2,801 | 6.67 | 6 | **1.11x** | 2.45 | 9.00 | 63.50 | 74.20 | 37.82 | 41.09 | 66.00 |

