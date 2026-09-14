# What limits `jacobi_apply`

## refined_yp82_rect_jacobi  (60.56 Mcell)

| kernel | launches | us | doubles/cell | min | vs min | ld sect/req | st sect/req | L2 hit % | DRAM %peak | SM %peak | occupancy % | registers |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `nvkernel_step_momentum__F1L153_6` | 3 | 16,168 | 16.87 | 16 | **1.05x** | 2.68 | 8.67 | 68.40 | 32.50 | 40.61 | 23.96 | 128.00 |
| `nvkernel_step_momentum__F1L366_8` | 3 | 2,800 | 6.67 | 6 | **1.11x** | 2.45 | 9.00 | 63.50 | 74.22 | 37.81 | 41.09 | 66.00 |

