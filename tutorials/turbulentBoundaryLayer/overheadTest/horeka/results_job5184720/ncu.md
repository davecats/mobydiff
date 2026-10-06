# What limits `jacobi_apply`

## refined_yp82_rect_jacobi  (60.56 Mcell)

| kernel | launches | us | doubles/cell | min | vs min | ld sect/req | st sect/req | L2 hit % | DRAM %peak | SM %peak | occupancy % | registers |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `jacobi_apply__F1L707_16` | 17 | 5,308 | 10.31 | 10 | **1.03x** | 2.36 | 8.97 | 50.85 | 60.50 | 40.54 | 46.04 | 64.00 |
| `interface_correct__F1L814_22` | 17 | 24 | 0.00 | - | - | 1.00 | 0.00 | 98.43 | 0.06 | 54.36 | 48.30 | 54.00 |
| `interface_correct__F1L841_24` | 17 | 58 | 0.08 | - | - | 1.74 | 9.00 | 30.35 | 42.85 | 46.09 | 50.65 | 54.00 |
| `interface_correct__F1L868_26` | 17 | 30 | 0.00 | - | - | 1.00 | 0.00 | 94.48 | 0.05 | 59.33 | 48.84 | 54.00 |
| `jacobi_compute_phi__F1L620_10` | 16 | 2,733 | 5.30 | 5 | **1.06x** | 2.89 | 9.00 | 50.78 | 60.44 | 52.44 | 34.98 | 80.00 |
| `nvkernel_comm_copy_local_scalar_same_level__F1L2300_26` | 16 | 307 | 0.39 | - | - | 1.62 | 17.50 | 66.29 | 39.97 | 25.53 | 46.97 | 64.00 |
| `nvkernel_comm_copy_local_scalar_entries__F1L2224_22` | 16 | 81 | 0.03 | - | - | 1.37 | 10.20 | 56.33 | 11.66 | 37.47 | 22.95 | 126.00 |
| `nvkernel_boundary_apply_scalar_bc__F1L1009_10` | 16 | 18 | 0.02 | - | - | 1.91 | 11.23 | 47.92 | 38.45 | 36.96 | 67.18 | 31.00 |
| `jacobi_apply__F1L671_14` | 16 | 1,535 | 3.08 | 3 | **1.03x** | 1.90 | 9.00 | 49.09 | 62.45 | 66.22 | 45.45 | 59.00 |
| `nvkernel_comm_copy_local_div_entries__F1L1864_2` | 14 | 210 | 0.28 | - | - | 1.45 | 17.44 | 55.93 | 41.18 | 21.20 | 41.41 | 72.00 |
| `nvkernel_comm_copy_local_same_level__F1L2053_14` | 5 | 980 | 1.62 | - | - | 2.38 | 17.50 | 65.07 | 51.27 | 12.00 | 42.13 | 72.00 |
| `nvkernel_comm_copy_local_cross_level__F1L2113_18` | 5 | 345 | 0.11 | - | - | 1.22 | 11.01 | 70.32 | 9.74 | 37.19 | 23.79 | 124.00 |
| `nvkernel_boundary_apply_bc__F1L827_2` | 3 | 137 | 0.27 | - | - | 2.93 | 9.88 | 59.19 | 61.53 | 26.76 | 34.67 | 78.00 |
| `nvkernel_bodyforce_fill_trip_kernel__F1L540_10` | 2 | 124 | 0.01 | - | - | 1.20 | 7.23 | 100.47 | 3.02 | 57.75 | 34.72 | 74.00 |
| `nvkernel_step_momentum__F1L220_6` | 2 | 15,217 | 16.84 | 16 | **1.05x** | 3.62 | 8.67 | 68.17 | 34.48 | 40.36 | 23.77 | 128.00 |
| `nvkernel_step_add_bodyforce_correction__F1L729_22` | 2 | 6,949 | 18.20 | - | - | 2.11 | 8.50 | 47.30 | 81.58 | 34.35 | 53.92 | 56.00 |
| `nvkernel_step_predict_outlet_faces__F1L588_14` | 2 | 70 | 0.07 | - | - | 2.92 | 11.23 | 50.46 | 31.85 | 18.11 | 16.58 | 130.00 |
| `nvkernel_step_momentum__F1L457_8` | 2 | 2,657 | 6.43 | 6 | **1.07x** | 1.48 | 9.00 | 63.77 | 75.40 | 62.70 | 58.72 | 48.00 |
| `nvkernel_step_get_timestep_rates__F1L1008_30` | 1 | 3,004 | 3.38 | - | - | 2.11 | 0.00 | 26.62 | 35.00 | 62.89 | 35.52 | 74.00 |
| `nvkernel_boundarylayer_stats_divergence_reductions__F1L351_6` | 1 | 5,114 | 3.15 | - | - | 2.92 | 0.00 | 45.23 | 19.21 | 46.09 | 30.21 | 90.00 |

_ncu_rough_jacobi: no header row -- see ncu.log_

