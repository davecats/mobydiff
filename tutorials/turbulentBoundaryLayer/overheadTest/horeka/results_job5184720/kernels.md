# Device time per kernel, per step

## base_jacobi_r1N1  (25.0 steps in the window; device busy 546563.05 us/step = 95.2 % of the 574.09 ms step (nsys run))

| kernel@line | launches/step | us/step | share of device | us/launch |
|---|---|---|---|---|
| `pressure_solver_jacobi_apply@707` | 18.0 | 204,115.7 | 37.3 % | 11339.8 |
| `pressure_solver_jacobi_compute_phi@620` | 18.0 | 104,254.9 | 19.1 % | 5791.9 |
| `step_momentum@220` | 3.0 | 101,169.3 | 18.5 % | 33723.1 |
| `pressure_solver_jacobi_apply@671` | 18.0 | 54,824.1 | 10.0 % | 3045.8 |
| `step_add_bodyforce_correction@729` | 3.0 | 46,334.8 | 8.5 % | 15444.9 |
| `step_momentum@457` | 3.0 | 15,774.0 | 2.9 % | 5258.0 |
| `bodyforce_fill_trip_kernel@540` | 3.0 | 10,971.6 | 2.0 % | 3657.2 |
| `step_get_timestep_rates@1008` | 1.0 | 5,329.7 | 1.0 % | 5329.7 |
| `boundary_apply_scalar_bc@1009` | 18.0 | 760.5 | 0.1 % | 42.2 |
| `comm_copy_local_scalar_same_level@2300` | 18.0 | 662.1 | 0.1 % | 36.8 |
| `boundary_apply_bc@827` | 3.0 | 552.6 | 0.1 % | 184.2 |
| `step_predict_outlet_faces@588` | 3.0 | 507.0 | 0.1 % | 169.0 |
| `comm_copy_local_same_level@2053` | 6.0 | 497.3 | 0.1 % | 82.9 |
| `comm_copy_local_div_entries@1864` | 15.0 | 409.3 | 0.1 % | 27.3 |
| `boundarylayer_stats_divergence_reductions@351` | 0.0 | 400.3 | 0.1 % | 10006.5 |

## min_channel_r1N1  (28.7 steps in the window; device busy 11557.83 us/step = 54.6 % of the 21.18 ms step (nsys run))

| kernel@line | launches/step | us/step | share of device | us/launch |
|---|---|---|---|---|
| `pressure_solver_jacobi_apply@707` | 36.0 | 3,463.7 | 30.0 % | 96.2 |
| `pressure_solver_jacobi_compute_phi@620` | 36.0 | 1,664.7 | 14.4 % | 46.2 |
| `comm_copy_local_scalar_same_level@2300` | 36.0 | 1,189.5 | 10.3 % | 33.0 |
| `comm_copy_local_div_entries@1864` | 33.0 | 906.3 | 7.8 % | 27.5 |
| `pressure_solver_cheb_combine@432` | 36.0 | 857.1 | 7.4 % | 23.8 |
| `pressure_solver_jacobi_apply@671` | 36.0 | 803.3 | 7.0 % | 22.3 |
| `comm_copy_local_same_level@2053` | 6.0 | 669.4 | 5.8 % | 110.9 |
| `step_momentum@220` | 3.0 | 643.8 | 5.6 % | 214.6 |
| `comm_copy_local_scalar_entries@2224` | 36.0 | 306.1 | 2.6 % | 8.5 |
| `pressure_solver_interface_correct@841` | 36.0 | 303.6 | 2.6 % | 8.4 |
| `pressure_solver_interface_correct@814` | 36.0 | 186.1 | 1.6 % | 5.2 |
| `step_momentum@457` | 3.0 | 173.7 | 1.5 % | 57.9 |
| `pressure_solver_interface_correct@868` | 36.0 | 154.4 | 1.3 % | 4.3 |
| `comm_copy_local_cross_level@2113` | 6.0 | 118.2 | 1.0 % | 19.6 |
| `step_get_timestep_rates@1008` | 1.0 | 74.3 | 0.6 % | 73.5 |
| `boundary_apply_bc@827` | 3.0 | 32.3 | 0.3 % | 10.7 |
| `ibmm_ibm_body_blocks@1329` | 0.1 | 8.7 | 0.1 % | 124.9 |
| `pressure_solver_compute_rdenom@569` | 0.0 | 2.5 | 0.0 % | 72.6 |

## refined_yp82_rect_jacobi_r1N1  (25.7 steps in the window; device busy 242715.50 us/step = 93.2 % of the 260.49 ms step (nsys run))

| kernel@line | launches/step | us/step | share of device | us/launch |
|---|---|---|---|---|
| `pressure_solver_jacobi_apply@707` | 18.0 | 86,432.8 | 35.6 % | 4801.8 |
| `pressure_solver_jacobi_compute_phi@620` | 18.0 | 43,633.1 | 18.0 % | 2424.1 |
| `step_momentum@220` | 3.0 | 38,584.0 | 15.9 % | 13030.6 |
| `pressure_solver_jacobi_apply@671` | 18.0 | 23,877.1 | 9.8 % | 1326.5 |
| `step_add_bodyforce_correction@729` | 3.0 | 20,553.4 | 8.5 % | 6851.1 |
| `step_momentum@457` | 3.0 | 7,478.4 | 3.1 % | 2492.8 |
| `comm_copy_local_same_level@2053` | 6.0 | 5,444.6 | 2.2 % | 907.4 |
| `comm_copy_local_scalar_same_level@2300` | 18.0 | 4,969.8 | 2.0 % | 276.1 |
| `comm_copy_local_div_entries@1864` | 15.0 | 2,900.5 | 1.2 % | 193.4 |
| `step_get_timestep_rates@1008` | 1.0 | 2,576.0 | 1.1 % | 2543.0 |
| `comm_copy_local_cross_level@2113` | 6.0 | 1,755.5 | 0.7 % | 292.6 |
| `comm_copy_local_scalar_entries@2224` | 18.0 | 1,376.0 | 0.6 % | 76.4 |
| `pressure_solver_interface_correct@841` | 18.0 | 959.1 | 0.4 % | 53.3 |
| `pressure_solver_interface_correct@868` | 18.0 | 403.2 | 0.2 % | 22.4 |
| `boundary_apply_bc@827` | 3.0 | 388.8 | 0.2 % | 129.6 |
| `boundary_apply_scalar_bc@1009` | 18.0 | 354.1 | 0.1 % | 19.7 |
| `pressure_solver_interface_correct@814` | 18.0 | 348.5 | 0.1 % | 19.4 |
| `bodyforce_fill_trip_kernel@540` | 3.0 | 309.5 | 0.1 % | 104.5 |
| `step_predict_outlet_faces@588` | 3.0 | 203.9 | 0.1 % | 68.0 |
| `boundarylayer_stats_divergence_reductions@351` | 0.0 | 167.3 | 0.1 % | 4294.3 |

## refined_yp82_rect_redblack_r1N1  (26.3 steps in the window; device busy 213635.16 us/step = 95.0 % of the 224.84 ms step (nsys run))

| kernel@line | launches/step | us/step | share of device | us/launch |
|---|---|---|---|---|
| `pressure_solver_redblack_sweep@1116` | 18.0 | 97,692.0 | 45.7 % | 5427.3 |
| `step_momentum@220` | 3.0 | 38,099.6 | 17.8 % | 12835.5 |
| `step_add_bodyforce_correction@729` | 3.0 | 20,327.8 | 9.5 % | 6848.3 |
| `pressure_solver_zero_phi@1061` | 18.0 | 19,822.2 | 9.3 % | 1101.2 |
| `comm_copy_local_same_level@2053` | 21.0 | 16,905.2 | 7.9 % | 804.8 |
| `step_momentum@457` | 3.0 | 7,387.1 | 3.5 % | 2488.7 |
| `comm_copy_local_scalar_same_level@2300` | 18.0 | 5,059.3 | 2.4 % | 280.5 |
| `step_get_timestep_rates@1008` | 1.0 | 2,562.1 | 1.2 % | 2493.6 |
| `comm_copy_local_cross_level@2113` | 6.0 | 1,722.5 | 0.8 % | 288.3 |
| `comm_copy_local_scalar_entries@2224` | 18.0 | 1,359.5 | 0.6 % | 75.4 |
| `pressure_solver_interface_correct@841` | 18.0 | 904.4 | 0.4 % | 50.1 |
| `pressure_solver_interface_correct@868` | 18.0 | 394.7 | 0.2 % | 21.9 |
| `boundary_apply_bc@827` | 3.0 | 391.5 | 0.2 % | 130.2 |
| `pressure_solver_interface_correct@814` | 18.0 | 343.1 | 0.2 % | 19.0 |
| `bodyforce_fill_trip_kernel@540` | 3.0 | 304.0 | 0.1 % | 102.4 |
| `step_predict_outlet_faces@588` | 3.0 | 199.7 | 0.1 % | 67.3 |
| `boundarylayer_stats_divergence_reductions@351` | 0.0 | 160.3 | 0.1 % | 4213.0 |

## rough_jacobi_r1N1  (28.7 steps in the window; device busy 567985.79 us/step = 94.5 % of the 600.84 ms step (nsys run))

| kernel@line | launches/step | us/step | share of device | us/launch |
|---|---|---|---|---|
| `pressure_solver_jacobi_apply@707` | 18.0 | 198,555.4 | 35.0 % | 11030.9 |
| `pressure_solver_jacobi_compute_phi@620` | 18.0 | 99,646.9 | 17.5 % | 5535.9 |
| `step_momentum@220` | 3.0 | 88,508.5 | 15.6 % | 29849.9 |
| `pressure_solver_jacobi_apply@671` | 18.0 | 54,593.0 | 9.6 % | 3032.9 |
| `step_add_bodyforce_correction@729` | 3.0 | 46,543.9 | 8.2 % | 15697.1 |
| `step_momentum@457` | 3.0 | 17,142.1 | 3.0 % | 5714.0 |
| `comm_copy_local_same_level@2053` | 6.0 | 13,633.7 | 2.4 % | 2272.3 |
| `comm_copy_local_scalar_same_level@2300` | 18.0 | 12,017.8 | 2.1 % | 667.7 |
| `ibmm_update_ibm_mu@1273` | 3.0 | 9,670.7 | 1.7 % | 3261.5 |
| `comm_copy_local_div_entries@1864` | 15.0 | 7,020.7 | 1.2 % | 468.0 |
| `step_add_penalization_state_correction@677` | 3.0 | 6,503.6 | 1.1 % | 2193.4 |
| `step_get_timestep_rates@1008` | 1.0 | 5,905.7 | 1.0 % | 5837.8 |
| `pressure_solver_compute_rdenom@569` | 3.0 | 5,611.6 | 1.0 % | 1870.5 |
| `boundary_apply_scalar_bc@1009` | 18.0 | 754.1 | 0.1 % | 41.9 |
| `boundary_apply_bc@827` | 3.0 | 593.5 | 0.1 % | 197.8 |
| `step_predict_outlet_faces@588` | 3.0 | 480.0 | 0.1 % | 160.0 |
| `bodyforce_fill_trip_kernel@540` | 3.0 | 458.7 | 0.1 % | 154.7 |
| `boundarylayer_stats_divergence_reductions@351` | 0.0 | 346.2 | 0.1 % | 9923.0 |

