
## /hkfs/work/workspace/scratch/xt8786-roughBoundaryLayer/optimiseBlockRefinement/moby-2to1-run/results_ab/after/nsys_rect_jacobi_r4N1/rep_0.sqlite

Steady window = 29.4 steps

| kernel | launches/step | H2D copies/launch | bytes/launch | size histogram per launch |
|---|---|---|---|---|
| `pressure_solver_jacobi_apply` | 36.0 | 2.0 | 9 | 3 B x2, 8 B x0 |
| `boundary_apply_bc` | 21.0 | 1.0 | 12 | 12 B x1 |
| `comm_pack_entries` | 21.0 | 1.0 | 16 | 16 B x1 |
| `comm_copy_local_same_level` | 21.0 | 1.0 | 16 | 16 B x1 |
| `comm_unpack_entries` | 21.0 | 1.0 | 16 | 16 B x1 |
| `pressure_solver_jacobi_compute_phi` | 18.0 | 1.0 | 8 | 8 B x1 |
| `comm_pack_scalar_entries` | 18.0 | 0.0 | 0 |  |
| `comm_copy_local_scalar_same_level` | 18.0 | 0.0 | 0 |  |
| `comm_unpack_scalar_entries` | 18.0 | 0.0 | 0 |  |
| `boundary_apply_scalar_bc` | 18.0 | 3.0 | 84 | 12 B x1, 24 B x1, 48 B x1 |
| `step_momentum` | 6.0 | 2.3 | 27 | 8 B x2, 24 B x0 |
| `ibmm_update_ibm_mu` | 3.0 | 1.0 | 8 | 8 B x1 |
| `bodyforce_fill_trip_kernel` | 3.0 | 10.0 | 816 | 8 B x6, 192 B x4 |
| `step_add_bodyforce_correction` | 3.0 | 1.0 | 8 | 8 B x1 |
| `pressure_solver_compute_rdenom` | 3.0 | 3.0 | 9 | 3 B x3 |
| `step_get_timestep_rates` | 1.0 | 0.0 | 0 |  |

## /hkfs/work/workspace/scratch/xt8786-roughBoundaryLayer/optimiseBlockRefinement/moby-2to1-run/results_ab/after/nsys_refined_yp82_rect_jacobi_r4N1/rep_0.sqlite

Steady window = 30.0 steps

| kernel | launches/step | H2D copies/launch | bytes/launch | size histogram per launch |
|---|---|---|---|---|
| `pressure_solver_interface_correct` | 54.0 | 3.0 | 9 | 3 B x3 |
| `pressure_solver_jacobi_apply` | 36.0 | 2.0 | 8 | 3 B x2, 8 B x0 |
| `boundary_apply_bc` | 21.0 | 1.0 | 575,817 | 12 B x1, 344 B x0, 1,088 B x0, 363,331,584 B x0 |
| `comm_pack_entries` | 21.0 | 1.0 | 16 | 16 B x1 |
| `comm_copy_local_same_level` | 21.0 | 1.0 | 16 | 16 B x1 |
| `comm_unpack_entries` | 21.0 | 1.0 | 16 | 16 B x1 |
| `pressure_solver_jacobi_compute_phi` | 18.0 | 1.0 | 8 | 8 B x1 |
| `comm_pack_scalar_entries` | 18.0 | 0.0 | 0 |  |
| `comm_copy_local_scalar_same_level` | 18.0 | 0.0 | 0 |  |
| `comm_copy_local_scalar_entries` | 18.0 | 0.0 | 0 |  |
| `comm_unpack_scalar_entries` | 18.0 | 0.0 | 0 |  |
| `boundary_apply_scalar_bc` | 18.0 | 3.0 | 84 | 12 B x1, 24 B x1, 48 B x1 |
| `comm_copy_local_cross_level` | 6.0 | 1.0 | 16 | 16 B x1 |
| `step_momentum` | 6.0 | 2.3 | 27 | 8 B x2, 24 B x0 |
| `ibmm_update_ibm_mu` | 3.0 | 1.0 | 8 | 8 B x1 |
| `bodyforce_fill_trip_kernel` | 3.0 | 10.0 | 816 | 8 B x6, 192 B x4 |
| `step_add_bodyforce_correction` | 3.0 | 1.0 | 8 | 8 B x1 |
| `pressure_solver_compute_rdenom` | 3.0 | 3.0 | 3,022,522 | 3 B x3, 272 B x0, 136,012,800 B x0 |
| `step_get_timestep_rates` | 1.0 | 0.0 | 0 |  |
