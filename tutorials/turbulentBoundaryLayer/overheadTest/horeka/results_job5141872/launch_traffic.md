# Host-to-device traffic per kernel launch

Produced by `horeka/exchange/analyse_launch_traffic.py` from job 5141872's nsys
sqlite exports. The full `.nsys-rep`, `.sqlite` and `cuda_gpu_trace.csv` (16 MB
each) are NOT committed; they live in `$RUN_DIR/results_timeline/` and can be
regenerated from the same job spec. Everything the report quotes is below.

## /hkfs/work/workspace/scratch/xt8786-roughBoundaryLayer/optimiseBlockRefinement/moby-2to1-run/results_timeline/nsys_rect_jacobi_r4N1/rep_0.sqlite

Steady window = 28.7 steps

| kernel | launches/step | H2D copies/launch | bytes/launch | size histogram per launch |
|---|---|---|---|---|
| `pressure_solver_jacobi_apply` | 36.0 | 2.0 | 8 | 3 B x2, 8 B x0 |
| `boundary_apply_bc` | 21.0 | 1.0 | 12 | 12 B x1 |
| `comm_pack_entries` | 21.0 | 14.0 | 10,264 | 152 B x6, 200 B x7, 7,952 B x1 |
| `comm_copy_local_same_level` | 21.0 | 21.0 | 9,264 | 8 B x13, 152 B x4, 200 B x3, 7,952 B x1 |
| `comm_unpack_entries` | 21.0 | 21.0 | 10,176 | 8 B x7, 152 B x9, 200 B x4, 7,952 B x1 |
| `pressure_solver_jacobi_compute_phi` | 18.0 | 11.8 | 95 | 8 B x12 |
| `comm_pack_scalar_entries` | 18.0 | 14.0 | 10,264 | 152 B x6, 200 B x7, 7,952 B x1 |
| `comm_copy_local_scalar_same_level` | 18.0 | 21.0 | 9,264 | 8 B x13, 152 B x4, 200 B x3, 7,952 B x1 |
| `comm_unpack_scalar_entries` | 18.0 | 16.0 | 9,368 | 8 B x7, 152 B x5, 200 B x3, 7,952 B x1 |
| `boundary_apply_scalar_bc` | 18.0 | 11.0 | 148 | 8 B x8, 12 B x1, 24 B x1, 48 B x1 |
| `step_momentum` | 6.0 | 2.3 | 27 | 8 B x2, 24 B x0 |
| `ibmm_update_ibm_mu` | 3.0 | 9.8 | 78 | 8 B x10 |
| `bodyforce_fill_trip_kernel` | 3.0 | 10.0 | 816 | 8 B x6, 192 B x4 |
| `step_add_bodyforce_correction` | 3.0 | 1.0 | 8 | 8 B x1 |
| `pressure_solver_compute_rdenom` | 3.0 | 16.0 | 113 | 3 B x3, 8 B x13 |
| `step_get_timestep_rates` | 1.0 | 13.0 | 104 | 8 B x13 |

## /hkfs/work/workspace/scratch/xt8786-roughBoundaryLayer/optimiseBlockRefinement/moby-2to1-run/results_timeline/nsys_refined_yp82_rect_jacobi_r4N1/rep_0.sqlite

Steady window = 30.0 steps

| kernel | launches/step | H2D copies/launch | bytes/launch | size histogram per launch |
|---|---|---|---|---|
| `pressure_solver_interface_correct` | 54.0 | 3.0 | 9 | 3 B x3 |
| `pressure_solver_jacobi_apply` | 36.0 | 2.0 | 8 | 3 B x2, 8 B x0 |
| `boundary_apply_bc` | 21.0 | 1.0 | 12 | 12 B x1 |
| `comm_pack_entries` | 21.0 | 14.0 | 10,264 | 152 B x6, 200 B x7, 7,952 B x1 |
| `comm_copy_local_same_level` | 21.0 | 21.0 | 9,264 | 8 B x13, 152 B x4, 200 B x3, 7,952 B x1 |
| `comm_unpack_entries` | 21.0 | 23.0 | 10,192 | 8 B x9, 152 B x9, 200 B x4, 7,952 B x1 |
| `pressure_solver_jacobi_compute_phi` | 18.0 | 11.8 | 95 | 8 B x12 |
| `comm_pack_scalar_entries` | 18.0 | 14.0 | 10,264 | 152 B x6, 200 B x7, 7,952 B x1 |
| `comm_copy_local_scalar_same_level` | 18.0 | 21.0 | 9,264 | 8 B x13, 152 B x4, 200 B x3, 7,952 B x1 |
| `comm_copy_local_scalar_entries` | 18.0 | 19.0 | 9,968 | 8 B x7, 152 B x5, 200 B x6, 7,952 B x1 |
| `comm_unpack_scalar_entries` | 18.0 | 20.0 | 9,400 | 8 B x11, 152 B x5, 200 B x3, 7,952 B x1 |
| `boundary_apply_scalar_bc` | 18.0 | 11.0 | 148 | 8 B x8, 12 B x1, 24 B x1, 48 B x1 |
| `comm_copy_local_cross_level` | 6.0 | 22.0 | 10,472 | 8 B x7, 152 B x7, 200 B x7, 7,952 B x1 |
| `step_momentum` | 5.9 | 2.3 | 27 | 8 B x2, 24 B x0 |
| `pressure_solver_compute_rdenom` | 3.0 | 15.9 | 3,022,625 | 3 B x3, 8 B x13, 272 B x0, 136,012,800 B x0 |
| `ibmm_update_ibm_mu` | 3.0 | 9.8 | 78 | 8 B x10 |
| `bodyforce_fill_trip_kernel` | 3.0 | 10.0 | 816 | 8 B x6, 192 B x4 |
| `step_add_bodyforce_correction` | 3.0 | 1.0 | 8 | 8 B x1 |
| `step_get_timestep_rates` | 1.0 | 13.0 | 104 | 8 B x13 |
