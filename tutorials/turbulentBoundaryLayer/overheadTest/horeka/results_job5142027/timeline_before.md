# What the per-launch fixed cost is made of


## rect_jacobi_r4N1

Steady window 5279 ms = 23.1 steps; GPU busy 81.4 % of it.

### Per kernel: what the host pays, and what the device does

| kernel | launches/step | device us/launch | H2D copies before each | bucket | host us/launch (untraced) | host − device |
|---|---|---|---|---|---|---|
| `comm_pack_entries` | 21.0 | **18.9** | 14.0 | `pack` | 99.8 | 80.9 |
| `comm_pack_scalar_entries` | 18.0 | **9.6** | 14.0 | `pack` | 99.8 | 90.2 |
| `comm_unpack_entries` | 21.0 | **22.2** | 21.0 | `unpack` | 91.4 | 69.3 |
| `comm_unpack_scalar_entries` | 18.0 | **6.7** | 16.0 | `unpack` | 91.4 | 84.7 |
| `comm_copy_local_same_level` | 21.0 | **504.6** | 21.0 | `local_copy` | 425.0 | -79.6 |
| `comm_copy_local_scalar_same_level` | 18.0 | **198.7** | 21.0 | `local_copy` | 425.0 | 226.3 |
| `pressure_solver_jacobi_compute_phi` | 18.0 | **1567.1** | 11.8 | `sweep` | 1550.4 | -16.7 |
| `pressure_solver_jacobi_apply` | 36.0 | **2232.1** | 2.0 | `apply` | 4471.8 | 2239.8 |

(The `bucket` column repeats for the velocity and scalar variants that
share one profiler bucket; its host figure is their average, so compare
it against the two device figures together, not against either alone.)

### Host-to-device copies

- **2737 per step**, 0.34 kB each on average (0.94 MB/step in total)
- 1.99 us each on the device, 5.45 ms/step of device time

### CUDA API calls per step (counts exact; DURATIONS INFLATED BY TRACING)

| API | calls/step | avg us (inflated) | share of API time |
|---|---|---|---|
| `cuStreamSynchronize` | 299 | 743.75 | 70.6 % |
| `cuMemcpyHtoDAsync_v2` | 3569 | 13.51 | 15.3 % |
| `cuMemcpyDtoHAsync_v2` | 154 | 252.70 | 12.3 % |
| `cuLaunchKernel` | 298 | 4.76 | 0.5 % |
| `cuMemsetD32Async` | 691 | 1.95 | 0.4 % |
| `cuMemHostAlloc` | 0 | 21691.12 | 0.3 % |
| `cuEventQuery` | 1642 | 0.44 | 0.2 % |
| `cuMemAlloc_v2` | 5 | 57.16 | 0.1 % |
| `cuMemcpyDtoDAsync_v2` | 51 | 3.62 | 0.1 % |
| `cuMemRetainAllocationHandle` | 110 | 1.65 | 0.1 % |
| `cuEventRecord` | 51 | 2.48 | 0.0 % |
| `cudaGetDeviceProperties_v2_v12000` | 0 | 2806.23 | 0.0 % |


## refined_yp82_rect_jacobi_r4N1

Steady window 2705 ms = 22.3 steps; GPU busy 68.7 % of it.

### Per kernel: what the host pays, and what the device does

| kernel | launches/step | device us/launch | H2D copies before each | bucket | host us/launch (untraced) | host − device |
|---|---|---|---|---|---|---|
| `comm_pack_entries` | 21.0 | **15.3** | 14.0 | `pack` | 101.5 | 86.2 |
| `comm_pack_scalar_entries` | 18.0 | **9.6** | 14.0 | `pack` | 101.5 | 92.0 |
| `comm_unpack_entries` | 21.0 | **11.8** | 23.0 | `unpack` | 87.8 | 76.0 |
| `comm_unpack_scalar_entries` | 18.0 | **6.3** | 20.0 | `unpack` | 87.8 | 81.5 |
| `comm_copy_local_same_level` | 21.0 | **195.5** | 21.0 | `local_copy` | 199.5 | 4.0 |
| `comm_copy_local_scalar_same_level` | 18.0 | **67.9** | 21.0 | `local_copy` | 199.5 | 131.6 |
| `comm_copy_local_cross_level` | 6.0 | **74.6** | 22.0 | `copy_cross` | 116.9 | 42.3 |
| `comm_copy_local_scalar_entries` | 18.0 | **22.3** | 19.0 | `copy_cross` | 116.9 | 94.7 |
| `pressure_solver_jacobi_compute_phi` | 18.0 | **688.9** | 11.8 | `sweep` | 697.6 | 8.7 |
| `pressure_solver_jacobi_apply` | 36.0 | **977.3** | 2.0 | `apply` | 2095.0 | 1117.7 |
| `pressure_solver_interface_correct` | 54.0 | **10.8** | 3.0 | `apply` | 2095.0 | 2084.2 |

(The `bucket` column repeats for the velocity and scalar variants that
share one profiler bucket; its host figure is their average, so compare
it against the two device figures together, not against either alone.)

### Host-to-device copies

- **3484 per step**, 0.32 kB each on average (1.13 MB/step in total)
- 1.97 us each on the device, 6.85 ms/step of device time

### CUDA API calls per step (counts exact; DURATIONS INFLATED BY TRACING)

| API | calls/step | avg us (inflated) | share of API time |
|---|---|---|---|
| `cuStreamSynchronize` | 414 | 250.04 | 64.3 % |
| `cuMemcpyHtoDAsync_v2` | 4693 | 6.30 | 18.4 % |
| `cuMemcpyDtoHAsync_v2` | 191 | 113.28 | 13.5 % |
| `cuMemsetD32Async` | 1036 | 1.83 | 1.2 % |
| `cuLaunchKernel` | 413 | 4.17 | 1.1 % |
| `cuMemHostAlloc` | 0 | 20802.69 | 0.6 % |
| `cuEventQuery` | 1067 | 0.45 | 0.3 % |
| `cuMemAlloc_v2` | 5 | 61.22 | 0.2 % |
| `cuMemRetainAllocationHandle` | 115 | 1.51 | 0.1 % |
| `cuMemcpyDtoDAsync_v2` | 52 | 3.23 | 0.1 % |
| `cuModuleLoadDataEx` | 0 | 2862.81 | 0.1 % |
| `cuEventRecord` | 52 | 2.23 | 0.1 % |

