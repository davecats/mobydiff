# What the per-launch fixed cost is made of


## rect_jacobi_r4N1

Steady window 5054 ms = 23.4 steps; GPU busy 86.6 % of it.

### Per kernel: what the host pays, and what the device does

| kernel | launches/step | device us/launch | H2D copies before each | bucket | host us/launch (untraced) | host − device |
|---|---|---|---|---|---|---|
| `comm_pack_entries` | 21.0 | **19.5** | 1.0 | `pack` | 33.2 | 13.7 |
| `comm_pack_scalar_entries` | 18.0 | **10.4** | 0.0 | `pack` | 33.2 | 22.9 |
| `comm_unpack_entries` | 21.0 | **22.8** | 1.0 | `unpack` | 33.7 | 10.9 |
| `comm_unpack_scalar_entries` | 18.0 | **7.1** | 0.0 | `unpack` | 33.7 | 26.5 |
| `comm_copy_local_same_level` | 21.0 | **505.3** | 1.0 | `local_copy` | 379.5 | -125.8 |
| `comm_copy_local_scalar_same_level` | 18.0 | **199.5** | 0.0 | `local_copy` | 379.5 | 180.0 |
| `pressure_solver_jacobi_compute_phi` | 18.0 | **1571.3** | 1.0 | `sweep` | 1550.1 | -21.2 |
| `pressure_solver_jacobi_apply` | 36.0 | **2238.0** | 2.0 | `apply` | 4469.8 | 2231.8 |

(The `bucket` column repeats for the velocity and scalar variants that
share one profiler bucket; its host figure is their average, so compare
it against the two device figures together, not against either alone.)

### Host-to-device copies

- **287 per step**, 0.00 kB each on average (0.00 MB/step in total)
- 1.15 us each on the device, 0.33 ms/step of device time

### CUDA API calls per step (counts exact; DURATIONS INFLATED BY TRACING)

| API | calls/step | avg us (inflated) | share of API time |
|---|---|---|---|
| `cuStreamSynchronize` | 294 | 810.66 | 81.0 % |
| `cuMemcpyHtoDAsync_v2` | 378 | 102.01 | 13.1 % |
| `cuMemcpyDtoHAsync_v2` | 2 | 7962.25 | 4.2 % |
| `cuMemsetD32Async` | 680 | 1.96 | 0.5 % |
| `cuLaunchKernel` | 293 | 4.13 | 0.4 % |
| `cuMemHostAlloc` | 0 | 20801.54 | 0.3 % |
| `cuEventQuery` | 1665 | 0.41 | 0.2 % |
| `cuMemcpyDtoDAsync_v2` | 50 | 4.45 | 0.1 % |
| `cuMemAlloc_v2` | 4 | 46.52 | 0.1 % |
| `cuMemRetainAllocationHandle` | 108 | 1.63 | 0.1 % |
| `cuModuleLoadDataEx` | 0 | 2825.59 | 0.0 % |
| `cuEventRecord` | 50 | 2.32 | 0.0 % |


## refined_yp82_rect_jacobi_r4N1

Steady window 2397 ms = 23.4 steps; GPU busy 81.6 % of it.

### Per kernel: what the host pays, and what the device does

| kernel | launches/step | device us/launch | H2D copies before each | bucket | host us/launch (untraced) | host − device |
|---|---|---|---|---|---|---|
| `comm_pack_entries` | 21.0 | **15.9** | 1.0 | `pack` | 31.2 | 15.3 |
| `comm_pack_scalar_entries` | 18.0 | **10.0** | 0.0 | `pack` | 31.2 | 21.2 |
| `comm_unpack_entries` | 21.0 | **12.4** | 1.0 | `unpack` | 27.7 | 15.3 |
| `comm_unpack_scalar_entries` | 18.0 | **6.5** | 0.0 | `unpack` | 27.7 | 21.2 |
| `comm_copy_local_same_level` | 21.0 | **195.9** | 1.0 | `local_copy` | 151.7 | -44.2 |
| `comm_copy_local_scalar_same_level` | 18.0 | **68.3** | 0.0 | `local_copy` | 151.7 | 83.4 |
| `comm_copy_local_cross_level` | 6.0 | **75.0** | 1.0 | `copy_cross` | 53.8 | -21.2 |
| `comm_copy_local_scalar_entries` | 18.0 | **22.4** | 0.0 | `copy_cross` | 53.8 | 31.4 |
| `pressure_solver_jacobi_compute_phi` | 18.0 | **688.9** | 1.0 | `sweep` | 692.1 | 3.3 |
| `pressure_solver_jacobi_apply` | 36.1 | **977.4** | 2.0 | `apply` | 2078.8 | 1101.4 |
| `pressure_solver_interface_correct` | 54.1 | **11.0** | 3.0 | `apply` | 2078.8 | 2067.8 |

(The `bucket` column repeats for the velocity and scalar variants that
share one profiler bucket; its host figure is their average, so compare
it against the two device figures together, not against either alone.)

### Host-to-device copies

- **456 per step**, 0.00 kB each on average (0.00 MB/step in total)
- 1.15 us each on the device, 0.53 ms/step of device time

### CUDA API calls per step (counts exact; DURATIONS INFLATED BY TRACING)

| API | calls/step | avg us (inflated) | share of API time |
|---|---|---|---|
| `cuStreamSynchronize` | 395 | 272.77 | 80.0 % |
| `cuMemcpyHtoDAsync_v2` | 594 | 30.14 | 13.3 % |
| `cuMemcpyDtoHAsync_v2` | 2 | 2038.61 | 2.3 % |
| `cuMemsetD32Async` | 990 | 1.85 | 1.4 % |
| `cuLaunchKernel` | 394 | 3.79 | 1.1 % |
| `cuMemHostAlloc` | 0 | 20632.79 | 0.7 % |
| `cuEventQuery` | 1114 | 0.43 | 0.4 % |
| `cuMemRetainAllocationHandle` | 109 | 2.69 | 0.2 % |
| `cuMemcpyDtoDAsync_v2` | 50 | 4.10 | 0.2 % |
| `cuMemAlloc_v2` | 4 | 32.41 | 0.1 % |
| `cuModuleLoadDataEx` | 0 | 2902.25 | 0.1 % |
| `cudaGetDeviceProperties_v2_v12000` | 0 | 2589.82 | 0.1 % |

