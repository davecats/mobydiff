# What the per-launch fixed cost is made of

```
job    : 5141872
date   : 2026-09-11 17:18:05 CEST
nodes  : 1  (hkn0402)
commit : f7a597a88927151e36f25ef1ddb0c5507a4fde34
dirty  : 0
gpu    : NVIDIA A100-SXM4-40GB
nsteps : 30
```


## rect_jacobi_r4N1

Steady window 5292 ms = 23.1 steps; GPU busy 81.2 % of it.

### Per kernel: what the host pays, and what the device does

| kernel | launches/step | device us/launch | H2D copies before each | bucket | host us/launch (untraced) | host − device |
|---|---|---|---|---|---|---|
| `comm_pack_entries` | 21.0 | **18.9** | 14.0 | `pack` | 102.9 | 84.0 |
| `comm_pack_scalar_entries` | 18.0 | **9.6** | 14.0 | `pack` | 102.9 | 93.3 |
| `comm_unpack_entries` | 21.0 | **22.2** | 21.0 | `unpack` | 92.8 | 70.6 |
| `comm_unpack_scalar_entries` | 18.0 | **6.8** | 16.0 | `unpack` | 92.8 | 86.1 |
| `comm_copy_local_same_level` | 21.0 | **504.8** | 21.0 | `local_copy` | 426.3 | -78.5 |
| `comm_copy_local_scalar_same_level` | 18.0 | **198.8** | 21.0 | `local_copy` | 426.3 | 227.4 |
| `pressure_solver_jacobi_compute_phi` | 18.0 | **1568.5** | 11.8 | `sweep` | 1553.4 | -15.1 |
| `pressure_solver_jacobi_apply` | 36.0 | **2234.3** | 2.0 | `apply` | 4478.2 | 2243.9 |

(The `bucket` column repeats for the velocity and scalar variants that
share one profiler bucket; its host figure is their average, so compare
it against the two device figures together, not against either alone.)

### Host-to-device copies

- **2736 per step**, 0.34 kB each on average (0.94 MB/step in total)
- 1.99 us each on the device, 5.45 ms/step of device time

### CUDA API calls per step (counts exact; DURATIONS INFLATED BY TRACING)

| API | calls/step | avg us (inflated) | share of API time |
|---|---|---|---|
| `cuStreamSynchronize` | 299 | 744.78 | 70.7 % |
| `cuMemcpyHtoDAsync_v2` | 3569 | 13.27 | 15.0 % |
| `cuMemcpyDtoHAsync_v2` | 154 | 257.53 | 12.6 % |
| `cuMemsetD32Async` | 691 | 1.92 | 0.4 % |
| `cuLaunchKernel` | 298 | 4.43 | 0.4 % |
| `cuMemHostAlloc` | 0 | 22670.09 | 0.3 % |
| `cuEventQuery` | 1756 | 0.42 | 0.2 % |
| `cuMemAlloc_v2` | 5 | 43.95 | 0.1 % |
| `cuMemcpyDtoDAsync_v2` | 51 | 3.51 | 0.1 % |
| `cuMemRetainAllocationHandle` | 110 | 1.58 | 0.1 % |
| `cuModuleLoadDataEx` | 0 | 2696.48 | 0.0 % |
| `cuEventRecord` | 51 | 2.30 | 0.0 % |


## refined_yp82_rect_jacobi_r4N1

Steady window 2685 ms = 22.7 steps; GPU busy 70.3 % of it.

### Per kernel: what the host pays, and what the device does

| kernel | launches/step | device us/launch | H2D copies before each | bucket | host us/launch (untraced) | host − device |
|---|---|---|---|---|---|---|
| `comm_pack_entries` | 21.0 | **15.3** | 14.0 | `pack` | 100.3 | 85.0 |
| `comm_pack_scalar_entries` | 18.0 | **9.5** | 14.0 | `pack` | 100.3 | 90.8 |
| `comm_unpack_entries` | 21.0 | **11.8** | 23.0 | `unpack` | 87.0 | 75.2 |
| `comm_unpack_scalar_entries` | 18.0 | **6.3** | 20.0 | `unpack` | 87.0 | 80.7 |
| `comm_copy_local_same_level` | 21.0 | **195.5** | 21.0 | `local_copy` | 199.0 | 3.5 |
| `comm_copy_local_scalar_same_level` | 18.0 | **68.0** | 21.0 | `local_copy` | 199.0 | 131.0 |
| `comm_copy_local_cross_level` | 6.0 | **74.6** | 22.0 | `copy_cross` | 116.7 | 42.1 |
| `comm_copy_local_scalar_entries` | 18.0 | **22.2** | 19.0 | `copy_cross` | 116.7 | 94.5 |
| `pressure_solver_jacobi_compute_phi` | 18.0 | **688.8** | 11.8 | `sweep` | 698.5 | 9.7 |
| `pressure_solver_jacobi_apply` | 36.0 | **977.2** | 2.0 | `apply` | 2095.8 | 1118.6 |
| `pressure_solver_interface_correct` | 54.0 | **10.8** | 3.0 | `apply` | 2095.8 | 2085.0 |

(The `bucket` column repeats for the velocity and scalar variants that
share one profiler bucket; its host figure is their average, so compare
it against the two device figures together, not against either alone.)

### Host-to-device copies

- **3484 per step**, 0.32 kB each on average (1.13 MB/step in total)
- 1.96 us each on the device, 6.84 ms/step of device time

### CUDA API calls per step (counts exact; DURATIONS INFLATED BY TRACING)

| API | calls/step | avg us (inflated) | share of API time |
|---|---|---|---|
| `cuStreamSynchronize` | 408 | 248.92 | 64.2 % |
| `cuMemcpyHtoDAsync_v2` | 4624 | 6.32 | 18.5 % |
| `cuMemcpyDtoHAsync_v2` | 188 | 111.64 | 13.3 % |
| `cuMemsetD32Async` | 1021 | 1.88 | 1.2 % |
| `cuLaunchKernel` | 407 | 4.31 | 1.1 % |
| `cuMemHostAlloc` | 0 | 21992.23 | 0.6 % |
| `cuEventQuery` | 1217 | 0.43 | 0.3 % |
| `cuMemRetainAllocationHandle` | 113 | 1.73 | 0.1 % |
| `cuMemcpyDtoDAsync_v2` | 52 | 3.48 | 0.1 % |
| `cuMemAlloc_v2` | 5 | 36.68 | 0.1 % |
| `cuEventRecord` | 52 | 2.47 | 0.1 % |
| `cuModuleLoadDataEx` | 0 | 2753.20 | 0.1 % |

