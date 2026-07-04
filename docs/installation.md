# Installation

`mobydiff` is built with CMake through the `compile.sh` wrapper. MPI is always enabled; the
GPU path is a build-time option. The CPU build is the reference used for debugging and for
bit-exact verification of the GPU build.

## Dependencies

| Component | Requirement | Notes |
|-----------|-------------|-------|
| CMake | ≥ 3.20 | |
| CPU compiler | `gcc` / `gfortran` ≥ 13 | reference build |
| GPU compiler | NVIDIA HPC SDK (NVHPC), CUDA 12+ | OpenMP target offload |
| HDF5 | ≥ 1.14, **parallel (MPI) build** | field and restart I/O |
| MPI | any MPI; **GPU-aware** for GPU runs | 3D Cartesian decomposition |
| Python 3 | `numpy`, `h5py`, `matplotlib` | preprocessing + analysis tools |

Parallel HDF5 is required — the solver writes collective/independent hyperslabs from every
rank into shared datasets.

## Building

```bash
./compile.sh cpu     # builds build_cpu/ (gfortran, CPU reference)
./compile.sh gpu     # builds build_gpu/ (NVHPC, GPU offload)
```

Each invocation configures and builds into its own directory (`build_cpu/`, `build_gpu/`),
so the two toolchains never collide. Both directories contain `main` (the solver) and
`mobygrid` (the serial grid export tool).

### GPU build

Load the NVHPC toolchain before configuring so CMake picks up `nvfortran` and the CUDA
runtime:

```bash
module load /opt/nvidia/hpc_sdk/modulefiles/nvhpc-hpcx-cuda13/26.3
./compile.sh gpu
```

The GPU build requires a GPU-aware MPI stack so device buffers can be passed straight to
the halo-exchange messages.

### Pointing CMake at HDF5

If CMake cannot find a parallel HDF5, set `HDF5_ROOT`:

```bash
HDF5_ROOT=/path/to/parallel-hdf5 ./compile.sh cpu
```

## Bit-exact builds for verification

Pure refactors are verified by comparing output fields between two builds (see
[Validation](validation.md)). Floating-point contraction (FMA) makes arithmetically
identical source differ by 1–2 ulp, so for a bit-exact comparison build **both** sides with
FMA contraction disabled:

- CPU: `-Mnofma`
- GPU: `-Mnofma -gpu=nofma`

## Troubleshooting

- **HDF5 not found / serial HDF5 picked up** — ensure the HDF5 you point at is the
  MPI-parallel build, and set `HDF5_ROOT` explicitly.
- **GPU run fails in MPI halo exchange** — confirm the MPI stack is GPU-aware (the module
  above provides HPC-X with CUDA support).
- **WSL instability** — on WSL hosts, commands can hang or fail oddly after a system
  update; retry a hung command before assuming a code problem.
