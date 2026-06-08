# Cluster Profiling Runs

This directory contains the cluster profiling harness for the MPI + OpenMP-offload red-black solver. GPU builds assume GPU-aware MPI; host-staged MPI is intentionally not supported for performance runs.

## Required Modules

Use the NVIDIA HPC SDK module that includes the GPU-aware MPI stack:

```bash
module purge
module load toolkit/nvidia-hpc-sdk/25.3
export HDF5_ROOT=/hkfs/work/workspace/scratch/xt8786-mobydiff/mobydiff/lib
```

The helper script used by the SLURM jobs does the same:

```bash
source profiling/scripts/load_hpc_modules.sh
```

The profiling scripts use `mpirun` by default. On this cluster, direct `srun` launch of the HPC SDK OpenMPI stack failed during `MPI_Init`, while `mpirun` correctly starts the GPU-aware runtime.

If an MPI stack needs extra MCA/UCX flags, pass them with `MPI_EXTRA_ARGS`, for example:

```bash
MPI_EXTRA_ARGS="--mca pml ucx --mca coll_hcoll_enable 0" profiling/scripts/run_case.sh
```

## Build

From the repository root:

```bash
source profiling/scripts/load_hpc_modules.sh
./compile.sh gpu
```

The executable is written to:

```text
build_gpu/main
```

The CPU validation build, if needed, is:

```bash
./compile.sh cpu
```

The only build modes are now `cpu` and `gpu`; MPI is always enabled.

## Static Inputs

Two profiling case directories are provided:

```text
profiling/cases/700/input_1gpu.ini
profiling/cases/700/input_2gpu.ini
profiling/cases/700/input_4gpu.ini
profiling/cases/700/input_8gpu.ini
profiling/cases/1024/input_4gpu.ini
profiling/cases/1024/input_8gpu.ini
profiling/cases/1024/input_16gpu.ini
profiling/cases/1024/input_32gpu.ini
profiling/cases/1024/input_64gpu.ini
```

The generated SLURM runs usually create fresh inputs in their result directories, but these files are useful for manual runs and quick inspection.

## SLURM Runs

### 700^3 Strong Scaling

This case fits on one A100-40 and is intended for 1, 2, 4, and 8 GPUs:

```bash
sbatch profiling/slurm/strong_700_gpu.slurm
```

Override defaults when useful:

```bash
sbatch --nodes=2 --export=ALL,NSTEPS=50,NITER=3,RANKS_LIST="1 2 4 8" profiling/slurm/strong_700_gpu.slurm
```

### 1024^3 Strong Scaling

This case is intended for multi-GPU runs. The default script requests up to sixteen nodes and runs 4, 8, 16, 32, and 64 ranks:

```bash
sbatch profiling/slurm/strong_1024_gpu.slurm
```

Useful overrides:

```bash
sbatch --nodes=1 --export=ALL,RANKS_LIST="4",NSTEPS=20 profiling/slurm/strong_1024_gpu.slurm
sbatch --nodes=16 --export=ALL,RANKS_LIST="4 8 16 32 64",NSTEPS=20 profiling/slurm/strong_1024_gpu.slurm
```

### Nsight Systems

Profile one selected point after timing is stable:

```bash
sbatch --export=ALL,SIZE=700,RANKS=4,NSTEPS=5 profiling/slurm/nsys_gpu.slurm
sbatch --nodes=4 --export=ALL,SIZE=1024,RANKS=16,NSTEPS=5 profiling/slurm/nsys_gpu.slurm
```

## Interactive Examples

After an interactive allocation:

```bash
source profiling/scripts/load_hpc_modules.sh
./compile.sh gpu
```

The manual examples below use `mpirun` internally and map ranks as `ppr:TASKS_PER_NODE:node`.

One GPU, 700^3:

```bash
NX=700 NY=700 NZ=700 DIMS="1 1 1" NSTEPS=20 NITER=3 RANKS=1 NODES=1 TASKS_PER_NODE=1 CASE_NAME=manual_700_1gpu profiling/scripts/run_case.sh
```

Four GPUs, 700^3:

```bash
NX=700 NY=700 NZ=700 DIMS="2 2 1" NSTEPS=20 NITER=3 RANKS=4 NODES=1 TASKS_PER_NODE=4 CASE_NAME=manual_700_4gpu profiling/scripts/run_case.sh
```

Four GPUs, 1024^3:

```bash
NX=1024 NY=1024 NZ=1024 DIMS="2 2 1" NSTEPS=20 NITER=3 RANKS=4 NODES=1 TASKS_PER_NODE=4 CASE_NAME=manual_1024_4gpu profiling/scripts/run_case.sh
```

## Results

Each run creates:

```text
profiling/results/<job-name>-<job-id>/<case>/
  input.ini
  command.txt
  env.txt
  nvidia-smi.txt
  stdout.txt
  stderr.txt
  time.txt
  summary.csv
```

Collect all summaries under a result root with:

```bash
python3 profiling/scripts/collect_timings.py profiling/results/<job-name>-<job-id> > profiling/results/<job-name>-<job-id>/timings.csv
```

## Remaining Validation Item

Velocity CPU/GPU checks were previously roundoff-level, but pressure output still needs a dedicated comparison. When HDF5/restart is validated on the final cluster stack, compare pressure with and without subtracting the mean and inspect interior versus halo values separately.
