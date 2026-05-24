# Scaling and GPU Profiling Harness

This directory contains reproducible scripts for strong scaling, weak scaling, and first-pass GPU profiling on the accelerated SLURM partition.

The scripts assume one MPI rank per GPU. The target node described for these templates has:

- 4x NVIDIA A100 40 GB
- 2x Intel Xeon Platinum 8368, 38 cores each
- 512 GiB host memory
- SLURM partition: `accelerated`
- SLURM account: `hk-project-exasim`

## Module Setup

The batch scripts load:

```bash
module load toolkit/nvidia-hpc-sdk/25.3-nompi
module load mpi/openmpi/5.0
module load lib/hdf5/1.14
```

If the site module exposes HDF5 under a different variable, set `HDF5_ROOT` before calling `./compile.sh gpu`.

## Quick Start

From the repository root on the HPC system:

```bash
sbatch profiling/slurm/strong_scaling_gpu.slurm
sbatch profiling/slurm/strong_scaling_1024_gpu.slurm
sbatch profiling/slurm/weak_scaling_gpu.slurm
sbatch profiling/slurm/nsys_gpu.slurm
```

Collect timing summaries after jobs complete:

```bash
python3 profiling/scripts/collect_timings.py profiling/results > profiling/results/timings.csv
```

## Useful Overrides

All SLURM scripts can be configured with exported variables:

```bash
sbatch --export=ALL,SIZE=1024,RANKS_LIST="4 8 16",NSTEPS=20,NITER=3 profiling/slurm/strong_scaling_gpu.slurm
sbatch --export=ALL,LOCAL_N=256,RANKS_LIST="1 2 4 8 16",NSTEPS=20 profiling/slurm/weak_scaling_gpu.slurm
sbatch --export=ALL,SIZE=512,RANKS=4,NSTEPS=5 profiling/slurm/nsys_gpu.slurm
```

Common variables:

- `SIZE`: global cubic size for strong scaling, e.g. `512` or `1024`.
- `RANKS_LIST`: list of MPI ranks/GPU counts to test.
- `LOCAL_N`: local cells per direction and rank for weak scaling.
- `NSTEPS`: number of timed steps.
- `NITER`: red-black SOR iterations.
- `IBM_ENABLED`: `true` or `false`.
- `DO_BUILD`: `1` builds before running, `0` reuses the existing executable.
- `RESULTS_ROOT`: output directory, default `profiling/results/<job-name>-<job-id>`.

## 1024^3 Case

A static input file is provided at `profiling/cases/strong_1024.ini`. For scaling runs, the scripts generate rank-specific inputs from the same defaults and set `[mpi] dims` explicitly.

The current large arrays are roughly `q(4) + qs(3) + oldrhs(3) + ibm%coef(3) = 13` double fields, about 104 bytes per owned cell before halos and small metric arrays. A `1024^3` case is therefore about 112 GiB total for these major arrays, or about 28 GiB per GPU on 4 ranks. It should not be expected to fit on one A100-40 with the current data layout.

## Recommended Scaling Runs

For one A100 node:

```bash
sbatch --nodes=1 --export=ALL,SIZE=512,RANKS_LIST="1 2 4" profiling/slurm/strong_scaling_gpu.slurm
sbatch --nodes=1 --export=ALL,SIZE=1024,RANKS_LIST="4" profiling/slurm/strong_scaling_gpu.slurm
```

For up to four A100 nodes:

```bash
sbatch --nodes=4 --export=ALL,SIZE=1024,RANKS_LIST="4 8 16" profiling/slurm/strong_scaling_gpu.slurm
sbatch --nodes=4 --export=ALL,LOCAL_N=256,RANKS_LIST="1 2 4 8 16" profiling/slurm/weak_scaling_gpu.slurm
```

## Results Layout

Each run creates one directory per case:

```text
profiling/results/<job-name>-<job-id>/
  strong_1024_4r/
    input.ini
    stdout.txt
    stderr.txt
    time.txt
    env.txt
    nvidia-smi.txt
    summary.csv
```

The solver's own timing line is parsed from `stdout.txt` into `summary.csv`.

## Nsight Systems

`profiling/slurm/nsys_gpu.slurm` runs Nsight Systems once per MPI rank and writes one report per rank. Use this after timing has identified the important scaling point. Profiling every scaling point is usually expensive and noisy.
