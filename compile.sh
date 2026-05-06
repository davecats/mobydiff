#!/bin/bash
set -euo pipefail

backend="${1:-cpu}"

build_cpu() {
    cmake -S . -B build_noibm -DUSE_IBM=OFF -DUSE_OPENMP_OFFLOAD=OFF -DUSE_CUFFT=OFF
    cmake --build build_noibm

    cmake -S . -B build_ibm -DUSE_IBM=ON -DUSE_OPENMP_OFFLOAD=OFF -DUSE_CUFFT=OFF
    cmake --build build_ibm
}

build_gpu() {
    export FC="${FC:-nvfortran}"

    if ! command -v "$FC" >/dev/null 2>&1; then
        echo "GPU build requires nvfortran. Load the NVIDIA HPC SDK module first." >&2
        echo "Example: module load /opt/nvidia/hpc_sdk/modulefiles/nvhpc-hpcx-cuda13/26.3" >&2
        exit 1
    fi

    cmake -S . -B build_gpu_noibm -DUSE_IBM=OFF -DUSE_OPENMP_OFFLOAD=ON -DUSE_CUFFT=ON
    cmake --build build_gpu_noibm

    cmake -S . -B build_gpu_ibm -DUSE_IBM=ON -DUSE_OPENMP_OFFLOAD=ON -DUSE_CUFFT=ON
    cmake --build build_gpu_ibm
}

case "$backend" in
    cpu)
        build_cpu
        ;;
    gpu)
        build_gpu
        ;;
    all)
        build_cpu
        build_gpu
        ;;
    *)
        echo "Usage: $0 [cpu|gpu|all]" >&2
        exit 1
        ;;
esac
