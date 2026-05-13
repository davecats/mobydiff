#!/bin/bash
set -euo pipefail

backend="${1:-cpu}"

build_cpu() {
    cmake -S . -B build_ibm -DUSE_OPENMP_OFFLOAD=OFF -DUSE_CUFFT=OFF -DUSE_REDBLACK=OFF
    cmake --build build_ibm
}

build_gpu() {
    export FC="${FC:-nvfortran}"

    if ! command -v "$FC" >/dev/null 2>&1; then
        echo "GPU build requires nvfortran. Load the NVIDIA HPC SDK module first." >&2
        echo "Example: module load /opt/nvidia/hpc_sdk/modulefiles/nvhpc-hpcx-cuda13/26.3" >&2
        exit 1
    fi

    cmake -S . -B build_gpu_ibm -DUSE_OPENMP_OFFLOAD=ON -DUSE_CUFFT=ON -DUSE_REDBLACK=OFF
    cmake --build build_gpu_ibm
}

build_redblack() {
    cmake -S . -B build_cpu_redblack -DUSE_REDBLACK=ON -DUSE_OPENMP_OFFLOAD=OFF -DUSE_CUFFT=OFF
    cmake --build build_cpu_redblack
}

build_redblack_gpu() {
    export FC="${FC:-nvfortran}"

    if ! command -v "$FC" >/dev/null 2>&1; then
        echo "GPU red-black build requires nvfortran. Load the NVIDIA HPC SDK module first." >&2
        echo "Example: module load /opt/nvidia/hpc_sdk/modulefiles/nvhpc-hpcx-cuda13/26.3" >&2
        exit 1
    fi

    cmake -S . -B build_gpu_redblack -DUSE_REDBLACK=ON -DUSE_OPENMP_OFFLOAD=ON -DUSE_CUFFT=OFF
    cmake --build build_gpu_redblack
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
	build_redblack
	build_redblack_gpu
        ;;
    redblack)
        build_redblack
        ;;
    redblack-gpu)
        build_redblack_gpu
        ;;
    *)
        echo "Usage: $0 [cpu|gpu|redblack|redblack-gpu|all]" >&2
        exit 1
        ;;
esac
