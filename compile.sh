#!/bin/bash
set -euo pipefail

backend="${1:-cpu}"

cache_value() {
    local build_dir="$1"
    local key="$2"
    local value

    value="$(grep -E "^${key}(:[^=]*)?=" "${build_dir}/CMakeCache.txt" 2>/dev/null | head -n 1 | cut -d= -f2-)"
    printf '%s' "${value:-<unset>}"
}

print_build_summary() {
    local name="$1"
    local build_dir="$2"
    local pressure_backend

    if [[ "$(cache_value "$build_dir" USE_REDBLACK)" == "ON" ]]; then
        pressure_backend="red-black SOR"
    elif [[ "$(cache_value "$build_dir" USE_CUFFT)" == "ON" ]]; then
        pressure_backend="FFT/cuFFT"
    else
        pressure_backend="FFT/FFTW"
    fi

    echo
    echo "========================================"
    echo "Build summary: ${name}"
    echo "========================================"
    echo "Build directory       : ${build_dir}"
    echo "Executable            : ${build_dir}/main"
    echo "Fortran compiler      : $(cache_value "$build_dir" CMAKE_Fortran_COMPILER)"
    echo "C compiler            : $(cache_value "$build_dir" CMAKE_C_COMPILER)"
    echo "Pressure backend      : ${pressure_backend}"
    echo "USE_OPENMP_OFFLOAD    : $(cache_value "$build_dir" USE_OPENMP_OFFLOAD)"
    echo "USE_CUFFT             : $(cache_value "$build_dir" USE_CUFFT)"
    echo "USE_REDBLACK          : $(cache_value "$build_dir" USE_REDBLACK)"
    echo "USE_IBM_SECONDORDER   : $(cache_value "$build_dir" USE_IBM_SECONDORDER)"
    echo "OPENMP_OFFLOAD_FLAGS  : $(cache_value "$build_dir" OPENMP_OFFLOAD_FLAGS)"
    echo "========================================"
    echo
}

build_cpu() {
    cmake -S . -B build_ibm -DUSE_OPENMP_OFFLOAD=OFF -DUSE_CUFFT=OFF -DUSE_REDBLACK=OFF
    cmake --build build_ibm
    print_build_summary "cpu" "build_ibm"
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
    print_build_summary "gpu" "build_gpu_ibm"
}

build_redblack() {
    cmake -S . -B build_cpu_redblack -DUSE_REDBLACK=ON -DUSE_OPENMP_OFFLOAD=OFF -DUSE_CUFFT=OFF
    cmake --build build_cpu_redblack
    print_build_summary "redblack" "build_cpu_redblack"
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
    print_build_summary "redblack-gpu" "build_gpu_redblack"
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
