#!/bin/bash
set -euo pipefail

mode="${1:-cpu}"

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

    echo
    echo "========================================"
    echo "Build summary: ${name}"
    echo "========================================"
    echo "Build directory       : ${build_dir}"
    echo "Executable            : ${build_dir}/main"
    echo "Fortran compiler      : $(cache_value "$build_dir" CMAKE_Fortran_COMPILER)"
    echo "C compiler            : $(cache_value "$build_dir" CMAKE_C_COMPILER)"
    echo "Pressure solver       : red-black SOR"
    echo "MPI                   : ON"
    echo "USE_CMAKE_FIND_MPI    : $(cache_value "$build_dir" USE_CMAKE_FIND_MPI)"
    echo "MPI wrapper compile   : $(cache_value "$build_dir" MPI_WRAPPER_COMPILE_FLAGS)"
    echo "MPI wrapper link      : $(cache_value "$build_dir" MPI_WRAPPER_LINK_FLAGS)"
    echo "USE_OPENMP_OFFLOAD    : $(cache_value "$build_dir" USE_OPENMP_OFFLOAD)"
    echo "USE_IBM_SECONDORDER   : $(cache_value "$build_dir" USE_IBM_SECONDORDER)"
    echo "OPENMP_OFFLOAD_FLAGS  : $(cache_value "$build_dir" OPENMP_OFFLOAD_FLAGS)"
    echo "HDF5 include dir      : $(cache_value "$build_dir" HDF5_C_INCLUDE_DIR)"
    echo "HDF5 C library        : $(cache_value "$build_dir" HDF5_C_LIB)"
    echo "========================================"
    echo
}

require_nvfortran() {
    export FC="${FC:-nvfortran}"

    if ! command -v "$FC" >/dev/null 2>&1; then
        echo "GPU build requires nvfortran. Load the NVIDIA HPC SDK module first." >&2
        echo "Example: module load /opt/nvidia/hpc_sdk/modulefiles/nvhpc-hpcx-cuda13/26.3" >&2
        exit 1
    fi
}

mpi_wrapper_showme() {
    local wrapper="$1"
    local kind="$2"

    if ! command -v "$wrapper" >/dev/null 2>&1; then
        echo "MPI wrapper not found: $wrapper" >&2
        exit 1
    fi

    "$wrapper" --showme:"$kind" 2>/dev/null || "$wrapper" -showme:"$kind" 2>/dev/null || true
}

configure_and_build() {
    local name="$1"
    local build_dir="$2"
    shift 2

    local mpi_args=()
    if [ "${FDM_FAST_MPI:-0}" = "1" ]; then
        local mpi_fc_wrapper="${MPI_FC_WRAPPER:-mpifort}"
        local mpi_compile_flags
        local mpi_link_flags

        mpi_compile_flags="$(mpi_wrapper_showme "$mpi_fc_wrapper" compile)"
        mpi_link_flags="$(mpi_wrapper_showme "$mpi_fc_wrapper" link)"
        if [ -z "$mpi_link_flags" ]; then
            echo "Could not query MPI link flags from $mpi_fc_wrapper." >&2
            echo "Try setting MPI_FC_WRAPPER=/path/to/mpifort, or unset FDM_FAST_MPI." >&2
            exit 1
        fi

        mpi_args+=(
            -DUSE_CMAKE_FIND_MPI=OFF
            -DMPI_WRAPPER_COMPILE_FLAGS="$mpi_compile_flags"
            -DMPI_WRAPPER_LINK_FLAGS="$mpi_link_flags"
        )
    else
        mpi_args+=(-DUSE_CMAKE_FIND_MPI=ON)
    fi

    cmake -S . -B "$build_dir" \
        -UHDF5_C_INCLUDE_DIR -UHDF5_C_LIB \
        -DHDF5_ROOT="${HDF5_ROOT:-/usr/lib/x86_64-linux-gnu/hdf5/openmpi}" \
        "${mpi_args[@]}" \
        "$@"
    cmake --build "$build_dir" -j
    print_build_summary "$name" "$build_dir"
}
build_cpu() {
    build_cpu_mpi
}

build_gpu() {
    build_gpu_mpi
}

build_cpu_mpi() {
    configure_and_build "cpu-mpi" "build_cpu_mpi" \
        -DUSE_OPENMP_OFFLOAD=OFF
}

build_gpu_mpi() {
    require_nvfortran
    configure_and_build "gpu-mpi" "build_gpu_mpi" \
        -DUSE_OPENMP_OFFLOAD=ON
}

case "$mode" in
    cpu)
        build_cpu
        ;;
    gpu)
        build_gpu
        ;;
    cpu-mpi)
        build_cpu_mpi
        ;;
    gpu-mpi)
        build_gpu_mpi
        ;;
    cpu-hpc)
        FDM_FAST_MPI=1 build_cpu_mpi
        ;;
    gpu-hpc)
        FDM_FAST_MPI=1 build_gpu_mpi
        ;;
    all)
        build_cpu_mpi
        build_gpu_mpi
        ;;
    all-hpc)
        FDM_FAST_MPI=1 build_cpu_mpi
        FDM_FAST_MPI=1 build_gpu_mpi
        ;;
    *)
        echo "Usage: $0 [cpu|gpu|cpu-mpi|gpu-mpi|cpu-hpc|gpu-hpc|all|all-hpc]" >&2
        echo "Note: MPI is always enabled; cpu/gpu are aliases for cpu-mpi/gpu-mpi." >&2
        echo "      *-hpc modes skip CMake FindMPI and use mpifort --showme flags." >&2
        exit 1
        ;;
esac
