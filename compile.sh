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
    echo "MPI Fortran mod dir   : $(cache_value "$build_dir" MPI_FORTRAN_MODULE_DIR)"
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

include_dirs_from_flags() {
    local next_is_include=0
    local token

    for token in "$@"; do
        if [ "$next_is_include" = "1" ]; then
            printf '%s
' "$token"
            next_is_include=0
            continue
        fi

        case "$token" in
            -I)
                next_is_include=1
                ;;
            -I*)
                printf '%s
' "${token#-I}"
                ;;
        esac
    done
}

find_mpi_fortran_module_dir() {
    local wrapper="$1"
    local wrapper_path
    local compile_flags
    local incdirs
    local libdirs
    local root
    local dir
    local found
    local candidates

    wrapper_path="$(command -v "$wrapper" 2>/dev/null || true)"
    if [ -z "$wrapper_path" ]; then
        return 1
    fi

    compile_flags="$(mpi_wrapper_showme "$wrapper" compile)"
    incdirs="$(mpi_wrapper_showme "$wrapper" incdirs)"
    libdirs="$(mpi_wrapper_showme "$wrapper" libdirs)"
    root="$(cd "$(dirname "$wrapper_path")/.." 2>/dev/null && pwd -P || true)"

    candidates="$(include_dirs_from_flags $compile_flags) $incdirs $libdirs"
    if [ -n "$root" ]; then
        candidates="$candidates $root/include $root/lib $root/lib64"
    fi

    for dir in $candidates; do
        if [ -f "$dir/mpi_f08.mod" ]; then
            printf '%s' "$dir"
            return 0
        fi
    done

    if [ -n "$root" ] && [ -d "$root" ]; then
        found="$(find "$root" -type f -iname 'mpi_f08.mod' -print -quit 2>/dev/null || true)"
        if [ -n "$found" ]; then
            dirname "$found"
            return 0
        fi
    fi

    return 1
}

configure_and_build() {
    local name="$1"
    local build_dir="$2"
    shift 2

    local mpi_args=()
    local use_cmake_find_mpi="${FDM_USE_CMAKE_FIND_MPI:-0}"
    local mpi_fc_wrapper="${MPI_FC_WRAPPER:-mpifort}"
    local mpi_c_wrapper="${MPI_C_WRAPPER:-mpicc}"
    local mpi_fortran_module_dir=""

    # Prefer wrapper flags on HPC systems. They carry the compiler-specific
    # include path for mpi_f08.mod, which CMake FindMPI can miss with NVHPC.
    if [ "${FDM_FAST_MPI:-}" = "1" ]; then
        use_cmake_find_mpi=0
    elif [ "${FDM_FAST_MPI:-}" = "0" ]; then
        use_cmake_find_mpi=1
    fi

    if command -v "$mpi_fc_wrapper" >/dev/null 2>&1; then
        mpi_fortran_module_dir="$(find_mpi_fortran_module_dir "$mpi_fc_wrapper" || true)"
        if [ -n "$mpi_fortran_module_dir" ]; then
            echo "MPI Fortran module dir: $mpi_fortran_module_dir"
            mpi_args+=(-DMPI_FORTRAN_MODULE_DIR="$mpi_fortran_module_dir")
        else
            echo "Warning: could not locate mpi_f08.mod from $mpi_fc_wrapper." >&2
            echo "         If compilation still fails, set MPI_FORTRAN_MODULE_DIR=/path/containing/mpi_f08.mod." >&2
        fi
    fi

    if [ -n "${MPI_FORTRAN_MODULE_DIR:-}" ]; then
        mpi_fortran_module_dir="$MPI_FORTRAN_MODULE_DIR"
        mpi_args+=(-DMPI_FORTRAN_MODULE_DIR="$mpi_fortran_module_dir")
    fi

    if [ "$use_cmake_find_mpi" = "1" ]; then
        mpi_args+=(-DUSE_CMAKE_FIND_MPI=ON)
    else
        local mpi_fc_compile_flags
        local mpi_c_compile_flags
        local mpi_compile_flags
        local mpi_link_flags

        mpi_fc_compile_flags="$(mpi_wrapper_showme "$mpi_fc_wrapper" compile)"
        mpi_c_compile_flags="$(mpi_wrapper_showme "$mpi_c_wrapper" compile)"
        mpi_compile_flags="${mpi_fc_compile_flags} ${mpi_c_compile_flags}"
        if [ -n "$mpi_fortran_module_dir" ]; then
            mpi_compile_flags="$mpi_compile_flags -I$mpi_fortran_module_dir"
        fi
        mpi_link_flags="$(mpi_wrapper_showme "$mpi_fc_wrapper" link)"
        if [ -z "$mpi_link_flags" ]; then
            echo "Could not query MPI link flags from $mpi_fc_wrapper." >&2
            echo "Try setting MPI_FC_WRAPPER=/path/to/mpifort." >&2
            exit 1
        fi

        mpi_args+=(
            -DUSE_CMAKE_FIND_MPI=OFF
            -DMPI_WRAPPER_COMPILE_FLAGS="$mpi_compile_flags"
            -DMPI_WRAPPER_LINK_FLAGS="$mpi_link_flags"
        )
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
        build_cpu_mpi
        ;;
    gpu-hpc)
        build_gpu_mpi
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
        echo "      MPI wrapper flags are used by default; set FDM_USE_CMAKE_FIND_MPI=1 to use CMake FindMPI." >&2
        exit 1
        ;;
esac
