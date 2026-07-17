#!/usr/bin/env bash
set -euo pipefail

mode="${1:-gpu}"

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
    echo "Executables           : ${build_dir}/moby_solve (+ main symlink), ${build_dir}/moby_prepare"
    echo "Fortran compiler      : $(cache_value "$build_dir" CMAKE_Fortran_COMPILER)"
    echo "C compiler            : $(cache_value "$build_dir" CMAKE_C_COMPILER)"
    echo "Pressure solver       : red-black SOR"
    echo "MPI                   : ON"
    echo "MPI wrapper compile   : $(cache_value "$build_dir" MPI_WRAPPER_COMPILE_FLAGS)"
    echo "MPI wrapper link      : $(cache_value "$build_dir" MPI_WRAPPER_LINK_FLAGS)"
    echo "MPI Fortran mod dir   : $(cache_value "$build_dir" MPI_FORTRAN_MODULE_DIR)"
    echo "USE_OPENMP_OFFLOAD    : $(cache_value "$build_dir" USE_OPENMP_OFFLOAD)"
    echo "USE_IBM_SECONDORDER   : $(cache_value "$build_dir" USE_IBM_SECONDORDER)"
    echo "OPENMP_OFFLOAD_FLAGS  : $(cache_value "$build_dir" OPENMP_OFFLOAD_FLAGS)"
    echo "HDF5_ROOT             : $(cache_value "$build_dir" HDF5_ROOT)"
    echo "HDF5 include dir      : $(cache_value "$build_dir" HDF5_C_INCLUDE_DIR)"
    echo "HDF5 C library        : $(cache_value "$build_dir" HDF5_C_LIB)"
    echo "========================================"
    echo
}

require_command() {
    local exe="$1"
    if ! command -v "$exe" >/dev/null 2>&1; then
        echo "Required command not found: $exe" >&2
        exit 1
    fi
}

mpi_wrapper_showme() {
    local wrapper="$1"
    local kind="$2"

    require_command "$wrapper"
    "$wrapper" --showme:"$kind" 2>/dev/null || "$wrapper" -showme:"$kind" 2>/dev/null || true
}

include_dirs_from_flags() {
    local next_is_include=0
    local token

    for token in "$@"; do
        if [ "$next_is_include" = "1" ]; then
            printf '%s\n' "$token"
            next_is_include=0
            continue
        fi

        case "$token" in
            -I)
                next_is_include=1
                ;;
            -I*)
                printf '%s\n' "${token#-I}"
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
    local root_dir
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
    root_dir="$(cd "$(dirname "$wrapper_path")/.." 2>/dev/null && pwd -P || true)"

    candidates="$(include_dirs_from_flags $compile_flags) $incdirs $libdirs"
    if [ -n "$root_dir" ]; then
        candidates="$candidates $root_dir/include $root_dir/lib $root_dir/lib64"
    fi

    for dir in $candidates; do
        if [ -f "$dir/mpi_f08.mod" ]; then
            printf '%s' "$dir"
            return 0
        fi
    done

    if [ -n "$root_dir" ] && [ -d "$root_dir" ]; then
        found="$(find "$root_dir" -type f -iname 'mpi_f08.mod' -print -quit 2>/dev/null || true)"
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
    local offload="$3"

    export FC="${FC:-${MPI_FC_WRAPPER:-mpifort}}"
    export CC="${CC:-${MPI_C_WRAPPER:-mpicc}}"

    require_command "$FC"
    require_command "$CC"
    if [ "$offload" = "ON" ]; then
        require_command nvfortran
    fi

    local mpi_fc_wrapper="${MPI_FC_WRAPPER:-mpifort}"
    local mpi_c_wrapper="${MPI_C_WRAPPER:-mpicc}"
    local mpi_fortran_module_dir="${MPI_FORTRAN_MODULE_DIR:-}"
    local mpi_fc_compile_flags
    local mpi_c_compile_flags
    local mpi_compile_flags
    local mpi_link_flags
    local hdf5_root
    local -a cmake_args

    if [ -z "$mpi_fortran_module_dir" ]; then
        mpi_fortran_module_dir="$(find_mpi_fortran_module_dir "$mpi_fc_wrapper" || true)"
    fi
    if [ -n "$mpi_fortran_module_dir" ]; then
        echo "MPI Fortran module dir: $mpi_fortran_module_dir"
    fi

    mpi_fc_compile_flags="$(mpi_wrapper_showme "$mpi_fc_wrapper" compile)"
    mpi_c_compile_flags="$(mpi_wrapper_showme "$mpi_c_wrapper" compile)"
    mpi_compile_flags="${mpi_fc_compile_flags} ${mpi_c_compile_flags}"
    if [ -n "$mpi_fortran_module_dir" ]; then
        mpi_compile_flags="$mpi_compile_flags -I$mpi_fortran_module_dir"
    fi

    mpi_link_flags="$(mpi_wrapper_showme "$mpi_fc_wrapper" link)"
    if [ -z "$mpi_link_flags" ]; then
        echo "Could not query MPI link flags from $mpi_fc_wrapper." >&2
        exit 1
    fi

    hdf5_root="${HDF5_ROOT:-}"
    if [ -n "$hdf5_root" ]; then
        echo "HDF5_ROOT: $hdf5_root"
    else
        echo "HDF5_ROOT: not set; using CMake/system search paths"
    fi

    cmake_args=(
        -S . -B "$build_dir"
        -UHDF5_C_INCLUDE_DIR -UHDF5_C_LIB
        -DHDF5_ROOT="$hdf5_root"
        -DMPI_WRAPPER_COMPILE_FLAGS="$mpi_compile_flags"
        -DMPI_WRAPPER_LINK_FLAGS="$mpi_link_flags"
        -DMPI_FORTRAN_MODULE_DIR="$mpi_fortran_module_dir"
        -DUSE_OPENMP_OFFLOAD="$offload"
    )

    cmake "${cmake_args[@]}"

    cmake --build "$build_dir" -j
    print_build_summary "$name" "$build_dir"
}

case "$mode" in
    cpu)
        configure_and_build "cpu" "build_cpu" OFF
        ;;
    gpu)
        configure_and_build "gpu" "build_gpu" ON
        ;;
    *)
        echo "Usage: $0 [cpu|gpu]" >&2
        echo "MPI is always enabled. GPU builds require a GPU-aware MPI stack." >&2
        echo "Set HDF5_ROOT=/path/to/parallel-hdf5 if CMake cannot find HDF5." >&2
        exit 1
        ;;
esac
