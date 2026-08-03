#!/usr/bin/env bash
# Bit-exactness builds: -Mnofma (CPU) / -Mnofma -gpu=nofma (GPU).
# Default FMA contraction makes the compiler introduce 1-2 ulp differences
# for arithmetically identical source, so every "max_abs 0" gate is run on
# nofma binaries (CLAUDE.md "Verification").
#
#   ./compile_nofma.sh cpu   -> build_cpu_nofma
#   ./compile_nofma.sh gpu   -> build_gpu_nofma
#
# Reuses the MPI/HDF5 discovery already done by ./compile.sh: the matching
# build_cpu / build_gpu cache must exist.
set -euo pipefail
cd "$(dirname "$0")/../.."

mode="${1:-cpu}"
case "$mode" in
    cpu) src=build_cpu; dst=build_cpu_nofma; offload=OFF; ompflags="" ;;
    gpu) src=build_gpu; dst=build_gpu_nofma; offload=ON;  ompflags="-mp=gpu -gpu=nofma" ;;
    *) echo "usage: $0 [cpu|gpu]" >&2; exit 1 ;;
esac

cache() { grep -E "^$2(:[^=]*)?=" "$1/CMakeCache.txt" | head -n1 | cut -d= -f2-; }
[ -f "$src/CMakeCache.txt" ] || { echo "run ./compile.sh $mode first" >&2; exit 1; }

cmake -S . -B "$dst" \
    -DCMAKE_Fortran_COMPILER="$(cache $src CMAKE_Fortran_COMPILER)" \
    -DCMAKE_C_COMPILER="$(cache $src CMAKE_C_COMPILER)" \
    -DCMAKE_Fortran_FLAGS="-Mnofma" \
    -DHDF5_ROOT="$(cache $src HDF5_ROOT)" \
    -DMPI_WRAPPER_COMPILE_FLAGS="$(cache $src MPI_WRAPPER_COMPILE_FLAGS)" \
    -DMPI_WRAPPER_LINK_FLAGS="$(cache $src MPI_WRAPPER_LINK_FLAGS)" \
    -DMPI_FORTRAN_MODULE_DIR="$(cache $src MPI_FORTRAN_MODULE_DIR)" \
    -DUSE_OPENMP_OFFLOAD="$offload" \
    -DOPENMP_OFFLOAD_FLAGS="$ompflags"
cmake --build "$dst" -j
