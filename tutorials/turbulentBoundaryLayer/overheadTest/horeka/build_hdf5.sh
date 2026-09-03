#!/bin/bash
# Download and build PARALLEL HDF5 (C library) with the NVIDIA HPC compilers,
# using HDF5's CMake build (HDF5 2.x is CMake-only).
#
#     build_hdf5.sh [install_dir]        # default install_dir = $HOME/hdf5
#
# HoreKa only ships a SERIAL prebuilt HDF5 module, but mobydiff's field I/O uses
# MPI-IO, so we need a parallel build. This compiles HDF5 with CC=mpicc (the MPI C
# wrapper from toolkit/nvidia-hpc-sdk/25.3, i.e. nvc under the bundled OpenMPI), so
# it is ABI-compatible with the solver's field_hdf5.c built by the same toolchain.
# mobydiff uses only the core C API (H5Dopen2/H5Pset_dxpl_mpio/... -- unchanged in
# HDF5 2.x), so the Fortran/C++/HL/tools components are left off.
#
# PREREQUISITES (loaded by the caller): `module load toolkit/nvidia-hpc-sdk/25.3`
# for mpicc, plus a `cmake` (the mobydiff solver build needs cmake too). Idempotent:
# skips the build if the library already exists.

set -euo pipefail

HDF5_DIR="${1:-${HDF5_DIR:-$HOME/hdf5}}"
# HDF5_VERSION is HARDCODED (not read from the environment): a loaded serial-HDF5
# module may export a stray HDF5_VERSION that would otherwise poison the paths.
# To use a different release, override HDF5_URL directly.
HDF5_VERSION=2.2.0
# GitHub tag archive (contains CMakeLists.txt; CMake needs no pre-generated
# configure). See https://github.com/HDFGroup/hdf5/releases for tag names.
HDF5_URL="${HDF5_URL:-https://github.com/HDFGroup/hdf5/archive/refs/tags/${HDF5_VERSION}.tar.gz}"

if [ -e "$HDF5_DIR/lib/libhdf5.so" ] || [ -e "$HDF5_DIR/lib/libhdf5.a" ]; then
    echo "HDF5 already installed at $HDF5_DIR -- skipping build."
    exit 0
fi

command -v mpicc >/dev/null || { echo "ERROR: mpicc not on PATH. module load toolkit/nvidia-hpc-sdk/25.3 first." >&2; exit 1; }
command -v cmake >/dev/null || { echo "ERROR: cmake not on PATH. Load a cmake module (also needed by the solver build)." >&2; exit 1; }
echo "=== building parallel HDF5 (CMake) -> $HDF5_DIR ==="
echo "    source : $HDF5_URL"
echo "    CC = mpicc -> $(mpicc -show 2>/dev/null | head -1 || echo nvc)"

WORK="$HDF5_DIR/_build"
mkdir -p "$WORK"
cd "$WORK"

# Version-independent tarball name so a mismatched env can't desync the paths.
TARBALL="hdf5-download.tar.gz"
if [ ! -f "$TARBALL" ]; then
    echo "=== downloading $HDF5_URL ==="
    if command -v wget >/dev/null; then wget -O "$TARBALL" "$HDF5_URL";
    else curl -fL -o "$TARBALL" "$HDF5_URL"; fi
fi
# Detect the archive's actual top-level directory instead of assuming its name
# (the GitHub tag archive extracts to hdf5-<tag>/, whatever the tag is).
SRCDIR="$(tar tzf "$TARBALL" | sed -n '1p' | cut -d/ -f1)"
[ -n "$SRCDIR" ] || { echo "ERROR: could not read the HDF5 source dir from $TARBALL" >&2; exit 1; }
rm -rf "$SRCDIR"
tar xzf "$TARBALL"
[ -d "$SRCDIR" ] || { echo "ERROR: extracted HDF5 source dir $SRCDIR not found" >&2; exit 1; }
echo "=== HDF5 source: $SRCDIR ==="

# Parallel, C-only, shared. HDF5_ENABLE_PARALLEL needs an MPI C compiler; passing
# CC=mpicc directly is the simplest route. Fortran/C++/HL/tools/tests/examples off.
cmake -S "$SRCDIR" -B cmbuild \
    -DCMAKE_C_COMPILER=mpicc \
    -DCMAKE_INSTALL_PREFIX="$HDF5_DIR" \
    -DCMAKE_BUILD_TYPE=Release \
    -DHDF5_ENABLE_PARALLEL=ON \
    -DBUILD_SHARED_LIBS=ON \
    -DONLY_SHARED_LIBS=ON \
    -DBUILD_TESTING=OFF \
    -DHDF5_BUILD_TOOLS=OFF \
    -DHDF5_BUILD_EXAMPLES=OFF \
    -DHDF5_BUILD_UTILS=OFF \
    -DHDF5_BUILD_FORTRAN=OFF \
    -DHDF5_BUILD_CPP_LIB=OFF \
    -DHDF5_BUILD_HL_LIB=OFF

cmake --build cmbuild -j"$(nproc)"
cmake --install cmbuild

echo "=== HDF5 installed under $HDF5_DIR ==="
ls -la "$HDF5_DIR/lib/"libhdf5.* 2>/dev/null | head
# HDF5 CMake may install headers under include/ -- confirm hdf5.h is where the
# solver's CMake search (HDF5_ROOT + include/) expects it.
[ -e "$HDF5_DIR/include/hdf5.h" ] || echo "NOTE: hdf5.h not at $HDF5_DIR/include -- check the install layout."
