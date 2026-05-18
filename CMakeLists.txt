cmake_minimum_required(VERSION 3.20)

project(3D_NSE_SOLVER LANGUAGES Fortran C)

# =========================
# Options
# =========================
option(USE_OPENMP_OFFLOAD "Enable OpenMP target offload for accelerator execution" OFF)
option(USE_CUFFT "Use cuFFT instead of FFTW for the pressure FFT backend" ${USE_OPENMP_OFFLOAD})
option(USE_REDBLACK "Use the red-black SOR pressure solver backend" OFF)
option(USE_IBM_SECONDORDER "Use second-order IBM interpolation/correction" ON)
set(OPENMP_OFFLOAD_FLAGS "" CACHE STRING "Additional Fortran compiler/linker flags for OpenMP target offload")

if(USE_REDBLACK)
    set(USE_CUFFT OFF CACHE BOOL "Use cuFFT instead of FFTW for the pressure FFT backend" FORCE)
endif()

if(USE_OPENMP_OFFLOAD AND NOT USE_REDBLACK AND NOT USE_CUFFT)
    message(FATAL_ERROR "USE_OPENMP_OFFLOAD requires USE_CUFFT=ON for the GPU pressure FFT backend.")
endif()

if(USE_CUFFT AND NOT USE_OPENMP_OFFLOAD)
    message(FATAL_ERROR "USE_CUFFT currently requires USE_OPENMP_OFFLOAD=ON so cuFFT receives OpenMP device pointers.")
endif()

# =========================
# Sources
# =========================
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(MODULE_DIR ${SRC_DIR}/modules)
set(PRESSURE_DIR ${MODULE_DIR}/pressure)

set(SOURCES
    ${SRC_DIR}/main.f90
    ${MODULE_DIR}/init.f90
    ${MODULE_DIR}/config.f90
    ${MODULE_DIR}/gpu_runtime.f90
    ${MODULE_DIR}/boundary.f90
    ${MODULE_DIR}/ibm.f90
    ${MODULE_DIR}/field_hdf5.c
    ${MODULE_DIR}/io.f90
    ${MODULE_DIR}/step.f90
    ${MODULE_DIR}/pressure_workspace.f90
    ${MODULE_DIR}/pressure_solver.f90
)

if(USE_REDBLACK)
    list(APPEND SOURCES
        ${PRESSURE_DIR}/pressure_redblack.f90
    )
elseif(USE_CUFFT)
    list(APPEND SOURCES
        ${PRESSURE_DIR}/cufft_bindings.f90
        ${PRESSURE_DIR}/pressure_fft.f90
        ${PRESSURE_DIR}/pressure_backend_cufft.f90
    )
else()
    list(APPEND SOURCES
        ${PRESSURE_DIR}/pressure_fft.f90
        ${PRESSURE_DIR}/pressure_backend_fftw.f90
    )
endif()

add_executable(main ${SOURCES})

# =========================
# Module directory
# =========================
set_target_properties(main PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod
)

# =========================
# Compiler flags
# =========================
target_compile_options(main PRIVATE
    $<$<COMPILE_LANGUAGE:Fortran>:-O2>
    $<$<COMPILE_LANGUAGE:C>:-O2>
)

# =========================
# HDF5 field output
# =========================
set(HDF5_ROOT "" CACHE PATH "Root directory of HDF5 installation")

find_path(HDF5_C_INCLUDE_DIR
    NAMES hdf5.h
    HINTS
        ${HDF5_ROOT}/include
        ${HDF5_ROOT}
        /usr/include/hdf5/serial
    PATH_SUFFIXES include
)

find_library(HDF5_C_LIB
    NAMES hdf5
    HINTS
        ${HDF5_ROOT}/lib
        ${HDF5_ROOT}/lib64
        ${HDF5_ROOT}
        /usr/lib/x86_64-linux-gnu/hdf5/serial
    PATH_SUFFIXES lib lib64
)

if(NOT HDF5_C_INCLUDE_DIR)
    message(FATAL_ERROR "Could not find hdf5.h. Set -DHDF5_ROOT=/path/to/hdf5 or -DHDF5_C_INCLUDE_DIR=/path/to/include")
endif()

if(NOT HDF5_C_LIB)
    message(FATAL_ERROR "Could not find HDF5 C library. Set -DHDF5_ROOT=/path/to/hdf5 or -DHDF5_C_LIB=/path/to/libhdf5.so")
endif()

message(STATUS "HDF5 include dir: ${HDF5_C_INCLUDE_DIR}")
message(STATUS "HDF5 C library: ${HDF5_C_LIB}")

target_include_directories(main PRIVATE ${HDF5_C_INCLUDE_DIR})
target_link_libraries(main PRIVATE ${HDF5_C_LIB})

# Enable preprocessing depending on compiler
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(main PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-cpp>)
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    target_compile_options(main PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-fpp>)
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
    target_compile_options(main PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-Mpreprocess>)
endif()

# =========================
# Define macros
# =========================
if(USE_OPENMP_OFFLOAD)
    target_compile_definitions(main PRIVATE USE_OPENMP_OFFLOAD)

    find_package(OpenMP REQUIRED COMPONENTS Fortran)
    target_link_libraries(main PRIVATE OpenMP::OpenMP_Fortran)

    if(NOT OPENMP_OFFLOAD_FLAGS)
        if(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
            set(OPENMP_OFFLOAD_FLAGS "-mp=gpu")
        elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
            set(OPENMP_OFFLOAD_FLAGS "-fopenmp -foffload=nvptx-none")
        elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
            set(OPENMP_OFFLOAD_FLAGS "-fiopenmp -fopenmp-targets=spir64")
        endif()
    endif()

    if(OPENMP_OFFLOAD_FLAGS)
        separate_arguments(OPENMP_OFFLOAD_FLAGS_LIST NATIVE_COMMAND "${OPENMP_OFFLOAD_FLAGS}")
        foreach(flag IN LISTS OPENMP_OFFLOAD_FLAGS_LIST)
            target_compile_options(main PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${flag}>)
        endforeach()
        target_link_options(main PRIVATE ${OPENMP_OFFLOAD_FLAGS_LIST})
    endif()
endif()

if(USE_REDBLACK)
    target_compile_definitions(main PRIVATE USE_REDBLACK)
endif()

if(USE_IBM_SECONDORDER)
    target_compile_definitions(main PRIVATE USE_IBM_SECONDORDER)
endif()


# =========================
# FFT backend
# =========================

if(USE_REDBLACK)
    if(USE_OPENMP_OFFLOAD)
        message(STATUS "Pressure backend: red-black SOR OpenMP offload")
    else()
        message(STATUS "Pressure backend: red-black SOR CPU")
    endif()
elseif(USE_CUFFT)
    target_compile_definitions(main PRIVATE USE_CUFFT)

    find_package(CUDAToolkit REQUIRED)

    message(STATUS "cuFFT library: CUDA::cufft")
    target_include_directories(main PRIVATE ${CUDAToolkit_INCLUDE_DIRS})
    target_link_libraries(main PRIVATE CUDA::cufft CUDA::cudart)
else()
    # Optional user-provided FFTW root directory, e.g.
    # cmake -B build -DFFTW3_ROOT=/path/to/fftw
    set(FFTW3_ROOT "" CACHE PATH "Root directory of FFTW3 installation")

    find_path(FFTW3_INCLUDE_DIR
        NAMES fftw3.f03
        HINTS
            ${FFTW3_ROOT}/include
            ${FFTW3_ROOT}
        PATH_SUFFIXES include
    )

    find_library(FFTW3_LIB
        NAMES fftw3
        HINTS
            ${FFTW3_ROOT}/lib
            ${FFTW3_ROOT}/lib64
            ${FFTW3_ROOT}
        PATH_SUFFIXES lib lib64
    )

    if(NOT FFTW3_INCLUDE_DIR)
        message(FATAL_ERROR "Could not find fftw3.f03. Set -DFFTW3_ROOT=/path/to/fftw or -DFFTW3_INCLUDE_DIR=/path/to/include")
    endif()

    if(NOT FFTW3_LIB)
        message(FATAL_ERROR "Could not find FFTW3 library. Set -DFFTW3_ROOT=/path/to/fftw or -DFFTW3_LIB=/path/to/libfftw3.so")
    endif()

    message(STATUS "FFTW3 include dir: ${FFTW3_INCLUDE_DIR}")
    message(STATUS "FFTW3 library: ${FFTW3_LIB}")

    target_include_directories(main PRIVATE ${FFTW3_INCLUDE_DIR})
    target_link_libraries(main PRIVATE ${FFTW3_LIB})
endif()


# =========================
# Optional: include dirs (only if really needed)
# =========================
# target_include_directories(main PRIVATE /usr/include)
