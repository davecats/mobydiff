cmake_minimum_required(VERSION 3.20)

project(3D_NSE_SOLVER LANGUAGES Fortran)

# =========================
# Options
# =========================
option(USE_IBM "Enable immersed boundary method (IBM)" OFF)
option(USE_OPENMP_OFFLOAD "Enable OpenMP target offload for accelerator execution" OFF)
option(USE_CUFFT "Use cuFFT instead of FFTW for the Poisson FFT backend" ${USE_OPENMP_OFFLOAD})
set(OPENMP_OFFLOAD_FLAGS "" CACHE STRING "Additional Fortran compiler/linker flags for OpenMP target offload")

if(USE_OPENMP_OFFLOAD AND NOT USE_CUFFT)
    message(FATAL_ERROR "USE_OPENMP_OFFLOAD requires USE_CUFFT=ON for the GPU Poisson FFT backend.")
endif()

if(USE_CUFFT AND NOT USE_OPENMP_OFFLOAD)
    message(FATAL_ERROR "USE_CUFFT currently requires USE_OPENMP_OFFLOAD=ON so cuFFT receives OpenMP device pointers.")
endif()

# =========================
# Sources
# =========================
set(SOURCES
    main.f90
    init.f90
    config.f90
    gpu_runtime.f90
    boundary.f90
    io.f90
    step.f90
    poisson_workspace.f90
    poisson_common.f90
    poisson_solver.f90
)

if(USE_CUFFT)
    list(APPEND SOURCES
        cufft_bindings.f90
        poisson_cufft.f90
    )
else()
    list(APPEND SOURCES
        poisson_fftw.f90
    )
endif()

# Add IBM sources only if enabled
if(USE_IBM)
    list(APPEND SOURCES
        ibm.f90
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
target_compile_options(main PRIVATE -O2)

# Enable preprocessing depending on compiler
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(main PRIVATE -cpp)
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    target_compile_options(main PRIVATE -fpp)
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
    target_compile_options(main PRIVATE -Mpreprocess)
endif()

# =========================
# Define macros
# =========================
if(USE_IBM)
    target_compile_definitions(main PRIVATE USE_IBM)
endif()

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
        target_compile_options(main PRIVATE ${OPENMP_OFFLOAD_FLAGS_LIST})
        target_link_options(main PRIVATE ${OPENMP_OFFLOAD_FLAGS_LIST})
    endif()
endif()


# =========================
# FFT backend
# =========================

if(USE_CUFFT)
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
