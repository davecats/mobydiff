cmake_minimum_required(VERSION 3.20)

project(3D_NSE_SOLVER LANGUAGES Fortran C)

# =========================
# Options
# =========================
option(USE_OPENMP_OFFLOAD "Enable OpenMP target offload for accelerator execution" OFF)
option(USE_IBM_SECONDORDER "Use second-order IBM interpolation/correction" ON)
set(OPENMP_OFFLOAD_FLAGS "" CACHE STRING "Additional Fortran compiler/linker flags for OpenMP target offload")
set(MPI_WRAPPER_COMPILE_FLAGS "" CACHE STRING "MPI wrapper compile flags")
set(MPI_WRAPPER_LINK_FLAGS "" CACHE STRING "MPI wrapper link flags")
set(MPI_FORTRAN_MODULE_DIR "" CACHE PATH "Directory containing mpi_f08.mod")

# =========================
# Sources
# =========================
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(MODULE_DIR ${SRC_DIR}/modules)

set(SOURCES
    ${SRC_DIR}/main.f90
    ${MODULE_DIR}/init.f90
    ${MODULE_DIR}/chron.f90
    ${MODULE_DIR}/config.f90
    ${MODULE_DIR}/gpu_runtime.f90
    ${MODULE_DIR}/boundary.f90
    ${MODULE_DIR}/ibm.f90
    ${MODULE_DIR}/field_hdf5.c
    ${MODULE_DIR}/io.f90
    ${MODULE_DIR}/comm.f90
    ${MODULE_DIR}/step.f90
    ${MODULE_DIR}/pressure_solver.f90
    ${MODULE_DIR}/flow_case_base.f90
    ${MODULE_DIR}/flow/case_config_helpers.f90
    ${MODULE_DIR}/flow/generic_flow.f90
    ${MODULE_DIR}/flow/channel/channel_profile.f90
    ${MODULE_DIR}/flow/channel/channel_stats.f90
    ${MODULE_DIR}/flow/channel/channel_flow.f90
    ${MODULE_DIR}/flow_case.f90
)

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

# Enable preprocessing depending on compiler.
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(main PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-cpp>)
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    target_compile_options(main PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-fpp>)
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
    target_compile_options(main PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-Mpreprocess>)
endif()

# =========================
# Parallel HDF5 field output
# =========================
set(HDF5_ROOT "" CACHE PATH "Root directory of HDF5 installation")

find_path(HDF5_C_INCLUDE_DIR
    NAMES hdf5.h
    HINTS
        ${HDF5_ROOT}
        /usr/include/hdf5/openmpi
        /usr/lib/x86_64-linux-gnu/hdf5/openmpi
    PATH_SUFFIXES include
)

find_library(HDF5_C_LIB
    NAMES hdf5
    HINTS
        ${HDF5_ROOT}
        /usr/lib/x86_64-linux-gnu/hdf5/openmpi
    PATH_SUFFIXES lib lib64
)

if(NOT HDF5_C_INCLUDE_DIR)
    message(FATAL_ERROR "Could not find parallel HDF5 hdf5.h. Set -DHDF5_ROOT=/path/to/parallel-hdf5")
endif()

if(NOT HDF5_C_LIB)
    message(FATAL_ERROR "Could not find parallel HDF5 C library. Set -DHDF5_ROOT=/path/to/parallel-hdf5")
endif()

message(STATUS "HDF5 include dir: ${HDF5_C_INCLUDE_DIR}")
message(STATUS "HDF5 C library: ${HDF5_C_LIB}")
message(STATUS "Pressure solver: red-black SOR")
message(STATUS "MPI support: required")
message(STATUS "MPI wrapper compile flags: ${MPI_WRAPPER_COMPILE_FLAGS}")
message(STATUS "MPI wrapper link flags: ${MPI_WRAPPER_LINK_FLAGS}")
if(MPI_FORTRAN_MODULE_DIR)
    message(STATUS "MPI Fortran module dir: ${MPI_FORTRAN_MODULE_DIR}")
endif()

if(NOT MPI_WRAPPER_LINK_FLAGS)
    message(FATAL_ERROR "MPI_WRAPPER_LINK_FLAGS is required")
endif()

target_include_directories(main PRIVATE ${HDF5_C_INCLUDE_DIR})
target_link_libraries(main PRIVATE ${HDF5_C_LIB})

if(MPI_FORTRAN_MODULE_DIR)
    target_include_directories(main PRIVATE ${MPI_FORTRAN_MODULE_DIR})
endif()

separate_arguments(MPI_WRAPPER_COMPILE_FLAGS_LIST NATIVE_COMMAND "${MPI_WRAPPER_COMPILE_FLAGS}")
separate_arguments(MPI_WRAPPER_LINK_FLAGS_LIST NATIVE_COMMAND "${MPI_WRAPPER_LINK_FLAGS}")

foreach(flag IN LISTS MPI_WRAPPER_COMPILE_FLAGS_LIST)
    target_compile_options(main PRIVATE
        $<$<COMPILE_LANGUAGE:Fortran>:${flag}>
        $<$<COMPILE_LANGUAGE:C>:${flag}>
    )
endforeach()

target_link_libraries(main PRIVATE ${MPI_WRAPPER_LINK_FLAGS_LIST})

# =========================
# Optional features
# =========================
if(USE_OPENMP_OFFLOAD)
    target_compile_definitions(main PRIVATE USE_OPENMP_OFFLOAD)

    find_package(OpenMP REQUIRED COMPONENTS Fortran)
    target_link_libraries(main PRIVATE OpenMP::OpenMP_Fortran)

    if(NOT OPENMP_OFFLOAD_FLAGS)
        if(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
            set(OPENMP_OFFLOAD_FLAGS "-mp=gpu" CACHE STRING "Additional Fortran compiler/linker flags for OpenMP target offload" FORCE)
        elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
            set(OPENMP_OFFLOAD_FLAGS "-fopenmp -foffload=nvptx-none" CACHE STRING "Additional Fortran compiler/linker flags for OpenMP target offload" FORCE)
        elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
            set(OPENMP_OFFLOAD_FLAGS "-fiopenmp -fopenmp-targets=spir64" CACHE STRING "Additional Fortran compiler/linker flags for OpenMP target offload" FORCE)
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

if(USE_IBM_SECONDORDER)
    target_compile_definitions(main PRIVATE USE_IBM_SECONDORDER)
endif()
