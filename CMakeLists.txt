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

set(MODULE_SOURCES
    ${MODULE_DIR}/init.f90
    ${MODULE_DIR}/blocks.f90
    ${MODULE_DIR}/chron.f90
    ${MODULE_DIR}/boundary.f90
    ${MODULE_DIR}/ibm.f90
    ${MODULE_DIR}/geometry_stl.f90
    ${MODULE_DIR}/walldist.f90
    ${MODULE_DIR}/turbulence.f90
    ${MODULE_DIR}/les.f90
    ${MODULE_DIR}/config.f90
    ${MODULE_DIR}/gpu_runtime.f90
    ${MODULE_DIR}/field_hdf5.c
    ${MODULE_DIR}/io.f90
    ${MODULE_DIR}/rans.f90
    ${MODULE_DIR}/bodyforce.f90
    ${MODULE_DIR}/comm.f90
    ${MODULE_DIR}/step.f90
    ${MODULE_DIR}/pressure_solver.f90
    ${MODULE_DIR}/flow_case_base.f90
    ${MODULE_DIR}/flow/case_config_helpers.f90
    ${MODULE_DIR}/flow/generic_flow.f90
    ${MODULE_DIR}/flow/channel/channel_profile.f90
    ${MODULE_DIR}/flow/channel/channel_stats.f90
    ${MODULE_DIR}/flow/channel/channel_flow.f90
    ${MODULE_DIR}/flow/airfoil/airfoil_flow.f90
    ${MODULE_DIR}/flow/boundaryLayer/boundarylayer_profile.f90
    ${MODULE_DIR}/flow/boundaryLayer/boundarylayer_stats.f90
    ${MODULE_DIR}/flow/boundaryLayer/boundarylayer_flow.f90
    ${MODULE_DIR}/flow_case.f90
)

set(SOLVER_TARGETS moby_solve moby_prepare walldist_test transition_test leaftable_test)

# The prepare/solve pair (docs/prepare_solve_strategy.md): moby_prepare
# writes the block-table case file, moby_solve runs it ([ibm] coeff_file).
# mobygrid was absorbed in P3 -- the case file carries the grid datasets.
add_executable(moby_solve ${SRC_DIR}/moby_solve.f90 ${MODULE_SOURCES})
add_executable(moby_prepare ${SRC_DIR}/moby_prepare.f90 ${MODULE_SOURCES})
# Compatibility name: every committed script/tutorial launches "main"
# (rm first -- older build dirs hold a regular file of that name).
add_custom_command(TARGET moby_solve POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E rm -f ${CMAKE_BINARY_DIR}/main
    COMMAND ${CMAKE_COMMAND} -E create_symlink moby_solve ${CMAKE_BINARY_DIR}/main)
# Gate driver for the geometry-agnostic wall distance (IDDES T1b gate b).
add_executable(walldist_test ${SRC_DIR}/test_walldist.f90 ${MODULE_SOURCES})
# Host-side unit test for the T4 transition correlations (rans.f90).
add_executable(transition_test ${SRC_DIR}/test_transition.f90 ${MODULE_SOURCES})
# Leaf-table builder dump vs mobygeom's Python mirror (R2D-0 gate b).
add_executable(leaftable_test ${SRC_DIR}/test_leaftable.f90 ${MODULE_SOURCES})

# =========================
# Module directory
# =========================
foreach(target IN LISTS SOLVER_TARGETS)
    set_target_properties(${target} PROPERTIES
        Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod/${target}
    )
endforeach()

# =========================
# Compiler flags
# =========================
foreach(target IN LISTS SOLVER_TARGETS)
    target_compile_options(${target} PRIVATE
        $<$<COMPILE_LANGUAGE:Fortran>:-O2>
        $<$<COMPILE_LANGUAGE:C>:-O2>
    )
endforeach()

# Enable preprocessing depending on compiler.
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    foreach(target IN LISTS SOLVER_TARGETS)
        target_compile_options(${target} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-cpp>)
    endforeach()
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    foreach(target IN LISTS SOLVER_TARGETS)
        target_compile_options(${target} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-fpp>)
    endforeach()
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
    foreach(target IN LISTS SOLVER_TARGETS)
        target_compile_options(${target} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-Mpreprocess>)
    endforeach()
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

foreach(target IN LISTS SOLVER_TARGETS)
    target_include_directories(${target} PRIVATE ${HDF5_C_INCLUDE_DIR})
    target_link_libraries(${target} PRIVATE ${HDF5_C_LIB})
endforeach()

if(MPI_FORTRAN_MODULE_DIR)
    foreach(target IN LISTS SOLVER_TARGETS)
        target_include_directories(${target} PRIVATE ${MPI_FORTRAN_MODULE_DIR})
    endforeach()
endif()

separate_arguments(MPI_WRAPPER_COMPILE_FLAGS_LIST NATIVE_COMMAND "${MPI_WRAPPER_COMPILE_FLAGS}")
separate_arguments(MPI_WRAPPER_LINK_FLAGS_LIST NATIVE_COMMAND "${MPI_WRAPPER_LINK_FLAGS}")

foreach(flag IN LISTS MPI_WRAPPER_COMPILE_FLAGS_LIST)
    foreach(target IN LISTS SOLVER_TARGETS)
        target_compile_options(${target} PRIVATE
            $<$<COMPILE_LANGUAGE:Fortran>:${flag}>
            $<$<COMPILE_LANGUAGE:C>:${flag}>
        )
    endforeach()
endforeach()

if(MPI_WRAPPER_LINK_FLAGS_LIST)
    foreach(target IN LISTS SOLVER_TARGETS)
        target_link_libraries(${target} PRIVATE ${MPI_WRAPPER_LINK_FLAGS_LIST})
    endforeach()
endif()

# =========================
# Optional features
# =========================
if(USE_OPENMP_OFFLOAD)
    foreach(target IN LISTS SOLVER_TARGETS)
        target_compile_definitions(${target} PRIVATE USE_OPENMP_OFFLOAD)
    endforeach()

    find_package(OpenMP REQUIRED COMPONENTS Fortran)
    foreach(target IN LISTS SOLVER_TARGETS)
        target_link_libraries(${target} PRIVATE OpenMP::OpenMP_Fortran)
    endforeach()

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
            foreach(target IN LISTS SOLVER_TARGETS)
                target_compile_options(${target} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:${flag}>)
            endforeach()
        endforeach()
        foreach(target IN LISTS SOLVER_TARGETS)
            target_link_options(${target} PRIVATE ${OPENMP_OFFLOAD_FLAGS_LIST})
        endforeach()
    endif()
endif()

if(USE_IBM_SECONDORDER)
    foreach(target IN LISTS SOLVER_TARGETS)
        target_compile_definitions(${target} PRIVATE USE_IBM_SECONDORDER)
    endforeach()
endif()
