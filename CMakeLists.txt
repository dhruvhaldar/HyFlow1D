cmake_minimum_required(VERSION 3.10)
project(HyFlow1D LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Compiler Hardening / Warnings
if(MSVC)
    add_compile_options(/W4 /WX)
else()
    add_compile_options(-Wall -Wextra -Wpedantic -Werror)

    # Security Hardening Flags
    add_compile_options(
        -fstack-protector-strong
        -Wformat
        -Wformat-security
    )

    # _FORTIFY_SOURCE requires optimization to work correctly.
    add_compile_definitions($<$<CONFIG:Release>:_FORTIFY_SOURCE=2>)
    add_compile_definitions($<$<CONFIG:RelWithDebInfo>:_FORTIFY_SOURCE=2>)
    add_compile_definitions($<$<CONFIG:MinSizeRel>:_FORTIFY_SOURCE=2>)

    # Enforce Position Independent Code/Executable (PIE/PIC)
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)

    # Linker Hardening: Full RELRO
    if(CMAKE_VERSION VERSION_LESS 3.13)
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-z,relro,-z,now")
    else()
        add_link_options(-Wl,-z,relro,-z,now)
    endif()
endif()

# OpenMP Multithreading Support
option(ENABLE_OPENMP "Enable OpenMP shared-memory multithreading" ON)
if(ENABLE_OPENMP)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        message(STATUS "OpenMP found: Enabling multithreading support")
        add_compile_definitions(HYFLOW_ENABLE_OPENMP)
        link_libraries(OpenMP::OpenMP_CXX)
    endif()
endif()

# MPI Distributed-Memory Support
option(ENABLE_MPI "Enable MPI distributed parallelism" ON)
if(ENABLE_MPI)
    find_package(MPI)
    if(MPI_CXX_FOUND)
        message(STATUS "MPI found: Enabling distributed domain decomposition")
        add_compile_definitions(HYFLOW_ENABLE_MPI)
        include_directories(SYSTEM ${MPI_CXX_INCLUDE_PATH})
        link_libraries(MPI::MPI_CXX)
    endif()
endif()

# Add include directory
include_directories(include)

# HyFlow1D core sources
set(HYFLOW_CORE_SOURCES
    src/numerics.cpp
    src/fv_solver.cpp
    src/dg_solver.cpp
    src/hybrid_coupling.cpp
)

if(MPI_CXX_FOUND)
    list(APPEND HYFLOW_CORE_SOURCES src/mpi_hybrid_coupling.cpp)
endif()

# Main executable
add_executable(hyflow1d src/main.cpp ${HYFLOW_CORE_SOURCES})

# Enable testing
enable_testing()

# Existing tests
add_executable(test_numerics tests/test_numerics.cpp src/numerics.cpp)
add_test(NAME numerics_test COMMAND test_numerics)

add_executable(test_solvers tests/test_solvers.cpp ${HYFLOW_CORE_SOURCES})
add_test(NAME solvers_test COMMAND test_solvers)

add_executable(test_validation tests/test_validation.cpp ${HYFLOW_CORE_SOURCES})
add_test(NAME validation_test COMMAND test_validation)

add_executable(test_security tests/test_security.cpp ${HYFLOW_CORE_SOURCES})
add_test(NAME security_test COMMAND test_security)

add_executable(test_hybrid_security tests/test_hybrid_security.cpp ${HYFLOW_CORE_SOURCES})
add_test(NAME hybrid_security_test COMMAND test_hybrid_security)

# New High-Order & HPC tests
add_executable(test_ssprk3 tests/test_ssprk3.cpp ${HYFLOW_CORE_SOURCES})
add_test(NAME ssprk3_test COMMAND test_ssprk3)

add_executable(test_openmp tests/test_openmp.cpp ${HYFLOW_CORE_SOURCES})
add_test(NAME openmp_test COMMAND test_openmp)

if(MPI_CXX_FOUND)
    add_executable(test_mpi tests/test_mpi.cpp ${HYFLOW_CORE_SOURCES})
    add_test(NAME mpi_test COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 $<TARGET_FILE:test_mpi>)
endif()

# Benchmarks
add_executable(benchmark_dg tests/benchmark_dg.cpp src/dg_solver.cpp src/numerics.cpp)
add_executable(benchmark_update tests/benchmark_update.cpp src/dg_solver.cpp src/numerics.cpp)
add_executable(benchmark_step tests/benchmark_step.cpp src/dg_solver.cpp src/numerics.cpp)
add_executable(benchmark_plot tests/benchmark_plot.cpp src/dg_solver.cpp src/numerics.cpp)
add_executable(benchmark_fv tests/benchmark_fv.cpp src/fv_solver.cpp)
