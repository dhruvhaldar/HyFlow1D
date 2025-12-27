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
    # It provides runtime buffer overflow checking.
    # We only enable it for Release builds to avoid warnings/errors in Debug (-O0).
    add_compile_definitions($<$<CONFIG:Release>:_FORTIFY_SOURCE=2>)
    add_compile_definitions($<$<CONFIG:RelWithDebInfo>:_FORTIFY_SOURCE=2>)
    add_compile_definitions($<$<CONFIG:MinSizeRel>:_FORTIFY_SOURCE=2>)

    # Enforce Position Independent Code/Executable (PIE/PIC)
    # This ensures ASLR effectiveness.
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)

    # Linker Hardening
    # -z relro -z now: Full RELRO (Read-Only Relocations)
    if(CMAKE_VERSION VERSION_LESS 3.13)
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-z,relro,-z,now")
    else()
        add_link_options(-Wl,-z,relro,-z,now)
    endif()
endif()

# Add include directory
include_directories(include)

# Create the executable
add_executable(hyflow1d 
    src/main.cpp 
    src/numerics.cpp 
    src/fv_solver.cpp 
    src/dg_solver.cpp 
    src/hybrid_coupling.cpp
)

# Enable testing
enable_testing()

# Add tests
add_executable(test_numerics tests/test_numerics.cpp src/numerics.cpp)
add_test(NAME numerics_test COMMAND test_numerics)

add_executable(test_solvers tests/test_solvers.cpp src/fv_solver.cpp src/dg_solver.cpp src/numerics.cpp)
add_test(NAME solvers_test COMMAND test_solvers)

add_executable(test_validation tests/test_validation.cpp src/fv_solver.cpp src/dg_solver.cpp src/numerics.cpp)
add_test(NAME validation_test COMMAND test_validation)
