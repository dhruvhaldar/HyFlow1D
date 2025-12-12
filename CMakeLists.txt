cmake_minimum_required(VERSION 3.10)
project(HyFlow1D LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

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
