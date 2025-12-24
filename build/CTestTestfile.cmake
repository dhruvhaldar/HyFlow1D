# CMake generated Testfile for
# Source directory: /app
# Build directory: /app/build
#
# This file includes the relevant testing commands required for
# testing this directory and lists subdirectories to be tested as well.
add_test(numerics_test "/app/build/test_numerics")
set_tests_properties(numerics_test PROPERTIES  _BACKTRACE_TRIPLES "/app/CMakeLists.txt;25;add_test;/app/CMakeLists.txt;0;")
add_test(solvers_test "/app/build/test_solvers")
set_tests_properties(solvers_test PROPERTIES  _BACKTRACE_TRIPLES "/app/CMakeLists.txt;28;add_test;/app/CMakeLists.txt;0;")
