cmake_minimum_required(VERSION 3.10)
project(Zephyrgebra C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Create a single static library for the project
add_library(Zephyrgebra STATIC src/Vector.c src/Matrix.c src/Polynomials.c
        "tests/Testing Suite.c")

# Expose public headers to users of the library
target_include_directories(Zephyrgebra PUBLIC ${PROJECT_SOURCE_DIR}/include)

# Enable testing
enable_testing()

# Create test executables
add_executable(overall_test "tests/Testing Suite.c")

# Link the library to the test executables
target_link_libraries(overall_test Zephyrgebra m)

# Register the tests with CTest
add_test(NAME zephyr_test COMMAND overall_test)

set_source_files_properties(src/Polynomials.c PROPERTIES LANGUAGE C)
