cmake_minimum_required(VERSION 3.10)
project(Zephyrgebra C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Include directories
include_directories(include)

# Source files
set(ZEPHYRGEBRA_SOURCES
    src/Matrix.c
    src/Vector.c
    src/Polynomials.c
)

# Add the library
add_library(Zephyrgebra ${ZEPHYRGEBRA_SOURCES})

# Link math library
target_link_libraries(Zephyrgebra PRIVATE m)

# Add test executable
add_executable(matrix_test_basic test/Matrix/matrix_test_basic.c)

# Link Zephyrgebra library to the test executable
target_link_libraries(matrix_test_basic PRIVATE Zephyrgebra)