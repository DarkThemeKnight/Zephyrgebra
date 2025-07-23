# Zephyr Vector Library Documentation

## Overview

The Zephyr Vector Library is a lightweight C library for performing vector operations in Euclidean space. It supports dynamic-size vectors and provides utility functions for algebraic computations, broadcasting, and geometric operations.

## Data Structure

```c
typedef struct zephyr_vector {
    double *data;     // Dynamic array of vector components
    size_t size;      // Dimension of the vector
} zephyr_vector;
```

---

## Memory Management

* `zephyr_vector *create_vector(size_t size);`

    * Allocates and returns a new vector of given size.

* `void destroy_vector(zephyr_vector *vector);`

    * Frees memory allocated for a vector.

* `zephyr_vector *duplicate_vector(const zephyr_vector *vector);`

    * Creates a deep copy of a vector.

* `zephyr_vector *resize(const zephyr_vector *old, size_t size, double fill_value, bool last_element_only);`

    * Creates a resized copy of a vector, optionally filling remaining elements.

---

## Initialization

* `zephyr_vector *vector_zero(size_t size);`

    * Returns a zero-filled vector of given size.

* `void vector_fill(zephyr_vector *v, double value);`

    * Fills all components of a vector with the given value.

---

## Arithmetic Operations

* `zephyr_vector *vector_sum(const zephyr_vector *a, const zephyr_vector *b);`
* `zephyr_vector *vector_subtract(const zephyr_vector *a, const zephyr_vector *b);`
* `zephyr_vector *scalar_multiply(const zephyr_vector *v, double scalar);`

> These functions perform element-wise addition, subtraction, or scalar multiplication. Vectors must be of equal size.

---

## Broadcasting Support

* `zephyr_vector *vector_scalar_sum_broadcasting(const zephyr_vector *v, double scalar);`
* `zephyr_vector *vector_scalar_subtract_broadcasting(const zephyr_vector *v, double scalar);`

> Adds/subtracts a scalar as if it were a vector of the same size.

---

## Algebraic Computations

* `double magnitude(const zephyr_vector *v);`

    * Returns the Euclidean norm.

* `zephyr_vector *unit_vector(const zephyr_vector *v);`

    * Returns the normalized vector.

* `double dot_product(const zephyr_vector *a, const zephyr_vector *b);`

    * Returns the dot product of two vectors.

* `double angle(const zephyr_vector *a, const zephyr_vector *b);`

    * Returns the angle in radians.

* `double angle_degrees(const zephyr_vector *a, const zephyr_vector *b);`

    * Returns the angle in degrees.

* `zephyr_vector *projection_of_v1_to_v2(const zephyr_vector *a, const zephyr_vector *b);`

    * Returns the projection of `a` onto `b`.

---

## Utility Functions

* `void print_vector(const zephyr_vector *v, int decimal_places);`

    * Prints a vector with configurable precision.

* `bool vector_equals(const zephyr_vector *a, const zephyr_vector *b, double epsilon);`

    * Checks if vectors are equal within a given tolerance.

* `bool is_right_angle_triangle(const zephyr_vector *A, const zephyr_vector *B, const zephyr_vector *C);`

    * Verifies if three points form a right-angled triangle.

---

## Notes

* Uses dynamic memory. Always `destroy_vector` to avoid leaks.
* Angles are computed with clamped cosine values to avoid domain errors.
* Math functions require linking with `-lm`.

---

## Build & Link

### CMake Example:

```cmake
add_library(zephyr_vector STATIC src/Vector.c)
target_include_directories(zephyr_vector PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(your_target zephyr_vector m)
```

---

## License

MIT or your preferred open license.

---

## Author

Omotola David Ayanfeoluwa (2025)

---
