#ifndef ZEPHYR_VECTOR_H
#define ZEPHYR_VECTOR_H

#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

    typedef struct zephyr_vector {
        double *data;
        size_t size;
    } zephyr_vector;

    // Memory management
    zephyr_vector *create_vector(size_t size);
    void destroy_vector(zephyr_vector *vector);
    zephyr_vector *duplicate_vector(const zephyr_vector *vector);
    zephyr_vector *resize(const zephyr_vector *old, size_t size, double fill_with_value, bool last_element_only);
    bool normalize_vector(const struct zephyr_vector *v);
    // Initialization
    void vector_fill(zephyr_vector *v, double value);
    zephyr_vector *vector_zero(size_t size);

    // Arithmetic operations
    zephyr_vector *vector_sum(const zephyr_vector *a, const zephyr_vector *b);
    zephyr_vector *vector_subtract(const zephyr_vector *a, const zephyr_vector *b);
    zephyr_vector *scalar_multiply(const zephyr_vector *v, double scalar);

    // Broadcasting support
    zephyr_vector *vector_scalar_sum_broadcasting(const zephyr_vector *v, double scalar);
    zephyr_vector *vector_scalar_subtract_broadcasting(const zephyr_vector *v, double scalar);

    // Algebraic computations
    double magnitude(const zephyr_vector *v);
    zephyr_vector *unit_vector(const zephyr_vector *v);
    double dot_product(const zephyr_vector *a, const zephyr_vector *b);
    double angle(const zephyr_vector *a, const zephyr_vector *b);
    double angle_degrees(const zephyr_vector *a, const zephyr_vector *b);
    zephyr_vector *projection_of_v1_to_v2(const zephyr_vector *a, const zephyr_vector *b);

    // Utility
    void print_vector(const zephyr_vector *v, int decimal_places);
    bool vector_equals(const zephyr_vector *a, const zephyr_vector *b, double epsilon);
    bool is_right_angle_triangle(const zephyr_vector *A, const zephyr_vector *B, const zephyr_vector *C);

#ifdef __cplusplus
}
#endif

#endif // ZEPHYR_VECTOR_H
