#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define RAD_TO_DEGREE (180.0 / M_PI)

struct zephyr_vector {
    double * data;
    size_t size;
};

// Memory Management
struct zephyr_vector * create_vector(const size_t size){
    struct zephyr_vector * vector = malloc(sizeof(struct zephyr_vector));
    if(vector == NULL) return NULL;
    vector->data = malloc(sizeof(double) * size);
    if(vector->data == NULL) {
        free(vector);
        return NULL;
    }
    vector->size = size;
    return vector;
}

void destroy_vector(struct zephyr_vector * vector) {
    if (vector != NULL) {
        free(vector->data);
        free(vector);
    }
}

struct zephyr_vector * duplicate_vector(const struct zephyr_vector * vector) {
    if (vector == NULL) return NULL;
    struct zephyr_vector * result = create_vector(vector->size);
    if (result == NULL) return NULL;
    for (size_t i = 0; i < vector->size; i++) {
        result->data[i] = vector->data[i];
    }
    return result;
}

struct zephyr_vector * resize(const struct zephyr_vector * old, const size_t size, const double fill_with_value, const bool last_element_only) {
    if (!old) return NULL;
    struct zephyr_vector * v = create_vector(size);
    if (!v) return NULL;
    size_t limit = old->size < size ? old->size : size;
    for (size_t i = 0; i < limit; i++) {
        v->data[i] = old->data[i];
    }
    if (size > old->size) {
        if (last_element_only) {
            v->data[size - 1] = fill_with_value;
        } else {
            for (size_t i = old->size; i < size; i++) {
                v->data[i] = fill_with_value;
            }
        }
    }
    return v;
}

// Initialization
void vector_fill(struct zephyr_vector *v, const double value) {
    if (!v) return;
    for (size_t i = 0; i < v->size; i++) {
        v->data[i] = value;
    }
}

struct zephyr_vector * vector_zero(const size_t size) {
    struct zephyr_vector * v = create_vector(size);
    if (!v) return NULL;
    vector_fill(v, 0.0);
    return v;
}

struct zephyr_vector *vector_from_array(double data[],size_t size) {
    struct zephyr_vector * v = create_vector(size);
    if (!v) return NULL;
    memcpy(v->data, data, size * sizeof(double));
    return v;
}

// Arithmetic Operations
struct zephyr_vector * vector_sum(const struct zephyr_vector * a, const struct zephyr_vector * b) {
    if (!a || !b || a->size != b->size) return NULL;
    struct zephyr_vector * result = create_vector(a->size);
    if (!result) return NULL;
    for (size_t i = 0; i < a->size; i++) {
        result->data[i] = a->data[i] + b->data[i];
    }
    return result;
}

struct zephyr_vector * vector_subtract(const struct zephyr_vector * a, const struct zephyr_vector * b) {
    if (!a || !b || a->size != b->size) return NULL;
    struct zephyr_vector * result = create_vector(a->size);
    if (!result) return NULL;
    for (size_t i = 0; i < a->size; i++) {
        result->data[i] = a->data[i] - b->data[i];
    }
    return result;
}

struct zephyr_vector * scalar_multiply(const struct zephyr_vector * v, const double scalar) {
    if (!v) return NULL;
    struct zephyr_vector * result = create_vector(v->size);
    if (!result) return NULL;
    for (size_t i = 0; i < v->size; i++) {
        result->data[i] = v->data[i] * scalar;
    }
    return result;
}

// Broadcasting
struct zephyr_vector * vector_scalar_sum_broadcasting(const struct zephyr_vector * v, const double scalar) {
    if (!v) return NULL;
    struct zephyr_vector * temp = create_vector(v->size);
    if (!temp) return NULL;
    for (size_t i = 0; i < v->size; i++) temp->data[i] = scalar;
    struct zephyr_vector * result = vector_sum(temp, v);
    destroy_vector(temp);
    return result;
}

struct zephyr_vector * vector_scalar_subtract_broadcasting(const struct zephyr_vector * v, const double scalar) {
    if (!v) return NULL;
    struct zephyr_vector * temp = create_vector(v->size);
    if (!temp) return NULL;
    for (size_t i = 0; i < v->size; i++) temp->data[i] = scalar;
    struct zephyr_vector * result = vector_subtract(temp, v);
    destroy_vector(temp);
    return result;
}

// Algebra

double magnitude(const struct zephyr_vector * v) {
    double sum = 0;
    for (size_t i = 0; i < v->size; i++) {
        sum += v->data[i] * v->data[i];
    }
    return sqrt(sum);
}

struct zephyr_vector * unit_vector(const struct zephyr_vector * v) {
    double mag = magnitude(v);
    if (mag == 0) return NULL;
    return scalar_multiply(v, 1.0 / mag);
}

double dot_product(const struct zephyr_vector * a, const struct zephyr_vector * b) {
    if (!a || !b || a->size != b->size) return NAN;
    double sum = 0;
    for (size_t i = 0; i < a->size; i++) sum += a->data[i] * b->data[i];
    return sum;
}

bool normalize_vector(const struct zephyr_vector *v) {
    if (!v) return false;
    double norm = sqrt(dot_product(v, v));
    if (norm < 1e-12) return false;
    for (size_t i = 0; i < v->size; i++) {
        v->data[i] /= norm;
    }
    return true;
}



double angle(const struct zephyr_vector * a, const struct zephyr_vector * b) {
    if (!a || !b) return NAN;
    double mag_a = magnitude(a);
    double mag_b = magnitude(b);
    if (mag_a == 0 || mag_b == 0) return NAN;
    double cosine = dot_product(a, b) / (mag_a * mag_b);
    return acos(fmax(-1.0, fmin(1.0, cosine)));
}

double angle_degrees(const struct zephyr_vector * a, const struct zephyr_vector * b) {
    return angle(a, b) * RAD_TO_DEGREE;
}

struct zephyr_vector * projection_of_v1_to_v2(const struct zephyr_vector * a, const struct zephyr_vector * b) {
    if (!a || !b) return NULL;
    double mag_b = magnitude(b);
    if (mag_b == 0) return NULL;
    double scale = dot_product(a, b) / (mag_b * mag_b);
    return scalar_multiply(b, scale);
}

// Utilities
void print_vector(const struct zephyr_vector *v, const int dp) {
    if (!v) {
        printf("NULL vector\n");
        return;
    }
    printf("[");
    for (size_t i = 0; i < v->size; i++) {
        char format[10];
        snprintf(format, sizeof(format), "%%.%df", dp);
        printf(format, v->data[i]);
        if (i < v->size - 1) printf(", ");
    }
    printf("]\n");
}

bool vector_equals(const struct zephyr_vector * a, const struct zephyr_vector * b, const double epsilon) {
    if (!a || !b || a->size != b->size) return false;
    for (size_t i = 0; i < a->size; i++) {
        if (fabs(a->data[i] - b->data[i]) > epsilon) return false;
    }
    return true;
}

bool is_right_angle_triangle(const struct zephyr_vector *A,
    const struct zephyr_vector *B,
    const struct zephyr_vector *C)
{
    struct zephyr_vector *AB = vector_subtract(B, A);
    struct zephyr_vector *BC = vector_subtract(C, B);
    struct zephyr_vector *CA = vector_subtract(A, C);
    const bool right =
        fabs(dot_product(AB, BC)) < 1e-6 ||
        fabs(dot_product(BC, CA)) < 1e-6 ||
        fabs(dot_product(CA, AB)) < 1e-6;

    destroy_vector(AB);
    destroy_vector(BC);
    destroy_vector(CA);
    return right;
}
