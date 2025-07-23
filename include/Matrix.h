#ifndef MATRIX_H
#define MATRIX_H

#include <stddef.h>
#include "Vector.h"

typedef struct zephyr_matrix {
    size_t m;
    size_t n;
    double * data;
} zephyr_matrix;

typedef enum {
    ROTATE_X,
    ROTATE_Y,
    ROTATE_Z
} rotation_axis;


typedef struct {
    zephyr_vector *values;       // Eigenvalues (size n)
    zephyr_matrix *vectors;      // Eigenvectors (n x n, column-major)
} zephyr_eigen_result;

zephyr_matrix * create_matrix(const size_t m, const size_t n);
zephyr_matrix * create_identity_matrix(const size_t m);
zephyr_matrix * copy_matrix(const zephyr_matrix * matrix);
zephyr_matrix * matrix_set(const size_t row, const size_t col, const double value, zephyr_matrix * matrix);
zephyr_matrix * add_matrices(const zephyr_matrix * matricxA, const zephyr_matrix * matricxB);
zephyr_matrix * subtract_matrices(const zephyr_matrix * matricxA, const zephyr_matrix * matricxB);
zephyr_matrix * scalar_multiply_matrix(const double scale, const zephyr_matrix * matrix);
zephyr_matrix * matrix_multiplication(const zephyr_matrix * a, const zephyr_matrix * b);
zephyr_matrix * transpose(const zephyr_matrix * matrix);
void destroy_matrix(zephyr_matrix * matrix);
zephyr_matrix * matrix_col_vector_mult(const zephyr_vector * zephyr_vector, const zephyr_matrix * matrix);
zephyr_matrix * matrix_row_vector_mult(const zephyr_vector * row_vector, const zephyr_matrix * matrix);
void print_matrix(const zephyr_matrix *matrix, const int decimal_places);
zephyr_matrix * matrix_get_row(const size_t row, const zephyr_matrix * matrix);
zephyr_matrix * matrix_get_col(const size_t col, const zephyr_matrix * matrix);
zephyr_matrix * slice(const size_t row_start, const size_t row_end, const size_t col_start, const size_t col_end, const zephyr_matrix * matrix);
zephyr_matrix * matrix_from_array(const double * data, const size_t m, const size_t n);
zephyr_matrix * matrix_from_vector(const zephyr_vector * vector);
zephyr_vector * vector_from_matrix(const zephyr_matrix * matrix);
zephyr_matrix * matrix_from_array_stride(const double * data, const size_t m, const size_t n, const size_t stride);
zephyr_matrix *minor_matrix(const zephyr_matrix *mat, const size_t skip_row, const size_t skip_col);
double determinant(const zephyr_matrix * matrix);
zephyr_matrix * cofactor_matrix(const zephyr_matrix *mat);
zephyr_matrix * adjoint_matrix(const zephyr_matrix *mat);
zephyr_matrix * matrix_inverse(const zephyr_matrix *mat);
zephyr_vector * simultaneous_eq(const zephyr_matrix * coefficients, const zephyr_matrix * solution);
zephyr_matrix *vandermonde_matrix(const double *x_values, size_t n);
bool eigen_solve(const zephyr_matrix *A, zephyr_eigen_result *result);
void destroy_eigen_result(zephyr_eigen_result *result);
zephyr_matrix *rotation_matrix_2d(const zephyr_vector *to_rotate, double theta);
zephyr_matrix *rotate_3d(const zephyr_vector *v, double theta, rotation_axis axis);
zephyr_matrix *rotate_about_axis(const zephyr_vector *v, const zephyr_vector *axis, double theta);
#endif //MATRIX_H