//
// Created by omotola-david-ayanfeoluwa on 7/21/25.
//

#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "../include/Vector.h"
#include "../include/Matrix.h"

#define EPSILON 1e-6


static size_t idx(const size_t i, const size_t j, const zephyr_matrix *matrix) {
    assert(i < matrix->m && j < matrix->n);
    return i * matrix->n + j;
}

zephyr_matrix * create_matrix(const size_t m, const size_t n) {
    zephyr_matrix *matrix = malloc(sizeof(zephyr_matrix));
    if (matrix == NULL) {
        return NULL;
    }
    matrix->m = m;
    matrix->n = n;
    matrix->data = calloc(m * n, sizeof(double));
    if (matrix->data == NULL) {
        free(matrix);
        return NULL;
    }
    return matrix;
}

zephyr_matrix *create_identity_matrix(const size_t m) {
    zephyr_matrix *matrix = create_matrix(m, m);
    if (!matrix) return NULL;

    for (size_t i = 0; i < m; i++) {
        matrix->data[i * m + i] = 1.0;
    }
    return matrix;
}

zephyr_matrix *copy_matrix(const zephyr_matrix *matrix) {
    if (matrix == NULL) return NULL;
    zephyr_matrix *copy = create_matrix(matrix->m, matrix->n);
    if (copy == NULL) return NULL;
    memcpy(copy->data, matrix->data, matrix->m * matrix->n * sizeof(double));
    return copy;
}


zephyr_matrix * matrix_set(const size_t row, const size_t col, const double value, zephyr_matrix * matrix) {
    if (!matrix || row >= matrix->m || col >= matrix->n) {
        return NULL;
    }
    const size_t index = idx(row, col, matrix);
    matrix->data[index] = value;
    return matrix;
}

double matrix_get(const size_t row, const size_t col, const zephyr_matrix * matrix) {
    if (!matrix || row >= matrix->m || col >= matrix->n) {
        return NAN;
    }
    const size_t index = idx(row, col, matrix);
    return matrix->data[index];
}

zephyr_matrix * add_matrices(const zephyr_matrix * matricxA, const zephyr_matrix * matricxB) {
    if (matricxA->m != matricxB->m || matricxA->n != matricxB->n) return NULL;
    zephyr_matrix * sum = create_matrix(matricxA->m, matricxA->n);
    if (sum == NULL) return NULL;
    for (size_t i = 0; i < (matricxA->m * matricxA->n); i++) {
        sum->data[i] = matricxA->data[i] + matricxB->data[i];
    }
    return sum;
}

zephyr_matrix * subtract_matrices(const zephyr_matrix * matricxA, const zephyr_matrix * matricxB) {
    if (matricxA->m != matricxB->m || matricxA->n != matricxB->n) return NULL;
    zephyr_matrix * sum = create_matrix(matricxA->m, matricxA->n);
    if (sum == NULL) return NULL;
    for (size_t i = 0; i < (matricxA->m * matricxA->n); i++) {
        sum->data[i] = matricxA->data[i] - matricxB->data[i];
    }
    return sum;
}

zephyr_matrix * scalar_multiply_matrix(const double scale, const zephyr_matrix * matrix) {
    zephyr_matrix *sum = create_matrix(matrix->m, matrix->n);
    if (sum == NULL) return NULL;
    for (size_t i = 0; i < (matrix->m * matrix->n); i++) {
        sum->data[i] = scale * matrix->data[i];
    }
    return sum;
}


static double get_sum_of_product(const zephyr_matrix * a, const zephyr_matrix * b, const size_t j, const size_t k) {
    double sum = 0.0;
    for (size_t r = 0; r < a->n; r++) {
        sum += a->data[idx(j, r, a)] * b->data[idx(r, k, b)];
    }
    return sum;
}

zephyr_matrix *matrix_multiplication(const zephyr_matrix *a, const zephyr_matrix *b) {
    if (a == NULL || b == NULL) return NULL;
    if (a->n != b->m) return NULL;
    zephyr_matrix *product = create_matrix(a->m, b->n);
    if (product == NULL) return NULL;
    for (size_t j = 0; j < product->m; j++) {
        for (size_t k = 0; k < product->n; k++) {
            const double sum = get_sum_of_product(a, b, j, k);

            product->data[idx(j, k, product)] = sum;
        }
    }
    return product;
}

zephyr_matrix *transpose(const zephyr_matrix *matrix) {
    if (matrix == NULL) return NULL;
    zephyr_matrix *transposed = create_matrix(matrix->n, matrix->m);
    if (transposed == NULL) return NULL;
    for (size_t j = 0; j < matrix->m; j++) {
        for (size_t k = 0; k < matrix->n; k++) {
            transposed->data[idx(k, j, transposed)] = matrix->data[idx(j, k, matrix)];
        }
    }
    return transposed;
}

void destroy_matrix(zephyr_matrix *matrix) {
    if (matrix == NULL) return;
    free(matrix->data);
    free(matrix);
}

zephyr_matrix *matrix_col_vector_mult(const zephyr_vector *zephyr_vector, const zephyr_matrix *matrix) {
    if (zephyr_vector == NULL) return NULL;
    if (matrix == NULL) return NULL;
    if (zephyr_vector->size != matrix->n) return NULL;
    zephyr_matrix *to_matrix = create_matrix(zephyr_vector->size, 1);
    memcpy(to_matrix->data, zephyr_vector->data, zephyr_vector->size * sizeof(double));
    zephyr_matrix *product = matrix_multiplication(matrix, to_matrix);
    destroy_matrix(to_matrix);
    return product;
}

zephyr_matrix *matrix_row_vector_mult(const zephyr_vector *row_vector, const zephyr_matrix *matrix) {
    if (row_vector == NULL || matrix == NULL) return NULL;
    if (row_vector->size != matrix->m) return NULL; // rows of matrix must match size of row vector

    // Convert row vector to 1 x n matrix
    zephyr_matrix *as_matrix = create_matrix(1, row_vector->size);
    if (as_matrix == NULL) return NULL;

    memcpy(as_matrix->data, row_vector->data, row_vector->size * sizeof(double));

    // Perform (1 x n) * (n x p) = (1 x p)
    zephyr_matrix *product = matrix_multiplication(as_matrix, matrix);
    destroy_matrix(as_matrix);

    return product;
}


void print_matrix(const zephyr_matrix *matrix, const int decimal_places) {
    if (!matrix || !matrix->data) {
        printf("NULL matrix\n");
        return;
    }
    char format[16];
    snprintf(format, sizeof(format), "%%.%df", decimal_places);
    for (size_t i = 0; i < matrix->m; i++) {
        printf("[ ");
        for (size_t j = 0; j < matrix->n; j++) {
            printf(format, matrix->data[i * matrix->n + j]);
            if (j < matrix->n - 1) {
                printf(", ");
            }
        }
        printf(" ]\n");
    }
}

zephyr_matrix * matrix_get_row(const size_t row, const zephyr_matrix * matrix) {
    if (!matrix || row >= matrix->m) return NULL;
    zephyr_matrix *row_vector = create_matrix(1, matrix->n);
    if (!row_vector) return NULL;
    for (size_t j = 0; j < matrix->n; j++) {
        row_vector->data[j] = matrix->data[idx(row, j, matrix)];
    }
    return row_vector;
}


zephyr_matrix * matrix_get_col(const size_t col, const zephyr_matrix * matrix) {
    if (!matrix || col >= matrix->n) return NULL;
    zephyr_matrix *col_vector = create_matrix(matrix->m, 1);
    if (!col_vector) return NULL;
    for (size_t i = 0; i < matrix->m; i++) {
        col_vector->data[i] = matrix->data[idx(i, col, matrix)];
    }
    return col_vector;
}


zephyr_matrix *slice(const size_t row_start, const size_t row_end, const size_t col_start, const size_t col_end,
                     const zephyr_matrix *matrix) {
    if (!matrix || row_start > row_end || col_start > col_end || row_end >= matrix->m || col_end >= matrix->n) {
        return NULL;
    }
    const size_t new_m = row_end - row_start + 1;
    const size_t new_n = col_end - col_start + 1;
    zephyr_matrix *sliced = create_matrix(new_m, new_n);
    if (!sliced) return NULL;
    for (size_t i = 0; i < new_m; i++) {
        for (size_t j = 0; j < new_n; j++) {
            sliced->data[idx(i, j, sliced)] = matrix->data[idx(row_start + i, col_start + j, matrix)];
        }
    }
    return sliced;
}


zephyr_matrix *matrix_from_array(const double *data, const size_t m, const size_t n) {
    if (!data || m == 0 || n == 0) return NULL;

    zephyr_matrix *matrix = create_matrix(m, n);
    if (!matrix) return NULL;

    memcpy(matrix->data, data, m * n * sizeof(double));

    return matrix;
}


zephyr_matrix *matrix_from_vector(const zephyr_vector *vector) {
    if (!vector) return NULL;
    zephyr_matrix *matrix = create_matrix(vector->size, 1);
    if (!matrix) return NULL;
    memcpy(matrix->data, vector->data, vector->size * sizeof(double));
    return matrix;
}


zephyr_vector *vector_from_matrix(const zephyr_matrix *matrix) {
    if (!matrix || (matrix->m != 1 && matrix->n != 1)) return NULL;
    const size_t size = matrix->m * matrix->n;
    zephyr_vector *vector = create_vector(size);
    if (!vector) return NULL;
    memcpy(vector->data, matrix->data, size * sizeof(double));
    return vector;
}

zephyr_matrix * matrix_from_array_stride(const double * data, const size_t m, const size_t n, const size_t stride) {
    if (!data || m == 0 || n == 0 || stride < n) return NULL;
    zephyr_matrix * matrix = create_matrix(m, n);
    if (!matrix) return NULL;
    for (size_t i = 0; i < m; i++) {
        memcpy(&matrix->data[i * n], &data[i * stride], n * sizeof(double));
    }
    return matrix;
}

zephyr_matrix *minor_matrix(const zephyr_matrix *mat, const size_t skip_row, const size_t skip_col) {
    if (!mat || mat->m != mat->n) return NULL;
    zephyr_matrix *minor = create_matrix(mat->m - 1, mat->n - 1);
    if (!minor) return NULL;
    size_t r = 0;
    for (size_t i = 0; i < mat->m; i++) {
        if (i == skip_row) continue;
        size_t c = 0;
        for (size_t j = 0; j < mat->n; j++) {
            if (j == skip_col) continue;
            minor->data[idx(r, c, minor)] = mat->data[idx(i, j, mat)];
            c++;
        }
        r++;
    }
    return minor;
}

double determinant(const zephyr_matrix *mat) {
    if (!mat || mat->m != mat->n) return NAN;
    const size_t n = mat->n;
    if (n == 1) return mat->data[0];
    if (n == 2) {
        return mat->data[0] * mat->data[3] - mat->data[1] * mat->data[2];
    }
    double det = 0.0;
    for (size_t col = 0; col < n; col++) {
        const double sign = (col % 2 == 0) ? 1.0 : -1.0;
        zephyr_matrix *minor = minor_matrix(mat, 0, col);
        if (!minor) return NAN;
        det += sign * mat->data[idx(0, col, mat)] * determinant(minor);
        destroy_matrix(minor);
    }
    return det;
}

zephyr_matrix *cofactor_matrix(const zephyr_matrix *mat) {
    if (!mat || mat->m != mat->n) return NULL;
    const size_t n = mat->n;
    zephyr_matrix *cof = create_matrix(n, n);
    if (!cof) return NULL;
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            zephyr_matrix *minor = minor_matrix(mat, i, j);
            if (!minor) {
                destroy_matrix(cof);
                return NULL;
            }
            else {
                double sign = ((i + j) % 2 == 0) ? 1.0 : -1.0;
                cof->data[idx(i, j, cof)] = sign * determinant(minor);
                destroy_matrix(minor);
            }
        }
    }
    return cof;
}

zephyr_matrix *adjoint_matrix(const zephyr_matrix *mat) {
    if (!mat || mat->m != mat->n) return NULL;
    zephyr_matrix *cof = cofactor_matrix(mat);
    if (!cof) return NULL;
    zephyr_matrix *adjoint = transpose(cof);
    destroy_matrix(cof);
    return adjoint;
}

zephyr_matrix *matrix_inverse(const zephyr_matrix *mat) {
    if (!mat || mat->m != mat->n) return NULL;
    const double det = determinant(mat);
    if (fabs(det) < 1e-10 || isnan(det)) return NULL;
    zephyr_matrix *adj = adjoint_matrix(mat);
    if (!adj) return NULL;
    zephyr_matrix *inv = scalar_multiply_matrix(1.0 / det, adj);
    destroy_matrix(adj);
    return inv;
}

zephyr_vector *simultaneous_eq(const zephyr_matrix *coefficients, const zephyr_matrix *solution) {
    if (!coefficients || !solution) return NULL;

    size_t n = coefficients->m;
    size_t m = coefficients->n;

    if (solution->m != n || solution->n != 1) return NULL;

    zephyr_matrix *aug = create_matrix(n, m + 1);
    if (!aug) return NULL;

    // Form the augmented matrix [A | b]
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < m; ++j) {
            aug->data[i * (m + 1) + j] = coefficients->data[i * m + j];
        }
        aug->data[i * (m + 1) + m] = solution->data[i];
    }

    // Gaussian elimination
    for (size_t k = 0; k < n; ++k) {
        // Find pivot
        size_t pivot_row = k;
        double pivot_val = fabs(aug->data[k * (m + 1) + k]);
        for (size_t i = k + 1; i < n; ++i) {
            double val = fabs(aug->data[i * (m + 1) + k]);
            if (val > pivot_val) {
                pivot_val = val;
                pivot_row = i;
            }
        }

        if (pivot_val < 1e-12) continue;  // Skip zero pivot (singular)

        // Swap rows if needed
        if (pivot_row != k) {
            for (size_t j = 0; j < m + 1; ++j) {
                double tmp = aug->data[k * (m + 1) + j];
                aug->data[k * (m + 1) + j] = aug->data[pivot_row * (m + 1) + j];
                aug->data[pivot_row * (m + 1) + j] = tmp;
            }
        }

        // Eliminate below
        for (size_t i = k + 1; i < n; ++i) {
            double factor = aug->data[i * (m + 1) + k] / aug->data[k * (m + 1) + k];
            for (size_t j = k; j < m + 1; ++j) {
                aug->data[i * (m + 1) + j] -= factor * aug->data[k * (m + 1) + j];
            }
        }
    }

    // Back-substitution
    zephyr_vector *x = create_vector(m);
    if (!x) {
        destroy_matrix(aug);
        return NULL;
    }

    for (ssize_t i = m - 1; i >= 0; --i) {
        double sum = aug->data[i * (m + 1) + m];
        for (size_t j = i + 1; j < m; ++j) {
            sum -= aug->data[i * (m + 1) + j] * x->data[j];
        }

        double coeff = aug->data[i * (m + 1) + i];
        x->data[i] = (fabs(coeff) < 1e-12) ? 1.0 : sum / coeff;  // Handle free vars
    }

    destroy_matrix(aug);
    return x;
}
bool power_iteration(const zephyr_matrix *A, zephyr_vector *eigenvector,
                     double *eigenvalue, int max_iter, double tol) {
    if (!A || A->m != A->n || !eigenvector || !eigenvalue || A->m != eigenvector->size)
        return false;

    zephyr_vector *b = create_vector(eigenvector->size);
    if (!b) return false;

    for (int i = 0; i < (int)b->size; i++) {
        b->data[i] = 1.0;  // Initial guess
    }
    normalize_vector(b);

    zephyr_matrix *b_mat = matrix_from_vector(b);
    for (int iter = 0; iter < max_iter; iter++) {
        zephyr_matrix *Ab = matrix_multiplication(A, b_mat);
        zephyr_vector *b_new = vector_from_matrix(Ab);
        destroy_matrix(Ab);

        if (!normalize_vector(b_new)) {
            destroy_vector(b);
            destroy_vector(b_new);
            destroy_matrix(b_mat);
            return false;
        }

        // Convergence check
        zephyr_vector *diff = create_vector(b->size);
        double diff_norm = 0.0;
        for (size_t i = 0; i < b->size; i++) {
            diff->data[i] = b_new->data[i] - b->data[i];
            diff_norm += diff->data[i] * diff->data[i];
        }
        diff_norm = sqrt(diff_norm);
        destroy_vector(diff);

        destroy_vector(b);
        b = b_new;
        destroy_matrix(b_mat);
        b_mat = matrix_from_vector(b);

        if (diff_norm < tol) break;
    }

    // Approximate eigenvalue using Rayleigh quotient
    zephyr_matrix *Av = matrix_multiplication(A, b_mat);
    zephyr_vector *Av_vec = vector_from_matrix(Av);
    *eigenvalue = dot_product(b, Av_vec);

    // Output vector
    memcpy(eigenvector->data, b->data, sizeof(double) * b->size);

    destroy_vector(Av_vec);
    destroy_matrix(Av);
    destroy_vector(b);
    destroy_matrix(b_mat);
    return true;
}
#include <stdbool.h>
#define MAX_ITER 1000
#define TOL 1e-9

static void householder_reflection(const zephyr_matrix *A, zephyr_matrix **Q, zephyr_matrix **R) {
    size_t n = A->m;
    *Q = create_identity_matrix(n);
    *R = copy_matrix(A);
    if (!*Q || !*R) return;

    for (size_t k = 0; k < n - 1; ++k) {
        double norm_x = 0.0;
        for (size_t i = k; i < n; ++i) {
            norm_x += (*R)->data[i * n + k] * (*R)->data[i * n + k];
        }
        norm_x = sqrt(norm_x);

        double sign = ((*R)->data[k * n + k] >= 0) ? 1.0 : -1.0;
        double u1 = (*R)->data[k * n + k] + sign * norm_x;

        zephyr_vector *v = create_vector(n);
        if (!v) return;
        for (size_t i = 0; i < n; ++i)
            v->data[i] = (i < k) ? 0.0 : (*R)->data[i * n + k];
        v->data[k] = u1;
        normalize_vector(v);

        zephyr_matrix *v_mat = matrix_from_vector(v);
        zephyr_matrix *v_mat_T = transpose(v_mat);
        zephyr_matrix *vvt = matrix_multiplication(v_mat, v_mat_T);
        zephyr_matrix *Pvvt = scalar_multiply_matrix(2.0, vvt);
        zephyr_matrix *X = create_identity_matrix(n);
        zephyr_matrix *P = subtract_matrices(X, Pvvt);

        zephyr_matrix *R_new = matrix_multiplication(P, *R);
        zephyr_matrix *Q_new = matrix_multiplication(*Q, transpose(P));

        destroy_matrix(*R);
        destroy_matrix(*Q);
        *R = R_new;
        *Q = Q_new;

        // Clean up
        destroy_vector(v);
        destroy_matrix(vvt);
        destroy_matrix(Pvvt);
        destroy_matrix(X);
        destroy_matrix(P);
        destroy_matrix(v_mat);
        destroy_matrix(v_mat_T);
    }
}
bool qr_algorithm(const zephyr_matrix *input, zephyr_vector **eigenvalues) {
    if (!input || input->m != input->n) return false;

    size_t n = input->m;
    zephyr_matrix *A = copy_matrix(input);
    if (!A) return false;

    for (int iter = 0; iter < MAX_ITER; ++iter) {
        zephyr_matrix *Q = NULL, *R = NULL;
        householder_reflection(A, &Q, &R);
        zephyr_matrix *A_next = matrix_multiplication(R, Q);

        destroy_matrix(A);
        A = A_next;

        destroy_matrix(Q);
        destroy_matrix(R);
    }

    *eigenvalues = create_vector(n);
    if (!*eigenvalues) {
        destroy_matrix(A);
        return false;
    }

    for (size_t i = 0; i < n; ++i) {
        (*eigenvalues)->data[i] = A->data[i * n + i];
    }

    destroy_matrix(A);
    return true;
}
#define SVD_EPSILON 1e-9
zephyr_matrix *svd_null_space(const zephyr_matrix *A) {
    if (!A || A->m < A->n) return NULL;  // Full SVD for m >= n

    // Compute A^T * A
    zephyr_matrix *At = transpose(A);
    if (!At) return NULL;
    zephyr_matrix *AtA = matrix_multiplication(At, A);
    destroy_matrix(At);
    if (!AtA) return NULL;

    // Eigen decomposition of AtA (symmetric positive semi-definite)
    zephyr_vector *singular_values_squared = NULL;
    if (!qr_algorithm(AtA, &singular_values_squared)) {
        destroy_matrix(AtA);
        return NULL;
    }

    // Sort eigenvalues and track indices
    size_t n = AtA->m;
    double *sigma = malloc(sizeof(double) * n);
    size_t *indices = malloc(sizeof(size_t) * n);
    for (size_t i = 0; i < n; ++i) {
        sigma[i] = sqrt(fmax(0.0, singular_values_squared->data[i]));
        indices[i] = i;
    }

    // Sort by ascending sigma
    for (size_t i = 0; i < n - 1; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            if (sigma[i] > sigma[j]) {
                double tmp = sigma[i]; sigma[i] = sigma[j]; sigma[j] = tmp;
                size_t ti = indices[i]; indices[i] = indices[j]; indices[j] = ti;
            }
        }
    }

    // Create null space matrix
    zephyr_matrix *null_space = create_matrix(n, 0);  // Start with 0 columns

    for (size_t k = 0; k < n; ++k) {
        if (sigma[k] > SVD_EPSILON) continue;  // Skip non-zero singular values

        double lambda = singular_values_squared->data[indices[k]];
        zephyr_matrix *shifted = create_matrix(n, n);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                double aij = AtA->data[i * n + j];
                if (i == j) aij -= lambda;
                shifted->data[i * n + j] = aij;
            }
        }
        zephyr_matrix *b_zero = create_matrix(n, 1);
        zephyr_vector *v = simultaneous_eq(shifted, b_zero);
        if (!v) {
            destroy_matrix(shifted);
            destroy_matrix(b_zero);
            continue;
        }
        normalize_vector(v);

        // Append column to null_space matrix
        zephyr_matrix *new_null = create_matrix(n, null_space->n + 1);
        for (size_t col = 0; col < null_space->n; ++col) {
            for (size_t row = 0; row < n; ++row) {
                new_null->data[row * new_null->n + col] = null_space->data[row * null_space->n + col];
            }
        }
        for (size_t row = 0; row < n; ++row) {
            new_null->data[row * new_null->n + null_space->n] = v->data[row];
        }
        destroy_matrix(null_space);
        null_space = new_null;

        destroy_vector(v);
        destroy_matrix(shifted);
        destroy_matrix(b_zero);
    }

    free(sigma);
    free(indices);
    destroy_vector(singular_values_squared);
    destroy_matrix(AtA);

    return null_space;  // Matrix with each column a null space basis vector
}
zephyr_matrix *eigenvectors_from_eigenvalues(const zephyr_matrix *A, const zephyr_vector *eigenvalues) {
    if (!A || !eigenvalues || A->m != A->n || A->m != eigenvalues->size) return NULL;

    const size_t n = A->m;
    zephyr_matrix *eigenvectors = create_matrix(n, n);
    if (!eigenvectors) return NULL;

    for (size_t k = 0; k < n; ++k) {
        double lambda = eigenvalues->data[k];
        zephyr_matrix *shifted = create_matrix(n, n);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                double aij = A->data[i * n + j];
                if (i == j) aij -= lambda;
                shifted->data[i * n + j] = aij;
            }
        }

        zephyr_matrix *b_zero = create_matrix(n, 1); // Ax = 0
        zephyr_vector *x = simultaneous_eq(shifted, b_zero);
        if (!x) {
            // fallback to SVD-like null-space approach
            zephyr_matrix *ns = svd_null_space(shifted);
            if (!ns || ns->n == 0) {
                destroy_matrix(ns);
                destroy_matrix(eigenvectors);
                return NULL;
            }
            x = create_vector(ns->m);
            for (size_t i = 0; i < ns->m; ++i) {
                x->data[i] = ns->data[i * ns->n + 0]; // pick first column
            }
            destroy_matrix(ns);
        }

        destroy_matrix(b_zero);
        destroy_matrix(shifted);

        if (!x) {
            destroy_matrix(eigenvectors);
            return NULL;
        }

        normalize_vector(x);
        for (size_t i = 0; i < n; ++i) {
            eigenvectors->data[i * n + k] = x->data[i];
        }
        destroy_vector(x);
    }

    return eigenvectors;
}
bool eigen_solve(const zephyr_matrix *A, zephyr_eigen_result *result) {
    if (!A || A->m != A->n || !result) return false;

    size_t n = A->m;
    zephyr_vector *eigenvalues = NULL;
    if (!qr_algorithm(A, &eigenvalues)) {
        return false;
    }
    zephyr_matrix *eigenvectors = eigenvectors_from_eigenvalues(A, eigenvalues);
    if (!eigenvectors) {
        destroy_vector(eigenvalues);
        return false;
    }

    result->values = eigenvalues;
    result->vectors = eigenvectors;
    return true;
}
void destroy_eigen_result(zephyr_eigen_result *result) {
    if (result) {
        destroy_vector(result->values);
        destroy_matrix(result->vectors);
    }
}


zephyr_matrix *rotation_matrix_2d(const zephyr_vector *to_rotate, double theta) {
    if (!to_rotate || to_rotate->size != 2) return NULL;

    zephyr_matrix *rotation_matrix_value = create_matrix(2, 2);
    if (!rotation_matrix_value) return NULL;

    double cos_theta = cos(theta);
    double sin_theta = sin(theta);

    rotation_matrix_value->data[0] = cos_theta;
    rotation_matrix_value->data[1] = -sin_theta;
    rotation_matrix_value->data[2] = sin_theta;
    rotation_matrix_value->data[3] = cos_theta;

    zephyr_matrix *result = matrix_col_vector_mult(to_rotate, rotation_matrix_value);
    destroy_matrix(rotation_matrix_value);
    return result;
}




zephyr_matrix *rotate_about_axis(const zephyr_vector *v, const zephyr_vector *axis, double theta) {
    if (!v || !axis || v->size != 3 || axis->size != 3) return NULL;

    double ux = axis->data[0], uy = axis->data[1], uz = axis->data[2];
    double c = cos(theta), s = sin(theta), t = 1 - c;

    // Normalize axis
    double norm = sqrt(ux * ux + uy * uy + uz * uz);
    if (norm < 1e-9) return NULL;
    ux /= norm; uy /= norm; uz /= norm;

    zephyr_matrix *R = create_matrix(3, 3);
    if (!R) return NULL;

    R->data[0] = t * ux * ux + c;
    R->data[1] = t * ux * uy - s * uz;
    R->data[2] = t * ux * uz + s * uy;

    R->data[3] = t * ux * uy + s * uz;
    R->data[4] = t * uy * uy + c;
    R->data[5] = t * uy * uz - s * ux;

    R->data[6] = t * ux * uz - s * uy;
    R->data[7] = t * uy * uz + s * ux;
    R->data[8] = t * uz * uz + c;

    zephyr_matrix *result = matrix_col_vector_mult(v, R);
    destroy_matrix(R);
    return result;
}

zephyr_vector *translate_2d_vector(const zephyr_vector *point, const zephyr_vector *translation) {
    if (!point || !translation || point->size != 2 || translation->size != 2)
        return NULL;

    zephyr_matrix *T = create_matrix(3, 3);
    if (!T) return NULL;

    // Construct 2D translation matrix
    T->data[0] = 1; T->data[1] = 0; T->data[2] = translation->data[0];  // tx
    T->data[3] = 0; T->data[4] = 1; T->data[5] = translation->data[1];  // ty
    T->data[6] = 0; T->data[7] = 0; T->data[8] = 1;

    // Convert point to 3x1 homogeneous vector
    zephyr_matrix *p = create_matrix(3, 1);
    if (!p) {
        destroy_matrix(T);
        return NULL;
    }

    p->data[0] = point->data[0];
    p->data[1] = point->data[1];
    p->data[2] = 1.0;

    zephyr_matrix *translated = matrix_multiplication(T, p);
    destroy_matrix(T);
    destroy_matrix(p);
    zephyr_vector * to_vector = vector_from_matrix(translated);
    destroy_matrix(translated);
    return to_vector;
}

zephyr_vector *translate_3d_vector(const zephyr_vector *point, const zephyr_vector *translation) {
    if (!point || !translation || point->size != 3 || translation->size != 3)
        return NULL;

    zephyr_matrix *T = create_matrix(4, 4);
    if (!T) return NULL;

    T->data[0] = 1; T->data[1] = 0; T->data[2] = 0; T->data[3] = translation->data[0];
    T->data[4] = 0; T->data[5] = 1; T->data[6] = 0; T->data[7] = translation->data[1];
    T->data[8] = 0; T->data[9] = 0; T->data[10] = 1; T->data[11] = translation->data[2];
    T->data[12] = 0; T->data[13] = 0; T->data[14] = 0; T->data[15] = 1;

    // Convert point to 4x1 homogeneous
    zephyr_matrix *p = create_matrix(4, 1);
    if (!p) {
        destroy_matrix(T);
        return NULL;
    }

    p->data[0] = point->data[0];
    p->data[1] = point->data[1];
    p->data[2] = point->data[2];
    p->data[3] = 1.0;

    zephyr_matrix *translated = matrix_multiplication(T, p);
    destroy_matrix(T);
    destroy_matrix(p);
    zephyr_vector * to_vector = vector_from_matrix(translated);
    destroy_matrix(translated);
    return to_vector;
}

zephyr_matrix *rotation_matrix_3d_xyz(double theta_x, double theta_y, double theta_z) {
    zephyr_matrix *Rx = create_matrix(3, 3);
    zephyr_matrix *Ry = create_matrix(3, 3);
    zephyr_matrix *Rz = create_matrix(3, 3);

    if (!Rx || !Ry || !Rz) return NULL;

    // Rotation about X-axis
    Rx->data[0] = 1;    Rx->data[1] = 0;            Rx->data[2] = 0;
    Rx->data[3] = 0;    Rx->data[4] = cos(theta_x); Rx->data[5] = -sin(theta_x);
    Rx->data[6] = 0;    Rx->data[7] = sin(theta_x); Rx->data[8] = cos(theta_x);

    // Rotation about Y-axis
    Ry->data[0] = cos(theta_y);  Ry->data[1] = 0; Ry->data[2] = sin(theta_y);
    Ry->data[3] = 0;             Ry->data[4] = 1; Ry->data[5] = 0;
    Ry->data[6] = -sin(theta_y); Ry->data[7] = 0; Ry->data[8] = cos(theta_y);

    // Rotation about Z-axis
    Rz->data[0] = cos(theta_z); Rz->data[1] = -sin(theta_z); Rz->data[2] = 0;
    Rz->data[3] = sin(theta_z); Rz->data[4] = cos(theta_z);  Rz->data[5] = 0;
    Rz->data[6] = 0;            Rz->data[7] = 0;             Rz->data[8] = 1;

    // Combine R = Rz * Ry * Rx
    zephyr_matrix *temp = matrix_multiplication(Rz, Ry);
    zephyr_matrix *R = matrix_multiplication(temp, Rx);

    destroy_matrix(Rx);
    destroy_matrix(Ry);
    destroy_matrix(Rz);
    destroy_matrix(temp);

    return R;
}

zephyr_vector *rotate_3d_per_axis(const zephyr_vector *v, double theta, rotation_axis axis) {
    if (!v || v->size != 3) return NULL;

    double tx = 0, ty = 0, tz = 0;

    switch (axis) {
        case ROTATE_X: tx = theta; break;
        case ROTATE_Y: ty = theta; break;
        case ROTATE_Z: tz = theta; break;
        default: return NULL;
    }

    zephyr_matrix *R = rotation_matrix_3d_xyz(tx, ty, tz);
    if (!R) return NULL;

    zephyr_matrix *result = matrix_col_vector_mult(v, R);
    destroy_matrix(R);
    zephyr_vector * vector = vector_from_matrix(result);
    destroy_matrix(result);
    return vector;
}

zephyr_vector *rotate_3d_vector(const zephyr_vector *v,
                                double theta_x, double theta_y, double theta_z) {
    if (!v || v->size != 3) return NULL;

    zephyr_matrix *R = rotation_matrix_3d_xyz(theta_x, theta_y, theta_z);
    if (!R) return NULL;

    zephyr_matrix *rotated_mat = matrix_col_vector_mult(v, R);
    destroy_matrix(R);
    if (!rotated_mat) return NULL;

    zephyr_vector *rotated_vector = vector_from_matrix(rotated_mat);
    destroy_matrix(rotated_mat);
    return rotated_vector;
}

zephyr_vector *transform_3d_with_matrix(const zephyr_vector *v, const zephyr_matrix *T) {
    if (!v || v->size != 3 || !T || T->m != 4 || T->n != 4)
        return NULL;

    // Homogenize input vector: [x, y, z, 1]^T
    zephyr_matrix *vec4 = create_matrix(4, 1);
    if (!vec4) return NULL;

    vec4->data[0] = v->data[0];
    vec4->data[1] = v->data[1];
    vec4->data[2] = v->data[2];
    vec4->data[3] = 1.0;

    // Multiply: T * vec4
    zephyr_matrix *result = matrix_multiplication(T, vec4);
    destroy_matrix(vec4);
    if (!result) return NULL;

    // Extract 3D result
    zephyr_vector *final = create_vector(3);
    for (size_t i = 0; i < 3; ++i) {
        final->data[i] = result->data[i];
    }

    destroy_matrix(result);
    return final;
}


zephyr_vector *transform_3d(const zephyr_vector *v,
                            double theta_x, double theta_y, double theta_z,
                            const zephyr_vector *translation) {
    if (!v || v->size != 3 || !translation || translation->size != 3)
        return NULL;

    zephyr_matrix *R3 = rotation_matrix_3d_xyz(theta_x, theta_y, theta_z);
    if (!R3) return NULL;

    zephyr_matrix *T = create_matrix(4, 4);
    if (!T) {
        destroy_matrix(R3);
        return NULL;
    }

    // Top-left 3x3 rotation part
    for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
            T->data[i * 4 + j] = R3->data[i * 3 + j];

    // Translation
    T->data[0 * 4 + 3] = translation->data[0];
    T->data[1 * 4 + 3] = translation->data[1];
    T->data[2 * 4 + 3] = translation->data[2];

    // Last row [0 0 0 1]
    T->data[3 * 4 + 3] = 1.0;

    zephyr_matrix *vec4 = create_matrix(4, 1);
    vec4->data[0] = v->data[0];
    vec4->data[1] = v->data[1];
    vec4->data[2] = v->data[2];
    vec4->data[3] = 1.0;

    zephyr_matrix *result = matrix_multiplication(T, vec4);

    // Extract only the x, y, z part as a vector
    zephyr_vector *final = create_vector(3);
    for (size_t i = 0; i < 3; ++i) {
        final->data[i] = result->data[i];
    }

    destroy_matrix(R3);
    destroy_matrix(T);
    destroy_matrix(vec4);
    destroy_matrix(result);

    return final;  // ✅ Now just a 3D vector
}

zephyr_vector * inverse_transformation_3d(const zephyr_vector *v,
                            double theta_x, double theta_y, double theta_z,
                            const zephyr_vector *translation) {
    if (!v || v->size != 3 || !translation || translation->size != 3)
        return NULL;

    zephyr_matrix *R3 = rotation_matrix_3d_xyz(theta_x, theta_y, theta_z);
    if (!R3) return NULL;

    zephyr_matrix *R3_transpose = transpose(R3);
    destroy_matrix(R3);
    if (!R3_transpose) return NULL;

    zephyr_matrix *T = create_matrix(4, 4);
    if (!T) {
        destroy_matrix(R3_transpose);
        return NULL;
    }

    // Fill top-left 3x3 with R^T
    for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
            T->data[i * 4 + j] = R3_transpose->data[i * 3 + j];

    // Multiply R^T with translation vector and negate result
    zephyr_matrix *t_mat = create_matrix(3, 1);
    for (size_t i = 0; i < 3; ++i)
        t_mat->data[i] = translation->data[i];

    zephyr_matrix *Rt_t = matrix_multiplication(R3_transpose, t_mat);
    destroy_matrix(t_mat);
    if (!Rt_t) {
        destroy_matrix(R3_transpose);
        destroy_matrix(T);
        return NULL;
    }

    for (size_t i = 0; i < 3; ++i)
        T->data[i * 4 + 3] = -Rt_t->data[i];

    // Last row [0 0 0 1]
    T->data[3 * 4 + 0] = 0.0;
    T->data[3 * 4 + 1] = 0.0;
    T->data[3 * 4 + 2] = 0.0;
    T->data[3 * 4 + 3] = 1.0;

    destroy_matrix(R3_transpose);
    destroy_matrix(Rt_t);

    // Convert v to homogeneous
    zephyr_matrix *v4 = create_matrix(4, 1);
    v4->data[0] = v->data[0];
    v4->data[1] = v->data[1];
    v4->data[2] = v->data[2];
    v4->data[3] = 1.0;

    zephyr_matrix *result = matrix_multiplication(T, v4);
    destroy_matrix(T);
    destroy_matrix(v4);

    // Return only xyz part as vector
    zephyr_vector *result_vec = create_vector(3);
    if (!result_vec) {
        destroy_matrix(result);
        return NULL;
    }
    for (size_t i = 0; i < 3; ++i)
        result_vec->data[i] = result->data[i];

    destroy_matrix(result);
    return result_vec;
}

zephyr_matrix *compose_transformations(const zephyr_matrix *R_AB, const zephyr_vector *P_AB,
                                       const zephyr_matrix *R_BC, const zephyr_vector *P_BC) {
    if (!R_AB || !P_AB || !R_BC || !P_BC)
        return NULL;

    // Validate sizes
    if (R_AB->m != 3 || R_AB->n != 3 || R_BC->m != 3 || R_BC->n != 3 ||
        P_AB->size != 3 || P_BC->size != 3)
        return NULL;

    // Step 1: R_AC = R_AB * R_BC
    zephyr_matrix *R_AC = matrix_multiplication(R_AB, R_BC);
    if (!R_AC) return NULL;

    // Step 2: P_AC = R_AB * P_BC + P_AB
    zephyr_matrix *P_BC_mat = create_matrix(3, 1);
    memcpy(P_BC_mat->data, P_BC->data, sizeof(double) * 3);

    zephyr_matrix *RAB_PBC = matrix_multiplication(R_AB, P_BC_mat);
    destroy_matrix(P_BC_mat);
    if (!RAB_PBC) {
        destroy_matrix(R_AC);
        return NULL;
    }

    zephyr_vector *P_AC = create_vector(3);
    for (size_t i = 0; i < 3; ++i)
        P_AC->data[i] = RAB_PBC->data[i] + P_AB->data[i];

    destroy_matrix(RAB_PBC);

    // Step 3: Construct full 4x4 homogeneous transform matrix
    zephyr_matrix *T_AC = create_matrix(4, 4);
    if (!T_AC) {
        destroy_matrix(R_AC);
        destroy_vector(P_AC);
        return NULL;
    }

    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            T_AC->data[i * 4 + j] = R_AC->data[i * 3 + j];
        }
        T_AC->data[i * 4 + 3] = P_AC->data[i];
    }
    T_AC->data[3 * 4 + 0] = 0.0;
    T_AC->data[3 * 4 + 1] = 0.0;
    T_AC->data[3 * 4 + 2] = 0.0;
    T_AC->data[3 * 4 + 3] = 1.0;

    destroy_matrix(R_AC);
    destroy_vector(P_AC);
    return T_AC;
}

zephyr_matrix *rigid_transform_about_axis(const zephyr_vector *axis_point,
                                          const zephyr_vector *axis_dir,
                                          double theta_rad) {
    if (!axis_point || !axis_dir || axis_point->size != 3 || axis_dir->size != 3)
        return NULL;
    // Normalize axis
    double ux = axis_dir->data[0], uy = axis_dir->data[1], uz = axis_dir->data[2];
    double norm = sqrt(ux * ux + uy * uy + uz * uz);
    if (norm < 1e-9) return NULL;
    ux /= norm; uy /= norm; uz /= norm;

    double c = cos(theta_rad);
    double s = sin(theta_rad);
    double t = 1 - c;

    zephyr_matrix *R = create_matrix(3, 3);
    R->data[0] = t * ux * ux + c;
    R->data[1] = t * ux * uy - s * uz;
    R->data[2] = t * ux * uz + s * uy;

    R->data[3] = t * ux * uy + s * uz;
    R->data[4] = t * uy * uy + c;
    R->data[5] = t * uy * uz - s * ux;

    R->data[6] = t * ux * uz - s * uy;
    R->data[7] = t * uy * uz + s * ux;
    R->data[8] = t * uz * uz + c;

    // Build 4x4 homogeneous transform
    zephyr_matrix *T = create_matrix(4, 4);
    if (!T) { destroy_matrix(R); return NULL; }

    // Fill rotation part
    for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
            T->data[i * 4 + j] = R->data[i * 3 + j];

    // Translation = -R * P + P
    for (size_t i = 0; i < 3; ++i) {
        double t_val = 0.0;
        for (size_t j = 0; j < 3; ++j)
            t_val += -R->data[i * 3 + j] * axis_point->data[j];
        t_val += axis_point->data[i];
        T->data[i * 4 + 3] = t_val;
    }

    // Final row
    T->data[3 * 4 + 0] = 0;
    T->data[3 * 4 + 1] = 0;
    T->data[3 * 4 + 2] = 0;
    T->data[3 * 4 + 3] = 1;

    destroy_matrix(R);
    return T;
}
