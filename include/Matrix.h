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

/**
 * @brief Creates a new matrix with specified dimensions and allocates memory.
 *
 * This function initializes a new `zephyr_matrix` structure with `m` rows and `n` columns.
 * It allocates memory for both the matrix structure and its internal data array.
 * All elements in the data array are initialized to 0 using `calloc`.
 *
 * @param m The number of rows in the matrix.
 * @param n The number of columns in the matrix.
 * @return A pointer to the newly created `zephyr_matrix` structure,
 *         or NULL if memory allocation fails.
 */
zephyr_matrix * create_matrix(const size_t m, const size_t n);
/**
 * @brief Creates a new identity matrix of size m x m.
 *
 * This function creates a square matrix with ones on the main diagonal and zeros elsewhere.
 *
 * @param m The number of rows and columns in the identity matrix.
 * @return A pointer to the newly created identity `zephyr_matrix` structure,
 *         or NULL if memory allocation fails.
 */
zephyr_matrix * create_identity_matrix(const size_t m);
/**
 * @brief Creates a deep copy of an existing matrix.
 *
 * This function allocates new memory for a new `zephyr_matrix` and copies
 * all elements from the source matrix to the new matrix.
 *
 * @param matrix A pointer to the source `zephyr_matrix` to be copied.
 * @return A pointer to the newly created `zephyr_matrix` copy,
 *         or NULL if the source matrix is NULL or memory allocation fails.
 */
zephyr_matrix * copy_matrix(const zephyr_matrix * matrix);
/**
 * @brief Sets the value of a specific element in a matrix.
 *
 * This function updates the element at the given row and column with the specified value.
 *
 * @param row The row index of the element to set (0-based).
 * @param col The column index of the element to set (0-based).
 * @param value The double value to set at the specified position.
 * @param matrix A pointer to the `zephyr_matrix` to modify.
 * @return A pointer to the modified `zephyr_matrix`,
 *         or NULL if the matrix is NULL or the row/column indices are out of bounds.
 */
zephyr_matrix * matrix_set(const size_t row, const size_t col, const double value, zephyr_matrix * matrix);
/**
 * @brief Retrieves the value of a specific element from a matrix.
 *
 * This function returns the value at the given row and column from the specified matrix.
 *
 * @param row The row index of the element to retrieve (0-based).
 * @param col The column index of the element to retrieve (0-based).
 * @param matrix A pointer to the `zephyr_matrix` to query.
 * @return The double value at the specified position,
 *         or NAN if the matrix is NULL or the row/column indices are out of bounds.
 */
double matrix_get(const size_t row, const size_t col, const zephyr_matrix * matrix);
/**
 * @brief Adds two matrices element-wise.
 *
 * This function creates a new matrix and stores the sum of the corresponding
 * elements of the two input matrices. Both matrices must have the same dimensions.
 *
 * @param matricxA A pointer to the first `zephyr_matrix`.
 * @param matricxB A pointer to the second `zephyr_matrix`.
 * @return A pointer to the newly created `zephyr_matrix` containing the sum,
 *         or NULL if the dimensions do not match or memory allocation fails.
 */
zephyr_matrix * add_matrices(const zephyr_matrix * matricxA, const zephyr_matrix * matricxB);
/**
 * @brief Subtracts one matrix from another element-wise.
 *
 * This function creates a new matrix and stores the result of subtracting
 * the elements of the second matrix from the corresponding elements of the first matrix.
 * Both matrices must have the same dimensions.
 *
 * @param matricxA A pointer to the first `zephyr_matrix`.
 * @param matricxB A pointer to the second `zephyr_matrix` to be subtracted.
 * @return A pointer to the newly created `zephyr_matrix` containing the difference,
 *         or NULL if the dimensions do not match or memory allocation fails.
 */
zephyr_matrix * subtract_matrices(const zephyr_matrix * matricxA, const zephyr_matrix * matricxB);
/**
 * @brief Multiplies a matrix by a scalar value.
 *
 * This function creates a new matrix where each element is the corresponding
 * element of the input matrix multiplied by the given scalar.
 *
 * @param scale The scalar value to multiply the matrix by.
 * @param matrix A pointer to the `zephyr_matrix` to be scaled.
 * @return A pointer to the newly created `zephyr_matrix` containing the scaled result,
 *         or NULL if memory allocation fails.
 */
zephyr_matrix * scalar_multiply_matrix(const double scale, const zephyr_matrix * matrix);
/**
 * @brief Multiplies two matrices.
 *
 * This function performs matrix multiplication of two input matrices, `a` and `b`.
 * The number of columns in matrix `a` must be equal to the number of rows in matrix `b`.
 *
 * @param a A pointer to the first `zephyr_matrix`.
 * @param b A pointer to the second `zephyr_matrix`.
 * @return A pointer to the newly created `zephyr_matrix` containing the product,
 *         or NULL if the dimensions are incompatible or memory allocation fails.
 */
zephyr_matrix * matrix_multiplication(const zephyr_matrix * a, const zephyr_matrix * b);
/**
 * @brief Computes the transpose of a matrix.
 *
 * This function creates a new matrix that is the transpose of the input matrix.
 * The rows of the original matrix become the columns of the new matrix, and vice-versa.
 *
 * @param matrix A pointer to the `zephyr_matrix` to be transposed.
 * @return A pointer to the newly created transposed `zephyr_matrix`,
 *         or NULL if the input matrix is NULL or memory allocation fails.
 */
zephyr_matrix * transpose(const zephyr_matrix * matrix);
/**
 * @brief Destroys a matrix and frees its allocated memory.
 *
 * This function deallocates the memory associated with the matrix data and the matrix structure itself.
 * It is safe to call this function with a NULL matrix pointer.
 *
 * @param matrix A pointer to the `zephyr_matrix` to be destroyed.
 */
void destroy_matrix(zephyr_matrix * matrix);
/**
 * @brief Multiplies a matrix by a column vector.
 *
 * This function treats the input `zephyr_vector` as a column matrix (n x 1)
 * and performs matrix multiplication with the given `zephyr_matrix`.
 * The number of columns in the matrix must match the size of the vector.
 *
 * @param zephyr_vector A pointer to the `zephyr_vector` (column vector) to multiply.
 * @param matrix A pointer to the `zephyr_matrix`.
 * @return A pointer to the newly created `zephyr_matrix` representing the product (a column vector),
 *         or NULL if inputs are invalid or dimensions are incompatible or memory allocation fails.
 */
zephyr_matrix * matrix_col_vector_mult(const zephyr_vector * zephyr_vector, const zephyr_matrix * matrix);
/**
 * @brief Multiplies a row vector by a matrix.
 *
 * This function treats the input `zephyr_vector` as a row matrix (1 x n)
 * and performs matrix multiplication with the given `zephyr_matrix`.
 * The size of the vector must match the number of rows in the matrix.
 *
 * @param row_vector A pointer to the `zephyr_vector` (row vector) to multiply.
 * @param matrix A pointer to the `zephyr_matrix`.
 * @return A pointer to the newly created `zephyr_matrix` representing the product (a row vector),
 *         or NULL if inputs are invalid or dimensions are incompatible or memory allocation fails.
 */
zephyr_matrix * matrix_row_vector_mult(const zephyr_vector * row_vector, const zephyr_matrix * matrix);
/**
 * @brief Prints the elements of a matrix to the console.
 *
 * This function iterates through the matrix and prints each element formatted to a specified number of decimal places.
 * If the matrix or its data is NULL, it prints "NULL matrix".
 *
 * @param matrix A pointer to the `zephyr_matrix` to be printed.
 * @param decimal_places The number of decimal places to display for each element.
 */
void print_matrix(const zephyr_matrix *matrix, const int decimal_places);
/**
 * @brief Retrieves a specific row from a matrix as a new matrix (row vector).
 *
 * This function creates a new 1 x n matrix containing the elements of the specified row
 * from the input matrix.
 *
 * @param row The row index to retrieve (0-based).
 * @param matrix A pointer to the `zephyr_matrix`.
 * @return A pointer to the newly created `zephyr_matrix` (row vector),
 *         or NULL if the matrix is NULL or the row index is out of bounds or memory allocation fails.
 */
zephyr_matrix * matrix_get_row(const size_t row, const zephyr_matrix * matrix);
/**
 * @brief Retrieves a specific column from a matrix as a new matrix (column vector).
 *
 * This function creates a new m x 1 matrix containing the elements of the specified column
 * from the input matrix.
 *
 * @param col The column index to retrieve (0-based).
 * @param matrix A pointer to the `zephyr_matrix`.
 * @return A pointer to the newly created `zephyr_matrix` (column vector),
 *         or NULL if the matrix is NULL or the column index is out of bounds or memory allocation fails.
 */
zephyr_matrix * matrix_get_col(const size_t col, const zephyr_matrix * matrix);
/**
 * @brief Extracts a sub-matrix (slice) from a given matrix.
 *
 * This function creates a new matrix containing a portion of the original matrix
 * defined by the start and end row and column indices.
 *
 * @param row_start The starting row index (0-based, inclusive).
 * @param row_end The ending row index (0-based, inclusive).
 * @param col_start The starting column index (0-based, inclusive).
 * @param col_end The ending column index (0-based, inclusive).
 * @param matrix A pointer to the source `zephyr_matrix`.
 * @return A pointer to the newly created sliced `zephyr_matrix`,
 *         or NULL if inputs are invalid, indices are out of bounds, or memory allocation fails.
 */
zephyr_matrix * slice(const size_t row_start, const size_t row_end, const size_t col_start, const size_t col_end, const zephyr_matrix * matrix);
/**
 * @brief Creates a new matrix from a 1D array of doubles.
 *
 * This function initializes a new `zephyr_matrix` with `m` rows and `n` columns,
 * populating its data from the provided 1D array. The array is assumed to be
 * row-major ordered.
 *
 * @param data A pointer to the 1D array of doubles containing the matrix elements.
 * @param m The number of rows for the new matrix.
 * @param n The number of columns for the new matrix.
 * @return A pointer to the newly created `zephyr_matrix` structure,
 *         or NULL if the data pointer is NULL, dimensions are zero, or memory allocation fails.
 */
zephyr_matrix * matrix_from_array(const double * data, const size_t m, const size_t n);
/**
 * @brief Creates a new matrix from a `zephyr_vector`.
 *
 * This function creates a new matrix with dimensions `vector->size` x 1,
 * populating its data from the provided `zephyr_vector`.
 *
 * @param vector A pointer to the `zephyr_vector` to convert.
 * @return A pointer to the newly created `zephyr_matrix` structure,
 *         or NULL if the vector is NULL or memory allocation fails.
 */
zephyr_matrix * matrix_from_vector(const zephyr_vector * vector);
/**
 * @brief Creates a new `zephyr_vector` from a matrix.
 *
 * This function converts a 1 x n or m x 1 matrix into a `zephyr_vector`.
 * It returns NULL if the matrix is not a row or column vector.
 *
 * @param matrix A pointer to the `zephyr_matrix` to convert.
 * @return A pointer to the newly created `zephyr_vector` structure,
 *         or NULL if the matrix is not a vector or memory allocation fails.
 */
zephyr_vector * vector_from_matrix(const zephyr_matrix * matrix);
/**
 * @brief Creates a new matrix from a 1D array of doubles with a specified stride.
 *
 * This function initializes a new `zephyr_matrix` with `m` rows and `n` columns,
 * populating its data from the provided 1D array. The `stride` parameter specifies
 * the number of elements to advance in the `data` array to get to the next row.
 * This is useful for extracting sub-matrices from larger arrays.
 *
 * @param data A pointer to the 1D array of doubles containing the matrix elements.
 * @param m The number of rows for the new matrix.
 * @param n The number of columns for the new matrix.
 * @param stride The number of elements in the `data` array between the start of consecutive rows.
 * @return A pointer to the newly created `zephyr_matrix` structure,
 *         or NULL if the data pointer is NULL, dimensions are zero, stride is less than n, or memory allocation fails.
 */
zephyr_matrix * matrix_from_array_stride(const double * data, const size_t m, const size_t n, const size_t stride);
/**
 * @brief Computes the minor matrix of a given matrix by skipping a specified row and column.
 *
 * This function creates a new matrix by removing the `skip_row` and `skip_col`
 * from the input matrix. This is typically used in determinant and inverse calculations.
 * The input matrix must be square.
 *
 * @param mat A pointer to the source `zephyr_matrix`.
 * @param skip_row The row index to skip (0-based).
 * @param skip_col The column index to skip (0-based).
 * @return A pointer to the newly created minor `zephyr_matrix`,
 *         or NULL if the input matrix is not square or memory allocation fails.
 */
/**
 * @brief Computes the minor matrix of a given matrix by skipping a specified row and column.
 *
 * This function creates a new matrix by removing the `skip_row` and `skip_col`
 * from the input matrix. This is typically used in determinant and inverse calculations.
 * The input matrix must be square.
 *
 * @param mat A pointer to the source `zephyr_matrix`.
 * @param skip_row The row index to skip (0-based).
 * @param skip_col The column index to skip (0-based).
 * @return A pointer to the newly created minor `zephyr_matrix`,
 *         or NULL if the input matrix is not square or memory allocation fails.
 */
zephyr_matrix *minor_matrix(const zephyr_matrix *mat, const size_t skip_row, const size_t skip_col);
/**
 * @brief Computes the determinant of a square matrix.
 *
 * This function calculates the determinant of the input square matrix using recursion
 * and cofactor expansion. It handles 1x1 and 2x2 matrices directly.
 *
 * @param matrix A pointer to the `zephyr_matrix` for which to compute the determinant.
 * @return The determinant of the matrix as a double,
 *         or NAN if the matrix is not square or NULL.
 */
double determinant(const zephyr_matrix * matrix);
/**
 * @brief Computes the cofactor matrix of a square matrix.
 *
 * This function calculates the cofactor of each element in the input square matrix
 * and returns a new matrix containing these cofactors. The cofactor of an element
 * is its minor multiplied by (-1)^(row+col).
 *
 * @param mat A pointer to the source `zephyr_matrix`.
 * @return A pointer to the newly created cofactor `zephyr_matrix`,
 *         or NULL if the input matrix is not square or memory allocation fails.
 */
/**
 * @brief Computes the cofactor matrix of a square matrix.
 *
 * This function calculates the cofactor of each element in the input square matrix
 * and returns a new matrix containing these cofactors. The cofactor of an element
 * is its minor multiplied by (-1)^(row+col).
 *
 * @param mat A pointer to the source `zephyr_matrix`.
 * @return A pointer to the newly created cofactor `zephyr_matrix`,
 *         or NULL if the input matrix is not square or memory allocation fails.
 */
zephyr_matrix * cofactor_matrix(const zephyr_matrix *mat);
/**
 * @brief Computes the adjoint matrix of a square matrix.
 *
 * The adjoint of a square matrix is the transpose of its cofactor matrix.
 *
 * @param mat A pointer to the source `zephyr_matrix`.
 * @return A pointer to the newly created adjoint `zephyr_matrix`,
 *         or NULL if the input matrix is not square or memory allocation fails.
 */
zephyr_matrix * adjoint_matrix(const zephyr_matrix *mat);
/**
 * @brief Computes the inverse of a square matrix.
 *
 * This function calculates the inverse of the input square matrix using the formula:
 * inverse(A) = (1/determinant(A)) * adjoint(A).
 * It returns NULL if the matrix is not square, its determinant is zero (or very close to zero),
 * or memory allocation fails.
 *
 * @param mat A pointer to the source `zephyr_matrix`.
 * @return A pointer to the newly created inverse `zephyr_matrix`,
 *         or NULL if the matrix is singular, not square, or memory allocation fails.
 */
zephyr_matrix * matrix_inverse(const zephyr_matrix *mat);
/**
 * @brief Solves a system of linear equations using Gaussian elimination.
 *
 * This function solves the system of linear equations Ax = b, where A is the
 * coefficients matrix and b is the solution vector (represented as a column matrix).
 * The function uses Gaussian elimination with partial pivoting and back-substitution.
 *
 * @param coefficients A pointer to the `zephyr_matrix` representing the coefficients matrix (A).
 * @param solution A pointer to the `zephyr_matrix` representing the solution vector (b), which must be a column matrix.
 * @return A pointer to the newly created `zephyr_vector` representing the solution vector x,
 *         or NULL if inputs are invalid, the system is singular, or memory allocation fails.
 */
zephyr_vector * simultaneous_eq(const zephyr_matrix * coefficients, const zephyr_matrix * solution);
/**
 * @brief Creates a Vandermonde matrix.
 *
 * A Vandermonde matrix is a matrix with the terms of a geometric progression in each row.
 * This function creates an n x n Vandermonde matrix from a given array of x values.
 *
 * @param x_values A pointer to an array of double values representing the x values.
 * @param n The size of the square Vandermonde matrix (n x n).
 * @return A pointer to the newly created `zephyr_matrix` (Vandermonde matrix),
 *         or NULL if memory allocation fails.
 */
/**
 * @brief Creates a Vandermonde matrix.
 *
 * A Vandermonde matrix is a matrix with the terms of a geometric progression in each row.
 * This function creates an n x n Vandermonde matrix from a given array of x values.
 *
 * @param x_values A pointer to an array of double values representing the x values.
 * @param n The size of the square Vandermonde matrix (n x n).
 * @return A pointer to the newly created `zephyr_matrix` (Vandermonde matrix),
 *         or NULL if memory allocation fails.
 */
zephyr_matrix *vandermonde_matrix(const double *x_values, size_t n);
/**
 * @brief Solves for the eigenvalues and eigenvectors of a square matrix.
 *
 * This function computes the eigenvalues and corresponding eigenvectors of a square matrix `A`.
 * The results are stored in the `zephyr_eigen_result` structure.
 *
 * @param A A pointer to the square `zephyr_matrix` for which to compute eigenvalues and eigenvectors.
 * @param result A pointer to a `zephyr_eigen_result` structure where the results will be stored.
 * @return True if the eigen-solution is successful, false otherwise (e.g., if the matrix is not square or memory allocation fails).
 */
bool eigen_solve(const zephyr_matrix *A, zephyr_eigen_result *result);
/**
 * @brief Destroys a `zephyr_eigen_result` structure and frees its allocated memory.
 *
 * This function deallocates the memory associated with the eigenvalues vector
 * and the eigenvectors matrix within the `zephyr_eigen_result` structure.
 * It is safe to call this function with a NULL result pointer.
 *
 * @param result A pointer to the `zephyr_eigen_result` to be destroyed.
 */
void destroy_eigen_result(zephyr_eigen_result *result);
/**
 * @brief Creates a 2D rotation matrix.
 *
 * This function generates a 2x2 rotation matrix for a given angle theta.
 *
 * @param to_rotate A pointer to the `zephyr_vector` to be rotated (expected to be 2D).
 * @param theta The angle of rotation in radians.
 * @return A pointer to the newly created 2x2 rotation `zephyr_matrix`,
 *         or NULL if the input vector is not 2D or memory allocation fails.
 */
zephyr_matrix *rotation_matrix_2d(const zephyr_vector *to_rotate, double theta);
/**
 * @brief Rotates a 3D vector about an arbitrary axis.
 *
 * This function applies a rotation to a 3D vector `v` around a specified `axis` by an angle `theta`.
 *
 * @param v A pointer to the `zephyr_vector` to be rotated (expected to be 3D).
 * @param axis A pointer to the `zephyr_vector` representing the axis of rotation (expected to be 3D).
 * @param theta The angle of rotation in radians.
 * @return A pointer to the newly created `zephyr_matrix` representing the rotated vector,
 *         or NULL if inputs are invalid or memory allocation fails.
 */
zephyr_matrix *rotate_about_axis(const zephyr_vector *v, const zephyr_vector *axis, double theta);
/**
 * @brief Translates a 2D vector by a given translation vector.
 *
 * This function applies a 2D translation to a point represented as a `zephyr_vector`.
 * It constructs a 3x3 homogeneous translation matrix and multiplies it with the
 * homogeneous representation of the input point.
 *
 * @param point A pointer to the `zephyr_vector` representing the 2D point to translate (expected to be 2D).
 * @param translation A pointer to the `zephyr_vector` representing the 2D translation to apply (expected to be 2D).
 * @return A pointer to the newly created `zephyr_vector` representing the translated point,
 *         or NULL if inputs are invalid or memory allocation fails.
 */
zephyr_vector *translate_2d_vector(const zephyr_vector *point, const zephyr_vector *translation);
/**
 * @brief Translates a 3D vector by a given translation vector.
 *
 * This function applies a 3D translation to a point represented as a `zephyr_vector`.
 * It constructs a 4x4 homogeneous translation matrix and multiplies it with the
 * homogeneous representation of the input point.
 *
 * @param point A pointer to the `zephyr_vector` representing the 3D point to translate (expected to be 3D).
 * @param translation A pointer to the `zephyr_vector` representing the 3D translation to apply (expected to be 3D).
 * @return A pointer to the newly created `zephyr_vector` representing the translated point,
 *         or NULL if inputs are invalid or memory allocation fails.
 */
zephyr_vector *translate_3d_vector(const zephyr_vector *point, const zephyr_vector *translation);
/**
 * @brief Rotates a 3D vector around a specified axis (X, Y, or Z).
 *
 * This function applies a rotation to a 3D vector `v` around one of the primary axes
 * (X, Y, or Z) by a given angle `theta`.
 *
 * @param v A pointer to the `zephyr_vector` to be rotated (expected to be 3D).
 * @param theta The angle of rotation in radians.
 * @param axis The axis of rotation (ROTATE_X, ROTATE_Y, or ROTATE_Z).
 * @return A pointer to the newly created `zephyr_vector` representing the rotated vector,
 *         or NULL if inputs are invalid or memory allocation fails.
 */
zephyr_vector *rotate_3d_per_axis(const zephyr_vector *v, double theta, rotation_axis axis);
/**
 * @brief Rotates a 3D vector by specified angles around the X, Y, and Z axes.
 *
 * This function applies a sequence of rotations (X, then Y, then Z) to a 3D vector `v`.
 *
 * @param v A pointer to the `zephyr_vector` to be rotated (expected to be 3D).
 * @param theta_x The angle of rotation around the X-axis in radians.
 * @param theta_y The angle of rotation around the Y-axis in radians.
 * @param theta_z The angle of rotation around the Z-axis in radians.
 * @return A pointer to the newly created `zephyr_vector` representing the rotated vector,
 *         or NULL if inputs are invalid or memory allocation fails.
 */
zephyr_vector *rotate_3d_vector(const zephyr_vector *v, double theta_x, double theta_y, double theta_z);
/**
 * @brief Creates a 3D rotation matrix from Euler angles (XYZ convention).
 *
 * This function constructs a 3x3 rotation matrix by composing individual rotations
 * around the X, Y, and Z axes in that order (R = Rz * Ry * Rx).
 *
 * @param theta_x The angle of rotation around the X-axis in radians.
 * @param theta_y The angle of rotation around the Y-axis in radians.
 * @param theta_z The angle of rotation around the Z-axis in radians.
 * @return A pointer to the newly created 3x3 rotation `zephyr_matrix`,
 *         or NULL if memory allocation fails.
 */
zephyr_matrix *rotation_matrix_3d_xyz(double theta_x, double theta_y, double theta_z);
/**
 * @brief Applies a 3D transformation (rotation and translation) to a vector.
 *
 * This function combines rotations around the X, Y, and Z axes with a translation
 * to transform a 3D vector. It constructs a 4x4 homogeneous transformation matrix
 * internally and applies it to the input vector.
 *
 * @param v A pointer to the `zephyr_vector` to be transformed (expected to be 3D).
 * @param theta_x The angle of rotation around the X-axis in radians.
 * @param theta_y The angle of rotation around the Y-axis in radians.
 * @param theta_z The angle of rotation around the Z-axis in radians.
 * @param translation A pointer to the `zephyr_vector` representing the 3D translation to apply.
 * @return A pointer to the newly created `zephyr_vector` representing the transformed result,
 *         or NULL if inputs are invalid or memory allocation fails.
 */
zephyr_vector *transform_3d(const zephyr_vector *v,
                            double theta_x, double theta_y, double theta_z,
                            const zephyr_vector *translation);
/**
 * @brief Applies the inverse of a 3D rigid body transformation to a vector.
 *
 * This function computes the inverse of a composite 3D transformation consisting of
 * rotation (around X, Y, and Z axes) and translation, and applies it to a 3D vector.
 *
 * The transformation is defined as:
 *     v' = Rᵀ * (v - t)
 * where R is the rotation matrix composed from Euler angles (theta_x, theta_y, theta_z),
 * Rᵀ is its transpose (i.e., inverse for orthogonal rotation matrices), and t is the translation vector.
 *
 * The function constructs a 4x4 inverse transformation matrix using Rᵀ and -Rᵀ * t,
 * applies it to the homogeneous form of the input vector, and returns the transformed 3D vector.
 *
 * @param v Pointer to the input 3D vector to be transformed.
 * @param theta_x Rotation angle (in radians) around the X-axis.
 * @param theta_y Rotation angle (in radians) around the Y-axis.
 * @param theta_z Rotation angle (in radians) around the Z-axis.
 * @param translation Pointer to the 3D translation vector.
 * @return A newly allocated 3D vector after applying the inverse transformation,
 *         or NULL on allocation or input error.
 */

zephyr_vector * inverse_transformation_3d(const zephyr_vector *v,
                            double theta_x, double theta_y, double theta_z,
                            const zephyr_vector *translation);

/**
 * @brief Composes two 3D rigid body transformations into a single transformation matrix.
 *
 * This function takes two transformations:
 *   - A → B: defined by rotation matrix R_AB and translation vector P_AB
 *   - B → C: defined by rotation matrix R_BC and translation vector P_BC
 *
 * It composes them into a single transformation A → C:
 *   - R_AC = R_AB * R_BC
 *   - P_AC = R_AB * P_BC + P_AB
 *
 * The result is returned as a 4x4 homogeneous transformation matrix combining
 * the rotation (upper-left 3x3) and translation (top-right 3x1), with the last
 * row being [0 0 0 1].
 *
 * @param R_AB Pointer to a 3x3 rotation matrix representing frame A to B.
 * @param P_AB Pointer to a 3D vector representing translation from A to B.
 * @param R_BC Pointer to a 3x3 rotation matrix representing frame B to C.
 * @param P_BC Pointer to a 3D vector representing translation from B to C.
 * @return A newly allocated 4x4 homogeneous transformation matrix from A to C,
 *         or NULL on invalid input or memory allocation failure.
 */

zephyr_matrix *compose_transformations(const zephyr_matrix *R_AB, const zephyr_vector *P_AB,
                                       const zephyr_matrix *R_BC, const zephyr_vector *P_BC);
/**
 * @brief Applies a 4x4 homogeneous transformation matrix to a 3D vector.
 *
 * This function takes a 3D vector and applies a 4x4 transformation matrix to it,
 * treating the vector as a point in homogeneous coordinates ([x, y, z, 1]^T).
 * It performs a matrix-vector multiplication and extracts the transformed 3D point.
 *
 * The transformation matrix can include rotation, translation, scaling, or any
 * affine transformation representable in homogeneous form.
 *
 * @param v Pointer to the input 3D vector.
 * @param T Pointer to a 4x4 transformation matrix.
 * @return A newly allocated transformed 3D vector, or NULL if input is invalid
 *         or memory allocation fails.
 */
zephyr_vector *transform_3d_with_matrix(const zephyr_vector *v, const zephyr_matrix *T);

/**
 * @brief Constructs a 4x4 homogeneous transformation matrix representing
 *        a rigid rotation about an arbitrary axis in 3D space.
 *
 * This function builds a transformation matrix that rotates a point around an arbitrary axis,
 * defined by a direction vector and a point on the axis, by an angle `theta_rad`.
 * The axis is first normalized. The rotation matrix is constructed using the Rodrigues'
 * rotation formula, and the translation is computed to ensure the axis point remains fixed.
 *
 * The final result is a 4x4 homogeneous matrix representing the rotation in space,
 * suitable for transforming 3D points with homogeneous coordinates.
 *
 * @param axis_point A 3D vector representing a point on the axis of rotation.
 * @param axis_dir A 3D vector representing the direction of the rotation axis.
 * @param theta_rad Rotation angle in radians.
 * @return A newly allocated 4x4 homogeneous transformation matrix,
 *         or NULL on invalid input or memory allocation failure.
 */

zephyr_matrix *rigid_transform_about_axis(const zephyr_vector *axis_point,
                                          const zephyr_vector *axis_dir,
                                          double theta_rad);

#endif //MATRIX_H