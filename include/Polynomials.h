#ifndef ZEPHYR_POLYNOMIAL_H
#define ZEPHYR_POLYNOMIAL_H

#include <stddef.h>
#include <stdbool.h>
#include <complex.h>
#include "Matrix.h"
#include "Vector.h"

#ifdef __cplusplus
extern "C" {
#endif


    typedef struct {
        double *coefficients; // Coefficients from lowest to highest degree
        size_t degree;        // Highest non-zero degree of the polynomial
    } zephyr_polynomial;

    // Memory Management
    zephyr_polynomial *create_polynomial(const double *coeffs, size_t degree);
    void destroy_polynomial(zephyr_polynomial *poly);
    zephyr_polynomial *copy_polynomial(const zephyr_polynomial *poly);

    // Arithmetic Operations
    zephyr_polynomial *add_polynomials(const zephyr_polynomial *a, const zephyr_polynomial *b);
    zephyr_polynomial *subtract_polynomials(const zephyr_polynomial *a, const zephyr_polynomial *b);
    zephyr_polynomial *fft_multiply(const zephyr_polynomial *a, const zephyr_polynomial *b);

    // Calculus
    zephyr_polynomial * derivative_polynomial(const zephyr_polynomial *polynomials, int nth);

    // Evaluation
    double evaluate_polynomial(const zephyr_polynomial *poly, double x);
    bool durand_kerner_solve(const zephyr_polynomial *poly, double complex *roots, int max_iterations, double tolerance);
    double hybrid_solve_polynomial(const zephyr_polynomial *poly,double a, double b,int max_iterations,double tolerance);
    double hybrid_secant_solve_polynomial(const zephyr_polynomial *poly,double x0, double x1,int max_iterations, double tolerance);
    // Utilities
    void print_polynomial(const zephyr_polynomial *poly, int decimal_places);
    bool equals_polynomial(const zephyr_polynomial *a, const zephyr_polynomial *b, double epsilon);
    zephyr_polynomial *deflate_polynomial(const zephyr_polynomial *poly, double root);
    zephyr_polynomial *expand_from_roots(const double *roots, size_t count);


#ifdef __cplusplus
}
#endif

#endif // ZEPHYR_POLYNOMIAL_H