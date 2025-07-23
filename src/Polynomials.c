//
// Created by omotola-david-ayanfeoluwa on 7/22/25.
//
#include "../include/Polynomials.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>

#include "Matrix.h"
#include <float.h>
#include <stdbool.h>

void normalize_polynomial(zephyr_polynomial *poly) {
    if (!poly || poly->degree == 0) return;
    size_t actual_degree = poly->degree;
    while (actual_degree > 0 && fabs(poly->coefficients[actual_degree]) < 1e-12) {
        actual_degree--;
    }
    if (actual_degree == poly->degree) return;
    double *trimmed = realloc(poly->coefficients, sizeof(double) * (actual_degree + 1));
    if (trimmed) {
        poly->coefficients = trimmed;
        poly->degree = actual_degree;
    }
}
int len(const zephyr_polynomial *poly) {
    if (!poly) return 0;
    const double n = log2((int) poly->degree);
    const int nAsInt = ceil(n);
    return nAsInt;
}
zephyr_polynomial *create_polynomial(const double *coeffs, const size_t degree) {
    if (coeffs == NULL || degree == 0) {
        return NULL;
    }
    zephyr_polynomial *poly = malloc(sizeof(zephyr_polynomial));
    if (poly == NULL) {
        return NULL;
    }
    poly->degree = degree;
    double *memCoefficient = malloc(sizeof(double) * (degree + 1));
    if (memCoefficient == NULL) {
        free(poly);
        return NULL;
    }
    poly->coefficients = memCoefficient;
    memcpy(poly->coefficients, coeffs, sizeof(double) * (degree + 1));
    normalize_polynomial(poly);
    return poly;
}


void destroy_polynomial(zephyr_polynomial *poly) {
    if (!poly) return;
    free(poly->coefficients);
    free(poly);
}

zephyr_polynomial *copy_polynomial(const zephyr_polynomial *poly) {
    if (!poly) return NULL;
    return create_polynomial(poly->coefficients, poly->degree);
}

void print_polynomial(const zephyr_polynomial *poly, int decimal_places) {
    if (!poly) {
        printf("NULL POLYNOMIAL\n");
        return;
    }

    char format[16];
    snprintf(format, sizeof(format), "%%.%df", decimal_places);

    for (size_t i = 0; i <= poly->degree; i++) {
        printf("xpow%zu=", i);
        printf(format, poly->coefficients[i]);
        if (i < poly->degree) {
            printf(", ");
        }
    }
    printf("\n");
}

zephyr_polynomial *add_polynomials(const zephyr_polynomial *a, const zephyr_polynomial *b) {
    if (!a || !b) return NULL;
    const size_t max_degree = (a->degree > b->degree) ? a->degree : b->degree;
    double *solution = calloc(max_degree + 1, sizeof(double));
    if (!solution) return NULL;
    for (size_t i = 0; i <= max_degree; i++) {
       const double coeff_a = (i <= a->degree) ? a->coefficients[i] : 0.0;
       const double coeff_b = (i <= b->degree) ? b->coefficients[i] : 0.0;
        solution[i] = coeff_a + coeff_b;
    }
    zephyr_polynomial *sum = create_polynomial(solution, max_degree);
    free(solution);
    return sum;
}

zephyr_polynomial *subtract_polynomials(const zephyr_polynomial *a, const zephyr_polynomial *b) {
    if (!a || !b) return NULL;
    const size_t max_degree = (a->degree > b->degree) ? a->degree : b->degree;
    double *solution = calloc(max_degree + 1, sizeof(double));
    if (!solution) return NULL;
    for (size_t i = 0; i <= max_degree; i++) {
        const double coeff_a = (i <= a->degree) ? a->coefficients[i] : 0.0;
        const double coeff_b = (i <= b->degree) ? b->coefficients[i] : 0.0;
        solution[i] = (coeff_a - coeff_b);
    }
    zephyr_polynomial *sum = create_polynomial(solution, max_degree);
    free(solution);
    return sum;
}

size_t next_power_of_two(size_t n) {
    size_t power = 1;
    while (power < n) power <<= 1;
    return power;
}

double complex *to_complex_array(const zephyr_polynomial *poly, size_t size) {
    double complex *arr = calloc(size, sizeof(double complex));
    for (size_t i = 0; i <= poly->degree; i++) {
        arr[i] = poly->coefficients[i];
    }
    return arr;
}

void FFT_recursive(double complex *a, size_t n, int invert) {
    if (n == 1) return;

    double complex *a0 = malloc(n / 2 * sizeof(double complex));
    double complex *a1 = malloc(n / 2 * sizeof(double complex));
    for (size_t i = 0; i < n / 2; i++) {
        a0[i] = a[i * 2];
        a1[i] = a[i * 2 + 1];
    }

    FFT_recursive(a0, n / 2, invert);
    FFT_recursive(a1, n / 2, invert);

    for (size_t k = 0; k < n / 2; k++) {
        double complex w = cexp((invert ? -2.0 : 2.0) * M_PI * I * k / n);
        a[k] = a0[k] + w * a1[k];
        a[k + n / 2] = a0[k] - w * a1[k];
        if (invert) {
            a[k] /= 2;
            a[k + n / 2] /= 2;
        }
    }

    free(a0);
    free(a1);
}

zephyr_polynomial *fft_multiply(const zephyr_polynomial *a, const zephyr_polynomial *b) {
    if (!a || !b) return NULL;

    size_t result_degree = a->degree + b->degree;
    size_t n = next_power_of_two(result_degree + 1);

    double complex *fa = to_complex_array(a, n);
    double complex *fb = to_complex_array(b, n);

    FFT_recursive(fa, n, 0);
    FFT_recursive(fb, n, 0);

    for (size_t i = 0; i < n; i++) {
        fa[i] *= fb[i];
    }

    FFT_recursive(fa, n, 1);

    double *result_coeffs = malloc(sizeof(double) * (result_degree + 1));
    for (size_t i = 0; i <= result_degree; i++) {
        result_coeffs[i] = creal(fa[i]);
    }

    double *reversed_coeffs = malloc(sizeof(double) * (result_degree + 1));
    if (!reversed_coeffs) {
        free(fa); free(fb); free(result_coeffs);
        return NULL;
    }

    zephyr_polynomial *result = create_polynomial(result_coeffs, result_degree);
    free(reversed_coeffs);
    free(fa);
    free(fb);
    free(result_coeffs);

    return result;
}

zephyr_polynomial *derivative_polynomial(const zephyr_polynomial *polynomials, const int nth) {
    if (!polynomials || nth < 0) return NULL;

    // Base case 0th derivative is a copy of the original
    if (nth == 0) return copy_polynomial(polynomials);

    // Base case: derivative of constant is zero
    if (polynomials->degree == 0) {
        double zero = 0.0;
        return create_polynomial(&zero, 0);
    }

    // Compute first derivative
    const size_t result_degree = polynomials->degree - 1;
    double result_coeffs[result_degree + 1];

    for (size_t i = 0; i <= result_degree; i++) {
        result_coeffs[i] = ((int)i + 1) * polynomials->coefficients[i + 1];
    }

    zephyr_polynomial *first_derivative = create_polynomial(result_coeffs, result_degree);

    // Recurse for (nth - 1) derivative
    zephyr_polynomial *nth_derivative = derivative_polynomial(first_derivative, nth - 1);
    destroy_polynomial(first_derivative);

    return nth_derivative;
}

double evaluate_polynomial(const zephyr_polynomial *poly, double x) {
    if (!poly) return 0.0;
    double result = poly->coefficients[poly->degree];
    for (int i = poly->degree - 1; i >= 0; i--) {
        result = result * x + poly->coefficients[i];
    }
    return result;
}

double hybrid_newton_solve_polynomial(const zephyr_polynomial *poly,
                                       double a, double b,
                                       double initial_guess,
                                       int max_iterations,
                                       double tolerance) {
    if (!poly || max_iterations <= 0) return NAN;

    double fa = evaluate_polynomial(poly, a);
    double fb = evaluate_polynomial(poly, b);

    // Ensure root is bracketed
    if (fa * fb > 0.0) return NAN;

    zephyr_polynomial *derivative = derivative_polynomial(poly, 1);
    if (!derivative) return NAN;

    double x = initial_guess;

    for (int i = 0; i < max_iterations; i++) {
        double fx = evaluate_polynomial(poly, x);
        double dfx = evaluate_polynomial(derivative, x);

        double step_valid = fabs(dfx) > DBL_EPSILON;
        double x_next;

        if (step_valid) {
            double newton_step = fx / dfx;
            x_next = x - newton_step;

            // Check if Newton step is out of bounds
            if (x_next < a || x_next > b) {
                x_next = (a + b) / 2.0;  // Fallback to bisection
            }
        } else {
            x_next = (a + b) / 2.0;  // Fallback to bisection
        }

        double f_next = evaluate_polynomial(poly, x_next);

        // Convergence check
        if (fabs(f_next) < tolerance || fabs(x_next - x) < tolerance) {
            destroy_polynomial(derivative);
            return x_next;
        }

        // Update bracket
        if (fa * f_next < 0.0) {
            b = x_next;
            fb = f_next;
        } else {
            a = x_next;
            fa = f_next;
        }

        x = x_next;
    }

    destroy_polynomial(derivative);
    return x;  // Best guess
}

double hybrid_secant_solve_polynomial(const zephyr_polynomial *poly,
                                       double x0, double x1,
                                       int max_iterations, double tolerance) {
    if (!poly || max_iterations <= 0) return NAN;

    double f0 = evaluate_polynomial(poly, x0);
    double f1 = evaluate_polynomial(poly, x1);

    // Optional: ensure bracketing
    int has_bracket = f0 * f1 < 0;

    for (int i = 0; i < max_iterations; i++) {
        double denominator = f1 - f0;

        double x2;
        int use_secant = fabs(denominator) > DBL_EPSILON;

        if (use_secant) {
            x2 = x1 - f1 * (x1 - x0) / denominator;
        } else {
            // Fallback midpoint if secant step is unstable
            if (has_bracket) {
                x2 = (x0 + x1) / 2.0;
            } else {
                return NAN;
            }
        }

        // If secant guess is outside [x0, x1] and bracketing is known, fallback
        if (has_bracket && (x2 < fmin(x0, x1) || x2 > fmax(x0, x1))) {
            x2 = (x0 + x1) / 2.0;
        }

        double f2 = evaluate_polynomial(poly, x2);

        if (fabs(f2) < tolerance || fabs(x2 - x1) < tolerance) {
            return x2;  // Converged
        }

        // Update history
        x0 = x1;
        f0 = f1;
        x1 = x2;
        f1 = f2;
    }

    return x1; // Return best estimate
}

double hybrid_solve_polynomial(const zephyr_polynomial *poly,
                                double a, double b,
                                int max_iterations,
                                double tolerance) {
    if (!poly) return NAN;

    double fa = evaluate_polynomial(poly, a);
    double fb = evaluate_polynomial(poly, b);

    // Ensure root is bracketed
    if (fa * fb > 0.0) return NAN;

    double prev = a;
    double curr = b;
    double fprev = fa;
    double fcurr = fb;
    double next = curr;

    for (int i = 0; i < max_iterations; i++) {
        // Secant step
        double denominator = fcurr - fprev;
        bool secant_valid = fabs(denominator) > DBL_EPSILON;

        if (secant_valid) {
            next = curr - fcurr * (curr - prev) / denominator;
        } else {
            secant_valid = false;
        }

        // If secant goes out of bounds or fails, fallback to bisection
        if (!secant_valid || next < a || next > b) {
            next = (a + b) / 2.0;
        }

        double fnext = evaluate_polynomial(poly, next);

        // Convergence check
        if (fabs(fnext) < tolerance || fabs(b - a) < tolerance) {
            return next;
        }

        // Update bracket
        if (fa * fnext < 0.0) {
            b = next;
            fb = fnext;
        } else {
            a = next;
            fa = fnext;
        }

        // Update for next secant step
        prev = curr;
        fprev = fcurr;
        curr = next;
        fcurr = fnext;
    }

    return curr; // Best guess
}

bool durand_kerner_solve(const zephyr_polynomial *poly, double complex *roots, int max_iterations, double tolerance) {
    if (!poly || !roots || poly->degree == 0) return false;

    int n = poly->degree;
    double complex *curr = malloc(sizeof(double complex) * n);
    if (!curr) return false;

    // Initialize guesses on the unit circle
    const double PI = acos(-1);
    for (int i = 0; i < n; ++i) {
        double angle = 2.0 * PI * i / n;
        curr[i] = cos(angle) + sin(angle) * I;
    }

    bool converged = false;

    for (int iter = 0; iter < max_iterations; ++iter) {
        converged = true;

        for (int i = 0; i < n; ++i) {
            double complex xi = curr[i];
            double complex fx = 0;

            // Evaluate P(xi) using Horner's method
            for (int j = n; j >= 0; --j) {
                fx = fx * xi + poly->coefficients[j];
            }

            double complex denom = 1.0;
            for (int j = 0; j < n; ++j) {
                if (j != i) {
                    denom *= (xi - curr[j]);
                }
            }

            if (cabs(denom) < 1e-12) {
                free(curr);
                return false;  // Risk of division by zero
            }

            double complex next = xi - fx / denom;

            if (cabs(next - xi) > tolerance) {
                converged = false;
            }

            curr[i] = next;
        }

        if (converged) break;
    }

    // Copy results to output array
    for (int i = 0; i < n; ++i) {
        roots[i] = curr[i];
    }

    free(curr);
    return converged;
}

zephyr_polynomial *deflate_polynomial(const zephyr_polynomial *poly, double root) {
    if (!poly || poly->degree == 0) return NULL;

    size_t new_degree = poly->degree - 1;
    double *new_coeffs = malloc(sizeof(double) * (new_degree + 1));
    if (!new_coeffs) return NULL;

    // Synthetic division
    new_coeffs[new_degree] = poly->coefficients[poly->degree];  // Start with leading coeff

    for (int i = (int)new_degree - 1; i >= 0; i--) {
        new_coeffs[i] = poly->coefficients[i + 1] + root * new_coeffs[i + 1];
    }

    double remainder = poly->coefficients[0] + root * new_coeffs[0];  // Optional

    // (Optional) Warn if root is not accurate
    if (fabs(remainder) > 1e-6) {
        fprintf(stderr, "⚠️ Warning: root %.10f is not accurate (remainder = %.10f)\n", root, remainder);
    }

    zephyr_polynomial *deflated = create_polynomial(new_coeffs, new_degree);
    free(new_coeffs);
    return deflated;
}

zephyr_polynomial *expand_from_roots(const double *roots, size_t count) {
    if (!roots || count == 0) return NULL;

    // Start with P(x) = 1
    double *coeffs = malloc(sizeof(double));
    coeffs[0] = 1;
    size_t degree = 0;

    for (size_t i = 0; i < count; i++) {
        double r = roots[i];

        // New degree will be degree + 1
        size_t new_degree = degree + 1;
        double *new_coeffs = calloc(new_degree + 1, sizeof(double));

        // Multiply (x - r) with current polynomial
        for (size_t j = 0; j <= degree; j++) {
            new_coeffs[j]     -= coeffs[j] * r;   // constant * -r
            new_coeffs[j + 1] += coeffs[j];       // constant * x
        }

        free(coeffs);
        coeffs = new_coeffs;
        degree = new_degree;
    }

    zephyr_polynomial *result = create_polynomial(coeffs, degree);
    free(coeffs);
    return result;
}
