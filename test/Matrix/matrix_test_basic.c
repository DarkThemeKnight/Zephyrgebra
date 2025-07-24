//
// Created by omotola-david-ayanfeoluwa on 7/24/25.
//
#include "../../include/Matrix.h"
#include <assert.h>

void test_create_destroy_matrix() {
    zephyr_matrix* mat = create_matrix(3, 4);
    assert(mat != NULL);
    assert(mat->m == 3);
    assert(mat->n == 4);
    destroy_matrix(mat);
}

void test_set_get_matrix() {
    zephyr_matrix* mat = create_matrix(3, 4);
    matrix_set(1, 2, 5.0, mat);
    assert(matrix_get(1, 2, mat) == 5.0);
    destroy_matrix(mat);
}

int main() {
    test_create_destroy_matrix();
    test_set_get_matrix();
    return 0;
}