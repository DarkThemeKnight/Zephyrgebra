//
// Created by omotola-david-ayanfeoluwa on 7/23/25.
//
#include "Matrix.h"
#include "Vector.h"
#include "Polynomials.h"

int main() {
    zephyr_vector * P1 = create_vector(3);
    P1->data[0] = 0; P1->data[1] = 2; P1->data[2] = 0;
    zephyr_vector * rotateAboutZ = rotate_3d(P1,30,ROTATE_Z);
    print_vector(rotateAboutZ,2);
}