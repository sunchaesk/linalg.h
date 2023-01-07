// TODO:
// function pointer

#ifndef LINALG_H
#define LINALG_H

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <time.h>
#include <math.h>
#include <string.h>


typedef struct {
    int size;
    double * data;
} Vec;

typedef struct {
    int rows;
    int cols;
    double ** data;
} Matrix;

// Vector
Vec * vector_add(Vec * v1, Vec * v2);
Vec * vector_subtract(Vec * v1, Vec * v2);
double vector_dot(Vec * v1, Vec * v2);
Vec * vector_cross(Vec * v1, Vec * v2);
Vec * vector_scalar_mult(Vec * v, int scalar);
void free_vector(Vec * v);
void print_vector(Vec * v);
bool vector_equal(Vec * v1, Vec * v2);
//////////////////////////////////////

// Matrix Operations
bool matrix_equal(Matrix * m1, Matrix * m2);
bool matrix_square(Matrix m);
void print_matrix(Matrix m);

Matrix * matrix_ones(size_t m, size_t n);
Matrix * matrix_zeros(size_t m, size_t n);
Matrix * matrix_identity(size_t dim);
Matrix * matrix_random(size_t m, size_t n, int lower_bound, int upper_bound);
Matrix * matrix_setval(size_t m, size_t n, double val);


Matrix * matrix_add(Matrix * m1, Matrix * m2);
Matrix * matrix_subtract(Matrix * m1, Matrix * m2);

//////////
Vec * matrix_mult_vec(Matrix * m, Vec * v);
Vec * vec_mult_matrix(Vec * v, Matrix * m);
Matrix * matrix_multiplication_method(Matrix * m1, Matrix * m2, char * method);
//#define add(a, b) _Generic(a, int: addi, char*: adds)(a, b)
Matrix * _matrix_multiplication_naive(Matrix * m1, Matrix * m2);
// TODO
Matrix * _matrix_multiplication_strassen(Matrix * m1, Matrix * m2);
//////////

// implementation of Frievald's algorithm
bool matrix_multiplication_guess(Matrix * m1, Matrix * m2, Matrix * mguess);


Matrix * matrix_transpose(Matrix * m);

Matrix * matrix_upper_triangle(Matrix * m1);
Matrix * matrix_gaussian_elimination(Matrix * m1);

// A = CR
int matrix_rank(Matrix * m);
int matrix_span(Matrix * m);



// add matrix inverse stuff
Matrix * matrix_inverse(Matrix * m);

double matrix_det(Matrix * m);

// Matrix Data

// Matrix Edit
void matrix_edit_remove_col(Matrix * m, size_t col);
void matrix_edit_remove_row(Matrix * m, size_t row);
void matrix_edit_change_cell(Matrix * m, size_t row, size_t col, int value);

// util
//void malloc_error();
void free_matrix(Matrix * m);

#endif
