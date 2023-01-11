
#include "linalg.h"
#include <limits.h>

//////////////////////////////////////////////////
//     Util
//////////////////////////////////////////////////
void check_malloc_error(void * data){
    if (data == NULL){
        fprintf(stderr, "Error: failed malloc\n");
    } else {
        // do nothing
    }
}

int randint(int lower, int upper){
    return (rand() % abs(upper - lower)) + lower;
}

Matrix * matrix_init(size_t m, size_t n){
    if ((m <= 0) || (n <= 0)){
        fprintf(stderr, "Error: Invalid matrix dimension\n");
        return NULL;
    }
    Matrix * ret_matrix = (Matrix *)malloc(sizeof(Matrix));
    check_malloc_error(ret_matrix);

    ret_matrix->rows = m;
    ret_matrix->cols = n;
    ret_matrix->data = calloc(ret_matrix->rows, sizeof(*ret_matrix->data));
    check_malloc_error(ret_matrix->data);

    for (int i = 0; i < ret_matrix->rows; i++){
        ret_matrix->data[i] = calloc(ret_matrix->cols, sizeof(**ret_matrix->data));
        check_malloc_error(ret_matrix->data[i]);
    }

    return ret_matrix;
}

//////////////////////////////////////////////////
//        Vector
//////////////////////////////////////////////////
Vec * vector_add(Vec * v1, Vec * v2){
    assert(v1->size == v2->size);

    Vec * ret_vec = (Vec *)malloc(sizeof(Vec));
    ret_vec->size = v1->size;
    check_malloc_error((void *)ret_vec);
    ret_vec->data = (double *)malloc(sizeof(double) * v1->size);
    check_malloc_error((void *)ret_vec);

    ret_vec->size = v1->size;
    for (int i = 0; i < ret_vec->size; i++){
        ret_vec->data[i] = v1->data[i] + v2->data[i];
    }

    return ret_vec;
}

Vec * vector_subtract(Vec * v1, Vec * v2){
    assert(v1->size == v2->size);

    Vec * ret_vec = (Vec *)malloc(sizeof(Vec));
    ret_vec->size = v1->size;
    check_malloc_error((void *)ret_vec);
    ret_vec->data = (double *)malloc(sizeof(double) * v1->size);
    check_malloc_error((void *)ret_vec->data);

    ret_vec->size = v1->size;
    for (int i = 0; i < ret_vec->size; i++){
        ret_vec->data[i] = v1->data[i] - v2->data[i];
    }

    return ret_vec;
}

double vector_dot(Vec * v1, Vec * v2){
    assert(v1->size == v2->size);

    double ret_val = 0;
    for (int i = 0; i < v2->size; i++){
        ret_val += v1->data[i] * v2->data[i];
    }

    return ret_val;
}

Vec * vector_cross(Vec * v1, Vec * v2){
    assert(v1->size == 3);
    assert(v2->size == 3);

    Vec * ret_vec = (Vec *)malloc(sizeof(Vec));
    ret_vec->size = 3;
    check_malloc_error(ret_vec);
    ret_vec->data = (double *)malloc(3 * sizeof(double));
    check_malloc_error(ret_vec->data);

    ret_vec->data[0] = v1->data[1] * v2->data[2] - v1->data[2] * v2->data[1];
    ret_vec->data[1] = v1->data[2] * v2->data[0] - v1->data[0] * v2->data[2];
    ret_vec->data[2] = v1->data[0] * v2->data[1] - v1->data[1] * v2->data[0];

    return ret_vec;
}

Vec * vector_scalar_mult(Vec * v, int scalar){
    Vec * ret_vec = (Vec *)malloc(sizeof(Vec));
    check_malloc_error(ret_vec);
    ret_vec->size = v->size;
    ret_vec->data = (double *)malloc(sizeof(double) * v->size);
    check_malloc_error(ret_vec->data);

    for (int i = 0; i < v->size; i++){
        ret_vec->data[i] = scalar * v->data[i];
        printf("%f\n", ret_vec->data[i]);
    }
    return ret_vec;
}

void free_vector(Vec * v){
    free(v->data);
    free(v);
}

void print_vector(Vec * v){
    printf("[");
    for(int i = 0; i < v->size; i++){
        printf(" %f", v->data[i]);
    }
    printf("]\n");
}

bool vector_equal(Vec * v1, Vec * v2){
    assert(v1->size == v2->size);

    for (int i = 0; i < v1->size; i++){
        if (v1->data[i] != v2->data[i]){
            return false;
        }
    }
    return true;
}

bool vector_zero(Vec * v){
    for (int i = 0; i < v->size; i++){
        if (v->data[i] != 0){
            return false;
        }
    }
    return true;
}

//////////////////////////////////////////////////
//        Matrix
//////////////////////////////////////////////////
void assert_dim_equal(Matrix * m1, Matrix * m2){
    assert(m1->rows == m2->rows);
    assert(m1->cols == m2->cols);
}

void assert_dim(Matrix * m){
    assert(m->rows > 0);
    assert(m->cols > 0);
}

void assert_matrix_mult(Matrix * m1, Matrix * m2){
    assert(m1->cols == m2->rows);
}

bool matrix_equal(Matrix * m1, Matrix * m2){
    //assert_dim_equal(m1, m2);
    if ((m1->rows != m2->rows) || (m1->cols != m2->rows)){
        return false;
    }
    for (int i = 0; i < m1->rows; i++){
        for (int j = 0; i < m1->cols; j++){
            if (m1->data[i][j] != m2->data[i][j]){
                return false;
            }
        }
    }
    return true;
}

void print_matrix(Matrix m){
    for (int i = 0; i < m.rows; i++){
        printf("[ ");
        for (int j = 0; j < m.cols; j++){
            printf("%f ", m.data[i][j]);
        }
        printf("]");
        printf("\n");
    }
    printf("\n");
}

bool matrix_square(Matrix m){
    return (m.rows == m.cols);
}

Matrix * matrix_ones(size_t m, size_t n){
    if ((m <= 0) || (n <= 0)){
        fprintf(stderr, "Error: Invalid matrix dimension\n");
        return NULL;
    }
    Matrix * ret_matrix = (Matrix *)malloc(sizeof(Matrix));
    check_malloc_error(ret_matrix);

    ret_matrix->rows = m;
    ret_matrix->cols = n;
    ret_matrix->data = calloc(ret_matrix->rows, sizeof(*ret_matrix->data));
    check_malloc_error(ret_matrix->data);

    for (int i = 0; i < ret_matrix->rows; i++){
        ret_matrix->data[i] = calloc(ret_matrix->cols, sizeof(**ret_matrix->data));
        check_malloc_error(ret_matrix->data[i]);
    }

    for (int i = 0; i < ret_matrix->rows; i++){
        for (int j = 0; j < ret_matrix->cols; j++){
            ret_matrix->data[i][j] = 1;
        }
    }

    return ret_matrix;
}

Matrix * matrix_zeros(size_t m, size_t n){
    if ((m <= 0) || (n <= 0)){
        fprintf(stderr, "Error: Invalid matrix dimension\n");
        return NULL;
    }
    Matrix * ret_matrix = (Matrix *)malloc(sizeof(Matrix));
    check_malloc_error(ret_matrix);

    ret_matrix->rows = m;
    ret_matrix->cols = n;
    ret_matrix->data = calloc(ret_matrix->rows, sizeof(*ret_matrix->data));
    check_malloc_error(ret_matrix->data);

    for (int i = 0; i < ret_matrix->rows; i++){
        ret_matrix->data[i] = calloc(ret_matrix->cols, sizeof(**ret_matrix->data));
        check_malloc_error(ret_matrix->data[i]);
    }

    for (int i = 0; i < ret_matrix->rows; i++){
        for (int j = 0; j < ret_matrix->cols; j++){
            ret_matrix->data[i][j] = 0;
        }
    }

    return ret_matrix;
}

Matrix * matrix_identity(size_t dim){
    if (dim <= 0){
        fprintf(stderr, "Error: Invalid matrix dimension\n");
        return NULL;
    }
    Matrix * ret_matrix = (Matrix *)malloc(sizeof(Matrix));
    check_malloc_error(ret_matrix);

    ret_matrix->rows = dim;
    ret_matrix->cols = dim;
    ret_matrix->data = calloc(ret_matrix->rows, sizeof(*ret_matrix->data));
    check_malloc_error(ret_matrix->data);

    for (int i = 0; i < ret_matrix->rows; i++){
        ret_matrix->data[i] = calloc(ret_matrix->cols, sizeof(**ret_matrix->data));
        check_malloc_error(ret_matrix->data[i]);
    }

    for (size_t i = 0; i < dim; i++){
        ret_matrix->data[i][i] = 1;
    }
    /* for (int i = 0; i < ret_matrix->rows; i++){ */
    /*     for (int j = 0; j < ret_matrix->cols; j++){ */
    /*         if (i == j){} */
    /*     } */
    /* } */
    return ret_matrix;
}


Matrix * matrix_random(size_t m, size_t n, int lower_bound, int upper_bound){
    srand(time(NULL));

    if ((m <= 0) || (n <= 0)){
        fprintf(stderr, "Error: Invalid matrix dimension\n");
        return NULL;
    }
    Matrix * ret_matrix = (Matrix *)malloc(sizeof(Matrix));
    check_malloc_error(ret_matrix);

    ret_matrix->rows = m;
    ret_matrix->cols = n;
    ret_matrix->data = calloc(ret_matrix->rows, sizeof(*ret_matrix->data));
    check_malloc_error(ret_matrix->data);

    for (int i = 0; i < ret_matrix->rows; i++){
        ret_matrix->data[i] = calloc(ret_matrix->cols, sizeof(**ret_matrix->data));
        check_malloc_error(ret_matrix->data[i]);
    }

    for (int i = 0; i < ret_matrix->rows; i++){
        for (int j = 0; j < ret_matrix->cols; j++){
            int rr = randint(lower_bound, upper_bound); // just var to save random int
            ret_matrix->data[i][j] = rr;
        }
    }

    return ret_matrix;
}

Matrix * matrix_setval(size_t m, size_t n, double val){
    if ((m <= 0) || (n <= 0)){
        fprintf(stderr, "Error: Invalid matrix dimension\n");
        return NULL;
    }
    Matrix * ret_matrix = (Matrix *)malloc(sizeof(Matrix));
    check_malloc_error(ret_matrix);

    ret_matrix->rows = m;
    ret_matrix->cols = n;
    ret_matrix->data = calloc(ret_matrix->rows, sizeof(*ret_matrix->data));
    check_malloc_error(ret_matrix->data);

    for (int i = 0; i < ret_matrix->rows; i++){
        ret_matrix->data[i] = calloc(ret_matrix->cols, sizeof(**ret_matrix->data));
        check_malloc_error(ret_matrix->data[i]);
    }

    for (int i = 0; i < ret_matrix->rows; i++){
        for (int j = 0; j < ret_matrix->cols; j++){
            ret_matrix->data[i][j] = val;
        }
    }

    return ret_matrix;
}

Matrix * matrix_add(Matrix * m1, Matrix * m2){
    assert(m1->rows == m2->rows);
    assert(m1->cols == m2->cols);
    Matrix * ret_mat = matrix_init(m1->rows, m1->cols);
    ret_mat->rows = m1->rows;
    ret_mat->cols = m1->cols;

    for (int i = 0; i < ret_mat->rows; i++){
        for (int j =0; j < ret_mat->cols; j++){
            ret_mat->data[i][j] = m1->data[i][j] + m2->data[i][j];
        }
    }
    return ret_mat;
}


Matrix * matrix_subtract(Matrix * m1, Matrix * m2){
    assert(m1->rows == m2->rows);
    assert(m1->cols == m2->cols);
    Matrix * ret_mat = matrix_init(m1->rows, m1->cols);
    ret_mat->rows = m1->rows;
    ret_mat->cols = m1->cols;

    for (int i = 0; i < ret_mat->rows; i++){
        for (int j =0; j < ret_mat->cols; j++){
            ret_mat->data[i][j] = m1->data[i][j] - m2->data[i][j];
        }
    }
    return ret_mat;
}

/////// Matrix multiplication
Matrix * _matrix_multiplication_naive(Matrix * m1, Matrix * m2){
    //assert(m1->cols == m2->rows);
    assert_matrix_mult(m1, m2);
    Matrix * ret_mat = matrix_init(m1->rows, m2->cols);
    ret_mat->rows = m1->rows;
    ret_mat->cols = m2->cols;
    //printf("row: %d....col: %d", ret_mat->rows, ret_mat->cols);
    for (int row = 0; row < ret_mat->rows; row++){
        for (int col = 0; col < ret_mat->cols; col++){

            double sum = 0;
            for (int k = 0; k < m1->cols; k++){
                sum += m1->data[row][k] * m2->data[k][col];
            }
            ret_mat->data[row][col] = sum;

        }
    }
    return ret_mat;
}

Matrix * _matrix_multiplication_strassen(Matrix * m1, Matrix * m2){
    return NULL;
}

Matrix * matrix_multiplication_method(Matrix * m1, Matrix * m2, char * method){
    if (strncmp(method, "naive", 100) == 0){
        return _matrix_multiplication_naive(m1, m2);
    } else if (strncmp(method, "strassen", 100) == 0){
        return _matrix_multiplication_strassen(m1, m2);
    }
    return NULL;
}

Vec * matrix_mult_vec(Matrix * m, Vec * v){
    Vec * ret = (Vec *)malloc(sizeof(Vec));
    check_malloc_error(ret);
    ret->size = m->rows;
    ret->data = (double *)malloc(sizeof(double) * ret->size);
    check_malloc_error(ret->data);

    int sum;
    for (int i = 0; i < ret->size; i++){
        sum = 0;
        for (int k = 0; k < m->cols; k++){
            sum += v->data[k] * m->data[i][k];
        }
        ret->data[i] = sum;
    }

    return ret;
}


// P_vec = A x (Br) - Cr
bool matrix_multiplication_guess(Matrix * A, Matrix * B, Matrix * C){
    assert_dim(A);
    assert_dim(B);
    // prereq for frievald's algorithm
    assert(A->rows == B->rows);
    assert(A->cols == B->cols);
    Vec * r_vec = (Vec *)malloc(sizeof(Vec));
    check_malloc_error(r_vec);
    r_vec->size = A->rows;
    r_vec->data = (double *)malloc(A->rows);
    check_malloc_error(r_vec->data);
    int rand_i;
    for (int i = 0; i < r_vec->size; i++){
        srand(time(NULL));
        rand_i = rand() % 2;
        r_vec->data[i] = rand_i;
    }
    Vec * Br = matrix_mult_vec(B, r_vec);
    Br->size = A->rows;

    Vec * ABr = (Vec *)malloc(sizeof(Vec));
    check_malloc_error(ABr);
    ABr->size = A->rows;
    ABr = matrix_mult_vec(A, Br);

    Vec * Cr = (Vec *)malloc(sizeof(Vec));
    check_malloc_error(Cr);
    Cr->size = A->rows;
    Cr = matrix_mult_vec(C, r_vec);

    Vec * P = vector_subtract(ABr, Cr);
    if (vector_zero(P)){
        return true;
    } else {
        return false;
    }
    free_vector(r_vec);
    free_vector(Br);
    free_vector(ABr);
    free_vector(Cr);
    free_vector(P);
}

Matrix * matrix_transpose(Matrix * m){
    Matrix * ret_mat = matrix_init(m->cols, m->rows);
    for (int row = 0; row < m->rows; row++){
        for (int col = 0; col < m->cols; col++){
            ret_mat->data[col][row] = m->data[row][col];
        }
    }
    return ret_mat;
}

void copy_matrix(Matrix * dest, Matrix * source){
    dest = matrix_init(source->rows, source->cols);
    for (int i = 0; i < dest->rows; i++){
        for (int j = 0; j < dest->cols; j++){
            dest->data[i][j] = source->data[i][j];
        }
    }
}

// elementary row operations
void matrix_elem_row_switch(Matrix * m, unsigned short row1, unsigned short row2){
    assert(row1 < m->rows && row2 < m->rows);
    for (int col = 0; col < m->cols; col++){
        int temp_row1[m->cols];
        temp_row1[col] = m->data[row1][col];
        m->data[row1][col] = m->data[row2][col];
        m->data[row2][col] = temp_row1[col];
    }
}
void matrix_elem_constant(Matrix * m, unsigned short row, int constant){
    assert(row < m->rows);
    assert(constant != 0);
    for (int i = 0; i < m->cols; i++){
        m->data[row][i] *= constant;
    }
}
void matrix_elem_row_add(Matrix * m, unsigned short row, unsigned short row2, int constant){
    assert(row < m->rows);
    assert(row2 < m->rows);
    assert(constant != 0);
    for (int i = 0; i < m->cols; i++){
        m->data[row][i] += constant * m->data[row2][i];
    }
}

bool matrix_is_zero(Matrix * m){
    for (int i = 0; i < m->rows; i++){
        for (int j = 0; j < m->cols; j++){
            if (m->data[i][j] != 0){
                return false;
            }
        }
    }
    return true;
}

// find the first non-zero element on the col column, under the row row
int _pivot_upper_triangle(Matrix * m, unsigned short col, unsigned short row){
    int i;
    for (i = row; i < m->rows; i++){
        if (fabs(m->data[i][col]) > INT_MIN){
            return i;
        }
    }
    return -1;
}

// returns a new matrix instead of in-place modification
Matrix * matrix_upper_triangle(Matrix * m){
    // for now assume square matrix
    /* if (matrix_is_zero(m)){ */
    /*     Matrix * ret = matrix_init(m->rows, m->cols); */
    /*     copy_matrix(ret, m); */
    /* } */
    //int lead = 0;
    Matrix * ret = matrix_init(m->rows, m->cols);
    copy_matrix(ret, m);

    int i, j, k, pivot;
    j = 0; i = 0;

}
