
    /* Vec * v = malloc(sizeof(Vec)); */
    /* v->size = 3; */
    /* v->data = malloc(3 * sizeof(double)); */
    /* v->data[0] = 1; */
    /* v->data[1] = 2; */
    /* v->data[2] = 3; */
    /* print_vector(v); */

    /* Vec * v1 = malloc(sizeof(Vec)); */
    /* v1->size = 3; */
    /* v1->data = malloc(3 * sizeof(double)); */
    /* v1->data[0] = 2; */
    /* v1->data[1] = 4; */
    /* v1->data[2] = 6; */
    /* print_vector(v1); */

    /* Vec * ret; */
    /* ret = vector_add(v1, v); */
    /* print_vector(ret); */
#include "linalg.h"

int test(int a, int b){
    return a + b;
}


int main(){
    Matrix * mat = matrix_setval(2, 2, 10);
    Vec v = {.size = 5, .data = (double[]){3, 5}};
    Vec * vp = &v;
    Vec * vv = matrix_mult_vec(mat, vp);
    print_vector(vv);

    return 0;
}
