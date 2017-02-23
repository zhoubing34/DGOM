//
// Created by li12242 on 12/10/16.
//

#include "utility_test.h"
#include "MatrixMultiply_data.cc"

int MatrixMultiply_test(){
    // global
    extern double A[M*K], B[K*N], C[M*N];

    // local
    double Ct[M*N];
    int fail = 0;

    // call
    Matrix_multiply(M, K, N, A, B, Ct);

    // check Ct equals to C
    fail = vector_double_test("MatrixMultiply_test", Ct, C, M * N);

    return fail;
}
