//
// Created by li12242 on 12/10/16.
//

#include "LibUtilities_test.h"
#include "MatrixMultiply_data.cc"

int MatrixMultiply_test(){
    // global
    extern double A[M*K], B[K*N], C[M*N];

    // local
    double Ct[M*N];
    int fail = 0;
    clock_t clockT1, clockT2;

    // call
    clockT1 = clock();
    Matrix_multiply(M, K, N, A, B, Ct);
    clockT2 = clock();

    // check Ct equals to C
    fail = Vector_test("MatrixMultiply_test", Ct, C, M*N, (double)((clockT2-clockT1)/CLOCKS_PER_SEC) );

    return fail;
}
