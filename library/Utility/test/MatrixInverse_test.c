//
// Created by li12242 on 12/10/16.
//
#include "MatrixInverse_data.cc"
#include "utility_test.h"

int MatrixInverse_test(){

    // global variable
    extern double invCt[Nrow*Nrow], Ct[Nrow*Nrow];

    // local variable
    int fail = 0;
    double temp[Nrow*Nrow];
    int i;

    // assignment
    for(i = 0; i<Nrow*Nrow; i++){
        temp[i] = Ct[i];
    }

    // call
    Matrix_inverse(temp, Nrow);

    // check
    fail = Vector_test("Matrix_Inverse", temp, invCt, Nrow*Nrow);

    return fail;
}
