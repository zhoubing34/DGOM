//
// Created by li12242 on 12/10/16.
//
#include "MatrixInverse_data.cc"
#include "LibUtilities_test.h"

int MatrixInverse_test(){

    // global variable
    extern double invCt[Nrow*Nrow], Ct[Nrow*Nrow];

    // local variable
    int fail = 0;
    double temp[Nrow*Nrow];
    int i;
    clock_t clockT1, clockT2;

    // assignment
    for(i = 0; i<Nrow*Nrow; i++){
        temp[i] = Ct[i];
    }

    // call
    clockT1 = clock();
    Matrix_Inverse(temp, Nrow);
    clockT2 = clock();

    // check
    fail = Vector_test("Matrix_Inverse", temp, invCt, Nrow*Nrow, (double)((clockT2-clockT1)/CLOCKS_PER_SEC) );

    return fail;
}
