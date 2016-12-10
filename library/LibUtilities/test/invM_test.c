//
// Created by li12242 on 12/10/16.
//
#include "invM_data.cc"
#include "LibUtilities_test.h"

int invM_test(){

    // global variable
    extern double invCt[Nrow*Nrow], Ct[Nrow*Nrow];
    extern int N;

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
    invM(temp, Nrow);
    clockT2 = clock();

    // check
    fail = Vector_test("invM", temp, invCt, Nrow*Nrow, (double)((clockT2-clockT1)/CLOCKS_PER_SEC) );

    return fail;
}
