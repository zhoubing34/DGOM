//
// Created by li12242 on 12/10/16.
//
#include "LibUtilities_test.h"

int main(){
    int i, flag[2];

    printf(HEADSTART "Running 2 test from LibUtilities_Test\n");

    flag[0] = MatrixInverse_test();
    flag[1] = MatrixMultiply_test();

    int failNum = 0;
    for(i=0; i<2; i++){
        if(flag[i])
            failNum++;
    }
    if(failNum)
        printf(HEADFINISH "%d test faild from LibUtilities_Test\n", failNum);
    else
        printf(HEADFINISH "2 test passed from LibUtilities_Test\n");

    return 0;
}