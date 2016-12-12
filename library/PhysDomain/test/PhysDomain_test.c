//
// Created by li12242 on 12/11/16.
//

#include "PhysDomain_test.h"

int main(){
    int i, flag[2];
    int failNum=0;

    printf(HEADSTART "Running 1 test from PhysDomain_test\n");

    for(i=0; i<2; i++){
        if(flag[i])
            failNum++;
    }
    if(failNum)
        printf(HEADFINISH "%d test faild from PhysDomain_test\n", failNum);
    else
    printf(HEADFINISH "1 test passed from PhysDomain_test\n");

    return 0;
}