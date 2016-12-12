//
// Created by li12242 on 12/11/16.
//

#include "PhysDomain_test.h"

#define NTEST 2

int main(){
    int i, flag[NTEST];
    int failNum=0;

    printf(HEADSTART "Running %d test from PhysDomain_test\n", NTEST);

    for(i=0; i<NTEST; i++){
        if(flag[i])
            failNum++;
    }

    if(failNum)
        printf(HEADFINISH "%d test faild from PhysDomain_test\n", failNum);
    else
        printf(HEADFINISH "%d test passed from PhysDomain_test\n", NTEST);

    return 0;
}