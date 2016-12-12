//
// Created by li12242 on 12/11/16.
//

#include "MulitRegions_test.h"
#include "NodePair_test.h"
#include "FacePair_test.h"
#include "VertexSort_data.h"

#define TESTNUM 6

int VertexSort_test();

int main(int argc, char **argv){

    int i, flag[TESTNUM];
    int failNum=0;

    printf(HEADSTART "Running %d test from MulitRegions_test\n", TESTNUM);
    MPI_Init(&argc, &argv);

    flag[0] = MultiTriRegions_VertexSort_test();
    flag[1] = MultiQuadRegions_VertexSort_test();
    flag[2] = MultiTriRegions_NodePair_Test();
    flag[3] = MultiQuadRegions_NodePair_Test();
    flag[4] = MultiTriRegions_FacePair_test();
    flag[5] = MultiQuadRegions_FacePair_test();
    MPI_Finalize();
    for(i=0; i<TESTNUM; i++){
        if(flag[i])
            failNum++;
    }
    if(failNum)
        printf(HEADFINISH "%d test faild from MulitRegions_test\n", failNum);
    else
        printf(HEADFINISH "%d test passed from MulitRegions_test\n", TESTNUM);

    return 0;
}