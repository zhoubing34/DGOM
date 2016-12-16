//
// Created by li12242 on 12/11/16.
//

#include <MultiRegions/MultiRegions.h>
#include <MultiRegions/MultiRegBC/MultiRegBC2d.h>
#include "MulitRegions_test.h"
#include "VarMapPair_test.h"
#include "CellPair_test.h"
#include "VertexSort_test.h"
#include "MultiRegBC2d_test.h"

#define TESTNUM 8

int main(int argc, char **argv){

    // get settings
    int ishelp, isverbose;
    UTest_Command(argc, argv, &ishelp, &isverbose);

    // print help information
    if(ishelp){
        printf(HEADFINISH "DGOM\n" HEADLINE "Unit tests for MultiRegions library\n"
                       HEADLINE "Example usages:\n"
                       HEADLINE "   mpirun -n 2 -host localhost ./MulitRegions_Test\n"
                       HEADLINE "\n"
                       HEADLINE "Optional features:\n"
                       HEADLINE "   -help     print help information\n"
                       HEADFINISH "   -verbose  print variables to log files\n\n");
        return 0;
    }

    // local vairable
    int i,flag[TESTNUM];
    int failNum=0;

    MPI_Init(&argc, &argv);

    int procid, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    if(!procid)
        printf(HEADSTART "Running %d test from MulitRegions_test\n", TESTNUM);

    MultiRegBC2d *triSurf = setTriTestMesh();
    MultiRegBC2d *quadSurf = setQuadTestMesh();
    MultiReg2d *triMesh = triSurf->mesh;
    MultiReg2d *quadMesh = quadSurf->mesh;

    flag[0] = MultiTriRegions_VertexSort_test(triMesh, isverbose);
    flag[1] = MultiQuadRegions_VertexSort_test(quadMesh, isverbose);
    flag[2] = MultiTriRegions_VarMapPair_Test(triSurf, isverbose);
    flag[3] = MultiQuadRegions_VarMapPair_Test(quadSurf, isverbose);
    flag[4] = MultiTriRegions_CellPair_test(triSurf, isverbose);
    flag[5] = MultiQuadRegions_CellPair_test(quadSurf, isverbose);
    flag[6] = MultiTriRegions_MultiRegBC2d_test(triSurf, isverbose);
    flag[7] = MultiQuadRegions_MultiRegBC2d_test(quadSurf, isverbose);


    StdRegions2d_free(triMesh->stdcell);
    StdRegions2d_free(quadMesh->stdcell);
    MultiReg2d_free(triMesh);
    MultiReg2d_free(quadMesh);
    MultiRegBC2d_free(triSurf);
    MultiRegBC2d_free(quadSurf);

    MPI_Finalize();
    for(i=0; i<TESTNUM; i++){
        if(flag[i])
            failNum++;
    }

    if(!procid){
        if(failNum)
            printf(HEADFINISH "%d test faild from MulitRegions_test\n", failNum);
        else
            printf(HEADFINISH "%d test passed from MulitRegions_test\n", TESTNUM);
    }

    return 0;
}