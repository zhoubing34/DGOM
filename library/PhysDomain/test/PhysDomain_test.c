//
// Created by li12242 on 12/11/16.
//

#include "PhysDomain_test.h"
#include "MultiRegions/test/SetTestMultiRegions.h"
#include "NodeFetch_test.h"

#define TESTNUM 2

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
        printf(HEADSTART "Running %d test from PhysDomain_test\n", TESTNUM);

    int Nfields = 2;
    MultiRegBC2d *triSurf = setTriTestMesh();
    MultiRegBC2d *quadSurf = setQuadTestMesh();
    MultiReg2d *triMesh = triSurf->mesh;
    MultiReg2d *quadMesh = quadSurf->mesh;
    PhysDomain2d *triPhys = PhysDomain2d_create(triMesh, triSurf, Nfields);
    PhysDomain2d *quadPhys = PhysDomain2d_create(quadMesh, quadSurf, Nfields);

    flag[0] = nodeFetch_test(quadPhys, isverbose, "QuadNodeFetchTest", "PhysDomain_quad_NodeFetch_test");
    flag[1] = nodeFetch_test(triPhys, isverbose, "TriNodeFetchTest", "PhysDomain_tri_NodeFetch_test");

    for(i=0; i<TESTNUM; i++){
        if(flag[i])
            failNum++;
    }

    StdRegions2d_free(triMesh->stdcell);
    StdRegions2d_free(quadMesh->stdcell);
    MultiReg2d_free(triMesh);
    MultiReg2d_free(quadMesh);
    MultiRegBC2d_free(triSurf);
    MultiRegBC2d_free(quadSurf);
    PhysDomain2d_free(triPhys);
    PhysDomain2d_free(quadPhys);

    if(!procid){
        if(failNum)
            printf(HEADFINISH "%d test faild from MulitRegions_test\n", failNum);
        else
            printf(HEADFINISH "%d test passed from MulitRegions_test\n", TESTNUM);
    }

    MPI_Finalize();

    return 0;
}