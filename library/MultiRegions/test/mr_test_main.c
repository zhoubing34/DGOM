//
// Created by li12242 on 12/19/16.
//

#include "LibUtilities/UTest.h"
//#include "LibUtilities/LibUtilities.h"
#include "MultiRegions/mr_grid.h"
#include "mr_grid_test.h"
#include "mr_reg_test.h"
#include "mr_mesh_test.h"

#define NTEST 3

int Mx = 4;
int My = 2;

int main(int argc, char **argv){

    // get settings
    int ishelp, isverbose;
    UTest_Command(argc, argv, &ishelp, &isverbose);

    MPI_Init(&argc, &argv);
    int procid, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    // print help information
    if(ishelp && !procid){
        printf(HEADEND "DGOM\n" HEADLINE "Unit tests for MultiRegions library\n"
                       HEADLINE "Usages:\n"
                       HEADLINE "   mpirun -n 2 -host localhost ./MulitRegions_Test\n"
                       HEADLINE "\n"
                       HEADLINE "Optional features:\n"
                       HEADLINE "   -help     print help information\n"
                       HEADEND "   -verbose  print variables to log files\n");
        return 0;
    }

    if(!procid)
        printf(HEADSTART "Running %d test from MultiRegions\n", NTEST);

    // local variable
    int failNum = 0, err[NTEST], i=0;
    // test
    err[i++] = mr_grid_test(isverbose);
    err[i++] = mr_reg_test(isverbose);
    err[i++] = mr_mesh_test(isverbose);


    MPI_Finalize();
    /* summary the fail test number */
    for(i=0; i<NTEST; i++){
        if(err[i])
            failNum++;
    }

    if(!procid){
        if(failNum)
            printf(HEADEND "%d test faild from MulitRegions_test\n", failNum);
        else
            printf(HEADEND "%d test passed from MulitRegions_test\n", NTEST);
    }
    return failNum;
}

