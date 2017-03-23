//
// Created by li12242 on 12/10/16.
//
#include "utility_test.h"
#include "nc_library_test.h"

int main(int argc, char **argv){

    int ishelp, isverbose;
    test_command(argc, argv, &ishelp, &isverbose);
    MPI_Init(&argc, &argv);

    int procid, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    // print help information

    if(ishelp && !procid){
        printf(HEADEND "DGOM\n" HEADLINE "Unit tests for Utility library\n"
                       HEADLINE "Usages:\n"
                       HEADLINE "   mpirun -n 2 -host localhost ./utility_test -verbose\n"
                       HEADLINE "\n"
                       HEADLINE "Optional features:\n"
                       HEADLINE "   -help     print help information\n"
                       HEADEND  "   -verbose  print variables to log files\n");
        return 0;
    }

    int i, flag[2];

    if(!procid) printf(HEADSTART "Running 2 test from LibUtilities_Test\n");

    flag[0] = MatrixInverse_test();
    flag[1] = MatrixMultiply_test();
    flag[2] = nc_file_read_from_file_test(isverbose);

    int failNum = 0;
    for(i=0; i<2; i++){
        if(flag[i]){ failNum++; }
    }
    if(failNum) {
        if(!procid) printf(HEADEND "%d test faild from LibUtilities_Test\n", failNum);
    } else {
        if(!procid) printf(HEADEND "2 test passed from LibUtilities_Test\n");
    }

    MPI_Finalize();
    return 0;
}