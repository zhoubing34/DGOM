//
// Created by li12242 on 12/10/16.
//

#include "sc_test.h"
#include "sc_tri_test.h"
#define NTEST 5

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

    int fail = 0, err[NTEST];
    int N = 3;
    stdCell *tri = sc_create(N, TRIANGLE);
    stdCell *quad = sc_create(N, QUADRIL);

    printf(HEADSTART "Running %d test for SandCell test\n", NTEST);

    /* test case */
    err[0] = sc_triCoor_test(tri, isverbose);
    err[1] = sc_triVandMatrix_test(tri, isverbose);
    err[2] = sc_triMassMatrix_test(tri, isverbose);
    err[3] = sc_triDeriMatrix_test(tri, isverbose);
    err[4] = sc_triLIFT_test(tri, isverbose);

    sc_free(tri);
    sc_free(quad);
    return fail;
}