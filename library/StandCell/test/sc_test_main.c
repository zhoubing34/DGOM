//
// Created by li12242 on 12/10/16.
//

#include "sc_test_main.h"
#include "StandCell/test/triangle/sc_tri_test.h"
#include "StandCell/test/quadrilateral/sc_quad_test.h"

#define NTEST 6

int main(int argc, char **argv){

    // get settings
    int ishelp, isverbose;
    UTest_Command(argc, argv, &ishelp, &isverbose);

    // print help information
    if(ishelp){
        printf(HEADEND "DGOM\n" HEADLINE "Unit tests for MultiRegions library\n"
                       HEADLINE "Example usages:\n"
                       HEADLINE "   ./StandCell_Test\n"
                       HEADLINE "\n"
                       HEADLINE "Optional features:\n"
                       HEADLINE "   -help     print help information\n"
                       HEADEND "   -verbose  print variables to log files\n\n");
        return 0;
    }

    int fail = 0, err[NTEST];
    int N = 3,i;
    printf(HEADSTART "Running %d test for SandCell_triangle test\n", NTEST);
    stdCell *tri = sc_create(N, TRIANGLE);
    /* triangle test case */
    i=0;
    err[i++] = sc_triCoor_test(tri, isverbose);
    err[i++] = sc_triVandMatrix_test(tri, isverbose);
    err[i++] = sc_triMassMatrix_test(tri, isverbose);
    err[i++] = sc_triDeriMatrix_test(tri, isverbose);
    err[i++] = sc_triLIFT_test(tri, isverbose);
    err[i++] = sc_triVertProj_test(tri, isverbose);


    sc_free(tri);

    for(i=0;i<NTEST;i++)
        fail += err[i];

    if(fail)
        printf(HEADEND "%d test faild from SandCell_triangle test\n", fail);
    else
        printf(HEADEND "%d test passed from SandCell_triangle test\n", NTEST);

    /* quadrilateral test case */
    stdCell *quad = sc_create(N, QUADRIL);
    printf(HEADSTART "Running %d test for SandCell_quadrilateral test\n", NTEST);
    i = 0;
    err[i++] = sc_quadCoor_test(quad, isverbose);
    err[i++] = sc_quadVandMatrix_test(quad, isverbose);
    err[i++] = sc_quadMassMatrix_test(quad, isverbose);
    err[i++] = sc_quadDeriMatrix_test(quad, isverbose);
    err[i++] = sc_quadLIFT_test(quad, isverbose);
    err[i++] = sc_quadVertProj_test(quad, isverbose);
    sc_free(quad);

    for(i=0;i<NTEST;i++)
        fail += err[i];

    if(fail)
        printf(HEADEND "%d test faild from SandCell_quadrilateral test\n", fail);
    else
        printf(HEADEND "%d test passed from SandCell_quadrilateral test\n", NTEST);

    return fail;
}