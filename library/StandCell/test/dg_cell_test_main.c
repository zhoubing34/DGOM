//
// Created by li12242 on 12/10/16.
//

#include "dg_cell_test_main.h"
#include "triangle/dg_cell_tri_test.h"
#include "quadrilateral/dg_cell_quad_test.h"
#include "point/dg_cell_point_test.h"
#include "line/dg_cell_line_test.h"

#define NTEST 8

int main(int argc, char **argv){

    // get settings
    int ishelp, isverbose;
    test_command(argc, argv, &ishelp, &isverbose);

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
    printf(HEADSTART "Running %d test for dg_cell_point test, verbose=%d\n",
           NTEST, isverbose);
    dg_cell *point = dg_cell_creat(N, POINT);
    i = 0;
    err[i++] = dg_point_info_test(point, isverbose);
    err[i++] = dg_point_node_test(point, isverbose);
    err[i++] = dg_point_vand_matrix_test(point, isverbose);
    err[i++] = dg_point_mass_matrix_test(point, isverbose);
    err[i++] = dg_point_deri_matrix_test(point, isverbose);
    err[i++] = dg_point_Fmask_test(point, isverbose);
    err[i++] = dg_point_LIFT_test(point, isverbose);
    err[i++] = dg_point_vert_proj_test(point, isverbose);

    dg_cell_free(point);
    for(i=0;i<NTEST;i++) {fail += err[i];}
    if(fail){
        printf(HEADEND "%d test faild from dg_cell_point test\n", fail);
    } else {
        printf(HEADEND "%d test passed from dg_cell_point test\n", NTEST);
    }

    printf(HEADSTART "Running %d test for dg_cell_line test, verbose=%d\n",
           NTEST, isverbose);
    dg_cell *line = dg_cell_creat(N, LINE);
    i = 0;
    err[i++] = dg_line_info_test(line, isverbose);
    err[i++] = dg_line_node_test(line, isverbose);
    err[i++] = dg_line_vand_matrix_test(line, isverbose);
    err[i++] = dg_line_mass_matrix_test(line, isverbose);
    err[i++] = dg_line_deri_matrix_test(line, isverbose);
    err[i++] = dg_line_Fmask_test(line, isverbose);
    err[i++] = dg_line_LIFT_test(line, isverbose);
    err[i++] = dg_line_vert_proj_test(line, isverbose);

    dg_cell_free(line);
    for(i=0;i<NTEST;i++) {fail += err[i];}
    if(fail){
        printf(HEADEND "%d test faild from dg_cell_line test\n", fail);
    } else {
        printf(HEADEND "%d test passed from dg_cell_line test\n", NTEST);
    }

    printf(HEADSTART "Running %d test for dg_cell_tri test, verbose=%d\n",
           NTEST, isverbose);
    dg_cell *tri = dg_cell_creat(N, TRIANGLE);
    /* triangle test case */
    i=0;
    err[i++] = dg_tri_info_test(tri, isverbose);
    err[i++] = dg_tri_nood_test(tri, isverbose);
    err[i++] = dg_tri_vand_matrix_test(tri, isverbose);
    err[i++] = dg_tri_mass_matrix_test(tri, isverbose);
    err[i++] = dg_tri_deri_matrix_test(tri, isverbose);
    err[i++] = dg_tri_Fmask_test(tri, isverbose);
    err[i++] = dg_tri_LIFT_test(tri, isverbose);
    err[i++] = dg_tri_vert_proj_test(tri, isverbose);
    dg_cell_free(tri);

    for(i=0;i<NTEST;i++) {fail += err[i];}
    if(fail){
        printf(HEADEND "%d test faild from dg_cell_tri test\n", fail);
    } else {
        printf(HEADEND "%d test passed from dg_cell_tri test\n", NTEST);
    }

    /* quadrilateral test case */
    dg_cell *quad = dg_cell_creat(N, QUADRIL);
    printf(HEADSTART "Running %d test for dg_cell_quad test, verbose=%d\n",
           NTEST, isverbose);
    i = 0;
    err[i++] = dg_quad_info_test(quad, isverbose);
    err[i++] = dg_quad_nood_test(quad, isverbose);
    err[i++] = dg_quad_vand_matrix_test(quad, isverbose);
    err[i++] = dg_quad_mass_matrix_test(quad, isverbose);
    err[i++] = dg_quad_deri_matrix_test(quad, isverbose);
    err[i++] = dg_quad_Fmask_test(quad, isverbose);
    err[i++] = dg_quad_LIFT_test(quad, isverbose);
    err[i++] = dg_quad_vert_proj_test(quad, isverbose);
    dg_cell_free(quad);

    for(i=0;i<NTEST;i++) {fail += err[i];}
    if(fail) {printf(HEADEND "%d test faild from dg_cell_quad test\n", fail);}
    else {printf(HEADEND "%d test passed from dg_cell_quad test\n", NTEST);}

    return fail;
}