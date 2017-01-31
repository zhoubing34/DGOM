//
// Created by li12242 on 12/11/16.
//


#include "MultiRegions/mr_grid_uniformGrid.h"
#include "MultiRegions/mr_mesh_addBoundary.h"
#include "pf_test.h"
#include "pf_nodeFetch_test.h"
#include "pf_strong_volume_flux2d_test.h"
#include "pf_strong_surface_flux2d_test.h"
#include "pf_strong_viscosity_LDG_flux2d_test.h"
#include "pf_cellMean_test.h"
#include "pf_cellFetch_test.h"
#include "pf_limiter_test.h"

#define TESTNUM 14

int Mx = 4;
int My = 2;
int N = 3;
int Nfield = 2;

int main(int argc, char **argv){

    // get settings
    int ishelp, isverbose;
    UTest_Command(argc, argv, &ishelp, &isverbose);

    // print help information
    if(ishelp){
        printf(HEADEND "DGOM\n" HEADLINE "Unit tests for PhysField library\n"
                       HEADLINE "Example usages:\n"
                       HEADLINE "   mpirun -n 2 -host localhost ./MulitRegions_Test\n"
                       HEADLINE "\n"
                       HEADLINE "Optional features:\n"
                       HEADLINE "   -help     print help information\n"
                       HEADEND "   -verbose  print variables to log files\n\n");
        return 0;
    }

    // local vairable
    int i=0,err[TESTNUM];
    int failNum=0;

    MPI_Init(&argc, &argv);

    int procid, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    if(!procid)
        printf(HEADSTART "Running %d test from PhysDomain_test\n", TESTNUM);


    extern int N, Nfield, Mx, My;
//    int **SFToV = IntMatrix_create(1, 3);
//    SFToV[0][0] = 0;
//    SFToV[0][1] = 1;
//    SFToV[0][2] = 0;

    stdCell *tri = sc_create(N, TRIANGLE);
    geoGrid *tri_grid = mr_grid_createUniformGrid_tri(tri, Mx, My, -1, 1, -1, 1, 1);
//    printf("procid=%d, finish set up tri-grid\n",procid);
    multiReg *tri_region = mr_reg_create(tri_grid);
    parallMesh *tri_mesh = mr_mesh_create(tri_region);
    mr_mesh_addBoundary2d(tri_mesh, 0, NULL);
    physField *tri_phys = pf_create(Nfield, tri_mesh);
//    printf("procid=%d, finish set up tri-physField\n", procid);

    err[i++] = phys_cellMean_test(tri_phys, isverbose,
                                  "phys_tri_cellMean_test",
                                  "phys_tri_cellMean_test");
    err[i++] = phys_nodeFetch_test(tri_phys, isverbose,
                                   "phys_tri_nodeFetch_test",
                                   "phys_tri_nodeFetch_test");
    err[i++] = phys_cellFetch_test(tri_phys, isverbose,
                                   "phys_tri_cellFetch_test",
                                   "phys_tri_cellFetch_test");
    err[i++] = phys_strong_volume_flux2d_test(tri_phys, isverbose,
                                              "phys_tri_strong_volume_flux2d_test",
                                              "phys_tri_strong_volume_flux2d_test");
    err[i++] = phys_strong_surface_flux2d_test(tri_phys, isverbose,
                                               "phys_tri_strong_surface_flux2d_test",
                                               "phys_tri_strong_surface_flux2d_test");

    err[i++] = phys_strong_viscosity_LDG_flux2d_test(tri_phys, isverbose,
                                                     "phys_tri_strong_viscosity_LDG_flux2d_test",
                                                     "phys_tri_strong_viscosity_LDG_flux2d_test");
    err[i++] = phys_limiter_test(tri_phys, isverbose,
                                 "phys_tri_limiter_test",
                                 "phys_tri_limiter_test");

    /* free memory */
    sc_free(tri);
    mr_grid_free(tri_grid);
    mr_reg_free(tri_region);
    mr_mesh_free(tri_mesh);
    pf_free(tri_phys);
//    printf("procid=%d, free tri-physField\n", procid);


    stdCell *quad = sc_create(N, QUADRIL);
    geoGrid *quad_grid = mr_grid_createUniformGrid_quad(quad, Mx, My, -1, 1, -1, 1);
    multiReg *quad_region = mr_reg_create(quad_grid);
    parallMesh *quad_mesh = mr_mesh_create(quad_region);
    mr_mesh_addBoundary2d(quad_mesh, 0, NULL);
    physField *quad_phys = pf_create(Nfield, quad_mesh);

    err[i++] = phys_cellMean_test(quad_phys, isverbose,
                                  "phys_quad_cellMean_test",
                                  "phys_quad_cellMean_test");
    err[i++] = phys_nodeFetch_test(quad_phys, isverbose,
                                   "phys_quad_nodeFetch_test",
                                   "phys_quad_nodeFetch_test");
    err[i++] = phys_cellFetch_test(quad_phys, isverbose,
                                   "phys_quad_cellFetch_test",
                                   "phys_quad_cellFetch_test");
    err[i++] = phys_strong_volume_flux2d_test(quad_phys, isverbose,
                                              "phys_quad_strong_volume_flux2d_test",
                                              "phys_quad_strong_volume_flux2d_test");
    err[i++] = phys_strong_surface_flux2d_test(quad_phys, isverbose,
                                               "phys_quad_strong_surface_flux2d_test",
                                               "phys_quad_strong_surface_flux2d_test");
    err[i++] = phys_strong_viscosity_LDG_flux2d_test(quad_phys, isverbose,
                                                     "phys_quad_strong_viscosity_LDG_flux2d_test",
                                                     "phys_quad_strong_viscosity_LDG_flux2d_test");
    err[i++] = phys_limiter_test(quad_phys, isverbose,
                                 "phys_quad_limiter_test",
                                 "phys_quad_limiter_test");
    /* free memory */
    sc_free(quad);
    mr_grid_free(quad_grid);
    mr_reg_free(quad_region);
    mr_mesh_free(quad_mesh);
    pf_free(quad_phys);

    for(i=0; i<TESTNUM; i++){
        if(err[i])
            failNum++;
    }

    if(!procid){
        if(failNum)
            printf(HEADEND "%d test faild from PhysField_test\n", failNum);
        else
            printf(HEADEND "%d test passed from PhysField_test\n", TESTNUM);
    }

    MPI_Finalize();

    return 0;
}