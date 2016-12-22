//
// Created by li12242 on 12/11/16.
//

#include "phys_test.h"
#include "phys_nodeFetch_test.h"
#include "MultiRegions/mr_grid_uniformGrid.h"
#include "MultiRegions/mr_mesh_addBoundary.h"

#define TESTNUM 2

int Mx = 4;
int My = 2;
int N = 1;
int Nfield = 2;

int main(int argc, char **argv){

    // get settings
    int ishelp, isverbose;
    UTest_Command(argc, argv, &ishelp, &isverbose);

    // print help information
    if(ishelp){
        printf(HEADFINISH "DGOM\n" HEADLINE "Unit tests for PhysField library\n"
                       HEADLINE "Example usages:\n"
                       HEADLINE "   mpirun -n 2 -host localhost ./MulitRegions_Test\n"
                       HEADLINE "\n"
                       HEADLINE "Optional features:\n"
                       HEADLINE "   -help     print help information\n"
                       HEADFINISH "   -verbose  print variables to log files\n\n");
        return 0;
    }

    // local vairable
    int i,err[TESTNUM];
    int failNum=0;

    MPI_Init(&argc, &argv);

    int procid, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    if(!procid)
        printf(HEADSTART "Running %d test from PhysDomain_test\n", TESTNUM);


    extern int N, Nfield, Mx, My;

    stdCell *shape = sc_create(N, TRIANGLE);
    geoGrid *grid = mr_grid_createUniformGrid_tri(shape, Mx, My, -1, 1, -1, 1, 1);
    multiReg *region = mr_reg_create(grid);
    parallMesh *mesh = mr_mesh_create(region);
    mr_mesh_addBoundary2d(mesh, 0, NULL);
    physField *phys = phys_create(Nfield, mesh);

    err[0] = phys_nodeFetch_test(phys, isverbose, "phys_trinodeFetch_test", "phys_tri_nodeFetch_test");

    /* free memory */
    sc_free(shape);
    mr_grid_free(grid);
    mr_reg_free(region);
    mr_mesh_free(mesh);
    PhysDomain2d_free(phys);


    shape = sc_create(N, QUADRIL);
    grid = mr_grid_createUniformGrid_quad(shape, Mx, My, -1, 1, -1, 1);
    region = mr_reg_create(grid);
    mesh = mr_mesh_create(region);
    mr_mesh_addBoundary2d(mesh, 0, NULL);
    phys = phys_create(Nfield, mesh);

    err[1] = phys_nodeFetch_test(phys, isverbose, "phys_quadnodeFetch_test", "phys_quad_nodeFetch_test");

    for(i=0; i<TESTNUM; i++){
        if(err[i])
            failNum++;
    }

    if(!procid){
        if(failNum)
            printf(HEADFINISH "%d test faild from MulitRegions_test\n", failNum);
        else
            printf(HEADFINISH "%d test passed from MulitRegions_test\n", TESTNUM);
    }

    MPI_Finalize();

    return 0;
}