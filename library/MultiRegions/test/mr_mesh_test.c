//
// Created by li12242 on 12/21/16.
//

#include <MultiRegions/mr_mesh.h>
#include "MultiRegions/mr_grid_uniformGrid.h"
#include "LibUtilities/UTest.h"


static int mr_triMesh_cellConnect_test(parallMesh *mesh, double dt, int verbose);
static int mr_quadMesh_cellConnect_test(parallMesh *mesh, double dt, int verbose);

int mr_mesh_test(int verbose){
    int fail = 0;

    int N=2;

    extern int Mx,My;

    double clockT1, clockT2;

    /* triangle regions */
    stdCell *shape = sc_create(N, TRIANGLE);
    geoGrid *grid = mr_grid_createUniformGrid_tri(shape, Mx, My, -1, 1, -1, 1, 1);
    multiReg *region = mr_reg_create(grid);
    clockT1 = MPI_Wtime();
    parallMesh *mesh = mr_mesh_create(region);
    clockT2 = MPI_Wtime();

    fail = mr_triMesh_cellConnect_test(mesh, clockT2-clockT1, verbose);

    /* free memory */
    sc_free(shape);
    mr_grid_free(grid);
    mr_reg_free(region);

    /* quadrilateral regions */
    shape = sc_create(N, QUADRIL);
    grid = mr_grid_createUniformGrid_quad(shape, Mx, My, -1, 1, -1, 1);
    region = mr_reg_create(grid);
    clockT1 = MPI_Wtime();
    mesh = mr_mesh_create(region);
    clockT2 = MPI_Wtime();

    fail = mr_quadMesh_cellConnect_test(mesh, clockT2-clockT1, verbose);

    /* free memory */
    sc_free(shape);
    mr_grid_free(grid);
    mr_reg_free(region);

    return fail;
}


static int mr_triMesh_cellConnect_test(parallMesh *mesh, double dt, int verbose){
    int fail = 0;
    const int procid = mesh->procid;
    const int nprocs = mesh->nprocs;

    const int K = mesh->grid->K;
    const int Nfaces = mesh->cell->Nfaces;
    const int Nfp = mesh->cell->Nfp;

    if(verbose){
        /* gen log filename */
        char casename[32] = "mr_triMesh_cellConnect_test";
        FILE *fp = CreateLog(casename, procid, nprocs);
        PrintIntMatrix2File(fp, "mesh->EToE", mesh->EToE, K, Nfaces);
        PrintIntMatrix2File(fp, "mesh->EToF", mesh->EToF, K, Nfaces);
        PrintIntMatrix2File(fp, "mesh->EToP", mesh->EToP, K, Nfaces);
        fprintf(fp, "mesh->parallCellNum: %d\n", mesh->parallCellNum);
        PrintIntVector2File(fp, "mesh->cellIndexIn", mesh->cellIndexIn, mesh->parallCellNum);
        PrintIntVector2File(fp, "mesh->cellIndexOut", mesh->cellIndexOut, mesh->parallCellNum);
        PrintIntVector2File(fp, "mesh->vmapM", mesh->vmapM, K*Nfaces*Nfp);
        PrintIntVector2File(fp, "mesh->vmapP", mesh->vmapP, K*Nfaces*Nfp);
        fprintf(fp, "mesh->parallNodeNum: %d\n", mesh->parallNodeNum);
        PrintIntVector2File(fp, "mesh->nodeIndexOut", mesh->nodeIndexOut, mesh->parallNodeNum);

        fclose(fp);
    }

    return fail;
}


static int mr_quadMesh_cellConnect_test(parallMesh *mesh, double dt, int verbose){
    int fail = 0;
    const int procid = mesh->procid;
    const int nprocs = mesh->nprocs;

    const int K = mesh->grid->K;
    const int Nfaces = mesh->cell->Nfaces;
    const int Nfp = mesh->cell->Nfp;

    if(verbose){
        /* gen log filename */
        char casename[32] = "mr_quadMesh_cellConnect_test";
        FILE *fp = CreateLog(casename, procid, nprocs);
        PrintIntMatrix2File(fp, "mesh->EToE", mesh->EToE, K, Nfaces);
        PrintIntMatrix2File(fp, "mesh->EToF", mesh->EToF, K, Nfaces);
        PrintIntMatrix2File(fp, "mesh->EToP", mesh->EToP, K, Nfaces);
        fprintf(fp, "mesh->parallCellNum: %d\n", mesh->parallCellNum);
        PrintIntVector2File(fp, "mesh->cellIndexIn", mesh->cellIndexIn, mesh->parallCellNum);
        PrintIntVector2File(fp, "mesh->cellIndexOut", mesh->cellIndexOut, mesh->parallCellNum);
        PrintIntVector2File(fp, "mesh->vmapM", mesh->vmapM, K*Nfaces*Nfp);
        PrintIntVector2File(fp, "mesh->vmapP", mesh->vmapP, K*Nfaces*Nfp);
        fprintf(fp, "mesh->parallNodeNum: %d\n", mesh->parallNodeNum);
        PrintIntVector2File(fp, "mesh->nodeIndexOut", mesh->nodeIndexOut, mesh->parallNodeNum);
        fclose(fp);
    }

    return fail;
}
