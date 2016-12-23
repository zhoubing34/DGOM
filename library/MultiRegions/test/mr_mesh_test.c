//
// Created by li12242 on 12/21/16.
//

#include <MultiRegions/mr_mesh.h>
#include "MultiRegions/mr_grid_uniformGrid.h"
#include "LibUtilities/UTest.h"
#include "MultiRegions/mr_mesh_addBoundary.h"

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

static void genBoundary(int *indicator, int **SFToV){

    extern int Mx, My;
    int Nfaces[4] = {Mx, Mx, My, My};
    int start_face_Id[4] = {0, Mx, Mx+Mx, Mx*2+My};
    int start_vert_Id[4] = {0, (Mx+1)*My, Mx, 0};
    int vert_stride[4] = {1, 1, Mx+1, Mx+1};
    // loop over 4 boundaries
    int b,f;
    for(b=0;b<4;b++){
        int ind = indicator[b];
        int startid = start_vert_Id[b];
        int vstride = vert_stride[b];
        int startf = start_face_Id[b];
        for(f=0;f<Nfaces[b];f++){
            SFToV[f+startf][0] = startid; startid+=vstride;
            SFToV[f+startf][1] = startid;
            SFToV[f+startf][2] = ind;
        }
    }

    return;
}


static int mr_triMesh_cellConnect_test(parallMesh *mesh, double dt, int verbose){
    int fail = 0;
    const int procid = mesh->procid;
    const int nprocs = mesh->nprocs;

    const int K = mesh->grid->K;
    const int Nfaces = mesh->cell->Nfaces;
    const int Nfp = mesh->cell->Nfp;

    /* add boundaries */
    extern int Mx,My;

    int ind[4] = {2,3,4,5};
    int Nsurf = 2*(Mx+My);
    int **SFToV = IntMatrix_create(Nsurf, 3);
    genBoundary(ind, SFToV);
    mr_mesh_addBoundary2d(mesh, 0, NULL);
    IntMatrix_free(SFToV);

    if(verbose){
        /* gen log filename */
        char casename[32] = "mr_triMesh_cellConnect_test";
        FILE *fp = CreateLog(casename, procid, nprocs);
        PrintIntMatrix2File(fp, "mesh->EToE", mesh->EToE, K, Nfaces);
        PrintIntMatrix2File(fp, "mesh->EToF", mesh->EToF, K, Nfaces);
        PrintIntMatrix2File(fp, "mesh->EToP", mesh->EToP, K, Nfaces);
        PrintIntVector2File(fp, "mesh->Npar", mesh->Npar, mesh->nprocs);
        fprintf(fp, "mesh->parallCellNum: %d\n", mesh->parallCellNum);
        PrintIntVector2File(fp, "mesh->cellIndexIn", mesh->cellIndexIn, mesh->parallCellNum);
        PrintIntVector2File(fp, "mesh->cellIndexOut", mesh->cellIndexOut, mesh->parallCellNum);
        PrintIntVector2File(fp, "mesh->vmapM", mesh->vmapM, K*Nfaces*Nfp);
        PrintIntVector2File(fp, "mesh->vmapP", mesh->vmapP, K*Nfaces*Nfp);
        fprintf(fp, "mesh->parallNodeNum: %d\n", mesh->parallNodeNum);
        PrintIntVector2File(fp, "mesh->nodeIndexOut", mesh->nodeIndexOut, mesh->parallNodeNum);

        fprintf(fp, "mesh->Nbc = %d\n", mesh->Nbc);
        PrintIntMatrix2File(fp, "mesh->EToBS", mesh->EToBS, K, Nfaces);
        PrintIntVector2File(fp, "mesh->bcIndList", mesh->bcIndList, mesh->Nbc);

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

    /* add boundaries */
    extern int Mx,My;

    int ind[4] = {2,3,4,5};
    int Nsurf = 2*(Mx+My);
    int **SFToV = IntMatrix_create(Nsurf, 3);
    genBoundary(ind, SFToV);
    mr_mesh_addBoundary2d(mesh, 0, NULL);
    IntMatrix_free(SFToV);


    if(verbose){
        /* gen log filename */
        char casename[32] = "mr_quadMesh_cellConnect_test";
        FILE *fp = CreateLog(casename, procid, nprocs);
        PrintIntMatrix2File(fp, "mesh->EToE", mesh->EToE, K, Nfaces);
        PrintIntMatrix2File(fp, "mesh->EToF", mesh->EToF, K, Nfaces);
        PrintIntMatrix2File(fp, "mesh->EToP", mesh->EToP, K, Nfaces);
        PrintIntVector2File(fp, "mesh->Npar", mesh->Npar, mesh->nprocs);
        fprintf(fp, "mesh->parallCellNum: %d\n", mesh->parallCellNum);
        PrintIntVector2File(fp, "mesh->cellIndexIn", mesh->cellIndexIn, mesh->parallCellNum);
        PrintIntVector2File(fp, "mesh->cellIndexOut", mesh->cellIndexOut, mesh->parallCellNum);
        PrintIntVector2File(fp, "mesh->vmapM", mesh->vmapM, K*Nfaces*Nfp);
        PrintIntVector2File(fp, "mesh->vmapP", mesh->vmapP, K*Nfaces*Nfp);
        fprintf(fp, "mesh->parallNodeNum: %d\n", mesh->parallNodeNum);
        PrintIntVector2File(fp, "mesh->nodeIndexOut", mesh->nodeIndexOut, mesh->parallNodeNum);

        fprintf(fp, "mesh->Nbc = %d\n", mesh->Nbc);
        PrintIntMatrix2File(fp, "mesh->EToBS", mesh->EToBS, K, Nfaces);
        PrintIntVector2File(fp, "mesh->bcIndList", mesh->bcIndList, mesh->Nbc);
        fclose(fp);
    }

    return fail;
}
