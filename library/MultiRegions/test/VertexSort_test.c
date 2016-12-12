//
// Created by li12242 on 12/11/16.
//

#include "MulitRegions_test.h"
#include "MultiRegions/VertexSort.h"
#include "VertexSort_data.h"
#include "SetTestMultiRegions.h"
#include "LibUtilities/GenUniformMesh.h"


int MultiTriRegions_VertexSort_test(){
    // global variable
    UnstructMesh *grid = UniformTriMesh_create(2, 2, -1, 1, -1, 1, 1);
    MultiReg2d *mesh = SetTriParallelMultiRegions();
    StdRegions2d *shape = mesh->stdcell;

    // local
    int **exEToV = IntMatrix_create(mesh->K, shape->Nv);
    int fail = 0,k,i;
    clock_t clockT1, clockT2;

    // assignment
    for(k=0; k<mesh->K; k++)
        for(i=0; i<shape->Nv; i++){
            exEToV[k][i] = mesh->EToV[k][i];
        }

    for(k=0; k<mesh->K; k++){
        for(i=0; i<shape->Nv; i++){
            int t = (shape->Nv - i)%shape->Nv;
            mesh->EToV[k][i] = exEToV[k][t];
        }
    }

    // // check
    // char casename[32] = "MultiTriRegions_VertexSort_Test";
    // FILE *fp = CreateLog(casename, mesh->procid, mesh->nprocs);
    // PrintIntMatrix2File(fp, "EToV", mesh->EToV, mesh->K, shape->Nv);


    // call
    clockT1 = clock();
    for(k=0; k<mesh->K; k++){
        VertexSort(shape->Nv, grid->vx, grid->vy, mesh->EToV[k]);
    }
    clockT2 = clock();

    // PrintIntMatrix2File(fp, "newEToV", mesh->EToV, mesh->K, shape->Nv);
    // fclose(fp);

    // check
    fail = IntMatrix_test("VertexSort_Quad", mesh->EToV, exEToV, mesh->K, shape->Nv,
                          (double)((clockT2 - clockT1)/CLOCKS_PER_SEC));

    IntMatrix_free(exEToV);
    StdRegions2d_free(shape);
    MultiReg2d_free(mesh);
    UnstructMesh_free(grid);

    return fail;
}


int MultiQuadRegions_VertexSort_test(){
    // global variable
    UnstructMesh *grid = UniformQuadMesh_create(2, 2, -1, 1, -1, 1);
    MultiReg2d *mesh = SetQuadParallelMultiRegions();
    StdRegions2d *shape = mesh->stdcell;

    // local
    int **exEToV = IntMatrix_create(mesh->K, shape->Nv);
    int fail = 0,k,i;
    clock_t clockT1, clockT2;

    // assignment
    for(k=0; k<mesh->K; k++)
        for(i=0; i<shape->Nv; i++){
            exEToV[k][i] = mesh->EToV[k][i];
        }

    for(k=0; k<mesh->K; k++){
        for(i=0; i<shape->Nv; i++){
            int t = (shape->Nv - i)%shape->Nv;
            mesh->EToV[k][i] = exEToV[k][t];
        }
    }

    // call
    clockT1 = clock();
    for(k=0; k<mesh->K; k++){
        VertexSort(shape->Nv, grid->vx, grid->vy, mesh->EToV[k]);
    }
    clockT2 = clock();

    // check
    fail = IntMatrix_test("VertexSort_Quad", mesh->EToV, exEToV, mesh->K, shape->Nv,
                          (double)((clockT2 - clockT1)/CLOCKS_PER_SEC));

    IntMatrix_free(exEToV);
    StdRegions2d_free(shape);
    MultiReg2d_free(mesh);
    UnstructMesh_free(grid);

    return fail;
}