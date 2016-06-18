//
// Created by li12242 on 16/6/18.
//

#include "MultiRegions.h"

/* private functions */
void SetCoord(StdRegions2d *shape, int K, double *VX, double *VY, int **EToV,
              double **GX, double **GY, double **x, double **y);

/**
 * @brief
 * Generation of two dimension mesh
 *
 * @param [StdRegions2d*] shape standard element
 * @param [int]           k     number of element
 * @param [int]           Nv    number of vertex
 * @param [int**]         EToV  element to vertex list
 * @param [double*]       VX    vertex coordinate
 * @param [double*]       VY    vertex coordinate
 *
 * @return
 * name     | type     | description of value
 * -------- |----------|----------------------
 * mesh     | MultiReg2d* |
 *
 * @note
 * The vertex index in EToV is start from 0.
 */
MultiReg2d* GenMultiReg2d(StdRegions2d *shape, int K, int Nv,
                          int **EToV, double *VX, double *VY){

    MultiReg2d *mesh = (MultiReg2d *)calloc(1, sizeof(MultiReg2d));

    MPI_Comm_rank(MPI_COMM_WORLD, &mesh->procid);
    MPI_Comm_size(MPI_COMM_WORLD, &mesh->nprocs);

    printf("procid = %d, nprocs = %d\n", mesh->procid, mesh->nprocs);

    /* element info */
    mesh->K  = K;
    mesh->Nv = Nv;

    /* standard element */
    mesh->StdElement = shape;

    mesh->GX = BuildMatrix(K, shape->Nv);
    mesh->GY = BuildMatrix(K, shape->Nv);

    mesh->x  = BuildMatrix(K, shape->Np);
    mesh->y  = BuildMatrix(K, shape->Np);

    SetCoord(shape, K, VX, VY, EToV, mesh->GX, mesh->GY, mesh->x, mesh->y);


    return mesh;
};

void SetCoord(StdRegions2d *shape, int K, double *VX, double *VY, int **EToV,
              double **GX, double **GY, double **x, double **y){
    int n,k;
    /* vertex coordinate */
    for(k=0;k<K;k++){
        for(n=0;n<shape->Nv;n++){
            GX[k][n] = VX[EToV[k][n]];
            GY[k][n] = VY[EToV[k][n]];
        }
    }

    /* node coordinate */
    for(k=0;k<K;k++) {
        if (shape->Nv == 3) {
            MapTriCoor(shape, GX[k], GY[k], x[k], y[k]);
        } else if (shape->Nv == 4) {
            MapQuadCoor(shape, GX[k], GY[k], x[k], y[k]);
        } else {
            printf("fatal error: wrong number of vertex %d in StdRegions2d", shape->Nv);
        }
    }
}


void FreeMultiReg2d(MultiReg2d *mesh){}

void LoadBalance(MultiReg2d *mesh){}