#include "MultiRegions.h"
#include "VertexSort.h"

/* private functions */
void SetVetxCoord2d(StdRegions2d *shape, int K, double *VX, double *VY, int **EToV, double **GX, double **GY);
void LoadBalance2d(StdRegions2d *shape, int K, int **EToV, double **GX, double **GY,
                 int *newK, int ***newEToV, double ***newx, double ***newy);
void SetFacePair2d(StdRegions2d *shape, int Klocal,
                 int **EToV, int **EToE, int **EToF, int **EToP,
                 int *Npar);//, int ***newParK, int ***newParF);
void SetNodeCoord2d(StdRegions2d *shape, int K, double **GX, double **GY, double **x, double **y);

void SetVolumeGeo(StdRegions2d *shape, int K, double **x, double **y,
                  double *J, double *area, double *ciradius, real *vgeo);

void Resort_Vertex(int K, int Nvert, int **EToV, double *VX, double *VY);

/**
 * @brief
 * Generation of two dimensional region.
 *
 * @param [in] shape standard element pointer.
 * @param [in] grid  Unstruct mesh pointer.
 *
 * @return mesh  MultiReg2d type pointer.
 *
 * @note
 * 1. The index of vertex in EToV is start from 0;
 * 2. The input parameter `EToV`, `VX` and `VY` are copied to the corresponding fields of mesh object,
 * the user may deallocate them safely after calling this function.
 * @warning
 * This function will allocate and initialize a `MultiReg2d` type pointer, so the user should remember to
 * call `MultiReg2d_free` function manually to free it in case of memory leak.
 */
MultiReg2d* MultiReg2d_create(StdRegions2d *shape, UnstructMesh *grid){

    /* assignment */
    int K      = grid->ne;
    int Nv     = grid->nv;
    int **EToV = grid->EToV;
    double *VX  = grid->vx;
    double *VY  = grid->vy;

    MultiReg2d *mesh = (MultiReg2d *)calloc(1, sizeof(MultiReg2d));

    MPI_Comm_rank(MPI_COMM_WORLD, &mesh->procid);
    MPI_Comm_size(MPI_COMM_WORLD, &mesh->nprocs);

    /* standard element */
    mesh->stdcell = shape;
    mesh->Nv = Nv;

    /* reorder vertex in each element */
    Resort_Vertex(K, shape->Nv, EToV, VX, VY);

    double **GX = Matrix_create(K, shape->Nv);
    double **GY = Matrix_create(K, shape->Nv);

    SetVetxCoord2d(shape, K, VX, VY, EToV, GX, GY);

    /* Redistribute the elements */
    LoadBalance2d(shape, K, EToV, GX, GY, &(mesh->K), &(mesh->EToV), &(mesh->GX), &(mesh->GY));

    Matrix_free(GX);
    Matrix_free(GY);

    /* Setup element connection, EToE,EToP,EToF & Npar,parK,parF */
    mesh->Npar = IntVector_create(mesh->nprocs);

    mesh->EToE = IntMatrix_create(mesh->K, shape->Nfaces);
    mesh->EToF = IntMatrix_create(mesh->K, shape->Nfaces);
    mesh->EToP = IntMatrix_create(mesh->K, shape->Nfaces);

    SetFacePair2d(shape,mesh->K,mesh->EToV,mesh->EToE,mesh->EToF,mesh->EToP,mesh->Npar);
    /* Setup nodes coordinate */
    mesh->x  = Matrix_create(mesh->K, shape->Np);
    mesh->y  = Matrix_create(mesh->K, shape->Np);
    SetNodeCoord2d(shape, mesh->K, mesh->GX, mesh->GY, mesh->x, mesh->y);

    /* mesh geo */
    int Nvgeo = 4;
    mesh->vgeo = (real*) calloc(Nvgeo*mesh->K*shape->Np, sizeof(real));
    mesh->J    = Vector_create(mesh->K * shape->Np);
    mesh->area = Vector_create(mesh->K);
    mesh->ciradius = Vector_create(mesh->K);
    SetVolumeGeo(shape, mesh->K, mesh->x, mesh->y, mesh->J, mesh->area, mesh->ciradius, mesh->vgeo);

    return mesh;
};

/**
 * @brief permute vertex anticlockwise in each element
 * @param [in] K number of element
 * @param [in] Nvert number of vertex in each element
 * @param [in,out] EToV element to vertex list
 * @param [in] VX vertex coordinate
 * @param [in] VY vertex coordiante
 *
 */
void Resort_Vertex(int K, int Nvert, int **EToV, double *VX, double *VY){
    int k;
    for(k=0; k<K; k++){
        VertexSort(Nvert, VX, VY, EToV[k]);
    }
}

void MultiReg2d_free(MultiReg2d *mesh){
    /* mesh info */
    IntMatrix_free(mesh->EToV);

    /* coordinate */
    Matrix_free(mesh->GX);
    Matrix_free(mesh->GY);

    /* element connection */
    IntVector_free(mesh->Npar);
    IntMatrix_free(mesh->EToE);
    IntMatrix_free(mesh->EToF);
    IntMatrix_free(mesh->EToP);

    /* nodes */
    Matrix_free(mesh->x);
    Matrix_free(mesh->y);
}

void SetVolumeGeo(StdRegions2d *shape, int K, double **x, double **y,
                  double *J, double *area, double *ciradius, real *vgeo){
    int k,n;
    double *drdx, *dsdx, *drdy, *dsdy, *eJ;

    int sz = shape->Np*sizeof(double);
    drdx = (double*) malloc(sz);
    drdy = (double*) malloc(sz);
    dsdx = (double*) malloc(sz);
    dsdy = (double*) malloc(sz);
    eJ   = (double*) malloc(sz);

    int sj = 0, sk = 0;
    for(k=0;k<K;++k){
        GeoFactor2d(shape->Np, x[k], y[k], shape->Dr, shape->Ds, drdx, dsdx, drdy, dsdy, eJ);
        for(n=0;n<shape->Np;n++){
            vgeo[sk++] = (real) drdx[n];
            vgeo[sk++] = (real) drdy[n];
            vgeo[sk++] = (real) dsdx[n];
            vgeo[sk++] = (real) dsdy[n];
            J[sj++] = eJ[n];

            area[k] += eJ[n]*shape->wv[n];
        }
        ciradius[k] = sqrt(area[k]/M_PI);
    }
}

/**
 * @brief
 * Set the vertex coordinate based on EToV
 *
 * @param [in]   shape StdRegions2d pointer: standard element object
 * @param [in]   K     number of elements
 * @param [in]   GX    input coordinate of vertex in each element
 * @param [in]   GY    input coordinate of vertex in each element
 *
 * @return
 * return values:
 * name     | type     | description of value
 * -------- |----------|----------------------
 * x | double[K][shape->Np]  | node coordinate
 * x | double[K][shape->Np]  | node coordinate
 *
 */
void SetNodeCoord2d(StdRegions2d *shape, int K, double **GX, double **GY, double **x, double **y){
    int k;
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

/**
 * @brief
 * Set the vertex coordinate based on EToV
 *
 * @param [in]   shape StdRegions2d pointer: standard element object
 * @param [in]   K     number of elements
 * @param [in]   VX    input coordinate of vertex
 * @param [in]   VY    input coordinate of vertex
 * @param [in]   EToV  element to vertex list
 *
 * @return
 * return values:
 * name     | type     | description of value
 * -------- |----------|----------------------
 * GX | double[K][shape->Nv]  | vertex coordinate
 * GX | double[K][shape->Nv]  | vertex coordinate
 *
 */
void SetVetxCoord2d(StdRegions2d *shape, int K, double *VX, double *VY, int **EToV, double **GX, double **GY){
    int n,k;
    /* vertex coordinate */
    for(k=0;k<K;k++){
        for(n=0;n<shape->Nv;n++){
            GX[k][n] = VX[EToV[k][n]];
            GY[k][n] = VY[EToV[k][n]];
        }
    }
}