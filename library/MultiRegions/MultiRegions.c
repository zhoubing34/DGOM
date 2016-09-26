#include "MultiRegions.h"

/* private functions */
void SetVetxCoord2d(StdRegions2d *shape, int K, const double *VX, const double *VY, const int **EToV, double **GX, double **GY);
void LoadBalance2d(StdRegions2d *shape, int K, const int **EToV, const double **GX, const double **GY,
                 int *newK, int ***newEToV, double ***newx, double ***newy);
void SetFacePair2d(StdRegions2d *shape, int Klocal,
                 int **EToV, int **EToE, int **EToF, int **EToP,
                 int *Npar, int ***newParK, int ***newParF);
void SetNodeCoord2d(StdRegions2d *shape, int K, double **GX, double **GY, double **x, double **y);
void SetNodePair2d(StdRegions2d *shape, int K, double **GX, double **GY,
                 int **EToE, int **EToF, int **EToP, double **x, double **y,
                 int *Npar, int *Ntotalout, int **mapOUT,
                 int *vmapM, int *vmapP);
void SetVolumeGeo(StdRegions2d *shape, int K, double **x, double **y,
                  double *J, double *area, double *ciradius, real *vgeo);

void SetElementPair(StdRegions2d *shape, MultiReg2d *mesh, int *parEtotalout, int **mapOUT);
/**
 * @brief
 * Generation of two dimension mesh
 *
 * @param [in] shape StdRegions2d pointer: standard element
 * @param [in] K     number of element
 * @param [in] Nv    number of vertex
 * @param [in] EToV  element to vertex list
 * @param [in] VX    vertex coordinate
 * @param [in] VY    vertex coordinate
 *
 * @return
 * name     | type     | description of value
 * -------- |----------|----------------------
 * mesh     | MultiReg2d* | mesh object
 *
 * @note
 * 1. The index of vertex in EToV is start from 0;
 * 2. The input parameter `EToV`, `VX` and `VY` are copied to the corresponding fields of mesh object,
 * the user may deallocate them safely after calling this function.
 * @warning
 * This function will allocate and initialize a `MultiReg2d` type pointer, so the user should remember to
 * call `FreeMultiReg2d` function manually to free it in case of memory leak.
 */
MultiReg2d* GenMultiReg2d(StdRegions2d *shape, int K, int Nv,
                          const int **EToV, const double *VX, const double *VY){

    MultiReg2d *mesh = (MultiReg2d *)calloc(1, sizeof(MultiReg2d));

    MPI_Comm_rank(MPI_COMM_WORLD, &mesh->procid);
    MPI_Comm_size(MPI_COMM_WORLD, &mesh->nprocs);

    /* standard element */
    mesh->stdcell = shape;
    mesh->Nv = Nv;

    double **GX = BuildMatrix(K, shape->Nv);
    double **GY = BuildMatrix(K, shape->Nv);

    SetVetxCoord2d(shape, K, VX, VY, EToV, GX, GY);

    /* Redistribute the elements */
    LoadBalance2d(shape, K, EToV, GX, GY, &(mesh->K), &(mesh->EToV), &(mesh->GX), &(mesh->GY));

    DestroyMatrix(GX);
    DestroyMatrix(GY);

    /* Setup element connection, EToE,EToP,EToF & Npar,parK,parF */
    mesh->Npar = BuildIntVector(mesh->nprocs);

    mesh->EToE = BuildIntMatrix(mesh->K, shape->Nfaces);
    mesh->EToF = BuildIntMatrix(mesh->K, shape->Nfaces);
    mesh->EToP = BuildIntMatrix(mesh->K, shape->Nfaces);

    SetFacePair2d(shape, mesh->K, mesh->EToV, mesh->EToE, mesh->EToF, mesh->EToP,
                mesh->Npar, &mesh->parK, &mesh->parF);
    /* Setup nodes coordinate */
    mesh->x  = BuildMatrix(mesh->K, shape->Np);
    mesh->y  = BuildMatrix(mesh->K, shape->Np);
    SetNodeCoord2d(shape, mesh->K, mesh->GX, mesh->GY, mesh->x, mesh->y);

    /* Setup boundary nodes connection */
    mesh->vmapM = BuildIntVector(shape->Nfp*shape->Nfaces*mesh->K);
    mesh->vmapP = BuildIntVector(shape->Nfp*shape->Nfaces*mesh->K);
    SetNodePair2d(shape, mesh->K, mesh->GX, mesh->GY,
                mesh->EToE, mesh->EToF, mesh->EToP, mesh->x, mesh->y,
                mesh->Npar, &(mesh->parNtotalout), &(mesh->parmapOUT),
                mesh->vmapM, mesh->vmapP);

    SetElementPair(shape, mesh, &(mesh->parEtotalout), &(mesh->elemapOut));
    /* mesh geo */
    int Nfactor = 4;
    mesh->vgeo = (real*) calloc(Nfactor*mesh->K*shape->Np, sizeof(real));
    mesh->J    = BuildVector(mesh->K*shape->Np);
    mesh->area = BuildVector(mesh->K);
    mesh->ciradius = BuildVector(mesh->K);
    SetVolumeGeo(shape, mesh->K, mesh->x, mesh->y, mesh->J, mesh->area, mesh->ciradius, mesh->vgeo);

    return mesh;
};

void FreeMultiReg2d(MultiReg2d *mesh){
    /* mesh info */
    DestroyIntMatrix(mesh->EToV);

    /* coordinate */
    DestroyMatrix(mesh->GX);
    DestroyMatrix(mesh->GY);

    /* element connection */
    DestroyIntVector(mesh->Npar);
    DestroyIntMatrix(mesh->EToE);
    DestroyIntMatrix(mesh->EToF);
    DestroyIntMatrix(mesh->EToP);

    /* nodes */
    DestroyMatrix(mesh->x);
    DestroyMatrix(mesh->y);

    /* node connection */
    DestroyIntVector(mesh->vmapM);
    DestroyIntVector(mesh->vmapP);
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
        GeoFactor2d(shape->Np, x[k], y[k], shape->Dr, shape->Ds,
                    drdx, dsdx, drdy, dsdy, eJ);
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
void SetVetxCoord2d(StdRegions2d *shape, int K, const double *VX, const double *VY, const int **EToV, double **GX, double **GY){
    int n,k;
    /* vertex coordinate */
    for(k=0;k<K;k++){
        for(n=0;n<shape->Nv;n++){
            GX[k][n] = VX[EToV[k][n]];
            GY[k][n] = VY[EToV[k][n]];
        }
    }
}