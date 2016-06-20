#include "MultiRegions.h"

/* private functions */
void SetVetxCoord2d(StdRegions2d *shape, int K, double *VX, double *VY, int **EToV, double **GX, double **GY);
void LoadBalance2d(StdRegions2d *shape, int K, int **EToV, double **GX, double **GY,
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
                  double *J, double *area, double *ciradius, float *vgeo);

/**
 * @brief
 * Generation of two dimension mesh
 *
 * @param [StdRegions2d*] shape standard element
 * @param [int]           K     number of element
 * @param [int]           Nv    number of vertex
 * @param [int**]         EToV  element to vertex list
 * @param [double*]       VX    vertex coordinate
 * @param [double*]       VY    vertex coordinate
 *
 * @return
 * name     | type     | description of value
 * -------- |----------|----------------------
 * mesh     | MultiReg2d* | mesh object
 *
 * @note
 * The vertex index in EToV is start from 0.
 */
MultiReg2d* GenMultiReg2d(StdRegions2d *shape, int K, int Nv,
                          int **EToV, double *VX, double *VY){

    MultiReg2d *mesh = (MultiReg2d *)calloc(1, sizeof(MultiReg2d));

    MPI_Comm_rank(MPI_COMM_WORLD, &mesh->procid);
    MPI_Comm_size(MPI_COMM_WORLD, &mesh->nprocs);

    /* standard element */
    mesh->StdElement = shape;
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

//    printf("procid:%d, finish SetFacePair\n", mesh->procid);
    /* Setup nodes coordinate */
    mesh->x  = BuildMatrix(mesh->K, shape->Np);
    mesh->y  = BuildMatrix(mesh->K, shape->Np);
    SetNodeCoord2d(shape, mesh->K, mesh->GX, mesh->GY, mesh->x, mesh->y);
//    printf("procid:%d, finish SetNodeCoor\n", mesh->procid);

    /* Setup boundary nodes connection */
    mesh->vmapM = BuildIntVector(shape->Nfp*shape->Nfaces*mesh->K);
    mesh->vmapP = BuildIntVector(shape->Nfp*shape->Nfaces*mesh->K);
    SetNodePair2d(shape, mesh->K, mesh->GX, mesh->GY,
                mesh->EToE, mesh->EToF, mesh->EToP, mesh->x, mesh->y,
                mesh->Npar, &(mesh->parNtotalout), &(mesh->parmapOUT),
                mesh->vmapM, mesh->vmapP);

    /* mesh geo */
    mesh->vgeo = (float*) calloc(4*mesh->K*shape->Np, sizeof(float));
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
                  double *J, double *area, double *ciradius, float *vgeo){

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
        GeometricFactors(shape->Np, x[k], y[k], shape->Dr, shape->Ds,
                         drdx, dsdx, drdy, dsdy, eJ);
        for(n=0;n<shape->Np;n++){
            vgeo[sk++] = (float) drdx[n];
            vgeo[sk++] = (float) drdy[n];
            vgeo[sk++] = (float) dsdx[n];
            vgeo[sk++] = (float) dsdy[n];
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
 * @param [StdRegions2d*]   shape standard element object
 * @param [int]             K     number of elements
 * @param [double**]        GX    input coordinate of vertex in each element
 * @param [double**]        GY    input coordinate of vertex in each element
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
 * @param [StdRegions2d*]   shape standard element object
 * @param [int]             K     number of elements
 * @param [double*]         VX    input coordinate of vertex
 * @param [double*]         VY    input coordinate of vertex
 * @param [int**]           EToV  element to vertex list
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