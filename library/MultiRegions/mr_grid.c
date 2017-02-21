#include <StandCell/sc_stdcell.h>
#include "mr_grid.h"
#include "mr_grid_resortEToV.h"
#include "mr_grid_loadBalance2d.h"

/* allocate and assignment the vertex coordinate in 2d geoGrid object */
void mr_copyVertex2d(geoGrid *grid, double *vx, double *vy);

/* allocate and assignment the vertex coordinate in 3d geoGrid object */
void mr_copyVertex3d(geoGrid *grid, double *vx, double *vy, double *vz);

/* partition of the whole grid into each process */
void mr_grid_partition(geoGrid *grid);

/* copy the vertex vx,vy,vz into geometry grid object */
void mr_copyVertex(geoGrid *grid, double *vx, double *vy, double *vz);

/* resort the vertex list in EToV */
void mr_resortEToV(geoGrid *grid);


/**
 * @brief
 * Deallocate unstructured grid structure.
 */
void mr_grid_free(geoGrid *grid) {

    const int dim = grid->dim;
    matrix_int_free(grid->EToV);

    switch (dim){ // allocate vertex
        case 2:
            vector_double_free(grid->vx);
            vector_double_free(grid->vy);
            break;
        case 3:
            vector_double_free(grid->vx);
            vector_double_free(grid->vy);
            vector_double_free(grid->vz);
            break;
        default:
            printf("MultuRegions (mr_grid_free): Wrong dimension %d.", dim);
            exit(-1);
    }

    free(grid);
    return;
}

/**
 * @brief
 * Create and allocate geometry grid object.
 *
 */
geoGrid* mr_grid_create(stdCell *shape, int K, int Nv, double *vx, double *vy, double *vz, int **EToV){

    geoGrid *grid = (geoGrid*) calloc(1, sizeof(geoGrid));

    grid->dim = shape->dim;
    grid->cell = shape;
    grid->type = shape->type;
    grid->K  = K;
    grid->Nv  = Nv;

    MPI_Comm_rank(MPI_COMM_WORLD, &grid->procid);
    MPI_Comm_size(MPI_COMM_WORLD, &grid->nprocs);

    // allocation and assignment of EToV
    int k,i;
    grid->EToV = matrix_int_create(K, shape->Nv);
    for(k=0;k<K;k++){
        for(i=0;i<shape->Nv;i++)
            grid->EToV[k][i] = EToV[k][i];
    }

    /* allocation and assignment of vertex */
    mr_copyVertex(grid, vx, vy, vz);
    /* resort the vertex list */
    mr_resortEToV(grid);

    /* partition of the whole grid */
    mr_grid_partition(grid);
    /* redistribute the elements to achieve the load balance */
    mr_grid_loadBalance2d(grid);

    return grid;
}

/**
 * @brief resort the vertex list in EToV
 * @param[in,out] grid geometry grid object
 */
void mr_resortEToV(geoGrid *grid){
    const int K = grid->K;
    const int dim = grid->cell->dim;
    stdCell *shape = grid->cell;
    int k;
    switch (dim){
        case 2:
            for(k=0;k<K;k++){
                mr_resortEToV2d(shape->Nv, grid->vx, grid->vy, grid->EToV[k]);
            }
            break;
        case 3:
        default:
            printf("MultuRegions (mr_resortEToV): Wrong dimensions %d\n", dim);
            exit(-1);
    }
}


/**
 * @brief copy the vertex vx,vy,vz into geometry grid object.
 *
 * @param[in,out] grid multi-region object.
 * @param[in] vx coordinate x of vertex
 * @param[in] vy coordinate y of vertex
 * @param[in] vz coordinate z of vertex
 */
void mr_copyVertex(geoGrid *grid, double *vx, double *vy, double *vz){
    const int dim = grid->cell->dim;
    switch (dim){
        case 2:
            mr_copyVertex2d(grid, vx, vy); break;
        case 3:
            mr_copyVertex3d(grid, vx, vy, vz); break;
        default:
            printf("MultuRegions (mr_grid_create): Wrong dimensions %d\n", dim);
            exit(-1);
    }
}


/**
 * @brief partition of the whole grid into each process
 * @details
 * The elements will be distributed into each process according to their index.
 * @param[in,out] grid geometry grid object
 */
void mr_grid_partition(geoGrid *grid){

    const int nprocs = grid->nprocs;
    const int procid = grid->procid;
    const int Nv = grid->cell->Nv;

    int K       = grid->K;
    int **EToV  = grid->EToV;

    int Kprocs[nprocs];
    int Klocal = K/nprocs;

    /* number of elements in each process */
    int p;
    for(p=0; p<(nprocs-1); ++p){
        Kprocs[p] = Klocal;
    }
    /* the rest elements is distributed to the last process */
    Kprocs[nprocs-1] = K - (nprocs-1)*Klocal;

    /* number of element in local process */
    Klocal = Kprocs[procid];

    /* start index of element */
    int Kstart = 0;
    for(p=0;p<procid;++p)
        Kstart += Kprocs[p];

    int **newEToV = matrix_int_create(Klocal, grid->cell->Nv);
    /* Assignment of local mesh information */
    int sk=0,n,i;
    for (n=0; n<K; n++){
        if(n>=Kstart && n<Kstart+Klocal) {
            for(i=0;i<Nv;i++) {newEToV[sk][i] = EToV[n][i];}
            ++sk;
        }
    }
    /* free the original EToV and assignment */
    matrix_int_free(grid->EToV);
    grid->EToV = newEToV;
    grid->K = Klocal;
}

/**
 * @brief allocate and assignment the vertex coordinate in 2d geoGrid object
 * @param[in,out] grid geoGrid object
 * @param[in] vx coordinate of vertex
 * @param[in] vy coordinate of vertex
 */
void mr_copyVertex2d(geoGrid *grid, double *vx, double *vy){
    const int Nv = grid->Nv;
    int i;
    grid->vx  = vector_double_create(Nv);
    grid->vy  = vector_double_create(Nv);
    grid->vz  = NULL;
    // allocation
    for(i=0;i<Nv;i++){
        grid->vx[i] = vx[i];
        grid->vy[i] = vy[i];
    }
}

/**
 * @brief allocate and assignment the vertex coordinate in 3d geoGrid object
 * @param[in,out] grid geoGrid object
 * @param[in] vx coordinate of vertex
 * @param[in] vy coordinate of vertex
 * @param[in] vz coordinate of vertex
 */
void mr_copyVertex3d(geoGrid *grid, double *vx, double *vy, double *vz){
    const int Nv = grid->Nv;
    int i;
    grid->vx  = vector_double_create(Nv);
    grid->vy  = vector_double_create(Nv);
    grid->vz  = vector_double_create(Nv);
    // allocation
    for(i=0;i<Nv;i++){
        grid->vx[i] = vx[i];
        grid->vy[i] = vy[i];
        grid->vz[i] = vz[i];
    }
}