#include "dg_grid.h"
#include "dg_grid_retreatEToV.h"
#include "dg_grid_loadBalance.h"

/* allocate and assignment the vertex coordinate in 2d dg_grid object */
void dg_grid_copyVert2d(dg_grid *grid, double *vx, double *vy, double *vz);

/* allocate and assignment the vertex coordinate in 3d dg_grid object */
void dg_grid_copyVert3d(dg_grid *grid, double *vx, double *vy, double *vz);

/* partition of the whole grid into each process */
void dg_grid_partition(dg_grid *grid);
void dg_grid_free2d(dg_grid *grid);
void dg_grid_free3d(dg_grid *grid);

typedef struct dg_grid_creator{
    void (*copy_vert)(dg_grid *grid, double *vx, double *vy, double *vz);
    void (*partition)(dg_grid *grid);
    void (*retreatEToV)(dg_grid *grid);
    void (*loadBalance)(dg_grid *grid);
    void (*free_func)(dg_grid *grid);
} dg_grid_creator;

const dg_grid_creator grid_2d_creator = {
        dg_grid_copyVert2d,
        dg_grid_partition,
        dg_grid_retreatEToV2d,
        mr_grid_loadBalance2d,
        dg_grid_free2d,
};

const dg_grid_creator grid_3d_creator = {
        dg_grid_copyVert3d,
        dg_grid_partition,
        dg_grid_retreatEToV3d,
        mr_grid_loadBalance3d,
        dg_grid_free3d,
};

/**
 * @brief
 * create dg_grid structure with input grid geometry.
 * @param cell standard cell
 * @param K number of element
 * @param Nv number of vertex
 * @param vx,vy,vz coordinate of vertex
 * @param EToV element to vertex list
 * @return
 */
dg_grid* dg_grid_create(dg_cell *cell, int K, int Nv, double *vx, double *vy, double *vz, int **EToV){

    dg_grid *grid = (dg_grid*) calloc(1, sizeof(dg_grid));
    /* basic info */
    grid->cell = cell;
    grid->K = K;
    grid->Nv = Nv;
    MPI_Comm_rank(MPI_COMM_WORLD, &grid->procid);
    MPI_Comm_size(MPI_COMM_WORLD, &grid->nprocs);

    const dg_grid_creator *creator;
    switch (cell->type){
        case TRIANGLE:
            creator = &grid_2d_creator; break;
        case QUADRIL:
            creator = &grid_2d_creator; break;
        default:
            fprintf(stderr, "%s (%d): Unknown cell type %d\n", __FUNCTION__, __LINE__, cell->type);
            exit(-1);
    }

    /* assignment of EToV */
    int k,i;
    grid->EToV = matrix_int_create(K, cell->Nv);
    for(k=0;k<K;k++){
        for(i=0;i<cell->Nv;i++)
            grid->EToV[k][i] = EToV[k][i];
    }
    /* copy vertex */
    creator->copy_vert(grid, vx, vy, vz);
    creator->partition(grid);
    creator->retreatEToV(grid);
    creator->loadBalance(grid);

    grid->free_func = creator->free_func;
    return grid;
}


void dg_grid_free(dg_grid *grid){
    grid->free_func(grid);
    return;
}

static void dg_grid_free2d(dg_grid *grid){
    matrix_int_free(grid->EToV);

    vector_double_free(grid->vx);
    vector_double_free(grid->vy);
    free(grid);
    return;
}

static void dg_grid_free3d(dg_grid *grid){
    matrix_int_free(grid->EToV);

    vector_double_free(grid->vx);
    vector_double_free(grid->vy);
    vector_double_free(grid->vz);
    free(grid);
    return;
}

/**
 * @brief
 * partition of the whole grid equally for each process
 * @details
 * The elements will be distributed into each process according to their index.
 * @param[in,out] grid dg_grid structure
 */
static void dg_grid_partition(dg_grid *grid){

    const int nprocs = grid->nprocs;
    const int procid = grid->procid;
    const int Nv = grid->cell->Nv;

    int K = grid->K;
    int **EToV = grid->EToV;

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
    for(p=0;p<procid;++p) {Kstart += Kprocs[p];}

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
static void dg_grid_copyVert2d(dg_grid *grid, double *vx, double *vy, double *vz){
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
static void dg_grid_copyVert3d(dg_grid *grid, double *vx, double *vy, double *vz){
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