#include "dg_grid.h"
#include "dg_grid_retreatEToV.h"
#include "dg_grid_loadBalance.h"
#include "dg_grid_connect.h"
#include "dg_grid_BS.h"

/* allocate and assignment the vertex coordinate in 2d dg_grid object */
static void dg_grid_set_vert2d(dg_grid *grid, double *vx, double *vy, double *vz);
/* allocate and assignment the vertex coordinate in 3d dg_grid object */
static void dg_grid_set_vert3d(dg_grid *grid, double *vx, double *vy, double *vz);
/* partition of the whole grid into each process */
static void dg_grid_partition(dg_grid *grid);
static void dg_grid_proj_vert2node(dg_grid *grid, int Nfield, double *vertval, double *nodeval);

/**
 * @brief creator for generating dg_grid structure.
 */
typedef struct dg_grid_creator{
    void (*copy_vert)(dg_grid *grid, double *vx, double *vy, double *vz);
    void (*partition)(dg_grid *grid);
    void (*restructEToV)(dg_grid *grid);
    void (*load_balance)(dg_grid *grid);
    void (*connect)(dg_grid *grid);
    void (*init_BS)(dg_grid *grid);
    void (*add_BS)(dg_grid *grid, int Nsurface, int **SFToV);
    void (*add_BS_from_file)(struct dg_grid *grid, char *casename);
} dg_grid_creator;

/** const creator for 2d grid */
const dg_grid_creator grid_2d_creator = {
        dg_grid_set_vert2d,
        dg_grid_partition,
        dg_grid_retreatEToV2d,
        dg_grid_load_balance2d,
        dg_grid_connect2d,
        dg_grid_init_BS,
        dg_grid_add_BS2d,
        dg_grid_add_BS_file2d,
};
/** const creator for 3d grid */
const dg_grid_creator grid_3d_creator = {
        dg_grid_set_vert3d,
        dg_grid_partition,
        dg_grid_retreatEToV3d,
        dg_grid_load_balance3d,
        dg_grid_connect3d,
        dg_grid_init_BS,
        dg_grid_add_BS3d,
        dg_grid_add_BS_file3d,
};

/**
 * @brief
 * create dg_grid structure with input grid geometry.
 * @param cell pointer to the dg_cell structure;
 * @param K number of element;
 * @param Nv number of vertex;
 * @param vx,vy,vz coordinate of vertex;
 * @param EToV element to vertex list;
 * @return
 * pointer to a new dg_grid structure.
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
    switch ( dg_cell_celltype(cell) ){
        case TRIANGLE:
            creator = &grid_2d_creator; break;
        case QUADRIL:
            creator = &grid_2d_creator; break;
        default:
            fprintf(stderr, "%s (%d): Unknown cell type %d\n",
                    __FUNCTION__, __LINE__, dg_cell_celltype(cell));
            exit(-1);
    }

    /* assignment of EToV */
    int k,i;
    int nv = dg_cell_Nv(cell);
    grid->EToV = matrix_int_create(K, nv);
    for(k=0;k<K;k++){
        for(i=0;i<nv;i++) {
            grid->EToV[k][i] = EToV[k][i];
        }
    }
    /* assignment of EToR */
    grid->EToR = vector_int_create(K);

    /* copy vertex */
    creator->copy_vert(grid, vx, vy, vz);
    creator->partition(grid);
    creator->restructEToV(grid);
    creator->load_balance(grid);
    creator->connect(grid);
    creator->init_BS(grid);
    grid->add_BS = creator->add_BS;
    grid->add_BS_from_file = creator->add_BS_from_file;
    grid->proj_vert2node = dg_grid_proj_vert2node;
    return grid;
}

/**
 * @brief free the memory of dg_grid type.
 * @param grid dg_grid structure;
 */
void dg_grid_free(dg_grid *grid){
    matrix_int_free(grid->EToV);
    matrix_int_free(grid->EToE);
    matrix_int_free(grid->EToF);
    matrix_int_free(grid->EToP);
    matrix_int_free(grid->EToBS);
    vector_int_free(grid->EToR);

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
    const int Nv = dg_cell_Nv(grid->cell);

    int K = dg_grid_K(grid);
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

    int **newEToV = matrix_int_create(Klocal, dg_cell_Nv(grid->cell));
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
    return;
}

/**
 * @brief allocate and assignment the vertex coordinate in 2d geoGrid object
 * @param[in,out] grid geoGrid object
 * @param[in] vx coordinate of vertex
 * @param[in] vy coordinate of vertex
 */
static void dg_grid_set_vert2d(dg_grid *grid, double *vx, double *vy, double *vz){
    const int Nv = dg_grid_Nv(grid);
    int i;
    grid->vx  = vector_double_create(Nv);
    grid->vy  = vector_double_create(Nv);
    grid->vz  = vector_double_create(Nv);
    // assignment
    for(i=0;i<Nv;i++){
        grid->vx[i] = vx[i];
        grid->vy[i] = vy[i];
    }
}

/**
 * @brief allocate and assignment the vertex coordinate to grid;
 * @param[in,out] grid geoGrid object
 * @param[in] vx coordinate of vertex
 * @param[in] vy coordinate of vertex
 * @param[in] vz coordinate of vertex
 */
static void dg_grid_set_vert3d(dg_grid *grid, double *vx, double *vy, double *vz){
    const int Nv = dg_grid_Nv(grid);
    int i;
    grid->vx  = vector_double_create(Nv);
    grid->vy  = vector_double_create(Nv);
    grid->vz  = vector_double_create(Nv);
    // assignment
    for(i=0;i<Nv;i++){
        grid->vx[i] = vx[i];
        grid->vy[i] = vy[i];
        grid->vz[i] = vz[i];
    }
}
/**
 * @brief project the vertex values to nodes.
 * @param grid
 * @param Nfield
 * @param vertval
 * @param nodeval
 */
static void dg_grid_proj_vert2node(dg_grid *grid, int Nfield, double *vertval, double *nodeval){

    /* read element file */
    int **EToV = dg_grid_EToV(grid);

    int k,n,fld;
    /* assign to node fields */
    dg_cell *cell = dg_grid_cell(grid);
    const int Np = dg_cell_Np(cell);
    const int Nv = dg_cell_Nv(cell);
    const int K = dg_grid_K(grid);

    dg_real initial_vert_value[Nv*Nfield];
    dg_real initial_node_value[Np*Nfield];
    for(k=0;k<K;k++){
        for(n=0;n<Nv;n++){ // vertex initial values
            int vertID = EToV[k][n];
            for(fld=0;fld<Nfield;fld++){
                initial_vert_value[n*Nfield+fld] = vertval[vertID*Nfield+fld];
            }
        }

        cell->proj_vert2node(cell, Nfield, initial_vert_value, initial_node_value);
        int sk = k*Np*Nfield;
        for(n=0;n<Np*Nfield;n++){ // assign to node values
            nodeval[sk + n] = initial_node_value[n];
        }
    }
    return;
}