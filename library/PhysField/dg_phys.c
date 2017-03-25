/**
 * @file
 * physical fields
 * @brief
 *
 * @author li12242, Tianjin University, li12242@tju.edu.cn
 */

#include "dg_phys.h"

static int dg_phys_fetch_cell_buffer(dg_phys *phys, MPI_Request *send_requests, MPI_Request *recv_requests);
static int dg_phys_fetch_node_buffer(dg_phys *phys, MPI_Request *send_requests, MPI_Request *recv_requests);
void dg_phys_obc_add(dg_phys *phys, char *filename);
void dg_phys_obc_update(dg_phys *phys, double elapseTime);
static void dg_phys_init_file(dg_phys *phys, char *casename);
static void dg_phys_set_limiter(dg_phys *phys, Limiter_Type type);
static void dg_phys_limit(dg_phys *phys, double parameter);

#define DEBUG 0

/**
 * @brief create a pointer to a new dg_phys structure.
 * @param Nfields number of physical field;
 * @param edge pointer to dg_edge structure;
 * @return
 * phys pointer to a new dg_phys structure.
 */
dg_phys* dg_phys_create(int Nfields, dg_edge *edge){

    dg_phys *phys = (dg_phys *) calloc(1, sizeof(dg_phys));
    phys->info = dg_phys_info_create(Nfields, edge);
    phys->obc = dg_phys_obc_create(phys->info);
    phys->limiter = dg_phys_limiter_create();

    phys->init_file = dg_phys_init_file;
    phys->fetch_cell_buffer = dg_phys_fetch_cell_buffer;
    phys->fetch_node_buffer = dg_phys_fetch_node_buffer;
    phys->obc_add = dg_phys_obc_add;
    phys->obc_update = dg_phys_obc_update;
    phys->set_limiter = dg_phys_set_limiter;
    phys->limit = dg_phys_limit;
    return phys;
}
/**
 * @brief free the memory of pointer to dg_phys structure.
 * @param phys pointer to a dg_phys structure;
 */
void dg_phys_free(dg_phys *phys){
    dg_phys_info_free(phys->info);
    dg_phys_obc_free(phys->obc);
    free(phys);
    return;
}

static void dg_phys_set_limiter(dg_phys *phys, Limiter_Type type){
    phys->limiter->set_limiter(phys->limiter, type);
    return;
}

static void dg_phys_limit(dg_phys *phys, double parameter){
    phys->limiter->limit(phys->info, parameter);
    return;
}

static void dg_phys_obc_add(dg_phys *phys, char *filename){
    phys->obc->add_obc(phys->obc, filename);
    return;
}
static void dg_phys_obc_update(dg_phys *phys, double elapseTime){
    phys->obc->update_obc(phys->obc, elapseTime);
    return;
}

/**
 * @brief
 * @param phys
 * @param send_requests
 * @param recv_requests
 * @return
 */
static int dg_phys_fetch_cell_buffer(dg_phys *phys,
                                      MPI_Request *send_requests,
                                      MPI_Request *recv_requests){
    dg_phys_info *info = phys->info;
    int Nmess;
    Nmess = info->fetch_cell_buffer(info, send_requests, recv_requests);
    return Nmess;
}
/**
 * @brief
 * @param phys
 * @param send_requests
 * @param recv_requests
 * @return
 */
static int dg_phys_fetch_node_buffer(dg_phys *phys,
                                      MPI_Request *send_requests,
                                      MPI_Request *recv_requests){
    dg_phys_info *info = phys->info;
    int Nmess;
    Nmess = info->fetch_node_buffer(info, send_requests, recv_requests);
    return Nmess;
}

/***
 * @brief initialize the physical fields with initial condition file.
 * @param phys pointer to dg_phys structure;
 * @param casename name of case;
 * @note
 * The physical field is initialized with the vertex value
 */
static void dg_phys_init_file(dg_phys *phys, char *casename){

    char init_filename[MAX_NAME_LENGTH];
    strcpy(init_filename, casename);
    strcat(init_filename, ".init");

    dg_grid *grid = dg_phys_grid(phys);

    /* open the init file */
    FILE *fp;
    dg_fopen(fp, init_filename, "Unable to open initial condition file");

    int Nvert, Nfield;
    fscanf(fp, "%d %d", &Nvert, &Nfield);
    if( Nvert != dg_grid_Nv(grid) ){
        fprintf(stderr, "%s (%d): Wrong number of vertex in initial file %s.\n",
                __FILE__, __LINE__, init_filename);
    }
    if( Nfield != dg_phys_Nfield(phys) ){
        fprintf(stderr, "%s (%d): Wrong physical field number in initial file %s.\n",
                __FILE__, __LINE__, init_filename);
    }

    /* read vertex initial values */
    register int k,n,fld;
    dg_real **init_value = matrix_real_create(Nvert, Nfield);
    for(n=0;n<Nvert;n++){
        int tmp; // node index
        fscanf(fp, "%d", &tmp);
        for(fld=0;fld<Nfield;fld++){
            fscanf(fp, "%lf", init_value[n]+fld);
        }
    }

    /* read element file */
    int **EToV = dg_grid_EToV(grid);

    /* assign to node fields */
    dg_cell *cell = dg_phys_cell(phys);
    const int Np = dg_cell_Np(cell);
    const int Nv = dg_cell_Nv(cell);
    const int K = dg_grid_K(grid);

    dg_real initial_vert_value[Nv*Nfield];
    dg_real initial_node_value[Np*Nfield];
    dg_real *f_Q = dg_phys_f_Q(phys);
    for(k=0;k<K;k++){
        for(n=0;n<Nv;n++){ // vertex initial values
            int vertID = EToV[k][n];
            for(fld=0;fld<Nfield;fld++){
                initial_vert_value[n*Nfield+fld] = init_value[vertID][fld];
            }
        }

        cell->proj_vert2node(cell, Nfield, initial_vert_value, initial_node_value);
        int sk = (k*Np+n)*Nfield;
        for(n=0;n<Np*Nfield;n++){ // assign to node values
            f_Q[sk + n] = initial_node_value[n];
        }

    }

    matrix_real_free(init_value);
    matrix_int_free(EToV);
    return;
}
