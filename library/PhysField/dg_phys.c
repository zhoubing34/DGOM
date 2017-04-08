/**
 * @file
 * physical fields
 * @brief
 *
 * @author li12242, Tianjin University, li12242@tju.edu.cn
 */

#include "dg_phys.h"

#define DEBUG 0
#if DEBUG
#include "Utility/unit_test.h"
#endif

static void dg_phys_cell_mesn(dg_phys *phys);
static int dg_phys_fetch_cell_buffer(dg_phys *phys, MPI_Request *send_requests, MPI_Request *recv_requests);
static int dg_phys_fetch_node_buffer(dg_phys *phys, MPI_Request *send_requests, MPI_Request *recv_requests);
static void dg_phys_obc_add(dg_phys *phys, char *filename);
static void dg_phys_obc_update(dg_phys *phys, double elapseTime);
static void dg_phys_init_file(dg_phys *phys, char *casename);
static void dg_phys_set_limiter(dg_phys *phys, Limiter_Type type);
static void dg_phys_set_indicator(dg_phys *phys, Indicator_Type type);
static void dg_phys_limit(dg_phys *phys, double parameter);

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
    phys->limiter = dg_phys_limiter_create(phys->info);
    phys->ldg = dg_phys_LDG_create(phys->info);

    phys->cell_mean = dg_phys_cell_mesn;
    phys->initialize_from_file = dg_phys_init_file;
    phys->fetch_cell_buffer = dg_phys_fetch_cell_buffer;
    phys->fetch_node_buffer = dg_phys_fetch_node_buffer;
    phys->attach_obc_ncfile = dg_phys_obc_add;
    phys->update_obc_data = dg_phys_obc_update;
    phys->set_limiter = dg_phys_set_limiter;
    phys->set_indicator = dg_phys_set_indicator;
    phys->limit = dg_phys_limit;
    return phys;
}
/**
 * @brief free the memory of pointer to dg_phys structure.
 * @param phys pointer to a dg_phys structure;
 */
void dg_phys_free(dg_phys *phys){
    dg_phys_LDG_free(phys->ldg);
    dg_phys_info_free(phys->info);
    dg_phys_obc_free(phys->obc);
    dg_phys_limiter_free(phys->limiter);
    free(phys);
    return;
}

static void dg_phys_cell_mesn(dg_phys *phys){
    phys->info->cell_mean(phys->info);
    return;
}

static void dg_phys_set_limiter(dg_phys *phys, Limiter_Type type){
    phys->limiter->set_limiter(phys->limiter, type);
    return;
}

static void dg_phys_set_indicator(dg_phys *phys, Indicator_Type type){
    phys->limiter->set_indicator(phys->limiter, type);
    return;
}

static void dg_phys_limit(dg_phys *phys, double parameter){
    phys->limiter->limit_trouble_cell(phys->limiter ,phys->info, parameter);
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
static void dg_phys_init_file(dg_phys *phys, char *init_filename){

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
    register int n,fld;
    dg_real **init_value = matrix_real_create(Nvert, Nfield);
    for(n=0;n<Nvert;n++){
        int tmp; // node index
        fscanf(fp, "%d", &tmp);
        for(fld=0;fld<Nfield;fld++){
            fscanf(fp, "%lf", init_value[n]+fld);
        }
    }
    fclose(fp);

    /* project the vertex values to nodes */
    grid->proj_vert2node(grid, Nfield, init_value[0], dg_phys_f_Q(phys));
    matrix_real_free(init_value);
    return;
}
