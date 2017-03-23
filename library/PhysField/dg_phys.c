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

    phys->fetch_cell_buffer = dg_phys_fetch_cell_buffer;
    phys->fetch_node_buffer = dg_phys_fetch_node_buffer;
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

void dg_phys_obc_add(dg_phys *phys, char *filename){
    phys->obc->add_obc(phys->obc, filename);
    return;
}
void dg_phys_obc_update(dg_phys *phys, double elapseTime){
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
