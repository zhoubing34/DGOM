//
// Created by li12242 on 17/3/22.
//

#ifndef DGOM_DG_PHYS_INFO_H
#define DGOM_DG_PHYS_INFO_H

#include "MultiRegions/Edge/dg_edge.h"

typedef struct dg_phys_info{
    int Nfield;

    dg_edge *edge; ///< parallel edge structure;
    dg_mesh *mesh; ///< parallel mesh structure;
    dg_region *region; ///< multi-region structure;
    dg_grid *grid; ///< geometry grid structure;
    dg_cell *cell; ///< standard cell structure;

    dg_real *c_Q; ///< cell information;
    dg_real *c_recvQ; ///< recv buffers for cell information;

    dg_real *f_Q; ///< nodal information;
    dg_real *f_recvQ; ///< recv buffers for nodal information;

    dg_real *f_rhsQ; ///< right hand side data;
    dg_real *f_resQ; ///< residual data;

    /** function to fetch node buffer with other process; */
    int (*fetch_node_buffer)(struct dg_phys_info *phys, MPI_Request *send_requests, MPI_Request *recv_requests);
    /** function to fetch cell buffer with other process; */
    int (*fetch_cell_buffer)(struct dg_phys_info *phys, MPI_Request *send_requests, MPI_Request *recv_requests);

}dg_phys_info;

dg_phys_info* dg_phys_info_create(int Nfields, dg_edge *edge);
void dg_phys_info_free(dg_phys_info *phys_info);

#endif //DGOM_DG_PHYS_INFO_H
