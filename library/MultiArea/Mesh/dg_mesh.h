//
// Created by li12242 on 12/19/16.
//

#ifndef DGOM_MR_MESH_H
#define DGOM_MR_MESH_H

#include "../Region/dg_region.h"

/**
 * @brief mesh structure to communicate the infomation with other processes.
 */
typedef struct dg_mesh{

    dg_region *region; ///< multi-region object;

    int *Nface2procs; ///< number of faces adjacent to each process;
    int NfetchFace; ///< total number of parallel faces;
    int *CBFToK; ///< map of cell index for buffer;
    int *CBFToF; ///< map of local face (0,1,2..) index for buffer;

    int *Nfp2procs; ///< number of nodes adjacent to each process;
    int NfetchNode; ///< total number of nodes to send and recv;
    int *NBFToN; ///< index list of nodes to send out;

    /** fetch the node buffer with other processes */
    int (*fetch_node_buffer)(struct dg_mesh *mesh, int Nfield, dg_real *f_Q, dg_real *f_recvQ,
                             MPI_Request *send_requests,
                             MPI_Request *recv_requests);
    /** fetch the cell buffer with other processes */
    int (*fetch_cell_buffer)(struct dg_mesh *mesh, int Nfield, dg_real *f_Q, dg_real *f_recvQ,
                             MPI_Request *send_requests,
                             MPI_Request *recv_requests);
} dg_mesh;

dg_mesh* dg_mesh_create(dg_region *region);
void dg_mesh_free(dg_mesh *mesh);

#define dg_mesh_region(mesh) mesh->region
#define dg_mesh_grid(mesh) dg_region_grid(dg_mesh_region(mesh))
#define dg_mesh_cell(mesh) dg_region_cell(dg_mesh_region(mesh))
#define dg_mesh_NfetchFace(mesh) mesh->NfetchFace
#define dg_mesh_NfetchNode(mesh) mesh->NfetchNode

#endif //DGOM_MR_MESH_H
