//
// Created by li12242 on 12/19/16.
//

#ifndef DGOM_MR_MESH_H
#define DGOM_MR_MESH_H

#include "MultiRegions/Region/dg_region.h"

typedef struct dg_mesh{

    int procid; ///< process id
    int nprocs; ///< number of process

    dg_region *region; ///< multi-region object
    dg_grid *grid; ///< geometry grid object
    dg_cell *cell; ///< standard element object

    int *Nface2procs; ///< number of faces adjacent to each process
    int NfetchFace; ///< total number of parallel faces
    int *CBFToK; ///< map of cell index for buffer
    int *CBFToF; ///< map of local face (0,1,2..) index for buffer

    int *Nfp2procs; ///< number of nodes adjacent to each process
    int NfetchNode; ///< total number of nodes to send recv
    int *NBFToN; ///< index list of nodes to send out

    int (*fetch_node_buffer)(struct dg_mesh *mesh, int Nfield, dg_real *f_Q, dg_real *f_recvQ,
                             MPI_Request *send_requests,
                             MPI_Request *recv_requests);
    int (*fetch_cell_buffer)(struct dg_mesh *mesh, int Nfield, dg_real *f_Q, dg_real *f_recvQ,
                             MPI_Request *send_requests,
                             MPI_Request *recv_requests);
} dg_mesh;

dg_mesh* dg_mesh_create(dg_region *region);
void dg_mesh_free(dg_mesh *mesh);

#define dg_mesh_procid(mesh) mesh->procid
#define dg_mesh_nprocs(mesh) mesh->nprocs
#define dg_mesh_Nparf(mesh) mesh->NfetchFace
#define dg_mesh_Nparn(mesh) mesh->NfetchNode

#endif //DGOM_MR_MESH_H