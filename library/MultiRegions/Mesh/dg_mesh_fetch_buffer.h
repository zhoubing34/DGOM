//
// Created by li12242 on 17/3/9.
//

#ifndef DGOM_DG_MESH_FETCH_BUFFER_H
#define DGOM_DG_MESH_FETCH_BUFFER_H

void dg_mesh_fetch_cell_buffer(dg_mesh *mesh, int Nfield, dg_real *f_Q, dg_real *f_recvQ,
                               MPI_Request *send_requests,
                               MPI_Request *recv_requests,
                               int *Nmess);

void dg_mesh_fetch_node_buffer(dg_mesh *mesh, int Nfield, dg_real *f_Q, dg_real *f_recvQ,
                               MPI_Request *send_requests,
                               MPI_Request *recv_requests,
                               int *Nmess);

#endif //DGOM_DG_MESH_FETCH_BUFFER_H
