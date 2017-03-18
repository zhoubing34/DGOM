//
// Created by li12242 on 12/22/16.
//

#include "dg_phys_fetch_buffer.h"

/**
 * @brief Send/rece nodal value node information through buffers
 *
 * @details
 * The internal boundary values of `f_Q` in local process is arranged into `f_outQ` buffer
 * to send to other processes. While `mpi_send_requests` get the request of each
 * MPI_Isend function. The incoming boundary values is stored into the `f_inQ` buffer
 * and the request are stored in `mpi_recv_requests`.
 * `Nmess` gets the number of send and receive requests.
 *
 * Usages:
 *
 *     MPI_Request mpi_out_requests[nprocs];
 *     MPI_Request mpi_in_requests[nprocs];
 *     int Nmess;
 *     Nmess = phys_fetchNodeBuffer2d(phys, mpi_out_requests, mpi_in_requests);
 *
 *     MPI_Status instatus[nprocs];
 *     MPI_Waitall(Nmess, mpi_in_requests, instatus);
 *
 * @param[in] phys pointer to dg_phys structure;
 * @param[in,out] mpi_send_requests MPI send request
 * @param[in,out] mpi_recv_requests MPI receive request
 * @return
 * Nmess number of messages in `mpi_send_requests` and `mpi_send_requests`
 *
 */
int dg_phys_fetch_node_buffer(dg_phys *phys,
                              MPI_Request *mpi_send_requests,
                              MPI_Request *mpi_recv_requests) {
    dg_mesh *mesh = phys->mesh;
    const int Nfield = dg_phys_Nfield(phys);
    int Nmess;
    Nmess = mesh->fetch_node_buffer(mesh, Nfield, phys->f_Q, phys->f_recvQ,
                                    mpi_send_requests, mpi_recv_requests);
    return Nmess;
}

/**
 * @brief Send/rece elemental value `c_Q` through buffers
 *
 * @details
 * The internal boundary elemental values of `c_Q` in local process is arranged into
 * `c_outQ` buffer to send to other processes. While `mpi_send_requests` get
 * the request of each `MPI_Isend` function. The incoming elemental values is stored
 * into the `c_inQ` buffer and the request are stored in `mpi_recv_requests`.
 * `Nmess` gets the number of send and receive requests.
 *
 * Usages:
 *
 *     MPI_Request mpi_out_requests[nprocs];
 *     MPI_Request mpi_in_requests[nprocs];
 *     int Nmess;
 *     phys_fetchCellBuffer(phys, mpi_out_requests, mpi_in_requests, &Nmess);
 *
 *     MPI_Status instatus[nprocs];
 *     MPI_Waitall(Nmess, mpi_in_requests, instatus);
 *
 * @param[in]       phys    PhysDomain2d pointer
 * @param[in,out]    mpi_send_requests  MPI_Request pointer
 * @param[in,out]    mpi_recv_requests  MPI_Request pointer
 * @param[in,out]    Nmessage number of messages stored in `mpi_send_requests` and `mpi_send_requests`
 *
 * @note Nmess is initialize inside the function
 *
 */
int dg_phys_fetch_cell_buffer(dg_phys *phys,
                              MPI_Request *mpi_send_requests,
                              MPI_Request *mpi_recv_requests){

    dg_mesh *mesh = phys->mesh;
    const int Nfield = dg_phys_Nfield(phys);
    int Nmess;
    Nmess = mesh->fetch_cell_buffer(mesh, Nfield, phys->f_Q, phys->f_recvQ,
                                    mpi_send_requests, mpi_recv_requests);
    return Nmess;
}