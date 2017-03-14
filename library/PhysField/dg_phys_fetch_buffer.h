//
// Created by li12242 on 12/22/16.
//

#ifndef DGOM_PHYS_FETCHBUFFER_H
#define DGOM_PHYS_FETCHBUFFER_H

#include "dg_phys.h"

/*
 * Send/rece nodal value `f_Q` through buffers
 *
 * Usages:
 *
 *     MPI_Request mpi_out_requests[nprocs];
 *     MPI_Request mpi_in_requests[nprocs];
 *     int Nmess;
 *     phys_fetchNodeBuffer2d(phys, mpi_out_requests, mpi_in_requests, &Nmess);
 *
 *     MPI_Status instatus[nprocs];
 *     MPI_Waitall(Nmess, mpi_in_requests, instatus);
 * */
int dg_phys_fetch_node_buffer(dg_phys *phys,
                              MPI_Request *mpi_send_requests,
                              MPI_Request *mpi_recv_requests);

/*
 * Send/rece elemental value `c_Q` through buffers
 *
 * Usages:
 *
 *     MPI_Request mpi_out_requests[nprocs];
 *     MPI_Request mpi_in_requests[nprocs];
 *     int Nmess;
 *     phys_fetchNodeBuffer2d(phys, mpi_out_requests, mpi_in_requests, &Nmess);
 *
 *     MPI_Status instatus[nprocs];
 *     MPI_Waitall(Nmess, mpi_in_requests, instatus);
 * */
int dg_phys_fetch_cell_buffer(dg_phys *phys,
                              MPI_Request *mpi_send_requests,
                              MPI_Request *mpi_recv_requests);

#endif //DGOM_PHYS_FETCHBUFFER_H
