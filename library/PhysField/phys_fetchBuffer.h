//
// Created by li12242 on 12/22/16.
//

#ifndef DGOM_PHYS_FETCHBUFFER_H
#define DGOM_PHYS_FETCHBUFFER_H

#include "phs_physField.h"

/*
 * Send/rece nodal value `f_Q` through buffers
 *
 * Usages:
 *
 *     MPI_Request mpi_out_requests[nprocs];
 *     MPI_Request mpi_in_requests[nprocs];
 *     int Nmess;
 *     fetchNodeBuffer2d(phys, mpi_out_requests, mpi_in_requests, &Nmess);
 *
 *     MPI_Status instatus[nprocs];
 *     MPI_Waitall(Nmess, mpi_in_requests, instatus);
 * */
void fetchNodeBuffer2d(physField *phys,
                       MPI_Request *mpi_send_requests,
                       MPI_Request *mpi_recv_requests,
                       int *Nmessage);

/*
 * Send/rece elemental value `c_Q` through buffers
 *
 * Usages:
 *
 *     MPI_Request mpi_out_requests[nprocs];
 *     MPI_Request mpi_in_requests[nprocs];
 *     int Nmess;
 *     fetchNodeBuffer2d(phys, mpi_out_requests, mpi_in_requests, &Nmess);
 *
 *     MPI_Status instatus[nprocs];
 *     MPI_Waitall(Nmess, mpi_in_requests, instatus);
 * */
void fetchCellBuffer(physField *phys,
                     MPI_Request *mpi_send_requests,
                     MPI_Request *mpi_recv_requests,
                     int *Nmessage);

#endif //DGOM_PHYS_FETCHBUFFER_H
