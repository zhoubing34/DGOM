//
// Created by li12242 on 12/22/16.
//

#include "phys_fetchBuffer.h"

/**
 * @brief Send/rece nodal value `f_Q` through buffers
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
 *     fetchNodeBuffer2d(phys, mpi_out_requests, mpi_in_requests, &Nmess);
 *
 *     MPI_Status instatus[nprocs];
 *     MPI_Waitall(Nmess, mpi_in_requests, instatus);
 *
 * @param[in]       phys    PhysDomain2d pointer
 * @param[inout]    mpi_send_requests  MPI_Request pointer
 * @param[inout]    mpi_send_requests  MPI_Request pointer
 * @param[inout]    Nmessage number of messages stored in `mpi_send_requests` and `mpi_send_requests`
 *
 * @note Nmess is initialize inside the function
 *
 */
void fetchNodeBuffer2d(physField *phys,
                       MPI_Request *mpi_send_requests,
                       MPI_Request *mpi_recv_requests,
                       int *Nmessage) {
    int n;
    /* buffer outgoing node data */
    for(n=0;n<phys->parallNodeNum;++n)
        phys->f_outQ[n] = phys->f_Q[phys->nodeIndexOut[n]];

    parallMesh *mesh = phys->mesh;

    const int nprocs = mesh->nprocs;
    const int procid = mesh->procid;
    const int Nfield = phys->Nfield;
    const int Nfp = phys->cell->Nfp;

    /* do sends */
    int sk = 0, Nmess = 0;
    int p, Nout;
    for(p=0;p<nprocs;++p){
        if(p!=procid){
            Nout = mesh->Npar[p]*Nfield*Nfp; // # of variables send to process p
            if(Nout){
                /* symmetric communications (different ordering) */
                MPI_Isend(phys->f_outQ+sk, Nout, MPI_SIZE, p, 6666+p,
                          MPI_COMM_WORLD, mpi_send_requests +Nmess);
                MPI_Irecv(phys->f_inQ+sk,  Nout, MPI_SIZE, p, 6666+mesh->procid,
                          MPI_COMM_WORLD,  mpi_recv_requests +Nmess);
                sk+=Nout;
                ++Nmess;
            }
        }
    }
    *Nmessage = Nmess; /* number of messages */
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
 *     fetchNodeBuffer2d(phys, mpi_out_requests, mpi_in_requests, &Nmess);
 *
 *     MPI_Status instatus[nprocs];
 *     MPI_Waitall(Nmess, mpi_in_requests, instatus);
 *
 * @param[in]       phys    PhysDomain2d pointer
 * @param[inout]    mpi_send_requests  MPI_Request pointer
 * @param[inout]    mpi_send_requests  MPI_Request pointer
 * @param[inout]    Nmessage number of messages stored in `mpi_send_requests` and `mpi_send_requests`
 *
 * @note Nmess is initialize inside the function
 *
 */
void fetchCellBuffer(physField *phys,
                     MPI_Request *mpi_send_requests,
                     MPI_Request *mpi_recv_requests,
                     int *Nmessage){

    parallMesh *mesh = phys->mesh;

    const int nprocs = mesh->nprocs;
    const int procid = mesh->procid;

    /* buffer outgoing node data */
    int n;
    for(n=0;n<phys->parallCellNum;++n)
        phys->c_outQ[n] = phys->c_Q[phys->cellIndexOut[n]];

    /* do sends */
    int sk=0, Nmess=0, p;
    for(p=0;p<nprocs;++p){
        if(p!=procid){
            int Nout = mesh->Npar[p]; // # of variables send to process p
            if(Nout){
                /* symmetric communications (different ordering) */
                MPI_Isend(phys->c_outQ+sk, Nout, MPI_SIZE, p, 6666+p,
                          MPI_COMM_WORLD, mpi_send_requests +Nmess);
                MPI_Irecv(phys->c_inQ+sk,  Nout, MPI_SIZE, p, 6666+mesh->procid,
                          MPI_COMM_WORLD,  mpi_recv_requests +Nmess);
                sk+=Nout;
                ++Nmess;
            }
        }
    }
    *Nmessage = Nmess;
}