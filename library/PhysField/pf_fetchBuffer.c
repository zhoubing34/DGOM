//
// Created by li12242 on 12/22/16.
//

#include <MultiRegions/Mesh/dg_mesh.h>
#include "pf_fetchBuffer.h"
#include "pf_phys.h"

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
 *     phys_fetchNodeBuffer2d(phys, mpi_out_requests, mpi_in_requests, &Nmess);
 *
 *     MPI_Status instatus[nprocs];
 *     MPI_Waitall(Nmess, mpi_in_requests, instatus);
 *
 * @param[in]       phys    PhysDomain2d pointer
 * @param[in,out]    mpi_send_requests  MPI_Request send request
 * @param[in,out]    mpi_recv_requests  MPI_Request receive request
 * @param[in,out]    Nmessage number of messages stored in `mpi_send_requests` and `mpi_send_requests`
 *
 * @note Nmess is initialize inside the function
 *
 */
void pf_fetchNodeBuffer2d(physField *phys,
                          MPI_Request *mpi_send_requests,
                          MPI_Request *mpi_recv_requests,
                          int *Nmessage) {

    int n;
    /* buffer outgoing node data */
    for(n=0;n<phys->parallNodeNum;++n)
        phys->f_outQ[n] = phys->f_Q[phys->nodeIndexOut[n]];

    dg_mesh *mesh = phys->mesh;

    const int nprocs = mesh->nprocs;
    const int procid = mesh->procid;
    const int Nfield = phys->Nfield;
    const int Nfp = phys->cell->Nfp;

    int Nout[nprocs];
    for(n=0;n<nprocs;n++){
        Nout[n] = mesh->Parf[n]*Nfield*Nfp;
    }

    /* do sends and recv */
    pf_fetchBuffer(procid, nprocs, Nout, phys->f_outQ, phys->f_inQ,
                   mpi_send_requests, mpi_recv_requests, Nmessage);

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
void pf_fetchCellBuffer(physField *phys,
                        MPI_Request *mpi_send_requests,
                        MPI_Request *mpi_recv_requests,
                        int *Nmessage){

    dg_mesh *mesh = phys->mesh;

    const int nprocs = mesh->nprocs;
    const int procid = mesh->procid;
    const int Nfield = phys->Nfield;

    /* buffer outgoing node data */
    int n;
    for(n=0;n<phys->parallCellNum;++n)
        phys->c_outQ[n] = phys->c_Q[phys->cellIndexOut[n]];

    int Nout[nprocs];
    for(n=0;n<nprocs;n++){
        Nout[n] = mesh->Parf[n]*Nfield;
    }

    /* do sends and recv */
    pf_fetchBuffer(procid, nprocs, Nout, phys->c_outQ, phys->c_inQ,
                   mpi_send_requests, mpi_recv_requests, Nmessage);
}


/**
 * @brief send and receive buffers to/from other processes.
 * @param [in]     procid No. of local process
 * @param [in]     number of all the process
 * @param [in] pout number of data to send to each process
 * @param [in] send_buffer buffer for sending
 * @param [in] recv_buffer buffer for recving
 * @param [in,out] mpi_send_requests MPI request for MPI_Isend operation
 * @param [in,out] mpi_recv_requests MPI request for MPI_Irecv operation
 * @param [in,out] Nmessage number of message
 *
 */
void pf_fetchBuffer(int procid, int nprocs, int *pout,
                    dg_real *send_buffer, dg_real *recv_buffer,
                    MPI_Request *mpi_send_requests,
                    MPI_Request *mpi_recv_requests,
                    int *Nmessage){
    /* do sends */
    int sk=0, Nmess=0, p;
    for(p=0;p<nprocs;++p){
        if(p!=procid){
            const int Nout = pout[p]; // # of variables send to process p
            if(Nout){
                /* symmetric communications (different ordering) */
                MPI_Isend(send_buffer+sk, Nout, MPI_TYPE, p, 5666+p,
                          MPI_COMM_WORLD, mpi_send_requests +Nmess);
                MPI_Irecv(recv_buffer+sk,  Nout, MPI_TYPE, p, 5666+procid,
                          MPI_COMM_WORLD,  mpi_recv_requests +Nmess);
                sk+=Nout;
                ++Nmess;
            }
        }
    }
    *Nmessage = Nmess;
}