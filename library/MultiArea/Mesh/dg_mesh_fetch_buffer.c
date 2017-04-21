//
// Created by li12242 on 17/3/9.
//
#include "dg_mesh_fetch_buffer.h"

static void fetch_buffer(int procid, int nprocs, int *pout, dg_real *send_buffer, dg_real *recv_buffer,
                         MPI_Request *mpi_send_requests,
                         MPI_Request *mpi_recv_requests,
                         int *Nmessage);
/**
 * @brief fetch the node buffer with other processes.
 * @param mesh pointer to the dg_mesh structure;
 * @param Nfield number of files;
 * @param f_Q node information;
 * @param f_recvQ receive buffer;
 * @param send_requests, recv_requests MPI send and receive request;
 * @return
 * Nmess number of message to send and receive.
 */
int dg_mesh_fetch_node_buffer(dg_mesh *mesh, int Nfield, dg_real *f_Q, dg_real *f_recvQ,
                               MPI_Request *send_requests,
                               MPI_Request *recv_requests){
    const int Nparn = dg_mesh_NfetchNode(mesh);
    int procid,nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    dg_real *f_sendQ = vector_double_create(Nfield*Nparn);
    /* prepare for send data */
    int n,fld;
    for(n=0;n<Nparn;n++){
        const int sk = n*Nfield;
        const int nind = mesh->NBFToN[n]*Nfield;
        for(fld=0;fld<Nfield;fld++){
            f_sendQ[sk+fld] = f_Q[nind+fld];
        }
    }
    int cout[nprocs];
    for(n=0;n<nprocs;n++){
        cout[n] = Nfield * mesh->Nfp2procs[n];
    }
    /* send and recv data */
    int Nmess = 0;
    fetch_buffer(procid, nprocs, cout, f_sendQ, f_recvQ, send_requests, recv_requests, &Nmess);
    vector_real_free(f_sendQ);
    return Nmess;
}

/**
 *
 * @param mesh
 * @param Nfield duplicate number on each cell
 * @param f_Q
 * @param f_recv
 * @param Nmess
 */
int dg_mesh_fetch_cell_buffer(dg_mesh *mesh, int Nfield, dg_real *f_Q, dg_real *f_recvQ,
                               MPI_Request *send_requests,
                               MPI_Request *recv_requests){
    const int Nparf = dg_mesh_NfetchFace(mesh);
    int procid,nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    dg_real *f_sendQ = vector_real_create(Nfield*Nparf);
    /* prepare for send data */
    int n,fld;
    for(n=0;n<Nparf;n++){
        const int sk = n*Nfield;
        const int cind = mesh->CBFToK[n]*Nfield;
        for(fld=0;fld<Nfield;fld++){
            f_sendQ[sk+fld] = f_Q[cind+fld];
        }
    }
    int cout[nprocs];
    for(n=0;n<nprocs;n++){
        cout[n] = Nfield * mesh->Nface2procs[n];
    }
    /* send and recv data */
    int Nmess = 0;
    fetch_buffer(procid, nprocs, cout, f_sendQ, f_recvQ, send_requests, recv_requests, &Nmess);
    vector_real_free(f_sendQ);
    return Nmess;
}


static void fetch_buffer(int procid, int nprocs, int *pout, dg_real *send_buffer, dg_real *recv_buffer,
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
                MPI_Isend(send_buffer+sk, Nout, MPI_DG_REAL, p, 5666+p,
                          MPI_COMM_WORLD, mpi_send_requests +Nmess);
                MPI_Irecv(recv_buffer+sk,  Nout, MPI_DG_REAL, p, 5666+procid,
                          MPI_COMM_WORLD,  mpi_recv_requests +Nmess);
                sk+=Nout;
                ++Nmess;
            }
        }
    }
    *Nmessage = Nmess;
}