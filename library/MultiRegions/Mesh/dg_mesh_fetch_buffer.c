//
// Created by li12242 on 17/3/9.
//
#include "dg_mesh.h"
#include "dg_mesh_fetch_buffer.h"

static void fetch_buffer(int procid, int nprocs, int *pout, dg_real *send_buffer, dg_real *recv_buffer,
                         MPI_Request *mpi_send_requests,
                         MPI_Request *mpi_recv_requests,
                         int *Nmessage);

void dg_mesh_fetch_node_buffer(dg_mesh *mesh, int Nfield, dg_real *f_Q, dg_real *f_recvQ,
                               MPI_Request *send_requests,
                               MPI_Request *recv_requests,
                               int *Nmess){
    const int Nparn = dg_mesh_Nparn(mesh);
    const int nprocs = dg_mesh_nprocs(mesh);
    dg_real *f_sendQ = vector_double_create(Nfield*Nparn);
    /* prepare for send data */
    int n,fld;
    for(n=0;n<Nparn;n++){
        const int sk = n*Nfield;
        const int nind = mesh->parnode[n]*Nfield;
        for(fld=0;fld<Nfield;fld++){
            f_sendQ[sk+fld] = f_Q[nind+fld];
        }
    }
    int cout[nprocs];
    for(n=0;n<nprocs;n++){
        cout[n] = Nfield * mesh->parnodeNum[n];
    }
    /* send and recv data */
    fetch_buffer(mesh->procid, nprocs, cout, f_sendQ, f_recvQ, send_requests, recv_requests, Nmess);
    vector_real_free(f_sendQ);
    return;
}

/**
 *
 * @param mesh
 * @param Nfield duplicate number on each cell
 * @param f_Q
 * @param f_recv
 * @param Nmess
 */
void dg_mesh_fetch_cell_buffer(dg_mesh *mesh, int Nfield, dg_real *f_Q, dg_real *f_recvQ,
                               MPI_Request *send_requests,
                               MPI_Request *recv_requests,
                               int *Nmess){
    const int Nparf = dg_mesh_Nparf(mesh);
    const int nprocs = dg_mesh_nprocs(mesh);
    dg_real *f_sendQ = vector_real_create(Nfield*Nparf);
    /* prepare for send data */
    int n,fld;
    for(n=0;n<Nparf;n++){
        const int sk = n*Nfield;
        const int cind = mesh->parcell[n]*Nfield;
        for(fld=0;fld<Nfield;fld++){
            f_sendQ[sk+fld] = f_Q[cind+fld];
        }
    }
    int cout[nprocs];
    for(n=0;n<nprocs;n++){
        cout[n] = Nfield * mesh->parfaceNum[n];
    }
    /* send and recv data */
    fetch_buffer(mesh->procid, nprocs, cout, f_sendQ, f_recvQ, send_requests, recv_requests, Nmess);
    vector_real_free(f_sendQ);
    return;
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