//
// Created by li12242 on 12/21/16.
//
#include "dg_mesh_test.h"

int dg_mesh_cell_fetch_buffer_test(dg_mesh *mesh, int verbose){

    int procid,nprocs,fail = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    const dg_grid *grid = dg_mesh_grid(mesh);

    const int K = dg_grid_K(grid);
    const int Nfield = 1;

    const int Nparf = mesh->NfetchFace;
    int **EToE = grid->EToE;
    /* allocation and assignment */
    dg_real *c_Q = vector_real_create(Nfield*K);
    dg_real *c_recvQ = vector_real_create(Nfield*Nparf);
    dg_real *c_extQ = vector_real_create(Nfield*Nparf);
    int k;
    for(k=0;k<K;k++){
        c_Q[k] = (dg_real)k;
    }
    MPI_Request mpi_send_requests[nprocs], mpi_recv_requests[nprocs];
    int Nmess = dg_mesh_fetch_cell_buffer(mesh, Nfield, c_Q, c_recvQ, mpi_send_requests, mpi_recv_requests);

    int n;
    for(n=0;n<Nparf;n++){
        int k1 = mesh->CBFToK[n];
        int f1 = mesh->CBFToF[n];
        c_extQ[n] = EToE[k1][f1];
    }

    MPI_Status instatus[nprocs];
    MPI_Waitall(Nmess, mpi_recv_requests, instatus);
    MPI_Waitall(Nmess, mpi_send_requests, instatus);

    fail = vector_double_test(__FUNCTION__, c_recvQ, c_extQ, Nfield*Nparf);

    if(verbose){
        FILE *fp = create_log(__FUNCTION__, procid, nprocs);
        print_double_vector2file(fp, "c_recvQ", c_recvQ, Nparf*Nfield);
        print_double_vector2file(fp, "c_extQ", c_extQ, Nparf*Nfield);
        fclose(fp);
    }
    vector_real_free(c_Q);
    vector_real_free(c_recvQ);
    vector_real_free(c_extQ);

    if(!procid) printf(HEADPASS "1 test passed from %s\n", __FUNCTION__);
    return fail;
}

int dg_mesh_node_fetch_buffer_test(dg_mesh *mesh, int verbose){
    int procid,nprocs,fail = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);

    const int Np = dg_cell_Np(dg_mesh_cell(mesh));
    const int K = dg_grid_K(dg_mesh_grid(mesh));

    const int Nfield = 2;
    const int Nparn = mesh->NfetchNode;

    double **x = mesh->region->x;
    double **y = mesh->region->y;

    /* allocation and assignment */
    dg_real *f_Q = vector_real_create(Nfield*K*Np);
    dg_real *f_recvQ = vector_real_create(Nfield*Nparn);
    dg_real *f_extQ = vector_real_create(Nfield*Nparn);
    int k,n,sk=0;
    for(k=0;k<K;k++){
        for(n=0;n<Np;n++){
            f_Q[sk++] = x[k][n];
            f_Q[sk++] = y[k][n];
        }
    }
    MPI_Request mpi_send_requests[nprocs], mpi_recv_requests[nprocs];
    int Nmess = dg_mesh_fetch_node_buffer(mesh, Nfield, f_Q, f_recvQ, mpi_send_requests, mpi_recv_requests);
    sk = 0;
    for(n=0;n<Nparn;n++){
        int ind = mesh->NBFToN[n];
        f_extQ[sk++] = f_Q[ind*Nfield];
        f_extQ[sk++] = f_Q[ind*Nfield+1];
    }

    MPI_Status instatus[nprocs];
    MPI_Waitall(Nmess, mpi_recv_requests, instatus);
    MPI_Waitall(Nmess, mpi_send_requests, instatus);

    fail = vector_double_test(__FUNCTION__, f_recvQ, f_extQ, Nfield*Nparn);
    if(verbose){
        FILE *fp = create_log(__FUNCTION__, procid, nprocs);
        print_double_vector2file(fp, "f_recvQ", f_recvQ, Nparn*Nfield);
        print_double_vector2file(fp, "f_extQ", f_extQ, Nparn*Nfield);
        fclose(fp);
    }

    vector_real_free(f_Q);
    vector_real_free(f_recvQ);
    vector_real_free(f_extQ);
    if(!procid) printf(HEADPASS "1 test passed from %s\n", __FUNCTION__);
    return fail;
}

int dg_mesh_parallel_test(dg_mesh *mesh, int verbose){
    int procid,nprocs,fail = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    if(verbose){
        FILE *fp = create_log(__FUNCTION__, procid, nprocs);
        fprintf(fp, "mesh->NfetchFace: %d\n", mesh->NfetchFace);
        print_int_vector2file(fp, "mesh->Nface2procs", mesh->Nface2procs, nprocs);
        print_int_vector2file(fp, "mesh->CBFToK", mesh->CBFToK, mesh->NfetchFace);
        print_int_vector2file(fp, "mesh->CBFToF", mesh->CBFToF, mesh->NfetchFace);
        fprintf(fp, "mesh->NfetchNode: %d\n", mesh->NfetchNode);
        print_int_vector2file(fp, "mesh->Nfp2procs", mesh->Nfp2procs, nprocs);
        print_int_vector2file(fp, "mesh->NBFToN", mesh->NBFToN, mesh->NfetchNode);
        fclose(fp);
    }
    if(!procid) printf(HEADPASS "1 test passed from %s\n", __FUNCTION__);
    return fail;
}
