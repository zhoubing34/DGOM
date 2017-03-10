//
// Created by li12242 on 17/1/9.
//

#include "pf_cellFetch_test.h"

int phys_cellFetch_test(dg_phys *phys, int verbose){
    int fail = 0;

    dg_grid *grid = phys->grid;
    dg_mesh *mesh = phys->mesh;

    const int Nfield = phys->Nfield;
    const int K = phys->grid->K;
    const int Nfaces = phys->cell->Nfaces;
    const int procid = phys->mesh->procid;
    const int nprocs = phys->mesh->nprocs;

    int k,f,n,sk=0;

    double *vx = grid->vx;
    double *vy = grid->vy;
    double par_coor[phys->parallCellNum];
    dg_real *c_Q = phys->c_Q;

    for(k=0;k<K;k++){
        for(f=0;f<Nfaces;f++){
            if( mesh->EToP[k][f] != procid){
                int v1 = grid->EToV[k][f];
                int v2 = grid->EToV[k][(f+1)%Nfaces];
                c_Q[k*Nfield + 0] = (vx[v1] + vx[v2])/2;
                c_Q[k*Nfield + 1] = (vy[v1] + vy[v2])/2;
            }
        }
    }

    MPI_Request mpi_out_requests[nprocs], mpi_in_requests[nprocs];
    int Nmess;

    /* do sends and recv */
    pf_fetchCellBuffer(phys, mpi_out_requests, mpi_in_requests, &Nmess);

    MPI_Status instatus[nprocs];
    MPI_Waitall(Nmess, mpi_in_requests, instatus);
    MPI_Waitall(Nmess, mpi_out_requests, instatus);

    sk = 0;
    for(n=0;n<mesh->Nparf;n++){
        k=mesh->Pcid_recv[n];
        f=mesh->Parfaceid[n];

        int v1 = grid->EToV[k][f];
        int v2 = grid->EToV[k][(f+1)%Nfaces];
        par_coor[sk++] = (vx[v1] + vx[v2])/2;
        par_coor[sk++] = (vy[v1] + vy[v2])/2;
    }

    if(!procid)
        fail = vector_double_test(__FUNCTION__, par_coor, phys->c_inQ, phys->parallCellNum);

    if(verbose){
        FILE *fp = create_log(__FUNCTION__, mesh->procid, mesh->nprocs);
        fprintf(fp, "Nparf = %d\n", phys->parallCellNum);
        print_int_vector2file(fp, "mesh->Pcid_recv", mesh->Pcid_recv, mesh->Nparf);
        print_int_vector2file(fp, "mesh->Parfaceid", mesh->Parfaceid, mesh->Nparf);
        print_double_vector2file(fp, "c_Q", phys->c_Q, K * Nfield);
        print_int_vector2file(fp, "Parcellid", phys->cellIndexOut, phys->parallCellNum);
        print_double_vector2file(fp, "c_inQ", phys->c_inQ, phys->parallCellNum);
        print_double_vector2file(fp, "c_outQ", phys->c_outQ, phys->parallCellNum);
    }

    if(!procid) {
        if(!fail) printf(HEADPASS "1 test passed from %s\n", __FUNCTION__);
    }
    return fail;
}
