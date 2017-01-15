//
// Created by li12242 on 17/1/9.
//

#include <StandCell/sc_stdcell.h>
#include <MultiRegions/mr_mesh.h>
#include "pf_cellFetch_test.h"
#include "PhysField/pf_fetchBuffer.h"


int phys_cellFetch_test(physField *phys, int verbose, char *message, char *filename){
    int fail = 0;

    geoGrid *grid = phys->grid;
    parallMesh *mesh = phys->mesh;

    const int Nfield = phys->Nfield;
    const int K = phys->grid->K;
    const int Nfaces = phys->cell->Nfaces;
    const int procid = phys->mesh->procid;
    const int nprocs = phys->mesh->nprocs;

    int k,f,n,sk=0;

    double *vx = grid->vx;
    double *vy = grid->vy;
    double par_coor[phys->parallCellNum];
    real *c_Q = phys->c_Q;

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
    for(n=0;n<mesh->parallCellNum;n++){
        k=mesh->cellIndexIn[n];
        f=mesh->faceIndexIn[n];

        int v1 = grid->EToV[k][f];
        int v2 = grid->EToV[k][(f+1)%Nfaces];
        par_coor[sk++] = (vx[v1] + vx[v2])/2;
        par_coor[sk++] = (vy[v1] + vy[v2])/2;
    }

    if(!procid)
        fail = Vector_test(message, par_coor, phys->c_inQ, phys->parallCellNum, 0);

    if(verbose){
        FILE *fp = CreateLog(filename, mesh->procid, mesh->nprocs);
        fprintf(fp, "parallCellNum = %d\n", phys->parallCellNum);
        PrintIntVector2File(fp, "mesh->cellIndexIn", mesh->cellIndexIn, mesh->parallCellNum);
        PrintIntVector2File(fp, "mesh->faceIndexIn", mesh->faceIndexIn, mesh->parallCellNum);
        PrintVector2File(fp, "c_Q", phys->c_Q, K*Nfield);
        PrintIntVector2File(fp, "cellIndexOut", phys->cellIndexOut, phys->parallCellNum);
        PrintVector2File(fp, "c_inQ", phys->c_inQ, phys->parallCellNum);
        PrintVector2File(fp, "c_outQ", phys->c_outQ, phys->parallCellNum);
    }

    return fail;
}
