//
// Created by li12242 on 17/1/5.
//

#include <MultiRegions/mr_grid.h>
#include <MultiRegions/mr_mesh.h>
#include "pf_strong_viscosity_LDG_flux2d_test.h"
#include "pf_test.h"
#include "PhysField/pf_add_LDG_solver.h"
#include "PhysField/pf_fetchBuffer.h"

int phys_strong_viscosity_LDG_flux2d_test(physField *phys, int verbose){
    int fail = 0;
    extern int Nfield;

    parallMesh *mesh = phys->mesh;
    multiReg *region = phys->region;
    stdCell *shape = phys->cell;

    const int K = phys->grid->K;
    const int Np = shape->Np;
    const int nprocs = phys->mesh->nprocs;

    int k,i;
    real miu = 0.5;
    real px_ext[Np*Nfield*K];
    real py_ext[Np*Nfield*K];
    real rhs_ext[Np*Nfield*K];

    pf_add_LDG_solver(phys);

    // assignment
    int sk = 0;
    for(k=0;k<K;k++){
        for(i=0;i<Np;i++){
            real u = region->x[k][i];
            real v = region->y[k][i];

            phys->f_rhsQ[sk] = 0;
            phys->f_Q[sk] = u*u; //
            phys->viscosity->vis_Q[sk] = miu;
            px_ext[sk] = 2*u*miu;
            py_ext[sk] = 0.0*miu;
            rhs_ext[sk++] = 2.0*miu*2;

            phys->f_rhsQ[sk] = 0;
            phys->f_Q[sk] = u*v+v*v; // field 1 of -(dEdx + dGdy)
            phys->viscosity->vis_Q[sk] = miu;
            px_ext[sk] = v*miu;
            py_ext[sk] = (u+2*v)*miu;
            rhs_ext[sk++] = 2.0*miu*2;
        }
    }

    MPI_Request mpi_out_requests[nprocs];
    MPI_Request mpi_in_requests[nprocs];
    int Nmess;
    pf_fetchNodeBuffer2d(phys, mpi_out_requests, mpi_in_requests, &Nmess);
    MPI_Status instatus[nprocs];
    MPI_Waitall(Nmess, mpi_in_requests, instatus);

    double clockT1 = MPI_Wtime();
    pf_strong_viscosity_LDG_flux2d(phys, NULL, NULL, 0, 0, 0);
    double clockT2 = MPI_Wtime();
    pf_strong_viscosity_LDG_flux2d(phys, NULL, NULL, 0, 0, 0);

    if(!phys->mesh->procid){
        Vector_test(__FUNCTION__, phys->viscosity->px_Q, px_ext, Np*Nfield*K, clockT2-clockT1);
        Vector_test(__FUNCTION__, phys->viscosity->py_Q, py_ext, Np*Nfield*K, clockT2-clockT1);
        Vector_test(__FUNCTION__, phys->f_rhsQ, rhs_ext, Np*Nfield*K, clockT2-clockT1);
    }

    if(verbose){
        FILE *fp = CreateLog(__FUNCTION__, mesh->procid, mesh->nprocs);
        fprintf(fp, "K = %d\n", phys->grid->K);
        fprintf(fp, "Nfield = %d\n", phys->Nfield);
        fprintf(fp, "Np = %d\n", phys->cell->Np);
        PrintIntMatrix2File(fp, "EToV", phys->grid->EToV, K, phys->cell->Nv);
        PrintVector2File(fp, "vx", phys->grid->vx, phys->grid->Nv);
        PrintVector2File(fp, "vy", phys->grid->vy, phys->grid->Nv);
        PrintVector2File(fp, "f_Q", phys->f_Q, Nfield*Np*k);
        PrintVector2File(fp, "px_Q", phys->viscosity->px_Q, Nfield*Np*K);
        PrintVector2File(fp, "py_Q", phys->viscosity->py_Q, Nfield*Np*K);
        PrintVector2File(fp, "f_rhsQ", phys->f_rhsQ, Nfield*Np*k);
        fclose(fp);
    }

    pf_delete_LDG_solver(phys);

    return fail;
}