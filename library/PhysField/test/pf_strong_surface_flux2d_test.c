//
// Created by li12242 on 12/27/16.
//

#include "pf_strong_surface_flux2d_test.h"
#include "PhysField/pf_strong_surface_flux2d.h"

static int nodal_flux(real *var, real *Eflux, real *Gflux){
    real u = var[0];
    real v = var[1];

    Eflux[0] = u;   Gflux[0] = v; // field 0
    Eflux[1] = u;   Gflux[1] = v; // field 1

    return 0;
}

static int numerical_flux(real nx, real ny, real *varM, real *varP, real *Fhs){
    real uM = varM[0], uP = varP[0];
    real vM = varM[1], vP = varP[1];

    Fhs[0] = 0;//uM*nx + vM*ny;
    Fhs[1] = 0;//uM*nx + vM*ny;

    return 0;
}

int phys_strong_surface_flux2d_test(physField *phys, int verbose, char *message, char *filename){
    int fail = 0;

    parallMesh *mesh = phys->mesh;
    multiReg *region = phys->region;
    stdCell *shape = phys->cell;

    const int K = phys->grid->K;
    const int Np = shape->Np;
    const int Nfield = phys->Nfield;

    int k,i;
    // initialize and assignment
    int sk = 0;
    for(k=0;k<K;k++){
        for(i=0;i<Np;i++){
            real u = region->x[k][i];
            real v = region->y[k][i];
            phys->f_rhsQ[sk] = 0.0;
            phys->f_Q[sk++] = u;
            phys->f_rhsQ[sk] = 0.0;
            phys->f_Q[sk++] = v;
        }
    }

    double clockT1 = MPI_Wtime();
    pf_strong_surface_flux2d(phys, NULL, NULL, nodal_flux, numerical_flux);
    double clockT2 = MPI_Wtime();

    if(verbose){
        FILE *fp = CreateLog(filename, mesh->procid, mesh->nprocs);
        fprintf(fp, "K = %d\n", phys->grid->K);
        fprintf(fp, "Nfield = %d\n", phys->Nfield);
        fprintf(fp, "Np = %d\n", phys->cell->Np);
        PrintVector2File(fp, "f_Q", phys->f_Q, Nfield*Np*k);
        PrintVector2File(fp, "f_rhsQ", phys->f_rhsQ, Nfield*Np*k);
    }

    return fail;
}