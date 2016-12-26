//
// Created by li12242 on 12/26/16.
//

#include <MultiRegions/mr_mesh.h>
#include "phys_strong_volume_flux2d_test.h"
#include "phys_test.h"
#include "PhysField/phys_strong_volume_flux2d.h"

int nodal_flux(real *var, real *Eflux, real *Gflux){
    real u = var[0];
    real v = var[1];

    Eflux[0] = u*v;   Gflux[0] = u-3*v; // field 0
    Eflux[1] = u+v;   Gflux[1] = -u*v; // field 1

    return 0;
}

int phys_strong_volume_flux2d_test(physField *phys, int verbose, char *message, char *filename){
    int fail = 0;
    extern int Nfield;

    parallMesh *mesh = phys->mesh;
    multiReg *region = phys->region;
    stdCell *shape = phys->cell;

    const int K = phys->grid->K;
    const int Np = shape->Np;

    int k,i;
    real rhs_ext[Np*Nfield*K];
    // assignment
    int sk = 0;
    for(k=0;k<K;k++){
        for(i=0;i<Np;i++){
            real u = region->x[k][i];
            real v = region->y[k][i];
            phys->f_Q[sk] = u; // field 0 of -(dEdx + dGdy)
            rhs_ext[sk++] = -(v - 3);
            phys->f_Q[sk] = v; // field 1 of -(dEdx + dGdy)
            rhs_ext[sk++] = -(1 - u);
        }
    }

    double clockT1 = MPI_Wtime();
    phys_strong_volume_flux2d(phys, nodal_flux);
    double clockT2 = MPI_Wtime();

    if(!phys->mesh->procid)
        Vector_test(message, phys->f_rhsQ, rhs_ext, Np*Nfield*K, clockT2-clockT1);

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