//
// Created by li12242 on 12/27/16.
//

#include "pf_strong_surface_flux2d_test.h"
#include "PhysField/pf_strong_surface_flux2d.h"

static int nodal_flux(dg_real *var, dg_real *Eflux, dg_real *Gflux){
    dg_real u = var[0];
    dg_real v = var[1];

    Eflux[0] = u;   Gflux[0] = v; // field 0
    Eflux[1] = u;   Gflux[1] = v; // field 1

    return 0;
}

static int numerical_flux(dg_real nx, dg_real ny, dg_real *varM, dg_real *varP, dg_real *Fhs){
    dg_real uM = varM[0], uP = varP[0];
    dg_real vM = varM[1], vP = varP[1];

    Fhs[0] = 0;//uM*nx + vM*ny;
    Fhs[1] = 0;//uM*nx + vM*ny;

    return 0;
}

static int wall_func(dg_real nx, dg_real ny, dg_real *varM, dg_real *varP){
    varP[0] = varM[0];
    varP[1] = varM[1];
    return 0;
}

int phys_strong_surface_flux2d_test(physField *phys, int verbose){
    int fail = 0;

    parallMesh *mesh = phys->mesh;
    multiReg *region = phys->region;
    dg_cell *shape = phys->cell;

    const int K = phys->grid->K;
    const int Np = shape->Np;
    const int Nfield = phys->Nfield;

    int k,i;
    // initialize and assignment
    int sk = 0;
    for(k=0;k<K;k++){
        for(i=0;i<Np;i++){
            dg_real u = region->x[k][i];
            dg_real v = region->y[k][i];
            phys->f_rhsQ[sk] = 0.0;
            phys->f_Q[sk++] = u;
            phys->f_rhsQ[sk] = 0.0;
            phys->f_Q[sk++] = v;
        }
    }

    pf_strong_surface_flux2d(phys, wall_func, wall_func, nodal_flux, numerical_flux);

    if(verbose){
        FILE *fp = create_log(__FUNCTION__, mesh->procid, mesh->nprocs);
        fprintf(fp, "K = %d\n", phys->grid->K);
        fprintf(fp, "Nfield = %d\n", phys->Nfield);
        fprintf(fp, "Np = %d\n", phys->cell->Np);
        print_double_vector2file(fp, "f_Q", phys->f_Q, Nfield * Np * k);
        print_double_vector2file(fp, "f_rhsQ", phys->f_rhsQ, Nfield * Np * k);
    }

    const int procid = region->procid;
    if(!procid) {
        if(!fail) printf(HEADPASS "1 test passed from %s\n", __FUNCTION__);
    }
    return fail;
}