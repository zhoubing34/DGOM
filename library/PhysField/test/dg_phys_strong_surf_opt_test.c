//
// Created by li12242 on 12/27/16.
//

#include "dg_phys_strong_surf_opt_test.h"

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

static int obc_func(dg_real nx, dg_real ny, dg_real *f_M, dg_real *f_ext, int obc_ind, dg_real *f_P){
    f_P[0] = f_M[0];
    f_P[1] = f_M[1];
    return 0;
}

int dg_phys_strong_surf_opt2d_test(dg_phys *phys, int verbose){
    int fail = 0;

    dg_mesh *mesh = dg_phys_mesh(phys);
    dg_region *region = dg_phys_region(phys);
    dg_grid *grid = dg_phys_grid(phys);
    dg_cell *cell = dg_phys_cell(phys);

    const int K = dg_grid_K(grid);
    const int Np = dg_cell_Np(cell);
    const int Nfield = dg_phys_Nfield(phys);

    dg_real *f_rhsQ = dg_phys_f_rhsQ(phys);
    dg_real *f_Q = dg_phys_f_Q(phys);

    int k,i;
    // initialize and assignment
    for(k=0;k<K;k++){
        for(i=0;i<Np;i++){
            dg_real u = region->x[k][i];
            dg_real v = region->y[k][i];
            int sk = k*Np*Nfield + i*Nfield;

            f_rhsQ[sk] = 0.0;
            f_Q[sk] = u;
            f_rhsQ[sk+1] = 0.0;
            f_Q[sk+1] = v;
        }
    }

    dg_phys_strong_surf_opt2d(phys, wall_func, wall_func, obc_func, nodal_flux, numerical_flux);

    const int procid = dg_phys_procid(phys);
    const int nprocs = dg_phys_nprocs(phys);
    if(verbose){
        FILE *fp = create_log(__FUNCTION__, procid, nprocs);
        fprintf(fp, "K = %d\n", K);
        fprintf(fp, "Nfield = %d\n", Nfield);
        fprintf(fp, "Np = %d\n", Np);
        print_double_vector2file(fp, "f_Q", f_Q, Nfield * Np * k);
        print_double_vector2file(fp, "f_rhsQ", f_rhsQ, Nfield * Np * k);
    }

    if(!procid) {
        if(!fail) printf(HEADPASS "1 test passed from %s\n", __FUNCTION__);
    }
    return fail;
}