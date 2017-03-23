//
// Created by li12242 on 12/26/16.
//

#include "dg_phys_strong_vol_opt_test.h"
#include "dg_phys_test.h"

dg_real t = 3.0f;
dg_real s = 4.0f;

static int nodal_flux(dg_real *var, dg_real *Eflux, dg_real *Gflux){
    dg_real u = var[0];
    dg_real v = var[1];

    Eflux[0] = u*v;   Gflux[0] = u+t*v; // field 0
    Eflux[1] = s*u+v;   Gflux[1] = -t*u*v; // field 1

    return 0;
}

int dg_phys_strong_vol_opt2d_test(dg_phys *phys, int verbose){
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
    dg_real rhs_ext[Np*Nfield*K];
    // assignment
    int sk = 0;
    for(k=0;k<K;k++){
        for(i=0;i<Np;i++){
            dg_real u = region->x[k][i];
            dg_real v = region->y[k][i];
            f_Q[sk] = u; // field 0 of -(dEdx + dGdy)
            rhs_ext[sk++] = -(v + t);
            f_Q[sk] = v; // field 1 of -(dEdx + dGdy)
            rhs_ext[sk++] = -(s - t*u);
        }
    }

    double clockT1 = MPI_Wtime();
    dg_phys_strong_vol_opt2d(phys, nodal_flux);
    double clockT2 = MPI_Wtime();

    fail = vector_double_test(__FUNCTION__, f_rhsQ, rhs_ext, Np * Nfield * K);

    if(verbose){
        FILE *fp = create_log(__FUNCTION__, mesh->procid, mesh->nprocs);
        fprintf(fp, "Nfield = %d\n", Nfield);
        fprintf(fp, "Np = %d\n", Np);
        print_double_vector2file(fp, "f_Q", f_Q, Nfield*Np*k);
        print_double_vector2file(fp, "f_rhsQ", f_rhsQ, Nfield*Np*k);
    }

    const int procid = region->procid;
    if(!procid) {
        if(!fail) printf(HEADPASS "1 test passed from %s\n", __FUNCTION__);
    }
    return fail;
}