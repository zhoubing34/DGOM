//
// Created by li12242 on 17/1/10.
//

#include "dg_phys_limiter_test.h"

int dg_phys_limiter_test(dg_phys *phys, int verbose){

    int fail = 0;

    dg_region *region = dg_phys_region(phys);
    const int K = dg_grid_K(dg_phys_grid(phys));
    const int Np = dg_cell_Np(dg_phys_cell(phys));
    const int Nfield = dg_phys_Nfield(phys);

    double *f_Q = dg_phys_f_Q(phys);

    int k,n;
    // assignment
    int sk = 0;
    for(k=0;k<K;k++){
        for(n=0;n<Np;n++){
            f_Q[sk++] = dg_region_x(region)[k][n];
            f_Q[sk++] = dg_region_y(region)[k][n];
        }
    }
    phys->limit(phys, 1.0);
    const int procid = dg_phys_procid(phys);
    const int nprocs = dg_phys_nprocs(phys);
    if(verbose) {
        FILE *fp = create_log(__FUNCTION__, procid, nprocs);
        print_double_vector2file(fp, "f_Q", f_Q, K*Np*Nfield);
        fclose(fp);
    }

    if(!procid) {
        if(!fail) printf(HEADPASS "1 test passed from %s\n", __FUNCTION__);
    }
    return fail;
}