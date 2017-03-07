//
// Created by li12242 on 17/1/10.
//

#include <MultiRegions/Region/mr_reg.h>
#include "pf_limiter_test.h"

int phys_limiter_test(physField *phys, int verbose){

    int fail = 0;

    dg_region *region = phys->region;
    const int K = phys->grid->K;
    const int Np = phys->cell->Np;

    int k,n;

    // assignment
    int sk = 0;
    for(k=0;k<K;k++){
        for(n=0;n<Np;n++){
            phys->f_Q[sk++] = region->x[k][n];
            phys->f_Q[sk++] = region->y[k][n];
        }
    }
    pf_limit_BJ2d(phys, 1.0);
    if(verbose) {
        FILE *fp = create_log(__FUNCTION__, phys->mesh->procid, phys->mesh->nprocs);
        print_double_vector2file(fp, "f_Q", phys->f_Q, K * Np * phys->Nfield);
        print_double_vector2file(fp, "c_Q", phys->c_Q, K * phys->Nfield);
        fclose(fp);
    }

    const int procid = region->procid;
    if(!procid) {
        if(!fail) printf(HEADPASS "1 test passed from %s\n", __FUNCTION__);
    }
    return fail;
}