//
// Created by li12242 on 16/12/16.
//

#include "pf_cellMean_test.h"
#include "PhysField/pf_cellMean.h"

int phys_cellMean_test(physField *phys, int verbose){

    // local variable
    int fail = 0;

    dg_region *region = phys->region;
    dg_mesh *mesh = phys->mesh;

    const int K = phys->grid->K;
    const int Np = phys->cell->Np;

    int k,i;
    // assignment
    int sk = 0;
    for(k=0;k<K;k++){
        for(i=0;i<Np;i++){
            phys->f_Q[sk++] = region->x[k][i];
            phys->f_Q[sk++] = region->y[k][i];
        }
    }
    pf_cellMean(phys);

    if(verbose) {
        FILE *fp = create_log(__FUNCTION__, mesh->procid, mesh->nprocs);
        print_double_vector2file(fp, "f_Q", phys->f_Q, K * Np * phys->Nfield);
        print_double_vector2file(fp, "c_Q", phys->c_Q, K * phys->Nfield);
        fclose(fp);
    }

    const int procid = mesh->procid;
    if(!procid) {
        if(!fail) printf(HEADPASS "1 test passed from %s\n", __FUNCTION__);
    }
    return fail;
}
