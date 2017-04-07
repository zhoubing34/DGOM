//
// Created by li12242 on 16/12/30.
//

#include "dg_phys_LDG.h"

dg_phys_LDG* dg_phys_LDG_create(dg_phys_info *info){
    dg_phys_LDG *ldg = (dg_phys_LDG *)calloc(1, sizeof(dg_phys_LDG));

    const int K = dg_grid_K(info->grid);
    const int Np = dg_cell_Np(info->cell);
    const int Nfield = info->Nfield;

    ldg->info = info;
    ldg->vissqrt = vector_real_create(Nfield*K*Np);
    ldg->px = vector_real_create(K*Np);
    ldg->py = vector_real_create(K*Np);
    ldg->pz = vector_real_create(K*Np);

    return ldg;
}

void dg_phys_LDG_free(dg_phys_LDG *ldg){

    vector_real_free(ldg->vissqrt);
    vector_real_free(ldg->px);
    vector_real_free(ldg->py);
    vector_real_free(ldg->pz);

    free(ldg);
    return;
}

static void ldg_set_vis(dg_phys_LDG *ldg, dg_real *vis){
    dg_phys_info *info = ldg->info;
    const int K = dg_grid_K(info->grid);
    const int Np = dg_cell_Np(info->cell);
    const int Nfield = info->Nfield;

    int k,n,fld;
    for(k=0;k<K;k++){
        for(n=0;n<Np;n++){
            for(fld=0;fld<Nfield;fld++){
                const int sk = (k*Np+n)*Nfield+fld;
                ldg->vissqrt[sk] = dg_sqrt(vis[sk]);
            }
        }
    }
    return;
}