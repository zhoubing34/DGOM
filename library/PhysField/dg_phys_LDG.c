//
// Created by li12242 on 16/12/30.
//

#include "dg_phys_LDG.h"

dg_phys_LDG* dg_phys_LDG_create(dg_phys_info *info){
    dg_phys_LDG *ldg = (dg_phys_LDG *)calloc(1, sizeof(dg_phys_LDG));

    const int K = dg_grid_K(info->grid);
    const int Np = dg_cell_Np(info->cell);
    const int Nfield = info->Nfield;
    const int NfetchNode = dg_mesh_NfetchNode(info->mesh);

    ldg->info = info;
    ldg->miux = vector_real_create(K*Np*Nfield);
    ldg->miuy = vector_real_create(K*Np*Nfield);
    ldg->miuz = vector_real_create(K*Np*Nfield);

    ldg->px = vector_real_create(K*Np*Nfield);
    ldg->py = vector_real_create(K*Np*Nfield);
    ldg->pz = vector_real_create(K*Np*Nfield);

    ldg->px_recv = vector_real_create(NfetchNode*Nfield);
    ldg->py_recv = vector_real_create(NfetchNode*Nfield);
    ldg->pz_recv = vector_real_create(NfetchNode*Nfield);

    return ldg;
}

void dg_phys_LDG_free(dg_phys_LDG *ldg){

    vector_real_free(ldg->miux);
    vector_real_free(ldg->miuy);
    vector_real_free(ldg->miuz);
    vector_real_free(ldg->px);
    vector_real_free(ldg->py);
    vector_real_free(ldg->pz);
    vector_real_free(ldg->px_recv);
    vector_real_free(ldg->py_recv);
    vector_real_free(ldg->pz_recv);

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
                ldg->miux[sk] = dg_sqrt(vis[sk]);
            }
        }
    }
    return;
}