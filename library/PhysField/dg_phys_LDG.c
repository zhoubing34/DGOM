//
// Created by li12242 on 16/12/30.
//

#include "dg_phys_LDG.h"

dg_phys_LDG* dg_phys_LDG_create(dg_phys_info *info){
    dg_phys_LDG *ldg = (dg_phys_LDG *)calloc(1, sizeof(dg_phys_LDG));

    const int K = dg_grid_K(dg_phys_info_grid(info));
    const int Np = dg_cell_Np(dg_phys_info_cell(info));
    const int Nfield = info->Nfield;
    const int NfetchNode = dg_mesh_NfetchNode(dg_phys_info_mesh(info));

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