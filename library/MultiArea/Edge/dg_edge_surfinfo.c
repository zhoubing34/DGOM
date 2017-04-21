//
// Created by li12242 on 17/3/10.
//

#include "dg_edge_surfinfo.h"

void dg_edge_surfinfo2d(dg_edge *edge){

    dg_cell *cell = dg_edge_cell(edge);
    dg_region *region = dg_edge_region(edge);
    const int Nnode = dg_edge_Nnode(edge);
    const int Nedge = dg_edge_Nedge(edge);
    dg_real *fsc = vector_real_create(Nnode);
    dg_real *nxe = vector_real_create(Nnode);
    dg_real *nye = vector_real_create(Nnode);

    int **Fmask = dg_cell_Fmask(cell);
    double **J = dg_region_J(region);
    double **sJ = dg_region_sJ(region);
    double **nx = dg_region_nx(region);
    double **ny = dg_region_ny(region);

    int f,n,sk=0;
    for(f=0;f<Nedge;f++){
        const int k1 = edge->varkM[f];
        const int f1 = edge->varfM[f];
        const double isJ = sJ[k1][f1];
        const double nxi = nx[k1][f1];
        const double nyi = ny[k1][f1];
        const int Nfp = dg_cell_Nfp(cell)[f1];
        for(n=0;n<Nfp;n++){
            nxe[sk] = nxi;
            nye[sk] = nyi;
            fsc[sk++] = isJ/J[k1][ Fmask[f1][n] ];
        }
    }

    /* assignment */
    edge->fsc = fsc;
    edge->nx = nxe;
    edge->ny = nye;
    edge->nz = NULL;
    return;
}