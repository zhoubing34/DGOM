//
// Created by li12242 on 17/3/7.
//

#include "dg_region_surfInfo.h"

/**
 * @brief calculate the outward normal vector and jacobi coefficient on faces
 * @param[in] region multi-regions
 */
void dg_reg_surfInfo2d(dg_region *region){
    int k,f;
    dg_cell *cell = dg_region_cell(region);
    const int K = dg_grid_K(region->grid);
    const int Nfaces = dg_cell_Nfaces(cell);

    double **x = region->x;
    double **y = region->y;
    int **Fmask = dg_cell_Fmask(cell);

    double **nx = matrix_double_create(K, Nfaces);
    double **ny = matrix_double_create(K, Nfaces);
    double **sJ = matrix_double_create(K, Nfaces);

    region->nx = nx;
    region->ny = ny;
    region->sJ = sJ;

    for(k=0;k<K;k++){
        for(f=0;f<Nfaces;f++){
            const int Nfp = dg_cell_Nfp(cell)[f];
            double x1 = x[k][Fmask[f][0]];
            double x2 = x[k][Fmask[f][Nfp-1]];
            double y1 = y[k][Fmask[f][0]];
            double y2 = y[k][Fmask[f][Nfp-1]];

            nx[k][f] =  (y2-y1);
            ny[k][f] = -(x2-x1);

            sJ[k][f] = sqrt(nx[k][f]*nx[k][f]+ny[k][f]*ny[k][f]);
            nx[k][f] /= sJ[k][f];
            ny[k][f] /= sJ[k][f];
            sJ[k][f] /= 2.;
        }
    }
    return;
}

void dg_reg_surfInfo3d(dg_region *region){
    return;
}