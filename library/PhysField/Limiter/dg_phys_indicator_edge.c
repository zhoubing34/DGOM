//
// Created by li12242 on 17/3/29.
//

#include "dg_phys_indicator_edge.h"

#define EXCEPTION(u, u1, u2) ( ((u-u1)>EPS)& ((u-u2)>EPS) )||( ((u1-u)>EPS)&((u2-u)>EPS))


static void dg_phys_local_face_mean(dg_phys_info *phys_info, dg_real *f_mean);

/**
 * @brief
 * @param phys
 * @param tind
 */
void dg_phys_indicator_edge(dg_phys_info *phys, int *tind){

    dg_grid *grid = phys->grid;
    const int K = dg_grid_K(grid);
    const int Nfield = phys->Nfield;
    const int Nfaces = dg_cell_Nfaces(phys->cell);
    const int procid = dg_grid_procid(grid);

    dg_real f_mean[K*Nfaces*Nfield];
    dg_phys_local_face_mean(phys, f_mean);

    register int k,n,f,fld;
    for(k=0;k<K;k++){
        dg_real *c_mean = phys->c_Q + k*Nfield; // mean value of k-th cell
        int **EToE = dg_grid_EToE(grid);
        int **EToP = dg_grid_EToP(grid);

        for(f=0;f<Nfaces;f++){
            int e = EToE[k][f];
            int p = EToP[k][f];
            // if adjacent cell is on other process, jump this cycle
            if(p!=procid) continue;

            int sf = k*Nfaces*Nfield + f*Nfield;
            for(fld=0;fld<Nfield;fld++){
                dg_real c_next = phys->c_Q[e*Nfield+fld];
                if( EXCEPTION(f_mean[sf+fld], c_next, c_mean[fld]) ) {
                    tind[k*Nfield + fld] = 1;
                }
            }
        }
    }

    /* parallel cell loop */
    dg_mesh *mesh = phys->mesh;
    const int Nfetchfaces = dg_mesh_NfetchFace(mesh);
    for(n=0;n<Nfetchfaces;n++){
        k = mesh->CBFToK[n];
        f = mesh->CBFToF[n];
        for(fld=0;fld<Nfield;fld++){
            dg_real c_mean = phys->c_Q[k*Nfield+fld];
            dg_real c_next = phys->c_recvQ[n*Nfield+fld];

            if( EXCEPTION(f_mean[k*Nfaces*Nfield+f*Nfield+fld], c_next, c_mean) ) {
                tind[k*Nfield + fld] = 1;
            }
        }
    }
    return;
}

/**
 * @brief
 * @param phys_info
 * @param f_mean
 */
static void dg_phys_local_face_mean(dg_phys_info *phys_info, dg_real *f_mean){
    dg_cell *cell = phys_info->cell;
    dg_region *region = phys_info->region;
    const int Nfield = phys_info->Nfield;
    const int Nfaces = dg_cell_Nfaces(cell);
    const int Np = dg_cell_Np(cell);
    const int K = dg_grid_K(phys_info->grid);

    dg_real *f_Q = phys_info->f_Q;

    register int k,f,fld,sk=0;
    for(k=0;k<K;k++){
        region->face_integral(region, Nfield, k, f_Q+k*Np*Nfield, f_mean+k*Nfaces*Nfield);
        for(f=0;f<Nfaces;f++){
            double ds = region->face_size[k][f];
            for(fld=0;fld<Nfield;fld++)
                f_mean[sk++] /= ds;
        }
    }
    return;
}