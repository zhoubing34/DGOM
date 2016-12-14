//
// Created by li12242 on 16/12/14.
//

#include <MultiRegions/MultiRegions.h>
#include <StdRegions/StdRegions.h>
#include "MultiRegBC2d.h"

/**
 * @brief create boundary type object
 * @details
 * @param [in] mesh
 * @param [in] SFoV surface to vertex
 *
 */

#define I_INNERLOC 0
#define I_INNERBS  1
#define I_OPENBS   2
#define I_SLIPWALL   3
#define I_NSLIPWALL  4


// sort numbers from small to large
int cmpVert(const void *a, const void *b){
    return *(int*)a-*(int*)b;
}

MultiRegBC2d* MultiRegBC2d_create(MultiReg2d *mesh, int Nsurf, int **SFToV){
    StdRegions2d *shape = mesh->stdcell;

    MultiRegBC2d *bc2d = (MultiRegBC2d *) malloc(sizeof(MultiRegBC2d));

    /* allocate EToBS */
    bc2d->EToBS = IntMatrix_create(mesh->K, shape->Nfaces);

    /* initialize EToBS */
    int k,f1,f2;
    int t[2], v[2];
    for(k=0;k<mesh->K;k++){
        for(f2=0; f2<shape->Nfaces; f2++){ // loop through all element surface
            /* the inner boundary surface */
            if(mesh->EToP[k][f2] != mesh->procid){
                bc2d->EToBS[k][f2] = INNERBS;
                continue;
            }

            /* get the usr-defined boundary */
            int n1 = f2;
            int n2 = (f2+1)%shape->Nv;
            t[0] = mesh->EToV[k][n1];
            t[1] = mesh->EToV[k][n2];
            qsort(t, 2, sizeof(int), cmpVert);

            int t_temp = t[0]*mesh->Nv + t[1];

            for(f1=0; f1<Nsurf; f1++){
                v[0] = SFToV[f1][0];
                v[1] = SFToV[f1][1];
                qsort(v, 2, sizeof(int), cmpVert);

                int v_temp = v[0]*mesh->Nv + v[1];

                if(t_temp == v_temp){
                    switch (SFToV[f1][2]){
                        case I_OPENBS: bc2d->EToBS[k][f2] = OPENBS; break;
                        case I_SLIPWALL: bc2d->EToBS[k][f2] = SLIPWALL; break;
                        case I_NSLIPWALL: bc2d->EToBS[k][f2] = NSLIPWALL; break;
                        default: printf("MultiRegBC2d_create: Unknown boundary type SFToV[%d][3] = %d\n",
                                        f1,SFToV[f1][2]);
                            printf("The boundary type must be one of\n");
                            printf("   %d - open boundary surface[default]\n", I_OPENBS);
                            printf("   %d - slip wall\n", I_SLIPWALL);
                            printf("   %d - non-slip wall\n", I_NSLIPWALL);
                    }
                    continue; // finish this loop
                }
            }

            /* inner surface */
            bc2d->EToBS[k][f2] = INNERLOC;
        }
    }

    return bc2d;
}

void MultiRegBC2d_free(MultiRegBC2d* bc2d){
    IntMatrix_free(bc2d->EToBS);
    free(bc2d);
}