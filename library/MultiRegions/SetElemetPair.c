#include "mr_reg.h"

/**
 * @brief
 * Build the parallel element connection.
 *
 * @details
 * EToE is modified inside the function. EToE contains the mapping from elements to incoming elements information,
 * which is marked as negative.
 *
 * @param[in] shape
 * @param[in] mesh
 *
 * @return
 * return values:
 * name     | type     | description of value
 * -------- |----------|----------------------
 * elemapOUT| int*     |
 *
 */
void SetElementPair(stdCell *shape, MultiReg2d *mesh, int *cellIndOut){

    int k, sk = 0;
    int p2, f1;
    /* build map from f_inE to element */
    for(p2=0;p2<mesh->nprocs;p2++){
        if(mesh->Npar[p2]){
            for(k=0;k<mesh->K;k++){
                for(f1=0;f1<shape->Nfaces;f1++){
                    if(mesh->EToP[k][f1]==p2 && p2!=mesh->procid)
                        cellIndOut[sk++] = k;
                }
            }
        }
    }

    int parcnt=-1;
    for(p2=0;p2<mesh->nprocs;p2++){
        for(k=0;k<mesh->K;k++){
            for(f1=0;f1<shape->Nfaces;f1++){
                if(mesh->EToP[k][f1]==p2 && p2!=mesh->procid) {
                    mesh->EToE[k][f1] = parcnt;
                    --parcnt;
                }
            }
        }
    }
}