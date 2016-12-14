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
#define I_SLIPWALL   2
#define I_NSLIPWALL  3
#define I_OPENBS   4

OBC2d* OBC2d_create(int Nsurf, int **SFToV, int typeid);

// sort numbers from small to large
int cmpVert(const void *a, const void *b){
    return *(int*)a-*(int*)b;
}

/**
 * @brief
 * @details
 * @param [in] mesh MultiReg2d mesh object
 * @param [in] Nsurf number of boundary surface in SFToV
 * @param [in] SFToV surface to vertex
 *
 */
MultiRegBC2d* MultiRegBC2d_create(MultiReg2d *mesh, int Nsurf, int **SFToV){
    StdRegions2d *shape = mesh->stdcell;

    MultiRegBC2d *bc2d = (MultiRegBC2d *) malloc(sizeof(MultiRegBC2d));

    int k,f1,f2;

    /* count obc number */
    int surfList[Nsurf];
    for(f1=0;f1<Nsurf;f1++){
        surfList[f1] = SFToV[f1][2];
        /* check SFToV indicator */
        if(SFToV[f1][2] == I_INNERLOC | SFToV[f1][2] == I_INNERBS ){
            printf("MultiRegBC2d_create: Error boundary type in SFToV[%d][3] = %d\n", f1,SFToV[f1][2]);
            printf("The boundary type indicator cannot be 0 or 1:\n");
            printf("   %d - local boundary surface[default]\n", I_INNERLOC);
            printf("   %d - parallel boundary surface\n", I_INNERBS);
        }
    }
    qsort(surfList, Nsurf, sizeof(int), cmpVert); // sort the boundary type
    bc2d->Nobc = 1;
    for(f1=0;f1<(Nsurf-1);f1++){
        if(surfList[f1+1]-surfList[f1] != 0)
            bc2d->Nobc++;
    }

    /* store the boundary type indicator (from smallest to largest) */
    bc2d->bcTypeList = IntVector_create(bc2d->Nobc);
    bc2d->bcTypeList[0] = surfList[0];
    int sk = 0;
    for(f1=0;f1<(Nsurf-1);f1++){
        if(surfList[f1+1]-surfList[f1] != 0)
            bc2d->bcTypeList[sk++] = surfList[f1+1];
    }

    /* allocate and initialize obc2d */
    bc2d->obc2d = (OBC2d**) malloc(bc2d->Nobc*sizeof(OBC2d*));
    for(k=0;k<bc2d->Nobc;k++){
        bc2d->obc2d[k] = OBC2d_create(Nsurf, SFToV, bc2d->bcTypeList[k]);
    }

    /* allocate and initialize EToBS */
    bc2d->EToBS = IntMatrix_create(mesh->K, shape->Nfaces);

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

                if(t_temp == v_temp){ // compare the face
                    bc2d->EToBS[k][f2] = SFToV[f1][2];
                    continue; // finish this loop
                }
            }

            /* inner surface */
            bc2d->EToBS[k][f2] = INNERLOC;
        }
    }

    return bc2d;
}


OBC2d* OBC2d_create(int Nsurf, int **SFToV, int typeid){

    OBC2d *obc2d = (OBC2d*) malloc(sizeof(OBC2d));
    /* count vertex number */
    int f1;
    int vertlist[Nsurf*2];
    for(f1=0; f1<Nsurf; f1++){
        vertlist[f1*2]   = 0;
        vertlist[f1*2+1] = 0;
        if( SFToV[f1][2] == typeid){
            vertlist[f1*2]   = SFToV[f1][0];
            vertlist[f1*2+1] = SFToV[f1][1];
        }
    }
    qsort(vertlist, Nsurf*2, sizeof(int), cmpVert); // sort the boundary type

    obc2d->Nv = 0;
    for(f1=0;f1<(Nsurf-1);f1++){
        if(vertlist[f1+1]-vertlist[f1] != 0)
            obc2d->Nv++;
    }

    obc2d->BVToV = IntVector_create(obc2d->Nv);
    int sk = 0;
    for(f1=0;f1<(Nsurf-1);f1++){
        if(vertlist[f1+1]-vertlist[f1] != 0)
            obc2d->BVToV[sk++] = vertlist[f1+1];
    }
    return obc2d;
}

void OBC2d_free(OBC2d *obc){
    IntVector_free(obc->BVToV);
}

void MultiRegBC2d_free(MultiRegBC2d* bc2d){
    IntMatrix_free(bc2d->EToBS);
    int i;
    for(i=0;i<bc2d->Nobc;i++){
        OBC2d_free(bc2d->obc2d[i]);
    }
    free(bc2d->obc2d);
    free(bc2d);
}