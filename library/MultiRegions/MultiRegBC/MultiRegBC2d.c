//
// Created by li12242 on 16/12/14.
//

#include <MultiRegions/MultiRegions.h>
#include "MultiRegBC2d.h"

/**
 * @brief create boundary type object
 * @details
 * @param [in] mesh
 * @param [in] SFoV surface to vertex
 *
 */
#define DEBUG 0
#define IPROC 0

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

#if DEBUG
    int procid, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    if(procid == IPROC)
        printf("Step into MultiRegBC2d_create\n");
#endif

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
    int sk = 1;
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
#if DEBUG
    if(procid == IPROC)
        printf("Ne=%d, Nfaces=%d\n", mesh->K, shape->Nfaces);
#endif
    for(k=0;k<mesh->K;k++){
        for(f2=0; f2<shape->Nfaces; f2++){ // loop through all element surface
            /* inner surface (default) */
            bc2d->EToBS[k][f2] = INNERLOC;

            /* the inner boundary surface */
            if(mesh->EToP[k][f2] != mesh->procid){
                bc2d->EToBS[k][f2] = INNERBS;
                continue; // finish this loop
            }

            /* get the usr-defined boundary */
            int n1 = f2;
            int n2 = (f2+1)%shape->Nv;
            t[0] = mesh->EToV[k][n1];
            t[1] = mesh->EToV[k][n2];
            qsort(t, 2, sizeof(int), cmpVert);

            int t_temp = t[0]*mesh->Nv + t[1];
#if DEBUG
            if(procid == IPROC)
                printf("k=%d, f=%d, n=[%d,%d], v1=%d, v2=%d, t_temp = %d\n", k, f2, n1, n2, t[0], t[1], t_temp);
#endif

            for(f1=0; f1<Nsurf; f1++){
                v[0] = SFToV[f1][0];
                v[1] = SFToV[f1][1];
                qsort(v, 2, sizeof(int), cmpVert);

                int v_temp = v[0]*mesh->Nv + v[1];
#if DEBUG
                if(procid == IPROC)
                    printf("surfid=%d, v1=%d, v2=%d, v_temp = %d, t_temp?=v_temp%d\n", f1, v[0], v[1], v_temp, t_temp==v_temp);
#endif
                if(t_temp == v_temp){ // compare the face
                    bc2d->EToBS[k][f2] = SFToV[f1][2];
                    break; // jump out loop
                }
            }
        }
    }

    /* allocate the vertex external data */
    bc2d->vert_ext = (real*) calloc(mesh->Nv, sizeof(real));

#if DEBUG
    if(procid == IPROC)
        for(k=0;k<mesh->K;k++){
            printf("EToBS[%d][:] = ", k);
            for(f1=0;f1<shape->Nv;f1++){
                printf("%d, ", bc2d->EToBS[k][f1]);
            }
            printf("\n");
        }
        printf("Step out MultiRegBC2d_create\n");
#endif

    return bc2d;
}


OBC2d* OBC2d_create(int Nsurf, int **SFToV, int typeid){

#if DEBUG
    int procid, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    if(procid == IPROC)
        printf("Step into OBC2d_create for typeId = %d\n", typeid);
#endif

    OBC2d *obc2d = (OBC2d*) malloc(sizeof(OBC2d));
    /* count vertex number */
    int f1, Nsurfv=0;
    int vertlist[Nsurf*2];
    for(f1=0; f1<Nsurf; f1++){
        if( SFToV[f1][2] == typeid){
            vertlist[Nsurfv++] = SFToV[f1][0];
            vertlist[Nsurfv++] = SFToV[f1][1];
        }
    }
    qsort(vertlist, Nsurfv, sizeof(int), cmpVert); // sort the vertex

#if DEBUG
    if(procid == IPROC){
        printf("vertex number = %d\n", Nsurfv);
        for(f1=0; f1<Nsurfv; f1++){
            printf("vertlist[%d] = %d, ", f1, vertlist[f1]);
        }
        printf("\n");
    }
#endif

    obc2d->Nv = 1;
    for(f1=0;f1<(Nsurfv-1);f1++){
        if(vertlist[f1+1]-vertlist[f1] != 0)
            obc2d->Nv++;
    }

#if DEBUG
    if(procid == IPROC){
        printf("different vertex number = %d\n", obc2d->Nv);
    }
#endif

    obc2d->BVToV = IntVector_create(obc2d->Nv);
    int sk = 0;
    obc2d->BVToV[sk++] = vertlist[0];
    for(f1=0;f1<(Nsurfv-1);f1++){
        if(vertlist[f1+1]-vertlist[f1] != 0)
            obc2d->BVToV[sk++] = vertlist[f1+1];
    }
#if DEBUG
    if(procid == IPROC)
        printf("Step out OBC2d_create for typeId = %d\n", typeid);
#endif

    return obc2d;
}

void OBC2d_free(OBC2d *obc){
    IntVector_free(obc->BVToV);
    free(obc);
}

void MultiRegBC2d_free(MultiRegBC2d* bc2d){
    IntMatrix_free(bc2d->EToBS);
    int i;
    for(i=0;i<bc2d->Nobc;i++){
        OBC2d_free(bc2d->obc2d[i]);
    }
    free(bc2d->obc2d);

    IntVector_free(bc2d->bcTypeList);
    free(bc2d->vert_ext);
    free(bc2d);
}