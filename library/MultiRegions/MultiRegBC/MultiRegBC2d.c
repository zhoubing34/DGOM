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
#define I_INNERBS 1

vertlist* OBC2d_create(int Nsurf, int **SFToV, int typeid);
// connect adjacent face nodes
void SetNodePair2d(stdCell *shape, int K, double **GX, double **GY,
                   int **EToE, int **EToF, int **EToP, double **x, double **y,
                   int *Npar, int *Ntotalout, int **mapOUT,
                   int *vmapM, int *vmapP);
void SetElementPair(stdCell *shape, MultiReg2d *mesh, int *cellIndOut);


// sort numbers from small to large
int cmpVert(const void *a, const void *b){
    return *(int*)a-*(int*)b;
}

/**
 * @brief count the number of uique elements in an array
 * @param [in]     len length of the array
 * @param [in,out] list array of integer
 * @param [out]    n number of unique elements
 */
int count_unique_integer(int len, int *list){
    qsort(list, len, sizeof(int), cmpVert); // sort the list
    int n = 1, i;
    for(i=0;i<(len-1);i++){
        if(list[i+1]-list[i] != 0)
            n++;
    }
    return n;
}

/**
 * @brief
 * @details
 * @param [in] mesh MultiReg2d mesh object
 * @param [in] Nsurf number of boundary surface in SFToV
 * @param [in] SFToV surface to vertex
 * @note
 * precondition:
 * 1. The index of vertex in SFToV is start from 0;
 * 2. The contents of SFToV is
 * postcondition:
 * 1. call MultiRegBC2d_free manually to free the MultiRegBC2d structure
 * in case of memory leak
 */
MultiRegBC2d* MultiRegBC2d_create(MultiReg2d *mesh, int Nsurf, int **SFToV){

#if DEBUG
    int procid, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    if(procid == IPROC)
        printf("Step into MultiRegBC2d_create\n");
#endif

    stdCell *shape = mesh->stdcell;
    MultiRegBC2d *surf2d = (MultiRegBC2d *) malloc(sizeof(MultiRegBC2d));
    surf2d->mesh = mesh;

    int k,f1,f2;

    /* count obc number */
    int surfList[Nsurf];
    for(f1=0;f1<Nsurf;f1++){
        surfList[f1] = SFToV[f1][2];
        /* check SFToV indicator */
        if(SFToV[f1][2] == I_INNERLOC | SFToV[f1][2] == I_INNERBS ){
            printf("MultiRegBC2d_create: Error boundary type in SFToV[%d][3] = %d\n", f1,SFToV[f1][2]);
            printf("The boundary type indicator cannot be %d or %d:\n", I_INNERLOC, I_INNERBS);
            printf("   %d - local boundary surface[default]\n", I_INNERLOC);
            printf("   %d - parallel boundary surface\n", I_INNERBS);
        }
    }

    surf2d->Nobc = count_unique_integer(Nsurf, surfList);

    /* store the boundary type indicator (from smallest to largest) */
    surf2d->bcTypeList = IntVector_create(surf2d->Nobc);
    surf2d->bcTypeList[0] = surfList[0];
    int sk = 1;
    for(f1=0;f1<(Nsurf-1);f1++){
        if(surfList[f1+1]-surfList[f1] != 0)
            surf2d->bcTypeList[sk++] = surfList[f1+1];
    }

    /* allocate and initialize oblist */
    surf2d->oblist = (vertlist**) malloc(surf2d->Nobc*sizeof(vertlist*));
    for(k=0;k<surf2d->Nobc;k++){
        surf2d->oblist[k] = OBC2d_create(Nsurf, SFToV, surf2d->bcTypeList[k]);
    }

    /* allocate and initialize EToE and EToBS */
    surf2d->EToE = mesh->EToE;
    surf2d->EToBS = IntMatrix_create(mesh->K, shape->Nfaces);

    int t[2], v[2];
    for(k=0;k<mesh->K;k++){
        for(f2=0; f2<shape->Nfaces; f2++){ // loop through all element surface
            /* inner surface (default) */
            surf2d->EToBS[k][f2] = INNERLOC;

            /* the inner boundary surface */
            if(mesh->EToP[k][f2] != mesh->procid){
                surf2d->EToBS[k][f2] = INNERBS;
                continue; // finish this loop
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
#if DEBUG
                if(procid == IPROC)
                    printf("surfid=%d, v1=%d, v2=%d, v_temp = %d, t_temp?=v_temp%d\n", f1, v[0], v[1], v_temp, t_temp==v_temp);
#endif
                if(t_temp == v_temp){ // compare the face
                    surf2d->EToBS[k][f2] = SFToV[f1][2];
                    break; // jump out loop of SFToV
                }
            }
        }
    }

    /* adjacent face node id */
    surf2d->vmapM = IntVector_create(shape->Nfp*shape->Nfaces*mesh->K);
    surf2d->vmapP = IntVector_create(shape->Nfp * shape->Nfaces * mesh->K);
    SetNodePair2d(shape, mesh->K, mesh->GX, mesh->GY,
                  mesh->EToE, mesh->EToF, mesh->EToP,
                  mesh->x, mesh->y,mesh->Npar,
                  &(surf2d->parNodeTotalOut), &(surf2d->nodeIndexOut),surf2d->vmapM, surf2d->vmapP);

    /* adjacent cell id */
    surf2d->parCellTotalOut = surf2d->parNodeTotalOut/shape->Nfp;
    surf2d->cellIndexOut = IntVector_create(surf2d->parCellTotalOut);
    SetElementPair(shape, mesh, surf2d->cellIndexOut);

#if DEBUG
    if(procid == IPROC)
        printf("Step out MultiRegBC2d_create\n");
#endif
    return surf2d;
}

/**
 * @brief create open boundary object
 *
 * @param [in]  Nsurf
 * @param [in]  SFToV
 * @param [in]  typeid
 *
 */
vertlist* OBC2d_create(int Nsurf, int **SFToV, int typeid){

#if DEBUG
    int procid, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    if(procid == IPROC)
        printf("Step into OBC2d_create for typeId = %d\n", typeid);
#endif

    vertlist *obc2d = (vertlist*) malloc(sizeof(vertlist));
    /* count vertex number */
    int f1, Nsurfv=0;
    int vertlist[Nsurf*2];
    for(f1=0; f1<Nsurf; f1++){
        if( SFToV[f1][2] == typeid){
            vertlist[Nsurfv++] = SFToV[f1][0];
            vertlist[Nsurfv++] = SFToV[f1][1];
        }
    }

    obc2d->Nv = count_unique_integer(Nsurfv, vertlist);

#if DEBUG
    if(procid == IPROC){
        printf("different vertex number = %d\n", oblist->Nv);
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

/**
 * @brief
 * @param [in] obc open boundary object
 *
 */
void OBC2d_free(vertlist *obc){
    IntVector_free(obc->BVToV);
    free(obc);
}

/**
 * @brief
 * @param [in] bc2d
 */
void MultiRegBC2d_free(MultiRegBC2d* bc2d){
    IntMatrix_free(bc2d->EToBS);
    int i;
    for(i=0;i<bc2d->Nobc;i++){
        OBC2d_free(bc2d->oblist[i]);
    }
    free(bc2d->oblist);

    IntVector_free(bc2d->nodeIndexOut);
    IntVector_free(bc2d->vmapM);
    IntVector_free(bc2d->vmapP);
    IntVector_free(bc2d->cellIndexOut);
    IntVector_free(bc2d->bcTypeList);
//    free(bc2d->vert_ext);
    free(bc2d);
}