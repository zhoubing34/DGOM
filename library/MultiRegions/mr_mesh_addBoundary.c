//
// Created by li12242 on 12/21/16.
//

#include "mr_mesh_addBoundary.h"

/* count the number of uique elements in an array */
static int count_unique_integer(int len, int *list);

/* sort numbers from small to large */
static int mr_mesh_cmp(const void *a, const void *b);

/* create open boundary vertex list */
static vertlist* mr_vertList2d_create(int Nsurf, int **SFToV, int typeid);

/* find the boundary type id, in bcTypeList */
#define _surfTypeId(mesh, ind, typeid) do{\
typeid = ind;\
}while(0)

/**
 * @brief add boundary conditions into the 2d mesh object
 *
 * @details
 * The `SFToV` matrix contains the vertex list and surface type indicator (integer) of each surface,
 * .e.g, [v1, v2, typeid].
 *
 * @param[in,out] mesh mesh object
 * @param[in] Nsurf number of surface
 * @param[in] SFToV surface to vertex list
 */
void mr_mesh_addBoundary2d(parallMesh *mesh, int Nsurf, int **SFToV){
    int k,f1,f2;

    const int K = mesh->grid->K;
    const int Nfaces = mesh->cell->Nfaces;
    const int Nv = mesh->cell->Nv; ///< number of vertex in each element
    const int Nvert = mesh->grid->Nv; ///< number of all the vertex

    int **EToV = mesh->grid->EToV;

    /* count obc number */
    int surfList[Nsurf];
    for(f1=0;f1<Nsurf;f1++){
        surfList[f1] = SFToV[f1][2];
        if(SFToV[f1][2] == INNERLOC | SFToV[f1][2] == INNERBS ){
            printf("MultiRegions (mr_mesh_addBoundary): Error boundary type in SFToV[%d][3] = %d\n", f1,SFToV[f1][2]);
            printf("The boundary type indicator cannot be %d or %d:\n", INNERLOC, INNERBS);
            printf("   %d - local boundary surface[default]\n", INNERLOC);
            printf("   %d - parallel boundary surface\n", INNERBS);
            printf("please use other integers for boundary ID:\n");
            printf("   %d - slip wall\n", SLIPWALL);
            printf("   %d - non-slip wall\n", NSLIPWALL);
            printf("   other ids - different open boundaries\n");
        }
    }
    int Nobc = count_unique_integer(Nsurf, surfList);
    mesh->Nbc = Nobc;

    /* store the boundary type indicator (from smallest to largest) */
    mesh->bcIndList = IntVector_create(Nobc);
    mesh->bcIndList[0] = surfList[0];
    int sk = 1;
    for(f1=0;f1<(Nsurf-1);f1++){
        if(surfList[f1+1]-surfList[f1] != 0)
            mesh->bcIndList[sk++] = surfList[f1+1];
    }

    /* allocate and initialize oblist */
    mesh->obvertlist = (vertlist**) malloc(Nobc*sizeof(vertlist*));
    for(k=0;k<Nobc;k++){
        mesh->obvertlist[k] = mr_vertList2d_create(Nsurf, SFToV, mesh->bcIndList[k]);
    }

    /* now allocate the element to surface type matrix */
    mesh->EToBS = IntMatrix_create(K, Nfaces);

    int t[2], v[2];
    for(k=0;k<K;k++){
        for(f2=0; f2<Nfaces; f2++){ // loop through all element surface
            /* inner surface (default) */
            mesh->EToBS[k][f2] = INNERLOC;

            /* the inner boundary surface */
            if(mesh->EToP[k][f2] != mesh->procid){
                mesh->EToBS[k][f2] = INNERBS;
                continue; // finish this loop
            }

            /* get the usr-defined boundary */
            int n1 = f2;
            int n2 = (f2+1)%Nv;
            t[0] = EToV[k][n1];
            t[1] = EToV[k][n2];
            qsort(t, 2, sizeof(int), mr_mesh_cmp);

            int t_temp = t[0]*Nvert + t[1];
#if 0
            if(mesh->procid == 0)
                printf("k=%d, f=%d, t1=%d, t2=%d, t_temp = %d, \n", k, f2, t[0], t[1], t_temp);
#endif
            for(f1=0; f1<Nsurf; f1++){
                v[0] = SFToV[f1][0];
                v[1] = SFToV[f1][1];
                qsort(v, 2, sizeof(int), mr_mesh_cmp);

                int v_temp = v[0]*Nvert + v[1];
#if 0
                if(mesh->procid == 0)
                    printf("surfid=%d, v1=%d, v2=%d, v_temp = %d\n", f1, v[0], v[1], v_temp);
#endif
                if(t_temp == v_temp){ // compare the face
                    _surfTypeId(mesh, SFToV[f1][2], mesh->EToBS[k][f2]);
                    break; // jump out loop of SFToV
                }
            }
        }
    }
}

/* sort numbers from small to large */
static int mr_mesh_cmp(const void *a, const void *b){
    return (* (int *)a) - (* (int *)b);
}

/**
 * @brief count the number of uique elements in an array
 * @param [in]     len length of the array
 * @param [in,out] list array of integer
 * @param [out]    n number of unique elements
 */
static int count_unique_integer(int len, int *list){

    if(len ==0 )
        return 0;

    qsort(list, len, sizeof(int), mr_mesh_cmp); // sort the list
    int n = 1, i;
    for(i=0;i<(len-1);i++){
        if(list[i+1]-list[i] != 0)
            n++;
    }
    return n;
}

/**
 * @brief create open boundary vertex list
 *
 * @param [in]  Nsurf number of surface in SFToV
 * @param [in]  SFToV surface to vertex list
 * @param [in]  typeid the boundary type
 */
static vertlist* mr_vertList2d_create(int Nsurf, int **SFToV, int typeid){

    vertlist *vertlist2d = (vertlist*) malloc(sizeof(vertlist));
    /* count vertex number */
    int f1, Nvert=0;
    int vertlist[Nsurf*2];
    for(f1=0; f1<Nsurf; f1++){
        if( SFToV[f1][2] == typeid){
            vertlist[Nvert++] = SFToV[f1][0];
            vertlist[Nvert++] = SFToV[f1][1];
        }
    }

    vertlist2d->Nv = count_unique_integer(Nvert, vertlist);
    vertlist2d->list = IntVector_create(vertlist2d->Nv);

    int sk = 0;
    vertlist2d->list[sk++] = vertlist[0];
    for(f1=0;f1<(Nvert-1);f1++){
        if(vertlist[f1+1]-vertlist[f1] != 0)
            vertlist2d->list[sk++] = vertlist[f1+1];
    }

    return vertlist2d;
}