//
// Created by li12242 on 17/3/10.
//

#include "dg_edge_face_map.h"

static int dg_mesh_face_cmp(const void *obj1, const void *obj2);

typedef struct face2d{
    int k1, k2; ///< cell id
    int f1, f2; ///< face id
    dg_grid_face_type ftype; ///< face type
    int tmp; ///< indicator for cell paris
}face2d;

void dg_edge_face_map3d(dg_edge *edge){

}

void dg_edge_face_map2d(dg_edge *edge){
    dg_grid *grid = edge->grid;
    const int K = dg_grid_K(grid);
    const int Nvert = dg_grid_Nv(grid);
    const int Nfaces = dg_cell_Nfaces(edge->cell);

    int **FToV = dg_cell_FToV(edge->cell);
    int **EToV = grid->EToV;
    int **EToE = grid->EToE;
    int **EToF = grid->EToF;
    int **EToBS = grid->EToBS;

    const int totalFace = K*Nfaces;
    face2d *my_face = (face2d *)calloc(totalFace, sizeof(face2d));
    int k,f,sk=0;
    for(k=0;k<K;k++){
        for(f=0;f<Nfaces;f++){
            int n1 = EToV[k][ FToV[f][0] ];
            int n2 = EToV[k][ FToV[f][1] ];
            my_face[sk].k1 = k;
            my_face[sk].k2 = EToE[k][f];
            my_face[sk].f1 = f;
            my_face[sk].f2 = EToF[k][f];
            my_face[sk].ftype = (dg_grid_face_type)EToBS[k][f];
            my_face[sk].tmp = min(n1, n2)*Nvert + max(n1, n2);
//            if(!edge->procid) printf("face[%d]: k1=%d, k2=%d, f1=%d, f2=%d, tmp=%d\n",
//                                     sk, my_face[sk].k1, my_face[sk].k2,
//                                     my_face[sk].f1, my_face[sk].f2, my_face[sk].tmp);
            sk++;
        }
    }

    qsort(my_face, (size_t)totalFace, sizeof(face2d), dg_mesh_face_cmp);
    /* count the number of edges */
    int ftmp[totalFace];
    for(f=0;f<totalFace;f++){ ftmp[f] = my_face[f].tmp; }
//    if(!edge->procid){
//        printf("ftmp=");
//        for(f=0;f<totalFace;f++){ printf("%d ", ftmp[f]); }
//        printf("\n");
//    }
    int Nedge = unique_int(totalFace, ftmp);

    int *varkM = vector_int_create(Nedge);
    int *varkP = vector_int_create(Nedge);
    int *varfM = vector_int_create(Nedge);
    int *varfP = vector_int_create(Nedge);
    int *ftype = vector_int_create(Nedge);

    varkM[0] = my_face[0].k1; varkP[0] = my_face[0].k2;
    varfM[0] = my_face[0].f1; varfP[0] = my_face[0].f2;
    ftype[0] = my_face[0].ftype;
    sk = 1;
    for(k=1;k<totalFace;k++){
        if(my_face[k].tmp != my_face[k-1].tmp ){
            varkM[sk] = my_face[k].k1;
            varkP[sk] = my_face[k].k2;
            varfM[sk] = my_face[k].f1;
            varfP[sk] = my_face[k].f2;
            ftype[sk] = my_face[k].ftype;
//            if(!edge->procid) printf("edge[%d]: k1=%d, k2=%d, f1=%d, f2=%d, ftype=%d\n",
//                                     sk, varkM[sk], varkP[sk], varfM[sk], varfP[sk], ftype[sk]);
            sk++;
        }
    }
    /* assignment */
    edge->Nedge = Nedge;
    edge->varkM = varkM;
    edge->varkP = varkP;
    edge->varfM = varfM;
    edge->varfP = varfP;
    edge->ftype = ftype;

    free(my_face);
}

/* sort numbers from small to large */
static int dg_mesh_face_cmp(const void *obj1, const void *obj2){
    face2d *e1 = (face2d*) obj1;
    face2d *e2 = (face2d*) obj2;
#if DEBUG
    int procid;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    printf("procid=%d, %s (%d), e1: k1=%ld, k2=%ld, f1=%d, f2=%d, tmp=%ld\n",
           procid, __FUNCTION__, __LINE__, e1->k1, e1->k2, e1->f1, e1->f2, e1->tmp);
    printf("procid=%d, %s (%d), e2: k1=%ld, k2=%ld, f1=%d, f2=%d, tmp=%ld\n",
           procid, __FUNCTION__, __LINE__, e2->k1, e2->k2, e2->f1, e2->f2, e2->tmp);
#endif
    if ( (e1->tmp) > (e2->tmp) ) {return 1;}
    else if( (e1->tmp) < (e2->tmp) ) {return -1;}
    return 0;
}