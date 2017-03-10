//
// Created by li12242 on 17/3/8.
//

#include "dg_grid.h"
#include "dg_grid_parallelPairs.h"

static int compare_pairs2d(const void *obj1, const void *obj2);
static int pairprocget2d(const void *obj1);
static int pairnumget2d(const void *obj1);
static void pairnumset2d(const void *obj1, int g);
static void pairmarry2d(const void *obj1, const void *obj2);

/* private functions and structures for FacePair */
typedef struct face {
    int p1, k1, f1;
    int p2, k2, f2;

    int va, vb; ///< index of vertex on face
    int g; ///< max index of vertex
}face2d;

void dg_grid_connect2d(dg_grid *grid){
    const int procid = dg_grid_procid(grid);
    const int Klocal = dg_grid_K(grid);
    const int Nfaces = dg_cell_Nfaces(grid->cell);
    const int Nv = dg_cell_Nv(grid->cell);

    int **EToV = grid->EToV;

    grid->EToE = matrix_int_create(Klocal, Nv);
    grid->EToF = matrix_int_create(Klocal, Nv);
    grid->EToP = matrix_int_create(Klocal, Nv);

    int n,k,e,sk=0;
    face2d *myfaces = (face2d*) calloc(Klocal*Nfaces, sizeof(face2d));

    int vind[Nfaces][2];
    // triangle      vind[3][2] = { {0,1}, {1,2}, {2,0} };
    // quadrilateral vind[4][2] = { {0,1}, {1,2}, {2,3}, {3,0}};
    for(n=0;n<Nfaces;n++){
        vind[n][0] = n%Nv;
        vind[n][1] = (n+1)%(Nv);
    }
    for(k=0;k<Klocal;++k){
        for(e=0;e<Nfaces;++e){
            int n1 = EToV[k][vind[e][0]];
            int n2 = EToV[k][vind[e][1]];

            myfaces[sk].p1 = procid; myfaces[sk].k1 = k; myfaces[sk].f1 = e;
            myfaces[sk].p2 = procid; myfaces[sk].k2 = k; myfaces[sk].f2 = e;

            myfaces[sk].va = max(n1,n2);
            myfaces[sk].vb = min(n1,n2);
            myfaces[sk].g  = max(n1,n2); /* marker for sorting into bins */
            ++sk;
        }
    }

    dg_grid_parallelPairs(myfaces, Klocal*Nfaces, sizeof(face2d),
                          pairnumget2d, pairnumset2d, pairprocget2d,
                          pairmarry2d, compare_pairs2d);

    for(n=0;n<Klocal*Nfaces;++n){
        int k1,k2,f1,f2,p1,p2;
        k1 = myfaces[n].k1; f1 = myfaces[n].f1; p1 = myfaces[n].p1;
        k2 = myfaces[n].k2; f2 = myfaces[n].f2; p2 = myfaces[n].p2;

        if(p1!=procid) {fprintf(stderr, "%s: %d\nWARNING WRONG proc\n", __FUNCTION__, __LINE__);}

        grid->EToE[k1][f1] = k2; /* adjacent element index */
        grid->EToF[k1][f1] = f2; /* local face2d index of adjacent element */
        grid->EToP[k1][f1] = p2; /* process id of adjacent element */
    }

    free(myfaces);
    return;
}

void dg_grid_connect3d(dg_grid *grid){
    return;
}

static int compare_pairs2d(const void *obj1, const void *obj2){
    face2d *e1 = (face2d*) obj1;
    face2d *e2 = (face2d*) obj2;
    int a1 = e1->va, b1 = e1->vb;
    int a2 = e2->va, b2 = e2->vb;
    int va1, vb1, va2, vb2;
    va1 = min(a1, b1); vb1 = max(a1, b1);
    va2 = min(a2, b2); vb2 = max(a2, b2);
    if(vb1<vb2)
        return -1;
    else if(vb1>vb2)
        return 1;
    else if(va1<va2)
        return -1;
    else if(va1>va2)
        return 1;

    return 0;
}


static int pairprocget2d(const void *obj1){
    face2d *e1 = (face2d*) obj1;
    return (e1->p1);
}

static int pairnumget2d(const void *obj1){
    face2d *e1 = (face2d*) obj1;
    return (e1->g);
}

static void pairnumset2d(const void *obj1, int g){
    face2d *e1 = (face2d*) obj1;
    e1->g = g;
}

static void pairmarry2d(const void *obj1, const void *obj2){
    face2d *e1 = (face2d*) obj1;
    face2d *e2 = (face2d*) obj2;
    e1->p2 = e2->p1;  e1->k2 = e2->k1;  e1->f2 = e2->f1;
    e2->p2 = e1->p1;  e2->k2 = e1->k1;  e2->f2 = e1->f1;
}