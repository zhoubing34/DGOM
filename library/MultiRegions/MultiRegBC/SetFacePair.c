#include "MultiRegions/MultiRegions.h"

/* ======================================================================== */
int compare_pairs(const void *obj1, const void *obj2);
int pairprocget(const void *obj1);
int pairnumget(const void *obj1);
void pairnumset(const void *obj1, int g);
void pairmarry(const void *obj1, const void *obj2);

/* private functions and structures for FacePair */
typedef struct face {
    int p1, k1, f1;
    int p2, k2, f2;
    /** max index of vertex on face */
    int va;
    /** min index of vertex on face */
    int vb;
    /** max index of vertex */
    int g;
}face;

int compare_pairs(const void *obj1, const void *obj2){

    face *e1 = (face*) obj1;
    face *e2 = (face*) obj2;

    int a1 = e1->va, b1 = e1->vb;
    int a2 = e2->va, b2 = e2->vb;

    int va1, vb1, va2, vb2;

    va1 = min(a1, b1);
    vb1 = max(a1, b1);

    va2 = min(a2, b2);
    vb2 = max(a2, b2);

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


int pairprocget(const void *obj1){
    face *e1 = (face*) obj1;
    return (e1->p1);
}


int pairnumget(const void *obj1){
    face *e1 = (face*) obj1;
    return (e1->g);
}

void pairnumset(const void *obj1, int g){
    face *e1 = (face*) obj1;
    e1->g = g;
}

void pairmarry(const void *obj1, const void *obj2){
    face *e1 = (face*) obj1;
    face *e2 = (face*) obj2;
    e1->p2 = e2->p1;  e1->k2 = e2->k1;  e1->f2 = e2->f1;
    e2->p2 = e1->p1;  e2->k2 = e1->k1;  e2->f2 = e1->f1;
}
/* ======================================================================== */

/**
 * @brief
 * Setup element connection
 * @details
 * Set EToE,EToP,EToF & Npar,parK,parF
 *
 * @param [stdCell*] shape standard element
 * @param [int]     Klocal  number of elements
 * @param [int**]   EToV    element to vertex list
 *
 * @return
 * name     | type     | description of value
 * -------- |----------|----------------------
 * EToE | int[Klocal][shape->Nv]   | element to element list
 * EToF | int[Klocal][shape->Nv]   | element to face index list
 * EToP | int[Klocal][shape->Nv]   | element to process list
 * Npar | int[nprocs] | number of faces adjacent to each process
 *
 * @note
 * EToE, EToF, EToP and Npar should be allocated before calling this function.
 */
void SetFacePair2d(stdCell *shape, int Klocal,
                 int **EToV, int **EToE, int **EToF, int **EToP,
                 int *Npar){

    int procid, nprocs;
    int Nfaces = shape->Nfaces;

    int n, k, e, sk;

    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

//    triangle      vnum[3][2] = { {0,1}, {1,2}, {2,0} };
//    quadrilateral vnum[4][2] = { {0,1}, {1,2}, {2,3}, {3,0}};
    int **vnum = IntMatrix_create(shape->Nfaces, 2);
    for(n=0;n<shape->Nfaces;n++){
        vnum[n][0] = n%(shape->Nv);
        vnum[n][1] = (n+1)%(shape->Nv);
    }

    face *myfaces = (face*) calloc(Klocal*Nfaces, sizeof(face));

    sk = 0;
    for(k=0;k<Klocal;++k){
        for(e=0;e<Nfaces;++e){
            int n1 = EToV[k][vnum[e][0]];
            int n2 = EToV[k][vnum[e][1]];

            myfaces[sk].p1 = procid; myfaces[sk].k1 = k; myfaces[sk].f1 = e;
            myfaces[sk].p2 = procid; myfaces[sk].k2 = k; myfaces[sk].f2 = e;

            myfaces[sk].va = max(n1,n2);
            myfaces[sk].vb = min(n1,n2);
            myfaces[sk].g  = max(n1,n2); /* marker for sorting into bins */
            ++sk;
        }
    }

    ParallelPairs(myfaces, Klocal*Nfaces, sizeof(face),
                  pairnumget, pairnumset, pairprocget,
                  pairmarry, compare_pairs);

    int k1, k2, f1, f2, p1, p2;

    for(n=0;n<Klocal*Nfaces;++n){

        k1 = myfaces[n].k1; f1 = myfaces[n].f1; p1 = myfaces[n].p1;
        k2 = myfaces[n].k2; f2 = myfaces[n].f2; p2 = myfaces[n].p2;

        if(p1!=procid)
            fprintf(stderr, "WARNING WRONG proc\n");

        EToE[k1][f1] = k2; /* adjacent element index */
        EToF[k1][f1] = f2; /* local face index of adjacent element */
        EToP[k1][f1] = p2; /* process id of adjacent element */

        if(p1!=p2){
            /* increment number of links */
            ++Npar[p2];
        }
    }

    free(myfaces);
    IntMatrix_free(vnum);

}
