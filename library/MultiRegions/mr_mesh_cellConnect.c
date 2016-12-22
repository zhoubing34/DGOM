#include "LibUtilities/LibUtilities.h"
#include "mr_mesh_parallelPairs.h"
#include "mr_mesh_cellConnect.h"

static int compare_pairs(const void *obj1, const void *obj2);
static int pairprocget(const void *obj1);
static int pairnumget(const void *obj1);
static void pairnumset(const void *obj1, int g);
static void pairmarry(const void *obj1, const void *obj2);

/* private functions and structures for FacePair */
typedef struct face {
    int p1, k1, f1;
    int p2, k2, f2;

    int va; ///< max index of vertex on face
    int vb; ///< min index of vertex on face
    int g; ///< max index of vertex
}face;

static int compare_pairs(const void *obj1, const void *obj2){

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


static int pairprocget(const void *obj1){
    face *e1 = (face*) obj1;
    return (e1->p1);
}


static int pairnumget(const void *obj1){
    face *e1 = (face*) obj1;
    return (e1->g);
}

static void pairnumset(const void *obj1, int g){
    face *e1 = (face*) obj1;
    e1->g = g;
}

static void pairmarry(const void *obj1, const void *obj2){
    face *e1 = (face*) obj1;
    face *e2 = (face*) obj2;
    e1->p2 = e2->p1;  e1->k2 = e2->k1;  e1->f2 = e2->f1;
    e2->p2 = e1->p1;  e2->k2 = e1->k1;  e2->f2 = e1->f1;
}

/**
 * @brief Setup element connection
 * @details
 * set properties of EToE,EToP,EToF, Npar, parCellTotalOut.
 * * EToE -  element to element index, includes the adjacent element in another process.
 * * EToF - element to face index.
 * * EToP - element to process id.
 * * Npar - number of faces adjacent to each process.
 * * parCellTotalOut - number of cell send buffer
 * * nodeIndexOut - index of element in the send buffer
 *
 * @param [in,out] shape standard element
 *
 */
void mr_mesh_cellConnect2d(parallMesh *mesh){

    const int procid = mesh->procid;
    const int nprocs = mesh->nprocs;

    const int Nfaces = mesh->cell->Nfaces;
    const int Nv = mesh->cell->Nv;
    const int Klocal = mesh->grid->K;

    int **EToV = mesh->grid->EToV;

    mesh->EToE = IntMatrix_create(Klocal, Nv);
    mesh->EToF = IntMatrix_create(Klocal, Nv);
    mesh->EToP = IntMatrix_create(Klocal, Nv);

    mesh->Npar = IntVector_create(nprocs);
    mesh->parallCellNum = 0;

    int *Npar = mesh->Npar;

    int n, k, e, sk;

    int vind[Nfaces][2];
    // triangle      vind[3][2] = { {0,1}, {1,2}, {2,0} };
    // quadrilateral vind[4][2] = { {0,1}, {1,2}, {2,3}, {3,0}};
    for(n=0;n<Nfaces;n++){
        vind[n][0] = n%Nv;
        vind[n][1] = (n+1)%(Nv);
    }

    face *myfaces = (face*) calloc(Klocal*Nfaces, sizeof(face));

    sk = 0;
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

    mr_mesh_parallelPairs(myfaces, Klocal*Nfaces, sizeof(face),
                          pairnumget, pairnumset, pairprocget,
                          pairmarry, compare_pairs);

    int k1, k2, f1, f2, p1, p2;

    for(n=0;n<Klocal*Nfaces;++n){

        k1 = myfaces[n].k1; f1 = myfaces[n].f1; p1 = myfaces[n].p1;
        k2 = myfaces[n].k2; f2 = myfaces[n].f2; p2 = myfaces[n].p2;

        if(p1!=procid)
            fprintf(stderr, "WARNING WRONG proc\n");

        mesh->EToE[k1][f1] = k2; /* adjacent element index */
        mesh->EToF[k1][f1] = f2; /* local face index of adjacent element */
        mesh->EToP[k1][f1] = p2; /* process id of adjacent element */

        if(p1!=p2){
            /* increment number of links */
            ++Npar[p2];
            mesh->parallCellNum++;
        }
    }

    free(myfaces);

    /* now build maps from incoming buffer to cells */
    mesh->cellIndexIn = IntVector_create(mesh->parallCellNum);

    int **Esend = (int**) calloc(nprocs, sizeof(int*));
    int **Erecv = (int**) calloc(nprocs, sizeof(int*));
    for(p2=0;p2<nprocs;++p2){
        if(Npar[p2]){
            Esend[p2] = IntVector_create(Npar[p2]);
            Erecv[p2] = IntVector_create(Npar[p2]);
        }
    }

    /* number of nodes adjacent to each process */
    int *skP = IntVector_create(nprocs);
    /* send coordinates in local order */
    for(k1=0;k1<Klocal;++k1){
        for(f1=0;f1<Nfaces;++f1){
            p2 = mesh->EToP[k1][f1];
            if(p2!=procid){
                Esend[p2][skP[p2]] = mesh->EToE[k1][f1];
                ++(skP[p2]);
            }
        }
    }

    MPI_Request Esendrequests[nprocs];
    MPI_Request Erecvrequests[nprocs];
    MPI_Status  status[nprocs];

    int cnt = 0;
    for(p2=0;p2<nprocs;++p2){
        if(p2!=procid && Npar[p2]!=0){
            int Nout = Npar[p2];

            MPI_Isend(Esend[p2], Nout, MPI_INT,    p2, 2666+p2, MPI_COMM_WORLD, Esendrequests+cnt);
            MPI_Irecv(Erecv[p2], Nout, MPI_INT,    p2, 2666+procid, MPI_COMM_WORLD, Erecvrequests+cnt);
            ++cnt;
        }
    }


    MPI_Waitall(cnt, Esendrequests, status);
    MPI_Waitall(cnt, Erecvrequests, status);

    /* now match up local cells with the requested (recv'ed cell) */
    sk = 0;
    for(p2=0;p2<nprocs;++p2){
        /* for each received face */
        for(n=0;n<skP[p2];++n){
            k1 = Erecv[p2][n]; /* adjacent element index */
            mesh->cellIndexIn[sk++] = k1;
        }
    }

    IntVector_free(skP);
    free(Esend);
    free(Erecv);

    /* now build maps from outgoing buffer to cells */
    mesh->cellIndexOut = IntVector_create(mesh->parallCellNum);

    sk = 0;
    for(p2=0;p2<mesh->nprocs;p2++){
        if(mesh->Npar[p2]){
            for(k=0;k<Klocal;k++){
                for(f1=0;f1<Nfaces;f1++){
                    if(mesh->EToP[k][f1]==p2 && p2!=mesh->procid)
                        mesh->cellIndexOut[sk++] = k;
                }
            }
        }
    }
}