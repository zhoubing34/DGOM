//
// Created by li12242 on 16/6/18.
//

#include "MultiRegions.h"

/* private functions */
void SetVetxCoord(StdRegions2d *shape, int K, double *VX, double *VY, int **EToV, double **GX, double **GY);
void LoadBalance(StdRegions2d *shape, int K, int **EToV, double **GX, double **GY,
                 int *newK, int ***newEToV, double ***newx, double ***newy);
void SetFacePair(StdRegions2d *shape, int Klocal,
                 int **EToV, int **EToE, int **EToF, int **EToP,
                 int *Npar, int ***newParK, int ***newParF);

/* private functions and structures for FacePair */
int compare_pairs(const void *obj1, const void *obj2);
int pairprocget(const void *obj1);
int pairnumget(const void *obj1);
void pairnumset(const void *obj1, int g);
void pairmarry(const void *obj1, const void *obj2);

/**
 * @brief
 * Generation of two dimension mesh
 *
 * @param [StdRegions2d*] shape standard element
 * @param [int]           K     number of element
 * @param [int]           Nv    number of vertex
 * @param [int**]         EToV  element to vertex list
 * @param [double*]       VX    vertex coordinate
 * @param [double*]       VY    vertex coordinate
 *
 * @return
 * name     | type     | description of value
 * -------- |----------|----------------------
 * mesh     | MultiReg2d* | mesh object
 *
 * @note
 * The vertex index in EToV is start from 0.
 */
MultiReg2d* GenMultiReg2d(StdRegions2d *shape, int K, int Nv,
                          int **EToV, double *VX, double *VY){

    MultiReg2d *mesh = (MultiReg2d *)calloc(1, sizeof(MultiReg2d));

    MPI_Comm_rank(MPI_COMM_WORLD, &mesh->procid);
    MPI_Comm_size(MPI_COMM_WORLD, &mesh->nprocs);

//    printf("procid = %d, nprocs = %d\n", mesh->procid, mesh->nprocs);

    /* standard element */
    mesh->StdElement = shape;

    double **GX = BuildMatrix(K, shape->Nv);
    double **GY = BuildMatrix(K, shape->Nv);

    SetVetxCoord(shape, K, VX, VY, EToV, GX, GY);

    /* Redistribute the elements */
    LoadBalance(shape, K, EToV, GX, GY, &(mesh->K), &(mesh->EToV), &(mesh->GX), &(mesh->GY));

    DestroyMatrix(GX);
    DestroyMatrix(GY);

    /* Setup element connection, EToE,EToP,EToF & Npar,parK,parF */
    mesh->Npar = BuildIntVector(mesh->nprocs);

    mesh->EToE = BuildIntMatrix(mesh->K, shape->Nfaces);
    mesh->EToF = BuildIntMatrix(mesh->K, shape->Nfaces);
    mesh->EToP = BuildIntMatrix(mesh->K, shape->Nfaces);

    SetFacePair(shape, mesh->K, mesh->EToV, mesh->EToE, mesh->EToF, mesh->EToP,
                mesh->Npar, &mesh->parK, &mesh->parF);

    /* Setup boundary nodes connection */

    /**/
//    double **x  = BuildMatrix(mesh->K, shape->Np);
//    double **y  = BuildMatrix(mesh->K, shape->Np);

    return mesh;
};


void FreeMultiReg2d(MultiReg2d *mesh){
    /* mesh info */
    DestroyIntMatrix(mesh->EToV);

    /* coordinate */
    DestroyMatrix(mesh->GX);
    DestroyMatrix(mesh->GY);

    /* element connection */
    DestroyIntVector(mesh->Npar);
    DestroyIntMatrix(mesh->EToE);
    DestroyIntMatrix(mesh->EToF);
    DestroyIntMatrix(mesh->EToP);
}

/**
 * @brief
 * Set the vertex coordinate based on EToV
 *
 * @param [StdRegions2d*]   shape standard element object
 * @param [int]             K     number of elements
 * @param [double*]         VX    input coordinate of vertex
 * @param [double*]         VY    input coordinate of vertex
 * @param [int**]           EToV  element to vertex list
 *
 * @return
 * return values:
 * name     | type     | description of value
 * -------- |----------|----------------------
 * GX | double[K][shape->Nv]  | vertex coordinate
 * GX | double[K][shape->Nv]  | vertex coordinate
 *
 */
void SetVetxCoord(StdRegions2d *shape, int K, double *VX, double *VY, int **EToV, double **GX, double **GY){
    int n,k;
    /* vertex coordinate */
    for(k=0;k<K;k++){
        for(n=0;n<shape->Nv;n++){
            GX[k][n] = VX[EToV[k][n]];
            GY[k][n] = VY[EToV[k][n]];
        }
    }

//    /* node coordinate */
//    for(k=0;k<K;k++) {
//        if (shape->Nv == 3) {
//            MapTriCoor(shape, GX[k], GY[k], x[k], y[k]);
//        } else if (shape->Nv == 4) {
//            MapQuadCoor(shape, GX[k], GY[k], x[k], y[k]);
//        } else {
//            printf("fatal error: wrong number of vertex %d in StdRegions2d", shape->Nv);
//        }
//    }
}


/**
 * @brief
 * Setup element connection
 * @details
 * Set EToE,EToP,EToF & Npar,parK,parF
 *
 * @param [StdRegions2d*] shape standard element
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
 * newParK | int[nprocs][*] | index of local element adjacent other process
 * newParF | int[nprocs][*] | index of local face adjacent other process
 *
 * @note
 * EToE, EToF, EToP and Npar should be allocated before calling this function, while
 * newParK and newParF is allocated inside of the function. Send the address of the
 * matrix pointer as parameters.
 */
void SetFacePair(StdRegions2d *shape, int Klocal,
                 int **EToV, int **EToE, int **EToF, int **EToP,
                 int *Npar, int ***newParK, int ***newParF){

    int procid, nprocs;
    int Nfaces = shape->Nfaces;

    int n, k, e, sk;

    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

//    triangle      vnum[3][2] = { {0,1}, {1,2}, {2,0} };
//    quadrilateral vnum[4][2] = { {0,1}, {1,2}, {2,3}, {3,0}};
    int **vnum = BuildIntMatrix(shape->Nfaces, 2);
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

        EToE[k1][f1] = k2;
        EToF[k1][f1] = f2;
        EToP[k1][f1] = p2;

        if(p1!=p2){
            /* increment number of links */
            ++Npar[p2];
        }
    }

    int parcnt=-1;
    for(p2=0;p2<nprocs;p2++){
        for(k=0;k<Klocal;k++){
            for(f1=0;f1<Nfaces;f1++){
                if(EToP[k][f1]==p2 && p2!=procid) {
                    EToE[k][f1] = parcnt;
                    --parcnt;
                }
            }
        }
    }

    int **parK = (int**) calloc(nprocs, sizeof(int*));
    int **parF = (int**) calloc(nprocs, sizeof(int*));
    for(p2=0;p2<nprocs;++p2){
        parK[p2] = BuildIntVector(Npar[p2]);
        parF[p2] = BuildIntVector(Npar[p2]);
        Npar[p2] = 0;
        for(n=0;n<Klocal*Nfaces;++n){
            if(myfaces[n].p2==p2 && p2!=procid){
                k1 = myfaces[n].k1, f1 = myfaces[n].f1;
                parK[p2][Npar[p2]  ] = k1;
                parF[p2][Npar[p2]++] = f1;
            }
        }
    }
    free(myfaces);

    /* assignment */
    *newParK = parK; *newParF = parF;
}


void SetNodePair(){}


/**
 * @brief
 * Redistribute the elements on each process
 *
 * @details
 * Call ParMetis library to repart the mesh
 *
 * @param [StdRegions2d*] shape     standard elements, triangle or quadrilateral
 * @param [int]         K           original number of elements
 * @param [int**]       EToV[K][shape->Nv] element to vertex list
 * @param [double**]    GX[K][shape->Nv]   vertex coordinate
 * @param [double**]    GY[K][shape->Nv]   vertex coordinate
 *
 * @return
 * return values:
 * name     | type     | description of value
 * -------- |----------|----------------------
 * newK   | int* | new number of elements after mesh redistribution
 * nEToV  | int*** | new element to vertex list after mesh redistribution
 * newGX  | double*** | new vertex coordinate after mesh redistribution
 * newGY  | double*** | new vertex coordinate after mesh redistribution
 *
 */
void LoadBalance(StdRegions2d *shape, int K, int **EToV, double **GX, double **GY,
                 int *newK, int ***nEToV, double ***newGX, double ***newGY){
    int n,p,k,v;

    int nprocs, procid;
//    int **EToV = mesh->EToV;
//    double **GX = mesh->GX;
//    double **GY = mesh->GY;

    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

//    printf("Entering LoadBlance\n");

    /* number of vertex in each element */
    int Nverts = shape->Nv;

    /* number of elements in each process */
    int *Kprocs = BuildIntVector(nprocs);

    /* local number of elements */
    int Klocal = K;

    /* find number of elements on all processors */
    MPI_Allgather(&Klocal, 1, MPI_INT, Kprocs, 1, MPI_INT, MPI_COMM_WORLD);

    /* element distribution -- cumulative element count on processes */
    idxtype *elmdist = idxmalloc(nprocs+1, "elmdist");

    elmdist[0] = 0;
    for(p=0;p<nprocs;++p)
        elmdist[p+1] = elmdist[p] + Kprocs[p];

    DestroyIntVector(Kprocs);

    /* list of element starts */
    idxtype *eptr = idxmalloc(Klocal+1, "eptr");

    eptr[0] = 0;
    for(k=0;k<Klocal;++k)
        eptr[k+1] = eptr[k] + Nverts;

    /* local element to vertex */
    idxtype *eind = idxmalloc(Nverts*Klocal, "eind");

    for(k=0;k<Klocal;++k)
        for(n=0;n<Nverts;++n)
            eind[k*Nverts+n] = EToV[k][n];

    /* weight per element */
    idxtype *elmwgt = idxmalloc(Klocal, "elmwgt");

    for(k=0;k<Klocal;++k)
        elmwgt[k] = 1;

    /* weight flag */
    int wgtflag = 0;

    /* number flag (1=fortran, 0=c) */
    int numflag = 0;

    /* ncon = 1 */
    int ncon = 1;

    /* nodes on element face */
    int ncommonnodes = 2;

    /* number of partitions */
    int nparts = nprocs;

    /* tpwgts */
    float *tpwgts = (float*) calloc(Klocal, sizeof(float));

    for(k=0;k<Klocal;++k)
        tpwgts[k] = (float)1./nprocs;

    float ubvec[MAXNCON];

    for (n=0; n<ncon; ++n)
        ubvec[n] = UNBALANCE_FRACTION;

    int options[10];

    options[0] = 1;
    options[PMV3_OPTION_DBGLVL] = 7;
    options[PMV3_OPTION_SEED] = 0;

    int edgecut;

    /** the process index of redistributing each element */
    idxtype *part = idxmalloc(Klocal, "part");

    MPI_Comm comm;
    MPI_Comm_dup(MPI_COMM_WORLD, &comm);

    /* call ParMetis lib to make mesh partition */
    ParMETIS_V3_PartMeshKway
            (elmdist,
             eptr,
             eind,
             elmwgt,
             &wgtflag,
             &numflag,
             &ncon,
             &ncommonnodes,
             &nparts,
             tpwgts,
             ubvec,
             options,
             &edgecut,
             part,
             &comm);

    free(tpwgts);

    /** number of elements send to each process */
    int *outK = (int*) calloc(nprocs, sizeof(int));
    /** number of elements to receive from each process */
    int *inK = (int*) calloc(nprocs, sizeof(int));


    MPI_Request *inrequests = (MPI_Request*) calloc(nprocs, sizeof(MPI_Request));
    MPI_Request *outrequests = (MPI_Request*) calloc(nprocs, sizeof(MPI_Request));
    MPI_Request *xinrequests = (MPI_Request*) calloc(nprocs, sizeof(MPI_Request));
    MPI_Request *xoutrequests = (MPI_Request*) calloc(nprocs, sizeof(MPI_Request));
    MPI_Request *yinrequests = (MPI_Request*) calloc(nprocs, sizeof(MPI_Request));
    MPI_Request *youtrequests = (MPI_Request*) calloc(nprocs, sizeof(MPI_Request));

    for(k=0;k<Klocal;++k)
        ++outK[part[k]];

    /* get count of incoming elements from each process in inK */
    MPI_Alltoall(outK, 1, MPI_INT, inK,  1, MPI_INT, MPI_COMM_WORLD);

    /** total element number receive from process */
    int totalinK = 0;
    for(p=0;p<nprocs;++p){
        totalinK += inK[p];
    }

    /* receive the element to vertex list and vertex coordinates from each process, including itself */
    int **newEToV = BuildIntMatrix(totalinK, Nverts);
    double **newx = BuildMatrix(totalinK, Nverts);
    double **newy = BuildMatrix(totalinK, Nverts);

    int cnt = 0;
    for(p=0; p<nprocs; ++p){
        MPI_Irecv(newEToV[cnt], Nverts*inK[p], MPI_INT, p, 666+p, MPI_COMM_WORLD, inrequests+p);
        MPI_Irecv(newx[cnt], Nverts*inK[p], MPI_DOUBLE, p, 1666+p, MPI_COMM_WORLD, xinrequests+p);
        MPI_Irecv(newy[cnt], Nverts*inK[p], MPI_DOUBLE, p, 2666+p, MPI_COMM_WORLD, yinrequests+p);
        cnt = cnt + inK[p];
    }

    int **outlist = (int**) calloc(nprocs, sizeof(int*));
    double **xoutlist = (double**) calloc(nprocs, sizeof(double*));
    double **youtlist = (double**) calloc(nprocs, sizeof(double*));

    /* send the element to vertex list and vertex coordinates to each process, including itself */
    for(p=0;p<nprocs;++p){
        int cnt = 0;
        outlist[p]  = BuildIntVector(Nverts*outK[p]);
        xoutlist[p]  = BuildVector(Nverts*outK[p]);
        youtlist[p]  = BuildVector(Nverts*outK[p]);

        for(k=0;k<Klocal;++k)
            if(part[k]==p){
                for(v=0;v<Nverts;++v){
                    outlist[p][cnt] = EToV[k][v];
                    xoutlist[p][cnt] = GX[k][v];
                    youtlist[p][cnt] = GY[k][v];
                    ++cnt;
                }
            }

        MPI_Isend(outlist[p], Nverts*outK[p], MPI_INT, p, 666+procid, MPI_COMM_WORLD, outrequests+p);
        MPI_Isend(xoutlist[p], Nverts*outK[p], MPI_DOUBLE, p, 1666+procid, MPI_COMM_WORLD, xoutrequests+p);
        MPI_Isend(youtlist[p], Nverts*outK[p], MPI_DOUBLE, p, 2666+procid, MPI_COMM_WORLD, youtrequests+p);

        /* deallocate mem */
        DestroyIntVector(outlist[p]);
        DestroyVector(xoutlist[p]);
        DestroyVector(youtlist[p]);
    }

    free(outlist); free(xoutlist); free(youtlist);

    /* waite for all send/recv messages */
    MPI_Status *instatus = (MPI_Status*) calloc(nprocs, sizeof(MPI_Status));
    MPI_Status *outstatus = (MPI_Status*) calloc(nprocs, sizeof(MPI_Status));

    MPI_Waitall(nprocs,  inrequests, instatus);
    MPI_Waitall(nprocs, xinrequests, instatus);
    MPI_Waitall(nprocs, yinrequests, instatus);

    MPI_Waitall(nprocs,  outrequests, outstatus);
    MPI_Waitall(nprocs, xoutrequests, outstatus);
    MPI_Waitall(nprocs, youtrequests, outstatus);

    /* assignment of the new element list and vertex coordinate */
    *nEToV = newEToV;
    *newGX = newx;
    *newGY = newy;
    *newK = totalinK;

    /* deallocate mem */

    free(inrequests);  free(xinrequests);  free(yinrequests);
    free(outrequests); free(xoutrequests); free(youtrequests);

//    printf("Finishing LoadBlance\n");
}


/* ======================================================================== */
/* private functions and structures for FacePair */
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