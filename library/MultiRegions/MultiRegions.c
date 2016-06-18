//
// Created by li12242 on 16/6/18.
//

#include "MultiRegions.h"

/* private functions */
void SetCoord(StdRegions2d *shape, int K, double *VX, double *VY, int **EToV,
              double **GX, double **GY, double **x, double **y);

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

    /* set element info from parameter */
    mesh->K  = K;
    mesh->Nv = Nv;
    mesh->EToV = BuildIntMatrix(K, shape->Nv);
    int k,n;
    for(k=0;k<K;k++){
        for(n=0;n<shape->Nv;n++)
            mesh->EToV[k][n] = EToV[k][n];
    }

    /* standard element */
    mesh->StdElement = shape;

    mesh->GX = BuildMatrix(K, shape->Nv);
    mesh->GY = BuildMatrix(K, shape->Nv);

    mesh->x  = BuildMatrix(K, shape->Np);
    mesh->y  = BuildMatrix(K, shape->Np);

    SetCoord(shape, K, VX, VY, EToV, mesh->GX, mesh->GY, mesh->x, mesh->y);

    /* Redistribute the elements */
    LoadBalance(mesh);

    return mesh;
};


void FreeMultiReg2d(MultiReg2d *mesh){
    /* mesh info */
    DestroyIntMatrix(mesh->EToV);

    /* coordinate */
    DestroyMatrix(mesh->x);
    DestroyMatrix(mesh->y);
    DestroyMatrix(mesh->GX);
    DestroyMatrix(mesh->GY);
}

/**
 * @brief
 * Set the vertex and node coordinate based on EToV
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
 * x  | double[K][shape->Np]  | node coordinate
 * y  | double[K][shape->Np]  | node coordinate
 *
 */
void SetCoord(StdRegions2d *shape, int K, double *VX, double *VY, int **EToV,
              double **GX, double **GY, double **x, double **y){
    int n,k;
    /* vertex coordinate */
    for(k=0;k<K;k++){
        for(n=0;n<shape->Nv;n++){
            GX[k][n] = VX[EToV[k][n]];
            GY[k][n] = VY[EToV[k][n]];
        }
    }

    /* node coordinate */
    for(k=0;k<K;k++) {
        if (shape->Nv == 3) {
            MapTriCoor(shape, GX[k], GY[k], x[k], y[k]);
        } else if (shape->Nv == 4) {
            MapQuadCoor(shape, GX[k], GY[k], x[k], y[k]);
        } else {
            printf("fatal error: wrong number of vertex %d in StdRegions2d", shape->Nv);
        }
    }
}


/**
 * @brief
 * Redistribute the elements on each process
 *
 * @details
 *
 * @param[in] beginPos
 * @param[in] order order>0: year/month/date;order=0: date/month/year
 *
 * @return
 * return values:
 * name     | type     | description of value
 * -------- |----------|----------------------
 * car_id   | int      |
 *
 */
void LoadBalance(MultiReg2d *mesh){
    int n,p,k,v;

    int nprocs = mesh->nprocs;
    int procid = mesh->procid;
    int **EToV = mesh->EToV;
    double **GX = mesh->GX;
    double **GY = mesh->GY;

//    printf("Entering LoadBlance\n");

    /* number of vertex in each element */
    int Nverts = mesh->StdElement->Nv;

    /* number of elements in each process */
    int *Kprocs = BuildIntVector(nprocs);

    /* local number of elements */
    int Klocal = mesh->K;

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
    double **newx   = BuildMatrix(totalinK, Nverts);
    double **newy   = BuildMatrix(totalinK, Nverts);

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

    /* deallocate mem before assignment */
    DestroyMatrix(mesh->GX);
    DestroyMatrix(mesh->GY);
    DestroyIntMatrix(mesh->EToV);

    /* assignment of the new element list and vertex coordinate */
    mesh->GX = newx;
    mesh->GY = newy;
    mesh->EToV = newEToV;
    mesh->K =  totalinK;

    /* deallocate mem */

    free(inrequests);  free(xinrequests);  free(yinrequests);
    free(outrequests); free(xoutrequests); free(youtrequests);

//    printf("Finishing LoadBlance\n");
}