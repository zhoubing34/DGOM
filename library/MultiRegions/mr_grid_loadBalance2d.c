#include <parmetisbin.h>
#include "mr_grid.h"

/**
 * @brief
 * Redistribute the elements on each process
 *
 * @details
 * Call ParMetis library to redistribute the mesh,
 * in order to achieve the load balance
 *
 * @param [in,out] grid geometry grid object
 *
 * @note
 * The properties of EToV and K is updated.
 */
void mr_grid_loadBalance2d(geoGrid *grid){

    int n,p,k,v;

    /* MPI process */
    const int nprocs = grid->nprocs;
    const int procid = grid->procid;

    /* number of vertex in each element */
    const int Nv = grid->cell->Nv;

    /* number of elements in each process */
    int Kprocs[nprocs];

    /* local number of elements */
    int Klocal = grid->K;
    int **EToV = grid->EToV;

    /* find number of elements on all processors */
    MPI_Allgather(&Klocal, 1, MPI_INT, Kprocs, 1, MPI_INT, MPI_COMM_WORLD);

    /* element distribution -- cumulative element count on processes */
    idxtype *elmdist = idxmalloc(nprocs+1, "elmdist");

    elmdist[0] = 0;
    for(p=0;p<nprocs;++p)
        elmdist[p+1] = elmdist[p] + Kprocs[p];

    /* list of element starts */
    idxtype *eptr = idxmalloc(Klocal+1, "eptr");

    eptr[0] = 0;
    for(k=0;k<Klocal;++k)
        eptr[k+1] = eptr[k] + Nv;

    /* local element to vertex */
    idxtype *eind = idxmalloc(Nv*Klocal, "eind");
    for(k=0;k<Klocal;++k)
        for(n=0;n<Nv;++n)
            eind[k*Nv+n] = EToV[k][n];

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
    int *inK = (int*) calloc(nprocs, sizeof(int));;

    for(k=0;k<Klocal;++k)
        ++outK[part[k]];

    /* get count of incoming elements from each process in inK */
    MPI_Alltoall(outK, 1, MPI_INT, inK,  1, MPI_INT, MPI_COMM_WORLD);

    /** total element number receive from process */
    int totalinK = 0;
    for(p=0;p<nprocs;++p){
        totalinK += inK[p];
    }

    /* receive the EToV from each process, including itself */
    int **newEToV = IntMatrix_create(totalinK, Nv);
    MPI_Request inRequests[nprocs], outRequests[nprocs];

    int cnt = 0;
    for(p=0; p<nprocs; ++p){
        MPI_Irecv(newEToV[cnt], Nv*inK[p], MPI_INT, p, 666+p, MPI_COMM_WORLD, inRequests+p);
        cnt = cnt + inK[p];
    }

    /* send the element to vertex list and vertex coordinates to each process, including itself */
    for(p=0;p<nprocs;++p){
        cnt = 0;
        int outlist[Nv*outK[p]];
        for(k=0;k<Klocal;++k)
            if(part[k]==p){
                for(v=0;v<Nv;++v){
                    outlist[cnt] = EToV[k][v];
                    ++cnt;
                }
            }

        MPI_Isend(outlist, Nv*outK[p], MPI_INT, p, 666+procid, MPI_COMM_WORLD, outRequests+p);
    }

    /* waite for all send/recv messages */
    MPI_Status instatus[nprocs], outstatus[nprocs];

    MPI_Waitall(nprocs, inRequests, instatus);
    MPI_Waitall(nprocs, outRequests, outstatus);

    /* assignment of the new element list and vertex coordinate */
    IntMatrix_free(grid->EToV);

    grid->EToV = newEToV;
    grid->K  = totalinK;

    free(outK);
    free(inK);
}