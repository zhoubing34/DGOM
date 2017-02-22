#include <StandCell/sc_stdcell.h>
#include "mr_mesh_nodeConnect.h"
#include "mr_mesh.h"

#define NODETOL 1e-10

/**
 * @brief Setup node connection
 * @details
 * Set fields of vmapM, vmapP
 *
 * @param [in,out] shape standard element
 *
 * @note
 * vmapP contains the mapping from nodes to incoming buffers, which is marked as negative.
 */
void mr_mesh_nodeConnect2d(parallMesh *mesh){

    stdCell *shape = mesh->cell;
    multiReg *region = mesh->region;

    /* constant variables */
    const int nprocs = mesh->nprocs;
    const int procid = mesh->procid;

    const int Nfaces = shape->Nfaces;
    const int Nfp    = shape->Nfp;
    const int Np     = shape->Np;
    const int K = mesh->grid->K;

    int m;
    int k1,f1, n1, id1, k2,f2,p2,n2;
    double x1, y1, x2, y2, d12;

    double **x = region->x;
    double **y = region->y;
    int **EToE = mesh->EToE;
    int **EToF = mesh->EToF;
    int **EToP = mesh->EToP;
    int **Fmask = shape->Fmask;
    int *Npar = mesh->Npar;

    /* allocation */
    int *vmapM = vector_int_create(K * Nfaces * Nfp);
    int *vmapP = vector_int_create(K * Nfaces * Nfp);

    mesh->vmapM = vmapM;
    mesh->vmapP = vmapP;
    mesh->parallNodeNum = mesh->parallCellNum*Nfp;

    int *parmapOUT = vector_int_create(mesh->parallNodeNum);
    mesh->nodeIndexOut = parmapOUT;


    /* first build local */
    for(k1=0;k1<K;++k1){

        /* get some information about the face geometries */
        for(f1=0;f1<Nfaces;++f1){

            /* volume -> face nodes */
            for(n1=0;n1<Nfp;++n1){
                id1 = n1+f1*Nfp+k1*Nfp*Nfaces; /* node index as face node */
                vmapM[id1] = Fmask[f1][n1] + k1*Np;
            }

            /* find neighbor */
            k2 = EToE[k1][f1]; /* adjacent element */
            f2 = EToF[k1][f1]; /* adjacent face index */
            p2 = EToP[k1][f1]; /* process id */

            if(k1==k2 || procid!=p2 ){
                /* adjacent to the boundary */
                for(n1=0;n1<Nfp;++n1){
                    id1 = n1+f1*Nfp+k1*Nfp*Nfaces;
                    /* set itself as the adjacent node index */
                    vmapP[id1] = k1*Np + Fmask[f1][n1];
                }
            }else{
                /* inner boundary */
                for(n1=0;n1<Nfp;++n1){
                    id1 = n1+f1*Nfp+k1*Nfp*Nfaces;
                    x1 = x[k1][Fmask[f1][n1]];
                    y1 = y[k1][Fmask[f1][n1]];

                    /* loop over of adjacent face node */
                    for(n2=0;n2<Nfp;++n2){

                        x2 = x[k2][Fmask[f2][n2]];
                        y2 = y[k2][Fmask[f2][n2]];

                        /* find normalized distance between these nodes */
                        /* [ use sJk as a measure of edge length (ignore factor of 2) ] */

                        d12 = ((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
                        /* judge adjacent node */
                        if(d12<NODETOL){
                            vmapP[id1] = k2*Np + Fmask[f2][n2];
                        }
                    }
                }
            }
        }
    }

    /* create map fron variables to incoming nodes */
    int parcnt = -1;
    for(p2=0;p2<nprocs;++p2){
        for(k1=0;k1<K;++k1){
            for(f1=0;f1<Nfaces;++f1){
                if(EToP[k1][f1]==p2 && p2!=procid){
                    for(n1=0;n1<Nfp;++n1){
                        id1 = n1+f1*Nfp+k1*Nfp*Nfaces;
                        vmapP[id1] = parcnt;
                        --parcnt;
                    }
                }
            }
        }
    }

    /* now build maps from incoming nodes to variables */
    double **xsend = (double**) calloc((size_t) nprocs, sizeof(double*));
    double **ysend = (double**) calloc((size_t) nprocs, sizeof(double*));
    double **xrecv = (double**) calloc((size_t) nprocs, sizeof(double*));
    double **yrecv = (double**) calloc((size_t) nprocs, sizeof(double*));

    int **Esend = (int**) calloc((size_t) nprocs, sizeof(int*));
    int **Fsend = (int**) calloc((size_t) nprocs, sizeof(int*));
    int **Erecv = (int**) calloc((size_t) nprocs, sizeof(int*));
    int **Frecv = (int**) calloc((size_t) nprocs, sizeof(int*));

    for(p2=0;p2<nprocs;++p2){
        if(Npar[p2]){
            xsend[p2] = vector_double_create(Npar[p2] * Nfp);
            ysend[p2] = vector_double_create(Npar[p2] * Nfp);
            Esend[p2] = vector_int_create(Npar[p2] * Nfp);
            Fsend[p2] = vector_int_create(Npar[p2] * Nfp);

            xrecv[p2] = vector_double_create(Npar[p2] * Nfp);
            yrecv[p2] = vector_double_create(Npar[p2] * Nfp);
            Erecv[p2] = vector_int_create(Npar[p2] * Nfp);
            Frecv[p2] = vector_int_create(Npar[p2] * Nfp);
        }
    }

    /* number of nodes adjacent to each process */
    int *skP = vector_int_create(nprocs);

    /* send coordinates in local order */
    int cnt = 0;
    for(k1=0;k1<K;++k1){
        for(f1=0;f1<Nfaces;++f1){
            p2 = EToP[k1][f1];
            if(p2!=procid){
                for(n1=0;n1<Nfp;++n1){
                    xsend[p2][skP[p2]] = x[k1][Fmask[f1][n1]];
                    ysend[p2][skP[p2]] = y[k1][Fmask[f1][n1]];
                    Esend[p2][skP[p2]] = EToE[k1][f1];
                    Fsend[p2][skP[p2]] = EToF[k1][f1];
                    ++(skP[p2]);
                }
            }
        }
    }

    MPI_Request *xsendrequests = (MPI_Request*) calloc((size_t) nprocs, sizeof(MPI_Request));
    MPI_Request *ysendrequests = (MPI_Request*) calloc((size_t) nprocs, sizeof(MPI_Request));
    MPI_Request *xrecvrequests = (MPI_Request*) calloc((size_t) nprocs, sizeof(MPI_Request));
    MPI_Request *yrecvrequests = (MPI_Request*) calloc((size_t) nprocs, sizeof(MPI_Request));
    MPI_Request *Esendrequests = (MPI_Request*) calloc((size_t) nprocs, sizeof(MPI_Request));
    MPI_Request *Fsendrequests = (MPI_Request*) calloc((size_t) nprocs, sizeof(MPI_Request));
    MPI_Request *Erecvrequests = (MPI_Request*) calloc((size_t) nprocs, sizeof(MPI_Request));
    MPI_Request *Frecvrequests = (MPI_Request*) calloc((size_t) nprocs, sizeof(MPI_Request));

    MPI_Status  *status = (MPI_Status*) calloc((size_t) nprocs, sizeof(MPI_Status));

    cnt = 0;
    for(p2=0;p2<nprocs;++p2){
        if(p2!=procid && Npar[p2]!=0){
            int Nout = Npar[p2]*Nfp;

            MPI_Isend(xsend[p2], Nout, MPI_DOUBLE, p2,  666+p2, MPI_COMM_WORLD, xsendrequests+cnt);
            MPI_Isend(ysend[p2], Nout, MPI_DOUBLE, p2, 1666+p2, MPI_COMM_WORLD, ysendrequests+cnt);
            MPI_Isend(Esend[p2], Nout, MPI_INT,    p2, 2666+p2, MPI_COMM_WORLD, Esendrequests+cnt);
            MPI_Isend(Fsend[p2], Nout, MPI_INT,    p2, 3666+p2, MPI_COMM_WORLD, Fsendrequests+cnt);

            MPI_Irecv(xrecv[p2], Nout, MPI_DOUBLE, p2,  666+procid, MPI_COMM_WORLD, xrecvrequests+cnt);
            MPI_Irecv(yrecv[p2], Nout, MPI_DOUBLE, p2, 1666+procid, MPI_COMM_WORLD, yrecvrequests+cnt);
            MPI_Irecv(Erecv[p2], Nout, MPI_INT,    p2, 2666+procid, MPI_COMM_WORLD, Erecvrequests+cnt);
            MPI_Irecv(Frecv[p2], Nout, MPI_INT,    p2, 3666+procid, MPI_COMM_WORLD, Frecvrequests+cnt);
            ++cnt;
        }
    }

    MPI_Waitall(cnt, xsendrequests, status);
    MPI_Waitall(cnt, ysendrequests, status);
    MPI_Waitall(cnt, Esendrequests, status);
    MPI_Waitall(cnt, Fsendrequests, status);
    MPI_Waitall(cnt, xrecvrequests, status);
    MPI_Waitall(cnt, yrecvrequests, status);
    MPI_Waitall(cnt, Erecvrequests, status);
    MPI_Waitall(cnt, Frecvrequests, status);

    /* add up the total number of outgoing/ingoing nodes */
    int parNtotalout = 0;
    for(p2=0;p2<nprocs;++p2) {
        parNtotalout += skP[p2];
    }

    if(mesh->parallNodeNum != parNtotalout){
        printf("%s (%s): %d\n"
                       "conflicts occurs between parallel node number "
                       "and parallel cell number\n",
               __FUNCTION__, __FILE__, __LINE__);
        exit(-1);
    }

    /* now match up local nodes with the requested (recv'ed nodes) */
    int sk = 0;
    for(p2=0;p2<nprocs;++p2){
        /* for each received face */
        for(m=0;m<skP[p2];++m){
            k1 = Erecv[p2][m]; /* adjacent element index */
            f1 = Frecv[p2][m]; /* adjacent face */
            x2 = xrecv[p2][m];
            y2 = yrecv[p2][m];

            for(n1=0;n1<Nfp;++n1){

                x1 = x[k1][Fmask[f1][n1]];
                y1 = y[k1][Fmask[f1][n1]];
                d12 = ((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
                if(d12<NODETOL){
                    /* sk and node index determine the map relationship */
                    parmapOUT[sk++] = (k1*Np+Fmask[f1][n1]);
                }
            }
        }
    }

    /* deallocate mem */
    vector_int_free(skP);

    free(xsendrequests); free(ysendrequests);
    free(xrecvrequests); free(yrecvrequests);
    free(Esendrequests); free(Fsendrequests);
    free(Erecvrequests); free(Frecvrequests);

    free(status);
}

