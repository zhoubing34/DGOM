#include "Convection2d/Convection2d.h"
#include <mpi.h>

/**
 * @brief
 * Check element parallel pair
 *
 * @details
 *
 * @author
 * li12242, Tianjin University, li12242@tju.edu.cn
 *
 */
void PrintBuilMapsLog(FILE *fig, Mesh *mesh){

    int parEtotalout = mesh->parNtotalout/p_Nfp;
    int n, f;
    for(n=0;n<parEtotalout;n++){
        fprintf(fig, "%d, ", mesh->elemapOUT[n]);
    }
    fprintf(fig, "\n");

    fprintf(fig, "EToE = \n");
    for(n=0;n<mesh->K;n++){
        for(f=0;f<p_Nfaces;f++)
            fprintf(fig, "%d, ", mesh->EToE[n][f]);

        fprintf(fig, "\n");
    }
}

/**
 * @brief
 * Build node pairing and collect outgoing node index
 *
 * @details
 * Build node pairing for vmapM & vmapP, and collect the outgoing node index for parmapOUT
 *
 * @author
 * li12242, Tianjin University, li12242@tju.edu.cn
 */
void BuildMaps(Mesh *mesh){

    int nprocs = mesh->nprocs;
    int procid = mesh->procid;

#if defined DEBUG
    printf("Procs %d: Entering BuildMaps\n", procid);
#endif

    int K = mesh->K;
    int Nfaces = mesh->Nfaces;

    mesh->vmapM = BuildIntVector(p_Nfp*p_Nfaces*K);
    mesh->vmapP = BuildIntVector(p_Nfp*p_Nfaces*K);

    int m;
    int k1,f1, n1, id1, k2,f2,p2,n2;

    double x1, y1, x2, y2, d12;

    double *nxk = BuildVector(Nfaces);
    double *nyk = BuildVector(Nfaces);
    double *sJk = BuildVector(Nfaces);

    /* first build local */
    for(k1=0;k1<K;++k1){

        /* get some information about the face geometries */
        Normals(mesh, k1, nxk, nyk, sJk);

        for(f1=0;f1<Nfaces;++f1){

            /* volume -> face nodes */
            for(n1=0;n1<p_Nfp;++n1){
                id1 = n1+f1*p_Nfp+k1*p_Nfp*p_Nfaces; /* node index as face node */
                mesh->vmapM[id1] = mesh->Fmask[f1][n1] + k1*p_Np;
            }

            /* find neighbor */
            k2 = mesh->EToE[k1][f1]; /* adjacent element */
            f2 = mesh->EToF[k1][f1]; /* adjacent face index */
            p2 = mesh->EToP[k1][f1]; /* process id */

            if(k1==k2 || procid!=p2 ){
                /* adjacent to the boundary */
                for(n1=0;n1<p_Nfp;++n1){
                    id1 = n1+f1*p_Nfp+k1*p_Nfp*p_Nfaces;
                    /* set itself as the adjacent node index */
                    mesh->vmapP[id1] = k1*p_Np + mesh->Fmask[f1][n1];
                }
            }else{
                /* inner boundary */
                for(n1=0;n1<p_Nfp;++n1){
                    id1 = n1+f1*p_Nfp+k1*p_Nfp*p_Nfaces;
                    x1 = mesh->x[k1][mesh->Fmask[f1][n1]];
                    y1 = mesh->y[k1][mesh->Fmask[f1][n1]];

                    /* loop over of adjacent face node */
                    for(n2=0;n2<p_Nfp;++n2){

                        x2 = mesh->x[k2][mesh->Fmask[f2][n2]];
                        y2 = mesh->y[k2][mesh->Fmask[f2][n2]];

                        /* find normalized distance between these nodes */
                        /* [ use sJk as a measure of edge length (ignore factor of 2) ] */

                        d12 = ((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2))/(sJk[f1]*sJk[f1]);
                        /* judge adjacent node */
                        if(d12<NODETOL){
                            mesh->vmapP[id1] = k2*p_Np + mesh->Fmask[f2][n2];
                        }
                    }
                }
            }
        }
    }

    /* now build parallel maps */
    double **xsend = (double**) calloc(nprocs, sizeof(double*));
    double **ysend = (double**) calloc(nprocs, sizeof(double*));
    double **xrecv = (double**) calloc(nprocs, sizeof(double*));
    double **yrecv = (double**) calloc(nprocs, sizeof(double*));

    int **Esend = (int**) calloc(nprocs, sizeof(int*));
    int **Fsend = (int**) calloc(nprocs, sizeof(int*));
    int **Erecv = (int**) calloc(nprocs, sizeof(int*));
    int **Frecv = (int**) calloc(nprocs, sizeof(int*));

    for(p2=0;p2<nprocs;++p2){
        if(mesh->Npar[p2]){
            xsend[p2] = BuildVector(mesh->Npar[p2]*p_Nfp);
            ysend[p2] = BuildVector(mesh->Npar[p2]*p_Nfp);
            Esend[p2] = BuildIntVector(mesh->Npar[p2]*p_Nfp);
            Fsend[p2] = BuildIntVector(mesh->Npar[p2]*p_Nfp);

            xrecv[p2] = BuildVector(mesh->Npar[p2]*p_Nfp);
            yrecv[p2] = BuildVector(mesh->Npar[p2]*p_Nfp);
            Erecv[p2] = BuildIntVector(mesh->Npar[p2]*p_Nfp);
            Frecv[p2] = BuildIntVector(mesh->Npar[p2]*p_Nfp);
        }
    }

    /* number of nodes adjacent to each process */
    int *skP = BuildIntVector(nprocs);

    /* send coordinates in local order */
    int cnt = 0;
    for(k1=0;k1<K;++k1){
        for(f1=0;f1<p_Nfaces;++f1){
            p2 = mesh->EToP[k1][f1];
            if(p2!=procid){
                for(n1=0;n1<p_Nfp;++n1){
                    xsend[p2][skP[p2]] = mesh->x[k1][mesh->Fmask[f1][n1]];
                    ysend[p2][skP[p2]] = mesh->y[k1][mesh->Fmask[f1][n1]];
                    Esend[p2][skP[p2]] = mesh->EToE[k1][f1];
                    Fsend[p2][skP[p2]] = mesh->EToF[k1][f1];
                    ++(skP[p2]);
                }
            }
        }
    }

    MPI_Request *xsendrequests = (MPI_Request*) calloc(nprocs, sizeof(MPI_Request));
    MPI_Request *ysendrequests = (MPI_Request*) calloc(nprocs, sizeof(MPI_Request));
    MPI_Request *xrecvrequests = (MPI_Request*) calloc(nprocs, sizeof(MPI_Request));
    MPI_Request *yrecvrequests = (MPI_Request*) calloc(nprocs, sizeof(MPI_Request));
    MPI_Request *Esendrequests = (MPI_Request*) calloc(nprocs, sizeof(MPI_Request));
    MPI_Request *Fsendrequests = (MPI_Request*) calloc(nprocs, sizeof(MPI_Request));
    MPI_Request *Erecvrequests = (MPI_Request*) calloc(nprocs, sizeof(MPI_Request));
    MPI_Request *Frecvrequests = (MPI_Request*) calloc(nprocs, sizeof(MPI_Request));

    MPI_Status  *status = (MPI_Status*) calloc(nprocs, sizeof(MPI_Status));

    cnt = 0;
    for(p2=0;p2<nprocs;++p2){
        if(p2!=procid && mesh->Npar[p2]!=0){
            int Nout = mesh->Npar[p2]*p_Nfp;

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
    mesh->parNtotalout = 0;
    for(p2=0;p2<nprocs;++p2)
        mesh->parNtotalout += skP[p2]*p_Nfields;

    mesh->parmapOUT = BuildIntVector(mesh->parNtotalout);

    /* now match up local nodes with the requested (recv'ed nodes) */
    int sk = 0, fld;
    for(p2=0;p2<nprocs;++p2){
        /* for each received face */
        for(m=0;m<skP[p2];++m){
            k1 = Erecv[p2][m]; /* adjacent element index */
            f1 = Frecv[p2][m]; /* adjacent face */
            x2 = xrecv[p2][m];
            y2 = yrecv[p2][m];

            Normals(mesh, k1, nxk, nyk, sJk);

            for(n1=0;n1<p_Nfp;++n1){

                x1 = mesh->x[k1][mesh->Fmask[f1][n1]];
                y1 = mesh->y[k1][mesh->Fmask[f1][n1]];

                d12 = ((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2))/(sJk[f1]*sJk[f1]);

                if(d12<NODETOL){
                    for(fld=0;fld<p_Nfields;++fld){
                        /* sk and node index determine the map relationship */
                        mesh->parmapOUT[sk++] = p_Nfields*(k1*p_Np+mesh->Fmask[f1][n1]) + fld;
                    }
                }
            }
        }
    }


    /* create incoming node map */
    int parcnt = -1;
    for(p2=0;p2<nprocs;++p2){
        for(k1=0;k1<K;++k1){
            for(f1=0;f1<p_Nfaces;++f1){
                if(mesh->EToP[k1][f1]==p2 && p2!=procid){
                    for(n1=0;n1<p_Nfp;++n1){
                        id1 = n1+f1*p_Nfp+k1*p_Nfp*p_Nfaces;
                        mesh->vmapP[id1] = parcnt;
                        --parcnt;
                    }
                }
            }
        }
    }

    /* buffers for communication */
    mesh->f_outQ = (float*) calloc(mesh->parNtotalout+1, sizeof(float));
    mesh->f_inQ  = (float*) calloc(mesh->parNtotalout+1, sizeof(float));


    /* build maps between element send buffer */
    int parEtotalout = mesh->parNtotalout/p_Nfp; /* total num of parallel faces*/
    mesh->elemapOUT = BuildIntVector(parEtotalout);

    int k;
    sk = 0;
    /* build map from f_inE to element */
    for(p2=0;p2<nprocs;p2++){
        for(m=0;m<mesh->Npar[p2];m++){
            for(k=0;k<K;k++){
                for(f1=0;f1<Nfaces;f1++){
                    if(mesh->EToP[k][f1]==p2 && p2!=procid)
                        mesh->elemapOUT[sk++] = k;
                }
            }
        }
    }

    parcnt=-1;
    for(p2=0;p2<nprocs;p2++){
        for(k=0;k<K;k++){
            for(f1=0;f1<Nfaces;f1++){
                if(mesh->EToP[k][f1]==p2 && p2!=procid) {
                    mesh->EToE[k][f1] = parcnt;
                    --parcnt;
                }
            }
        }
    }

    /* deallocate mem */
    DestroyVector(nxk);
    DestroyVector(nyk);
    DestroyVector(sJk);
    DestroyIntVector(skP);

    free(xsendrequests); free(ysendrequests);
    free(xrecvrequests); free(yrecvrequests);
    free(Esendrequests); free(Fsendrequests);
    free(Erecvrequests); free(Frecvrequests);

    free(status);

//    char funname[24] = "elepar";
//    FILE *fig = CreateLog(funname, nprocs, procid);
//    PrintBuilMapsLog(fig, mesh);
//    fclose(fig);

#if defined DEBUG
    printf("Procs %d: Leaving BuildMaps\n", procid);
#endif
}

