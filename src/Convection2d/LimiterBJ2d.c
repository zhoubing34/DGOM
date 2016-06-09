#include "Convection2d/Convection2d.h"

void PrintLimiterLog(FILE *fig, Mesh *mesh, float *Cm){
    int K = mesh->K;
    int k,n;
    fprintf(fig, "Mean Value:\n");
    for(k=0;k<K;k++){
        fprintf(fig, "%f, ", Cm[k]);
    }
    fprintf(fig, "\n\n");

    fprintf(fig, "\nElement info to send:\n");
    for(n=0;n<mesh->parFtotalout;n++)
        fprintf(fig, "%f, ", mesh->f_outE[n]);
    fprintf(fig, "\n\n");

    fprintf(fig, "\nElement info recv:\n");
    for(n=0;n<mesh->parFtotalout;n++)
        fprintf(fig, "%f, ", mesh->f_inE[n]);
    fprintf(fig, "\n\n");
}

void PrintCalculationLog(FILE *fig, int k, float minQ, float maxQ, float Cm, float alpha, float *fQ){
    int n;
    fprintf(fig, "Ele %d: minQ = %f, maxQ = %f, meanC = %f, alpha = %e \n", k, minQ, maxQ, Cm, alpha);
    fprintf(fig, "fQ = ");
    for(n=0;n<p_Np;n++){
        fprintf(fig, "%f, ", fQ[n]);
    }
    fprintf(fig, "\n");
}


void LimiterBJ2d(Mesh *mesh){

    int k1, k2, n, f, p, Nout;

    const int K = mesh->K;
    const int procid = mesh->procid;

    float *f_Q    = mesh->f_Q;
    float *f_inE  = mesh->f_inE;
    float *f_outE = mesh->f_outE;

    float *Cmean = (float*) calloc(K, sizeof(float));

    /* pointer to local element storage */
    float *qpt;

#if defined DEBUG
    char *funname = "limiter2d";
    FILE* fig = CreateLog(funname, mesh->nprocs, procid);
#endif

    int sj =0;
    /* mean value in element */
    for(k1=0;k1<K;k1++){
        qpt = f_Q + p_Nfields*p_Np*k1;
        for(n=0;n<p_Np;n++){
            /* volume integral */
            Cmean[k1] += mesh->wv[n]*qpt[n]*mesh->J[sj++];
        }
        Cmean[k1] /= mesh->area[k1];
    }

    /* send and recv ele info */
    /* mpi request buffer */
    MPI_Request *mpi_out_requests = (MPI_Request*) calloc(mesh->nprocs, sizeof(MPI_Request));
    MPI_Request *mpi_in_requests  = (MPI_Request*) calloc(mesh->nprocs, sizeof(MPI_Request));

    /* buffer outgoing node data */
    for(n=0;n<mesh->parFtotalout;++n)
        mesh->f_outE[n] = Cmean[mesh->elemapOUT[n]];

    /* do sends */
    int sk = 0, Nmess = 0, uk;
    for(p=0;p<mesh->nprocs;++p){
        if(p!=mesh->procid){
            Nout = mesh->Npar[p]; // # of variables send to process p
            if(Nout){
                /* symmetric communications (different ordering) */
                MPI_Isend(f_outE+sk, Nout, MPI_FLOAT, p, 6666+p,            MPI_COMM_WORLD, mpi_out_requests +Nmess);
                MPI_Irecv(f_inE+sk,  Nout, MPI_FLOAT, p, 6666+mesh->procid, MPI_COMM_WORLD,  mpi_in_requests +Nmess);
                sk+=Nout;
                ++Nmess;
            }
        }
    }

    /* DO RECV */
    MPI_Status *instatus  = (MPI_Status*) calloc(mesh->nprocs, sizeof(MPI_Status));
    MPI_Waitall(Nmess, mpi_in_requests, instatus);
    free(instatus);

#if defined DEBUG
    PrintLimiterLog(fig, mesh, Cmean);
#endif

    float maxQ, minQ, alpha;
    /* Limit */
    for(k1=0;k1<K;k1++){

        /* check if k1 is trouble ele */
        if(mesh->tcflag[k1]>DETECTOR) {

            const float Cm = Cmean[k1];
            maxQ = Cm, minQ = Cm;
            for (f = 0; f < p_Nfaces; f++) {
                k2 = mesh->EToE[k1][f];
                if (k2 < 0) {
                    k2 = -k2 - 1;
                    maxQ = max(maxQ, f_inE[k2]);
                    minQ = min(minQ, f_inE[k2]);
                } else {
                    maxQ = max(maxQ, Cmean[k2]);
                    minQ = min(minQ, Cmean[k2]);
                }
            }

            qpt = f_Q + p_Nfields * p_Np * k1;

            /* init slope limiter */
            alpha = 1.0f;

#if defined DEBUG
            PrintCalculationLog(fig, k1, minQ, maxQ, Cm ,alpha, qpt);
#endif
            /* calculate slope limit */
            for (n = 0; n < p_Np; n++) {
                if (qpt[n] > Cm) {
                    alpha = min(alpha, fabsf((maxQ - Cm) / (qpt[n] - Cm)));
                }
                if (qpt[n] < Cm) {
                    alpha = min(alpha, fabsf((minQ - Cm) / (qpt[n] - Cm)));
                }
            }
            /* Do limit */
            for (n = 0; n < p_Np; n++) {
                qpt[n] = alpha * (qpt[n] - Cm) + Cm;
            }


        }
#if defined DEBUG
        PrintCalculationLog(fig, k1, minQ, maxQ, Cm, alpha, qpt);
#endif
    }

    /* make sure all messages went out */
    MPI_Status *outstatus  = (MPI_Status*) calloc(mesh->nprocs, sizeof(MPI_Status));
    MPI_Waitall(Nmess, mpi_out_requests, outstatus);
    free(outstatus);

#if defined DEBUG
    PrintLimiterLog(fig, mesh, Cmean);
    fclose(fig);
#endif

    /* dealloc mem */
    free(Cmean);
    free(mpi_out_requests);
    free(mpi_in_requests);
}

