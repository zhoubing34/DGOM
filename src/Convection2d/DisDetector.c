#include "Convection2d/Convection2d.h"

/**
 * @brief
 * Check calculation with writing to log file
 *
 * @author
 * li12242, Tianjin University, li12242@tju.edu.cn
 *
 */
void PrintDisDetectorLog(FILE *fig, int k,
                         double len, double disQ,
                         double maxQ, float I, float flag){
    fprintf(fig, "ele: %d, len: %f, disQ: %f, maxQ: %e, detector: %f, flag: %f\n",
            k, len, disQ, maxQ, I, flag);
}

/**
 * @brief
 * Discontinuity detector of Krivodonova (2003)
 *
 * @details
 * The discontinuity detector \f$I_j\f$ of jth cell is obtained from
 * Krivodonova et. al (2003) and giving as
 * \f[ I_j = \frac{ \left| \int_{\partial \Omega_j^-} \left( Q_j - Q_{nb_j}
 * \right) ds \right| }{ r \left| \partial \Omega_j \right| \left\|Q_j \right\| } \f]
 *
 * where \f$ r \f$ is the radius of the circumscribed circle of jth cell and
 * \f$ \partial \Omega_j \f$ is the total boundary length.
 *
 * @attention
 * The radius \f$r \f$ is not powered by \f$(p+1)/2\f$ and the
 * \f$ \left| \partial \Omega_j \right| \f$ in the denominator is the total length
 * of all boundaries, which is different with original formular of Krivodonova et. al (2003).
 *
 * @author
 * li12242, Tianjin University, li12242@tju.edu.cn
 *
 */
void DisDetector(Mesh *mesh){

#if defined(DEBUG)
    char name[24] = "DisDetector";
    /* create log file */
    FILE *fig = CreateLog(name, mesh->nprocs, mesh->procid);
#endif

    /* discontinuity indicator */
    float I;

    /* mesh parameters */
    const int K = mesh->K;

    double inlen, maxQ;
    /* radius of the circumscribed circle */
    const float  *f_Q    = mesh->f_Q;
    const float  *f_s    = mesh->f_s;
    const float  *surfinfo = mesh->surfinfo;
    float  dC;

    float *f_inQ  = mesh->f_inQ;
    float *f_outQ  = mesh->f_outQ;

    register unsigned int k, f, m, n;

    /* mpi request buffer */
    MPI_Request *mpi_out_requests = (MPI_Request*) calloc(mesh->nprocs, sizeof(MPI_Request));
    MPI_Request *mpi_in_requests  = (MPI_Request*) calloc(mesh->nprocs, sizeof(MPI_Request));

    /* buffer outgoing node data */
    for(n=0;n<mesh->parNtotalout;++n)
        mesh->f_outQ[n] = f_Q[mesh->parmapOUT[n]];

    /* do sends */
    int p, Nout;
    int sk = 0, Nmess = 0;
    for(p=0;p<mesh->nprocs;++p){
        if(p!=mesh->procid){
            Nout = mesh->Npar[p]*p_Nfields*p_Nfp; // # of variables send to process p
            if(Nout){
                /* symmetric communications (different ordering) */
                MPI_Isend(f_outQ+sk, Nout, MPI_FLOAT, p, 6666+p,            MPI_COMM_WORLD, mpi_out_requests +Nmess);
                MPI_Irecv(f_inQ+sk,  Nout, MPI_FLOAT, p, 6666+mesh->procid, MPI_COMM_WORLD,  mpi_in_requests +Nmess);
                sk+=Nout;
                ++Nmess;
            }
        }
    }

    /* DO RECV */
    MPI_Status *instatus  = (MPI_Status*) calloc(mesh->nprocs, sizeof(MPI_Status));
    MPI_Waitall(Nmess, mpi_in_requests, instatus);
    free(instatus);

    /* difference & edge length on infow boundary */
    for(k=0; k<K; k++){
        mesh->tcflag[k] = 0;

        /* init difference and edge length */
        I = 0; inlen =0;

        /* NOTE: index into geometric factors */
        int surfid=k*6*p_Nfp*p_Nfaces;

        for(f=0; f<p_Nfaces; f++){
            for(m=0; m<p_Nfp; m++){

                int idM = surfinfo[surfid++];
                int idP = surfinfo[surfid++];
                const float FSc = surfinfo[surfid++];
                const float SJc = surfinfo[surfid++];
                const float NXf = surfinfo[surfid++];
                const float NYf = surfinfo[surfid++];

                const float uM = f_s[idM * 2];
                const float vM = f_s[idM * 2 + 1];
                const float un = (NXf * uM + NYf * vM);

//                if (un < 0){
                    if (idP < 0) {
                        idP = p_Nfields * (-1 - idP);
                        dC = fabsf(f_inQ[idP] - f_Q[idM]);
                    }
                    else{
                        dC = fabsf(f_Q[idP] - f_Q[idM]);
                    }
                    /* difference and edge length */
                    I     += mesh->w[m]*dC*SJc;
                    inlen += mesh->w[m]*SJc;
//                }
            }

        }

        /* NOTE: buffer element k into local storage */
        const float *qpt = f_Q + p_Nfields*p_Np*k;

        maxQ = 0;
        for(n=0;n<p_Nfields*p_Np;n++){
            maxQ = max(maxQ, fabsf(qpt[n]) );
        }

#if defined DEBUG
        double disQ = I;
        double mQ = maxQ;
#endif

//        maxQ *= powf(mesh->ciradius[k], (float)(p_N+1)/2 );
        maxQ *= mesh->ciradius[k];
        maxQ *= inlen;

        /* check if denominator greater than 0 */
        if( maxQ > NODETOL ){
            I /= maxQ;
        }else{
            I = 0;
        }

        /* assignment */
        mesh->tcflag[k] = I;

#if defined(DEBUG)
        /* Print log information */
        PrintDisDetectorLog(fig, k, inlen, disQ, mQ, I, mesh->tcflag[k]);
#endif
    }

#if defined(DEBUG)
    fclose(fig);
#endif

    /* make sure all messages went out */
    MPI_Status *outstatus  = (MPI_Status*) calloc(mesh->nprocs, sizeof(MPI_Status));
    MPI_Waitall(Nmess, mpi_out_requests, outstatus);
    free(outstatus);

    free(mpi_out_requests);
    free(mpi_in_requests);

}


