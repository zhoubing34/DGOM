#include "Convection2d/Convection2d.h"

/**
 * @brief
 * Check MPI send/recv data
 *
 * @author
 * li12242, Tianjin University, li12242@tju.edu.cn
 *
 */
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

/**
 * @brief
 * Check calculation with writing to log file
 *
 * @author
 * li12242, Tianjin University, li12242@tju.edu.cn
 *
 */
void PrintCalculationLog(FILE *fig, int k, float minQ, float maxQ, float Cm, float alpha, float *fQ){
    int n;
    fprintf(fig, "Ele %d: minQ = %f, maxQ = %f, meanC = %f, alpha = %e \n", k, minQ, maxQ, Cm, alpha);
    fprintf(fig, "fQ = ");
    for(n=0;n<p_Np;n++){
        fprintf(fig, "%f, ", fQ[n]);
    }
    fprintf(fig, "\n");
}

/**
 * @brief
 * Slope limiter from Barth and Jesperson
 *
 * @details
 * The BJ slope limiter \f$ \alpha_e \f$ is used to limit the slope of cell in the form
 * \f[
 * u_h(x) = u_c + \alpha_e \left( \nabla u \right)_e \cdot \left(x - x_c \right), \quad
 * 0\le \alpha_e \le 1, \quad x\in \Omega_e
 * \f]
 * In order to simplify the reconstruction procedure, the calculation of element's gradient
 * is avoided with another form giving as
 * \f[ u_i^{limit} = u_c + \alpha_e \left(u_i - u_c \right), \quad
 * 0\le \alpha_e \le 1, \quad i=1,2,\cdots N_p \f]
 *
 * where \f$ N_p \f$ is the total number of nodes in single element.
 *
 * The limiter \f$ \alpha_e \f$ in each element is calculated through cycling all the nodes
 * in the element \f$ \Omega_e \f$ by
 *
 * \f[
 * \alpha_e = min_i \left\{ \begin{array}{ll}
 * min \left\{1, \frac{u_c^{max} - u_c}{u_i - u_c} \right\} & \text{if} \quad u_i - u_c >0 \cr
 * 1, & \text{if} \quad u_i - u_c =0 \cr
 * min \left\{1, \frac{u_c^{min} - u_c}{u_i - u_c} \right\} & \text{if} \quad u_i - u_c <0 \cr
 * \end{array} \right.
 * \f]
 *
 * where \f$ u_c^{max} \f$ and \f$ u_c^{min} \f$ is the maximum and minimum averaged cell value
 * surrounding \f$ \Omega_e \f$ (including itself):
 *
 * \f[ u_c^{max} = max\{ u_e, u_p \}, \quad  p \in \{ \text{neighbour} \} \f]
 *
 * The limited distribution of \f$ u_h \f$ fulfill the maximum principle
 * \f$ u_c^{min} \le u_i^{limit} \le u_c^{max}, \quad i=1,2,\cdots N_p \f$
 *
 * @author
 * li12242, Tianjin University, li12242@tju.edu.cn
 *
 *
 * @return
 */
void LimiterBJ2d(Mesh *mesh){

    int k1, k2, n, f, p, Nout;

    const int K = mesh->K;
    const int procid = mesh->procid;
    const int nprocs = mesh->nprocs;

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
    MPI_Request *mpi_out_requests = (MPI_Request*) calloc(nprocs, sizeof(MPI_Request));
    MPI_Request *mpi_in_requests  = (MPI_Request*) calloc(nprocs, sizeof(MPI_Request));

    /* buffer outgoing node data */
    for(n=0;n<mesh->parFtotalout;++n)
        mesh->f_outE[n] = Cmean[mesh->elemapOUT[n]];

    /* do sends */
    int sk = 0, Nmess = 0, uk;
    for(p=0;p<nprocs;++p){
        if(p!=procid){
            Nout = mesh->Npar[p]; // # of variables send to process p
            if(Nout){
                /* symmetric communications (different ordering) */
                MPI_Isend(f_outE+sk, Nout, MPI_FLOAT, p, 6666+p,            MPI_COMM_WORLD, mpi_out_requests +Nmess);
                MPI_Irecv(f_inE+sk,  Nout, MPI_FLOAT, p, 6666+procid, MPI_COMM_WORLD,  mpi_in_requests +Nmess);
                sk+=Nout;
                ++Nmess;
            }
        }
    }

    /* DO RECV */
    MPI_Status *instatus  = (MPI_Status*) calloc(nprocs, sizeof(MPI_Status));
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
        fprintf(fig, "\n");
#endif
    }

    /* make sure all messages went out */
    MPI_Status *outstatus  = (MPI_Status*) calloc(nprocs, sizeof(MPI_Status));
    MPI_Waitall(Nmess, mpi_out_requests, outstatus);
    free(outstatus);

#if defined DEBUG
    fclose(fig);
#endif

    /* dealloc mem */
    free(Cmean);
    free(mpi_out_requests);
    free(mpi_in_requests);
}

