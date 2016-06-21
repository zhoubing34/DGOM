#include "ConvectionDriver2d.h"

/**
 * @brief
 * Post process of the 2d convection problem
 *
 * @details
 * Calculation of the \f$L_2\f$ and \f$L_{\infty}\f$ error. The formula is giving as
 * \f[ L_2 = \sqrt{\frac{1}{N} \sum_i^N \left(u_i - u(x_i) \right)^2}, \quad
 * L_{\infty} = max\{ \left| u_i - u(x_i) \right| \}, \quad i=1,2,\cdots,N_p \f]
 *
 * @author
 * li12242, Tianjin University, li12242@tju.edu.cn
 *
 * @warning
 * @attention
 * @note
 * @todo
 */
void Postprocess(PhysDomain2d *phys){

    double L2   = 0.0, Linf  = 0.0;
    double gL2  = 0.0, gLinf = 0.0;
    float *Qe, locerr, dis;
    int k, n;

    MultiReg2d *mesh = phys->mesh;
    StdRegions2d *shape = mesh->stdcell;

    const int K      = mesh->K;
    const int procid = mesh->procid;
//    const int nprocs = mesh->nprocs;
    float *f_Q = phys->f_Q;

    Qe = (float*) calloc(mesh->K*shape->Np, sizeof(float));
    /* get the exact solution */
    void GenExactSolution(MultiReg2d *, float *);
    GenExactSolution(mesh, Qe);

    /* calculate the L2 and Linf error */
    for(k=0;k<K;k++){
        locerr = 0.0f;
        float *qpt = f_Q + phys->Nfields*shape->Np*k;
        float *qet = Qe + phys->Nfields*shape->Np*k;
        for(n=0;n<shape->Np;n++){
            dis = qpt[n] - qet[n];
            locerr = (dis*dis);
            Linf = max(Linf, fabsf(dis));
        }
        L2 += locerr;
    }
    L2 = sqrt( L2/(K*shape->Np) );

    /* collect all process */
    MPI_Allreduce(&Linf, &gLinf, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&L2, &gL2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    /* print output */
    if(!procid) printf("L2:%le, Linf:%le\n", gL2, gLinf);

    /* finish */
    free(Qe);
}


void GenExactSolution(MultiReg2d *mesh, float *Qe){
    int k, n, sk=0;
    int K = mesh->K;
    double t, sigma = 125*1e3/(33*33);
    /* initial position */
    double xc = 0.0, yc = 0.6;

    StdRegions2d *shape = mesh->stdcell;

    /* initial scalar field */
    for(k=0;k<K;++k){
        for(n=0;n<shape->Np;++n){
            t = -sigma * ( ( mesh->x[k][n] - xc )*( mesh->x[k][n] - xc )
                           + ( mesh->y[k][n] - yc )*( mesh->y[k][n] - yc ) );
            Qe[sk++] = (float) exp(t);
        }
    }
}
