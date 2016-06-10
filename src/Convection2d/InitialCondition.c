#include "Convection2d/Convection2d.h"

/**
 * @brief
 * Set initial condition
 *
 * @details
 * Set initial scalar field distribution and the constant flow field
 * 1. scalar distribution
 * 2. constant velocity field
 *
 * @author
 * li12242, Tianjin University, li12242@tju.edu.cn
 *
 */
void InitData(Mesh * mesh){
    double sigma = 125*1e3/(33*33);
    double xc, yc, t;
    double w;

    const int K = mesh->K;
    int k, n, sk=0;

#if defined DEBUG
    int procid = mesh->procid;
    if(!procid) printf("Root: Entering InitData\n");
#endif

    /* initial position */
    xc = 0.0; yc = 0.6;

    /* initial scalar field */
    for(k=0;k<K;++k){
        for(n=0;n<p_Np;++n){

            t = -sigma * ( ( mesh->x[k][n] - xc )*( mesh->x[k][n] - xc )
                           + ( mesh->y[k][n] - yc )*( mesh->y[k][n] - yc ) );
            mesh->f_Q[sk++] = (float) exp(t);
        }
    }

    /* flow rate field */
    w = 5*M_PI/6;
    sk=0;
    for (k=0; k<K; ++k){
        for (n=0; n<p_Np; ++n){
            mesh->f_s[sk++] = (float)(-w * mesh->y[k][n]); // flow rate at x-coordinate
            mesh->f_s[sk++] = (float)(w * mesh->x[k][n]);  // flow rate at y-coordinate
        }
    }

#if defined DEBUG
    if(!procid) printf("Root: Leaving InitData\n");
#endif
}
