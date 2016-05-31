#include "Convection2d/Convection2d.h"

void InitData(Mesh * mesh){
    double sigma = 125*1e3/(33*33);
    double xc, yc, t;
    double w;

    const int K = mesh->K;
    int k, n, sk=0;

    // initial position
    xc = 0.0; yc = 0.6;

    // allocate memory
    mesh->f_Q    = (float*) calloc(mesh->K*BSIZE*p_Nfields, sizeof(float));
    mesh->f_rhsQ = (float*) calloc(mesh->K*BSIZE*p_Nfields, sizeof(float));
    mesh->f_resQ = (float*) calloc(mesh->K*BSIZE*p_Nfields, sizeof(float));

    mesh->f_s = (float*) calloc(mesh->K*BSIZE*2, sizeof(float));

    // initial scalar field
    for(k=0;k<K;++k){
        for(n=0;n<p_Np;++n){
//            mesh->f_Q[sk++] = mesh->x[k][n];
            t = -sigma * ( ( mesh->x[k][n] - xc )*( mesh->x[k][n] - xc )
                           + ( mesh->y[k][n] - yc )*( mesh->y[k][n] - yc ) );
            mesh->f_Q[sk++] = (float) exp(t);
        }
    }


    w = 5*M_PI/6;
    sk=0;
    // flow rate field
    for (k=0; k<K; ++k){
        for (n=0; n<p_Np; ++n){
//            mesh->f_s[sk++] = 1;
//            mesh->f_s[sk++] = 0;
            mesh->f_s[sk++] = (float)(-w * mesh->y[k][n]); // flow rate at x-coordinate
            mesh->f_s[sk++] = (float)(w * mesh->x[k][n]);  // flow rate at y-coordinate
        }
    }
}
