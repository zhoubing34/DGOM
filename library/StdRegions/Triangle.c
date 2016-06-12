#include <stdlib.h>

/**
 * @file
 * Triangle.c
 *
 * @brief
 * Functions related with standard triangle elements
 *
 * @author
 * li12242, Tianjin University, li12242@tju.edu.cn
 */


#include "StdRegions.h"

/* local functions */
void xytors(int Np, double *x, double *y, double *r, double *s);

/**
 * @brief
 * Generation of standard triangle element
 *
 * @details
 *
 * @param[in] N polynomial order
 *
 * @return
 * return values:
 * name     | type     | description of value
 * -------- |----------|----------------------
 * tri | StdRegions2d* |
 *
 */
StdRegions2d* GenStdTriEle(int N){
    StdRegions2d *tri = (StdRegions2d *) calloc(1, sizeof(StdRegions2d));

    /* basic info */
    tri->N = N;
    tri->Np = ((N+1)*(N+2)/2);
    tri->Nfaces = 3;
    tri->Nfp = N+1;
    tri->Nv = 3;

    /* coordinate */
    tri->r = BuildVector(tri->Np);
    tri->s = BuildVector(tri->Np);

    return tri;
}


void FreeStdTriEle(StdRegions2d * triangle){
    /* coordinate */
    free(triangle->r);
    free(triangle->s);
}

/**
 * @brief
 * Get the nature coordinate of triangle element
 *
 * @details
 *
 * @param[in] N polynomial order
 *
 * @return
 * return values:
 * name     | type     | description of value
 * -------- |----------|----------------------
 * r | double[Np] | coordinate r
 * s | double[Np] | coordinate s
 * where Np = (N+1)(N+2)/2
 *
 */
void GetTriCoord(int N, double *r, double *s){
    double alpopt[15] =  {0.0000, 0.0000, 1.4152, 0.1001, 0.2751, 0.9800, 1.0999, 1.2832,
                          1.3648, 1.4773, 1.4959, 1.5743, 1.5770, 1.6223, 1.6258};
    int i,j,sk=0,Np=(N+1)*(N+2)/2;
    double *L1, *L2, *L3, *dL;
    double *x, *y;
    double *warpf1;
    double alpha, temp;

    if(N<16){
        alpha = alpopt[N-1];
    }else{
        alpha = 5.0/3.0;
    }

    L1 = BuildVector(Np);
    L2 = BuildVector(Np);
    L3 = BuildVector(Np);
    dL = BuildVector(Np);
    x  = BuildVector(Np);
    y  = BuildVector(Np);

    warpf1 = BuildVector(Np);

    printf("Np = %d\n", Np);
    for(i=0;i<N+1;i++){
        for(j=0;j<(N-i)+1;j++){
            L1[sk] = (double)i/N;
            L3[sk] = (double)j/N;
            L2[sk] = 1.0 - L1[sk] - L3[sk];

            x[sk] = -L2[sk] + L3[sk];
            y[sk] = (-L2[sk] - L3[sk] + 2.0*L1[sk])/sqrt(3.0);
            sk++;
        }
    }

    /* warp1 */
    for(i=0;i<Np;i++){
        dL[i] = L3[i] - L2[i];
    }
    Warpfactor(N, dL, Np, warpf1);

    for(i=0;i<Np;i++){
        temp = alpha*L1[i];
        warpf1[i] *= 4.0*L2[i]*L3[i]*(1.0 + temp*temp);
        x[i] += 1.0*warpf1[i];
//        y[sk] += 0.0;
    }

    /* warp2 */
    for(i=0;i<Np;i++){
        dL[i] = L1[i] - L3[i];
    }
    Warpfactor(N, dL, Np, warpf1);

    for(i=0;i<Np;i++){
        temp = alpha*L2[i];
        warpf1[i] *= 4.0*L1[i]*L3[i]*(1.0 + temp*temp);
        x[i] += cos(2.0*M_PI/3.0)*warpf1[i];
        y[i] += sin(2.0*M_PI/3.0)*warpf1[i];
    }

    /* warp3 */
    for(i=0;i<Np;i++){
        dL[i] = L2[i] - L1[i];
    }
    Warpfactor(N, dL, Np, warpf1);

    for(i=0;i<Np;i++){
        temp = alpha*L3[i];
        warpf1[i] *= 4.0*L1[i]*L2[i]*(1.0 + temp*temp);
        x[i] += cos(4.0*M_PI/3.0)*warpf1[i];
        y[i] += sin(4.0*M_PI/3.0)*warpf1[i];
    }

    /* coordinate transfer */
    xytors(Np, x, y, r, s);

    /* deallocate mem */
    DestroyVector(L1); DestroyVector(L2); DestroyVector(L3);
    DestroyVector(x);
    DestroyVector(y);
    DestroyVector(warpf1);
}

/**
 * @brief
 * Warp factor to connnect the Legendre-Gauss-Lobatto and equidistant nodes
 *
 * @details
 * The warp factor w is used to connect the equidistant nodes \f$\{r_i^e\}\f$ and the
 * Legendre-Gauss-Lobatto nodes \f$\{r_i^{LGL}\}\f$ and calculated by following formula
 *
 * \f[ w(r) = \sum_{i=1}^{Np}\left( r_i^{LGL} - r_i^e \right)l_i^e(r) \f]
 *
 * where \f$l_i^e(r)\f$ are the Lagrange polynomials based on \f$r_i^e\f$.
 *
 * @param[in] N order of degree
 * @param[in] r input points
 * @param[in] Nr number of input points
 *
 * @return
 * return values:
 * name     | type     | description of value
 * -------- |----------|----------------------
 * w | double[Nr]  | waper factor
 *
 */
void Warpfactor(int N, double *r, int Nr, double *w){
    int i, j, Np = N+1;
    double *ye, *l, *re, *rlgl, *wlgl;
    double temp;

    ye   = BuildVector(Np);
    l    = BuildVector(Nr);
    re   = BuildVector(Np); /* equidistant nodes */
    rlgl = BuildVector(Np); /* Gauss-Lobatto-Jacobi nodes */
    wlgl = BuildVector(Np); /* Gauss-Lobatto-Jacobi weights */

    /* initialization */
    for(i=0;i<Nr;i++){
        w[i] = 0.0;
    }

    for(i=0;i<Np;i++){
        /* equidistant between [-1,1] */
        re[i] = (double)i/N*2.0 - 1.0;
    }

    /* get Gauss-Lobatto-Jacobi zeros and weights */
    zwglj(rlgl, wlgl, Np, 0, 0);

    for(i=0;i<Np;i++){
        /* get lagrange basis l at r */
        ye[i] = 1.0;
        laginterp(Np, re, ye, Nr, r, l);

        for(j=0;j<Nr;j++){
            w[j] += l[j]*(rlgl[i] - re[i]);
        }
        ye[i] = 0.0;
    }

    for(i=0;i<Nr;i++){
        temp = (1.0 - r[i]*r[i]);
        if(temp > 1.0e-10){
            w[i] /= temp;
        }
    }

    /* deallocate mem */
    DestroyVector(ye);
    DestroyVector(l);
    DestroyVector(re);
    DestroyVector(rlgl);
    DestroyVector(wlgl);
}

/**
 * @brief
 *
 *
 * @details
 * transfer coordinate (x,y) on equilateral triangle to natural coordinate (r,s)
 * on right triangle
 *
 * @param[in] Np number of points
 * @param[in] x coordinate
 * @param[in] y coordinate
 *
 * @return
 * return values:
 * name     | type     | description of value
 * -------- |----------|----------------------
 * r | double[Np] |
 * s | double[Np] |
 *
 */
void xytors(int Np, double *x, double *y, double *r, double *s){
    double L1, L2, L3;
    int i;

    for(i=0;i<Np;i++){
        L1 = (sqrt(3.0)*y[i] + 1)/3.0;
        L2 = (-3.0*x[i] - sqrt(3.0)*y[i] + 2.0)/6.0;
        L3 = ( 3.0*x[i] - sqrt(3.0)*y[i] + 2.0)/6.0;

        r[i] = -L2+L3-L1;
        s[i] = -L2-L3+L1;
    }
}
