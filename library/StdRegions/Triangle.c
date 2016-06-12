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
    x = BuildVector(Np); y = BuildVector(Np);

    warpf1 = BuildVector(Np);

    printf("Np = %d\n", Np);
    for(i=0;i<N+1;i++){
        for(j=0;j<(N-i)+1;j++){
            L1[sk] = (double)i/N;
            L3[sk] = (double)j/N;
            L2[sk] = 1.0 - L1[sk] - L3[sk];

            x[sk] = -L2[sk] + L3[sk];
            y[sk] = (-L2[sk] - L3[sk] + 2.0*L1[sk])/sqrt(3.0);
//            printf("L1[%d]=%lg, \tL2[%d]=%lg, \tL3[%d]=%lg, "
//                           "\tx[%d]=%lg, \ty[%d]=%lg\n",sk, L1[sk],
//            sk, L2[sk], sk, L3[sk], sk, x[sk], sk, y[sk]);
            sk++;
        }
    }
    printf("sk = %d\n", sk);

    /* warp1 */
    for(i=0;i<Np;i++){
        dL[i] = L3[i] - L2[i];
    }
    Warpfactor(N, dL, Np, warpf1);

    printf("\nwarpf1\n");
    for(i=0;i<Np;i++){
        printf("L3[%d]-L2[%d]=%12.4f, warpf[%d]=%12.4f\n",i, i, dL[i], i, warpf1[i]);
    }

    for(i=0;i<Np;i++){
        temp = alpha*L1[sk];
        warpf1[sk] *= 4.0*L2[i]*L3[i]*(1.0 + temp*temp);
        x[sk] += 1.0*warpf1[sk];
//        y[sk] += 0.0;
    }

    /* warp2 */
    for(i=0;i<Np;i++){
        dL[i] = L1[i] - L3[i];
    }
    Warpfactor(N, dL, Np, warpf1);

    printf("\nwarpf2\n");
    for(i=0;i<Np;i++){
        printf("L1[%d]-L3[%d]=%12.4f, warpf[%d]=%12.4f\n",i, i, dL[i], i, warpf1[i]);
    }

    for(i=0;i<Np;i++){
        temp = alpha*L2[sk];
        warpf1[sk] *= 4.0*L1[i]*L3[i]*(1.0 + temp*temp);
        x[sk] += cos(2.0*M_PI/3.0)*warpf1[sk];
        y[sk] += sin(2.0*M_PI/3.0)*warpf1[sk];
    }

    /* warp3 */
    for(i=0;i<Np;i++){
        dL[i] = L2[i] - L1[i];
    }
    Warpfactor(N, dL, Np, warpf1);

    printf("\nwarpf3\n");
    for(i=0;i<Np;i++){
        printf("L2[%d]-L1[%d]=%12.4f, warpf[%d]=%12.4f\n",i, i, dL[i], i, warpf1[i]);
    }

    for(i=0;i<Np;i++){
        temp = alpha*L3[sk];
        warpf1[sk] *= 4.0*L1[i]*L2[i]*(1.0 + temp*temp);
        x[sk] += cos(4.0*M_PI/3.0)*warpf1[sk];
        y[sk] += sin(4.0*M_PI/3.0)*warpf1[sk];
    }

    for(i=0;i<Np;i++){
        printf("x[%d]=%12.4f, \ty[%d] = %12.4f\n", i, x[i], i, y[i]);
    }

    printf("\nleaving GetTriCoord\n");
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
        printf("w[%d] = %12.4f\n", i, w[i]);
        temp = (1.0 - r[i]*r[i]);
        if(temp > 1.0e-10){
            w[i] /= temp;
        }
    }

    /**/
    DestroyVector(ye);
    DestroyVector(l);
    DestroyVector(re);
    DestroyVector(rlgl);
    DestroyVector(wlgl);
}
