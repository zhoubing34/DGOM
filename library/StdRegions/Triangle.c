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

/* private functions */
void xytors(int Np, double *x, double *y, double *r, double *s);
void rstoad(int Np, double *r, double *s, double *a, double *b);

/**
 * @brief
 * Free variables fields in StdRegions2d
 */
void FreeStdRegions2d(StdRegions2d *triangle){
    /* coordinate */
    DestroyVector(triangle->r);
    DestroyVector(triangle->s);
    /* vandermonde matrix */
    DestroyMatrix(triangle->V);
    /* mass matrix */
    DestroyMatrix(triangle->M);
}

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
StdRegions2d* GenStdTriEle(const unsigned N){
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

    GetTriCoord(N, tri->r, tri->s);

    /* vandermonde matrix */
    tri->V = BuildMatrix(tri->Np, tri->Np);
    GetTriV(N, tri->Np, tri->r, tri->s, tri->V);

    /* mass matrix */
    tri->M = BuildMatrix(tri->Np, tri->Np);
    GetTriM(tri->Np, tri->V, tri->M);

    /* Derivative Matrix */

    return tri;
}

/**
 * @brief
 * Generate the mass matrix
 * @details
 * The mass matrix is calculated with
 * \f[ \mathbf{M} = (\mathbf{V}^T)^{-1} \cdot \mathbf{V}^{-1} \f]
 *
 * @param[in] Np
 * @param[in] V Vandermonde matrix
 *
 * @return
 * return values:
 * name     | type     | description of value
 * -------- |----------|----------------------
 * M   | double[Np][Np] | Mass Matrix
 */
void GetTriM(unsigned Np, double **V, double **M){
    doublereal *temp = BuildVector(Np*Np);
    doublereal *invt = BuildVector(Np*Np);
    double *Mv = BuildVector(Np*Np);
    const unsigned n = Np;

    int i,j,sk=0;

    for(i=0;i<Np;i++){
        for(j=0;j<Np;j++){
            /* row counts first */
            temp[sk++] = V[i][j];
        }
    }

    invM(temp, Np);

    for(i=0;i<Np;i++){
        for(j=0;j<Np;j++){
            invt[j*Np + i] = temp[i*Np + j];
        }
    }

    dgemm_(n, n, n, n, invt, temp, Mv);

    for(i=0;i<Np;i++){
        for(j=0;j<Np;j++){
            M[i][j] = Mv[i*Np + j];
        }
    }
}

/**
 * @brief
 * Evaluate 2D orthonormal polynomial on simplex at (a,b) of order (i,j).
 *
 * @details
 *
 * @param[in] Np number of points
 * @param[in] a coordinate
 * @param[in] b coordinate
 * @param[in] i order
 * @param[in] j order
 *
 * @return
 * return values:
 * name     | type     | description of value
 * -------- |----------|----------------------
 * poly   | double[Np] |
 *
 */
void Simplex2dP(int Np, double *a, double *b, int i, int j, double *poly){
    double *h1 = BuildVector(Np);
    double *h2 = BuildVector(Np);
    int n;

    jacobiP(Np, a, h1, i, 0.0, 0.0);
    jacobiP(Np, b, h2, j, 2*i+1, 0.0);

    for(n=0;n<Np;n++){
        poly[n] = sqrt(2.0)*h1[n]*h2[n]*pow(1-b[n], i);
    }

    DestroyVector(h1);
    DestroyVector(h2);
}

/**
 * @brief
 * Generate the Vandermonde matrix of triangle
 *
 * @details
 *
 * @param[in] r
 * @param[in] s
 * @param[in] Nr number of r
 * @param[in] N order
 *
 * @return
 * return values:
 * name     | type     | description of value
 * -------- |----------|----------------------
 * V  | double[Np][Np] |
 *
 */
void GetTriV(int N, int Nr, double *r, double *s, double **V){
    int i,j,k,sk=0;

    double *temp = BuildVector(Nr);
    double *a=BuildVector(Nr);
    double *b=BuildVector(Nr);

    rstoad(Nr, r, s, a, b);

    for(i=0;i<N+1;i++){
        for(j=0;j<N-i+1;j++){
            Simplex2dP(Nr, a, b, i, j, temp);
            for(k=0;k<Nr;k++){
                V[k][sk] = temp[k];
            }
            sk++;
        }
    }

    DestroyVector(a);
    DestroyVector(b);
    DestroyVector(temp);
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
 * @note
 * r and s shuold be allocated before calling GetTriCoord
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
 * @note
 * w shuold be allocated before calling Warpfactor
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
 * @note
 * x,y,r and s shuold be allocated before calling xytors
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

/**
 * @brief
 *
 * @details 详细说明
 *
 *
 * @param[in] r
 * @param[in] s
 * @param[in] Np
 *
 * @return
 * return values:
 * name     | type     | description of value
 * -------- |----------|----------------------
 * a | int      |
 * b | object   |
 *
 * @note
 * r,s,a and b shuold be allocated before calling rstoad
 */

void rstoad(int Np, double *r, double *s, double *a, double *b){
    int i;
    for(i=0;i<Np;i++){
        if( fabs(s[i] - 1.0) > 1.0e-10){
            a[i] = 2.0*(1.0+r[i])/(1.0-s[i])-1;
        }else{
            a[i] = -1.0;
        }
        b[i] = s[i];
    }
}