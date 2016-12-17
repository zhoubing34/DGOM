/**
 * @file
 * Standard triangle element
 *
 * @brief
 * Functions related with standard triangle elements
 *
 * @author
 * li12242, Tianjin University, li12242@tju.edu.cn
 */

#include "sc_stdcell.h"

/* transform the index of orthogonal function to [ti,tj] */
void sc_triTransInd(int N, int ind, int *ti, int *tj);

/* transfer coordinate [x,y] on equilateral triangle to natural coordinate [r,s] */
void sc_trixytors(int Np, double *x, double *y, double *r, double *s);

/* transform natural coordinate to collapse coordinate */
void sc_rstoad(int Np, double *r, double *s, double *a, double *b);

/* get the gradient of the modal basis (id,jd) on the 2D simplex at (a,b). */
void sc_triGradSimplex2DP(int Np, double *a, double *b, int id, int jd, double *dmodedr, double *dmodeds);

/* build the nodes index matrix on each faces */
int** sc_triFmask(stdCell *tri);

/* Warp factor to connnect the Legendre-Gauss-Lobatto and equidistant nodes */
void sc_triWarpfactor(int, double *, int, double *);

/* get the coordinate of interpolations points on standard triangle element */
void sc_triCoord(stdCell *);

/* evaluate 2D orthonormal polynomial on simplex at (a,b) of order (i,j). */
void sc_triSimplex2dP(int Np, double *a, double *b, int i, int j, double *poly);

/* derivative of orthogonal function */
void sc_triDeriOrthogFunc(stdCell *tri, int ind, double *dr, double *ds);

/* get orthogonal function value at interpolation nodes */
void sc_triOrthogFunc(stdCell *cell, int ind, double *func);

/* create standard triangle element */
stdCell* sc_createTri(int N);



/**
 * @brief transform the index of orthogonal function to [ti,tj] for triangle elements
 * @details the index is arranged as
 * \f[ [i,j] = \left\{ \begin{array}{lllll}
 * (0,0) & (0,1) & \cdots & (0,N-1) & (0, N) \cr
 * (1,0) & (1,1) & \cdots & (1,N-1) \cr
 * \cdots \cr
 * (0,N) \end{array} \right\} \f]
 *
 * @param [in] N polynomial order
 * @param [in] ind orthogonal polynomial index
 * @param [out] ti
 * @param [out] tj
 */
void sc_triTransInd(int N, int ind, int *ti, int *tj){
    int i,j,sk=0;
    for(i=0;i<N+1;i++){
        for(j=0;j<N-i+1;j++){
            if(sk == ind){
                *ti = i;
                *tj = j;
            }
            sk++;
        }
    }
}

/**
 * @brief get orthogonal function value at interpolation nodes
 *
 * @param [in] cell triangle cell
 * @param [in] ind index of orthogonal function
 * @param [out] func value of orthogonal function
 */
void sc_triOrthogFunc(stdCell *cell, int ind, double *func){
    const int N = cell->N;
    const int Np = cell->Np;
    double *r = cell->r;
    double *s = cell->s;

    double a[Np], b[Np];
    sc_rstoad(Np, r, s, a, b);
    int i,j;
    sc_triTransInd(N, ind, &i, &j);
    sc_triSimplex2dP(Np, a, b, i, j, func);
}

/**
 * @brief create standard triangle element
 * @param[in] N polynomial order
 * @note
 * postcondition: the stdCell object should be free manually with @ref sc_free
 */
stdCell* sc_createTri(int N){
    stdCell *tri = (stdCell *) calloc(1, sizeof(stdCell));

    /* cell type */
    tri->type = TRIANGLE;

    const int Np = ((N+1)*(N+2)/2);
    const int Nv = 3;
    const int Nfaces = 3;
    const int Nfp = N+1;
    /* basic info */
    tri->N      = N;
    tri->Np     = Np;
    tri->Nv     = Nv;
    tri->Nfaces = Nfaces;
    tri->Nfp    = Nfp;

    /* nodes at faces, Fmask */
    tri->Fmask = sc_triFmask(tri);

    /* coordinate, r and s */
    sc_triCoord(tri);
    /* vandermonde matrix, V */
    tri->V = sc_VandMatrix2d(tri, sc_triOrthogFunc);
    /* mass matrix, M */
    tri->M = sc_massMatrix(tri);

    /* Derivative Matrix, Dr and Ds */
    sc_deriMatrix2d(tri, sc_triDeriOrthogFunc);

    /* suface LIFT matrix, LIFT */
    tri->LIFT = sc_liftMatrix2d(tri);

    /* integration coefficients, ws and wv */
    sc_GaussQuadrature2d(tri);

    /* float version */
    size_t sz = Np*Nfp*Nfaces*sizeof(real);
    tri->f_LIFT = (real *) malloc(sz);
    sz = Np*Np*sizeof(real);
    tri->f_Dr = (real*) malloc(sz);
    tri->f_Ds = (real*) malloc(sz);

    int sk = 0, n, m;
    for(n=0;n<tri->Np;++n){
        for(m=0;m<tri->Nfp*tri->Nfaces;++m){
            tri->f_LIFT[sk++] = (real) tri->LIFT[n][m];
        }
    }

    sk = 0;
    for(n=0;n<tri->Np;++n){
        for(m=0;m<tri->Np;++m){
            tri->f_Dr[sk  ] = (real) tri->Dr[n][m];
            tri->f_Ds[sk++] = (real) tri->Ds[n][m];
        }
    }
    return tri;
}

/**
 * @brief calculate the value of derivative function at interpolation points.
 * @details
 * @param [in] tri standard triangle element
 * @param [in] ind index of orthogonal function
 * @param [out] dr derivative function with r
 * @param [out] ds derivative function with s
 * @note
 * precondition: dr and ds should be allocated before calling @ref sc_triDeriOrthogFunc
 */
void sc_triDeriOrthogFunc(stdCell *tri, int ind, double *dr, double *ds){
    const int Np = tri->Np;
    double a[Np], b[Np];
    sc_rstoad(Np, tri->r, tri->s, a, b);
    int ti, tj;
    sc_triTransInd(tri->N, ind, &ti, &tj);
    sc_triGradSimplex2DP(Np, a, b, ti, tj, dr, ds);
}

/**
 * @brief
 * get the gradient of the modal basis (id,jd) on the 2D simplex at (a,b).
 * @details
 * Coordinate (a,b) is the collapse coordinate, given by
 * \f[ a = 2\frac{1+r}{1-s}, \quad b=s \f]
 * Therefore, the derivatives of (r,s) is obtained by
 * \f[\begin{array}{ll}
 * \frac{\partial }{\partial r} = \frac{\partial a}{\partial r}\frac{\partial }{\partial a} +
 * \frac{\partial b}{\partial r}\frac{\partial }{\partial b} =
 * \frac{2}{1-b}\frac{\partial }{\partial a} \cr
 * \frac{\partial }{\partial s} = \frac{\partial a}{\partial s}\frac{\partial }{\partial a} +
 * \frac{\partial b}{\partial s}\frac{\partial }{\partial b} =
 * \frac{(1+a)/2}{(1-b)/2}\frac{\partial }{\partial a} + \frac{\partial }{\partial b} \cr
 * \end{array}\f]
 *
 * @param [in] Np number of coordinate
 * @param [in] a collapse coordinate
 * @param [in] b collapse coordinate
 * @param [in] id index of
 * @param [in] jd
 * @param [out] dmodedr
 * @param [out] dmodeds
 * @note
 * precondition: dmodedr and dmodeds should be allocated before calling @ref sc_triGradSimplex2DP
 */
void sc_triGradSimplex2DP(int Np, double *a, double *b, int id, int jd, double *dmodedr, double *dmodeds){

    double fa[Np],dfa[Np],gb[Np],dgb[Np],temp[Np];
    int i;

    jacobiP(Np, a, fa, id, 0, 0);
    jacobiP(Np, b, gb, jd, 2*id+1, 0);
    GradjacobiP(Np, a, dfa, id, 0, 0);
    GradjacobiP(Np, b, dgb, jd, 2*id+1, 0);

    /* r-derivative */
    for(i=0;i<Np;i++)
        dmodedr[i] = dfa[i]*gb[i];
    if(id>0){
        for(i=0;i<Np;i++){
            dmodedr[i] *= pow(0.5*(1-b[i]), (id-1));
        }
    }
    /* s-derivative */
    for(i=0;i<Np;i++){
        dmodeds[i] = dfa[i]*gb[i]*0.5*(1.0+a[i]);
    }
    if(id>0){
        for(i=0;i<Np;i++){
            dmodeds[i] *= pow(0.5*(1.0-b[i]), (id-1));
        }
    }

    for(i=0;i<Np;i++){
        temp[i] = dgb[i]*pow(0.5*(1.0-b[i]), id);
    }
    if(id>0){
        for(i=0;i<Np;i++){
            temp[i] -= 0.5*(double)id*gb[i] * pow(0.5*(1.0-b[i]), (id-1));
        }
    }

    for(i=0;i<Np;i++) {
        dmodeds[i] += fa[i]*temp[i];
    }

    /* Normalize */
    for(i=0;i<Np;i++) {
        dmodedr[i] *= pow(2.0, id+0.5);
        dmodeds[i] *= pow(2.0, id+0.5);
    }

}

/**
 * @brief
 * evaluate 2D orthonormal polynomial on simplex at (a,b) of order (i,j).
 * @details
 * The orthogonal basis function \f$ \varphi(\mathbf{r}) \f$ is obtained by
 * \f[ \varphi_n(\mathbf{r}) = \sqrt{2}P_i(a)P_j^{(2i+1, 0)}(b)(1-b)^i \f]
 *
 * @param [in] Np number of points
 * @param [in] a coordinate
 * @param [in] b coordinate
 * @param [in] i order
 * @param [in] j order
 * @param [out] poly value of orthogonal basis function value
 *
 */
void sc_triSimplex2dP(int Np, double *a, double *b, int i, int j, double *poly){
    double *h1 = Vector_create(Np);
    double *h2 = Vector_create(Np);
    int n;

    jacobiP(Np, a, h1, i, 0.0, 0.0);
    jacobiP(Np, b, h2, j, 2*i+1, 0.0);

    for(n=0;n<Np;n++){
        poly[n] = sqrt(2.0)*h1[n]*h2[n]*pow(1-b[n], i);
    }
    Vector_free(h1);
    Vector_free(h2);
}

/**
 * @brief get the coordinate of interpolations points on standard triangle element
 * @param[in] N polynomial order
 * @param[in,out] r coordinate
 * @param[in,out] s coordinate
 *
 */
void sc_triCoord(stdCell *tri){

    // local
    double alpopt[15] =  {0.0000, 0.0000, 1.4152, 0.1001, 0.2751, 0.9800, 1.0999, 1.2832,
                          1.3648, 1.4773, 1.4959, 1.5743, 1.5770, 1.6223, 1.6258};
    const int N = tri->N;
    const int Np = tri->Np;
    int i,j,sk=0;
    double alpha, temp;

    if(N<16){
        alpha = alpopt[N-1];
    }else{
        alpha = 5.0/3.0;
    }

    // allocations
    tri->r = Vector_create(Np);
    tri->s = Vector_create(Np);

    double L1[Np], L2[Np], L3[Np];
    double x[Np], y[Np], dL[Np], warpf1[Np];

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
    sc_triWarpfactor(N, dL, Np, warpf1);

    for(i=0;i<Np;i++){
        temp = alpha*L1[i];
        warpf1[i] *= 4.0*L2[i]*L3[i]*(1.0 + temp*temp);
        x[i] += 1.0*warpf1[i];
    }

    /* warp2 */
    for(i=0;i<Np;i++){
        dL[i] = L1[i] - L3[i];
    }
    sc_triWarpfactor(N, dL, Np, warpf1);

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
    sc_triWarpfactor(N, dL, Np, warpf1);

    for(i=0;i<Np;i++){
        temp = alpha*L3[i];
        warpf1[i] *= 4.0*L1[i]*L2[i]*(1.0 + temp*temp);
        x[i] += cos(4.0*M_PI/3.0)*warpf1[i];
        y[i] += sin(4.0*M_PI/3.0)*warpf1[i];
    }

    /* coordinate transfer */
    sc_trixytors(Np, x, y, tri->r, tri->s);
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
 * @param [in] N order of degree
 * @param [in] r input coordinate
 * @param [in] Nr number of input points
 * @param [out] w waper factors
 *
 * @note
 * precondition: w shuold be allocated before calling sc_triWarpfactor
 */
void sc_triWarpfactor(int N, double *r, int Nr, double *w){
    int i, j, Np = N+1;
    double *ye, *l, *re, *rlgl, *wlgl;
    double temp;

    ye   = Vector_create(Np);
    l    = Vector_create(Nr);
    re   = Vector_create(Np); /* equidistant nodes */
    rlgl = Vector_create(Np); /* Gauss-Lobatto-Jacobi nodes */
    wlgl = Vector_create(Np); /* Gauss-Lobatto-Jacobi weights */

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
    Vector_free(ye);
    Vector_free(l);
    Vector_free(re);
    Vector_free(rlgl);
    Vector_free(wlgl);
}

/**
 * @brief build the nodes index matrix on each faces
 * @details
 * Three faces of the standard triangle element is
 * \f[ s=-1, \quad r+s=0, \quad r=-1 \f]
 *
 * @param [in] N order
 * @return int Fmask[Nfaces][Nfp]
 *
 */
int** sc_triFmask(stdCell *tri){

    const int Nfp = tri->Nfp;
    const int Nfaces = tri->Nfaces;
    int **Fmask = IntMatrix_create(Nfaces, Nfp);

    int temp[Nfp];
    int nrp[Nfp]; /* # of nodes from s=-1 to s=1 */
    int i;

    /* face 1, s=-1 */
    for(i=0;i<Nfp;i++)
        Fmask[0][i] = i;

    /* face 3, r=-1 */
    for(i=0;i<Nfp;i++)
        nrp[i] = Nfp - i;

    temp[0] = 1; /* node index on r=-1 from s=-1 to s=1 */
    for(i=1;i<Nfp;i++){
        temp[i] = temp[i-1] + nrp[i-1];
    }
    for(i=0;i<Nfp;i++){
        /* node index on r=-1 from s=1 to s=-1 */
        Fmask[2][i] = temp[Nfp-1-i] - 1;
    }
    /* face 2, r+s=0 */
    nrp[0] = Nfp;
    for(i=1;i<Nfp-1;i++)
        nrp[i] = temp[i+1]-1;
    nrp[Nfp-1] = temp[Nfp-1];
    for(i=0;i<Nfp;i++){
        Fmask[1][i] = nrp[i] - 1;
    }

    return Fmask;
}


/**
 * @brief
 * transfer coordinate [x,y] on equilateral triangle to natural coordinate [r,s]
 *
 * @param [in] Np number of points
 * @param [in] x coordinate of equilateral triangle
 * @param [in] y coordinate of equilateral triangle
 * @param [out] r coordinate of standard triangle
 * @param [out] s coordinate of standard triangle
 *
 * @note
 * precondition: [x, y] and [r, s] should be allocated before calling @ref sc_trixytors
 */
void sc_trixytors(int Np, double *x, double *y, double *r, double *s){
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
 * @brief transform natural coordinate to collapse coordinate
 *
 * @param [in] Np number of points
 * @param [in] r coordinate
 * @param [in] s coordinate
 * @param [out] a collapse coordinate
 * @param [out] b collapse coordinate
 *
 * @note
 * precondition: a and b should be allocated before calling @ref sc_rstoad
 */

void sc_rstoad(int Np, double *r, double *s, double *a, double *b){
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