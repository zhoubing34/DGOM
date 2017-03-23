/**
 * @file Standard triangle element
 * @brief Functions related with standard triangle elements
 * @author li12242, Tianjin University, li12242@tju.edu.cn
 */

#include "dg_cell_triangle.h"
#include "Polylib/polylib.h"

/* transform the index of orthogonal function to [ti,tj] */
static void dg_tri_transInd(int N, int ind, int *ti, int *tj);
/* transfer coordinate [x,y] on equilateral triangle to natural coordinate [r,s] */
static void dg_tri_xytors(int Np, double *x, double *y, double *r, double *s);
/* transform natural coordinate to collapse coordinate */
static void dg_tri_rstoad(int Np, double *r, double *s, double *a, double *b);
/* get the gradient of the modal basis (ncid,jd) on the 2D simplex at (a,b). */
static void dg_tri_gradSimplex2DP(int Np, double *a, double *b, int id, int jd,
                                  double *dmodedr, double *dmodeds);
/* Warp factor to connnect the Legendre-Gauss-Lobatto and equidistant nodes */
static void dg_tri_warpfactor(int, double *, int, double *);
/* evaluate 2D orthonormal polynomial on simplex at (a,b) of order (i,j). */
static void dg_tri_simplex2DP(int Np, double *a, double *b, int i, int j, double *poly);

dg_cell_info* dg_cell_tri_info(int N){
    dg_cell_info *info = (dg_cell_info *)calloc(1, sizeof(dg_cell_info));
    const int Nv = 3;
    const int Nfaces = 3;
    const int Nfv = 2;
    info->N = N;
    info->Nv = Nv;
    info->Nfaces = Nfaces;
    info->type = TRIANGLE;
    info->face_type = calloc(Nfaces, sizeof(dg_cell_type));
    info->FToV = matrix_int_create(Nfaces, Nfv);
    info->vr = (double *)calloc(Nv, sizeof(double));
    info->vs = (double *)calloc(Nv, sizeof(double));
    info->vt = (double *)calloc(Nv, sizeof(double));
    /* faces */
    info->face_type[0] = LINE; info->FToV[0][0] = 0; info->FToV[0][1] = 1;
    info->face_type[1] = LINE; info->FToV[1][0] = 1; info->FToV[1][1] = 2;
    info->face_type[2] = LINE; info->FToV[2][0] = 2; info->FToV[2][1] = 0;
    /* vertex */
    info->vr[0] = -1; info->vr[1] =  1; info->vr[2] = -1;
    info->vs[0] = -1; info->vs[1] = -1; info->vs[2] =  1;
    return info;
}

/**
 * @brief
 * transform the index of orthogonal function to [ti,tj] for triangle elements
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
static void dg_tri_transInd(int N, int ind, int *ti, int *tj){
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
 * @brief
 * get orthogonal function value at interpolation nodes
 * @param [in] cell triangle cell
 * @param [in] ind index of orthogonal function
 * @param [out] func value of orthogonal function
 */
void dg_cell_tri_orthog_func(int N, int ind, int Np, double *r, double *s, double *t, double *func){
    double a[Np], b[Np];
    dg_tri_rstoad(Np, r, s, a, b);
    int i,j;
    dg_tri_transInd(N, ind, &i, &j);
    dg_tri_simplex2DP(Np, a, b, i, j, func);
}

/**
 * @brief calculate the value of derivative function at interpolation points.
 * @details
 * @param [in] tri standard triangle element
 * @param [in] ind index of orthogonal function
 * @param [out] dr derivative function with r
 * @param [out] ds derivative function with s
 * @note
 * precondition: dr and ds should be allocated before calling @ref sc_deriOrthogFunc_tri
 */
void dg_cell_tri_deriorthog_func(int N, int ind, int Np, double *r, double *s, double *t,
                                 double *dr, double *ds, double *dt){
    double a[Np], b[Np];
    dg_tri_rstoad(Np, r, s, a, b);
    int ti, tj;
    dg_tri_transInd(N, ind, &ti, &tj);
    dg_tri_gradSimplex2DP(Np, a, b, ti, tj, dr, ds);
    return;
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
 * precondition: dmodedr and dmodeds should be allocated before
 * calling @ref sc_gradSimplex2DP_tri
 */
static void dg_tri_gradSimplex2DP(int Np, double *a, double *b,
                                  int id, int jd,
                                  double *dmodedr, double *dmodeds){

    double fa[Np],dfa[Np],gb[Np],dgb[Np],temp[Np];
    int i;

    jacobiP(Np, a, fa, id, 0, 0);
    jacobiP(Np, b, gb, jd, 2*id+1, 0);
    GradjacobiP(Np, a, dfa, id, 0, 0);
    GradjacobiP(Np, b, dgb, jd, 2*id+1, 0);

    /* r-derivative */
    for(i=0;i<Np;i++) {dmodedr[i] = dfa[i]*gb[i];}
    if(id>0){
        for(i=0;i<Np;i++){ dmodedr[i] *= pow(0.5*(1-b[i]), (id-1)); }
    }
    /* s-derivative */
    for(i=0;i<Np;i++){ dmodeds[i] = dfa[i]*gb[i]*0.5*(1.0+a[i]); }
    if(id>0){
        for(i=0;i<Np;i++){ dmodeds[i] *= pow(0.5*(1.0-b[i]), (id-1)); }
    }
    for(i=0;i<Np;i++){
        temp[i] = dgb[i]*pow(0.5*(1.0-b[i]), id);
    }
    if(id>0){
        for(i=0;i<Np;i++){ temp[i] -= 0.5*(double)id*gb[i] * pow(0.5*(1.0-b[i]), (id-1)); }
    }
    for(i=0;i<Np;i++) { dmodeds[i] += fa[i]*temp[i]; }
    /* Normalize */
    for(i=0;i<Np;i++) {
        dmodedr[i] *= pow(2.0, id+0.5);
        dmodeds[i] *= pow(2.0, id+0.5);
    }
    return;
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
static void dg_tri_simplex2DP(int Np, double *a, double *b, int i, int j, double *poly){
    double *h1 = vector_double_create(Np);
    double *h2 = vector_double_create(Np);
    int n;

    jacobiP(Np, a, h1, i, 0.0, 0.0);
    jacobiP(Np, b, h2, j, 2*i+1, 0.0);

    for(n=0;n<Np;n++){ poly[n] = sqrt(2.0)*h1[n]*h2[n]*pow(1-b[n], i); }
    vector_double_free(h1);
    vector_double_free(h2);
    return;
}

/**
 * @brief get the coordinate of interpolations points on standard triangle element
 * @param[in] N polynomial order
 * @param[in,out] r coordinate
 * @param[in,out] s coordinate
 *
 */
void dg_cell_tri_set_node(dg_cell *cell, int *Np, double **r, double **s, double **t){

    // local
    double alpopt[15] =  {0.0000, 0.0000, 1.4152, 0.1001, 0.2751, 0.9800, 1.0999, 1.2832,
                          1.3648, 1.4773, 1.4959, 1.5743, 1.5770, 1.6223, 1.6258};

    const int N = dg_cell_N(cell);
    const int Npt = ((N+1)*(N+2)/2);
    int i,j,sk=0;
    double alpha, temp;

    if(N<16){ alpha = alpopt[N-1]; }
    else{  alpha = 5.0/3.0; }

    // allocations
    *Np = Npt;
    *r = vector_double_create(Npt);
    *s = vector_double_create(Npt);
    *t = vector_double_create(Npt);

    double L1[Npt], L2[Npt], L3[Npt];
    double x[Npt], y[Npt], dL[Npt], warpf1[Npt];

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
    for(i=0;i<Npt;i++){ dL[i] = L3[i] - L2[i]; }
    dg_tri_warpfactor(N, dL, Npt, warpf1);

    for(i=0;i<Npt;i++){
        temp = alpha*L1[i];
        warpf1[i] *= 4.0*L2[i]*L3[i]*(1.0 + temp*temp);
        x[i] += 1.0*warpf1[i];
    }

    /* warp2 */
    for(i=0;i<Npt;i++){ dL[i] = L1[i] - L3[i]; }
    dg_tri_warpfactor(N, dL, Npt, warpf1);

    for(i=0;i<Npt;i++){
        temp = alpha*L2[i];
        warpf1[i] *= 4.0*L1[i]*L3[i]*(1.0 + temp*temp);
        x[i] += cos(2.0*M_PI/3.0)*warpf1[i];
        y[i] += sin(2.0*M_PI/3.0)*warpf1[i];
    }

    /* warp3 */
    for(i=0;i<Npt;i++){ dL[i] = L2[i] - L1[i]; }
    dg_tri_warpfactor(N, dL, Npt, warpf1);

    for(i=0;i<Npt;i++){
        temp = alpha*L3[i];
        warpf1[i] *= 4.0*L1[i]*L2[i]*(1.0 + temp*temp);
        x[i] += cos(4.0*M_PI/3.0)*warpf1[i];
        y[i] += sin(4.0*M_PI/3.0)*warpf1[i];
    }

    /* coordinate transfer */
    dg_tri_xytors(Npt, x, y, *r, *s);
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
 * precondition: w shuold be allocated before calling sc_warpfactor_tri
 */
static void dg_tri_warpfactor(int N, double *r, int Nr, double *w){
    int i, j, Np = N+1;
    double *ye, *l, *re, *rlgl, *wlgl;
    double temp;

    ye   = vector_double_create(Np);
    l    = vector_double_create(Nr);
    re   = vector_double_create(Np); /* equidistant nodes */
    rlgl = vector_double_create(Np); /* Gauss-Lobatto-Jacobi nodes */
    wlgl = vector_double_create(Np); /* Gauss-Lobatto-Jacobi weights */
    /* initialization */
    for(i=0;i<Nr;i++){ w[i] = 0.0; }
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

        for(j=0;j<Nr;j++){ w[j] += l[j]*(rlgl[i] - re[i]); }
        ye[i] = 0.0;
    }
    for(i=0;i<Nr;i++){
        temp = (1.0 - r[i]*r[i]);
        if(temp > 1.0e-10){ w[i] /= temp; }
    }
    /* deallocate mem */
    vector_double_free(ye);
    vector_double_free(l);
    vector_double_free(re);
    vector_double_free(rlgl);
    vector_double_free(wlgl);
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
 * precondition: [x, y] and [r, s] should be allocated before calling @ref sc_xytors_tri
 */
static void dg_tri_xytors(int Np, double *x, double *y, double *r, double *s){
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
static void dg_tri_rstoad(int Np, double *r, double *s, double *a, double *b){
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

/**
 * @brief
 * Project the triangular vertex value to interpolation nodes.
 * @param[in] cell pointer to dg_cell structure;
 * @param[in] Nfield number of field;
 * @param[in] vertVal value of vertex;
 * @param[out] nodeVal value of nodes;
 *
 */
void dg_cell_tri_proj(dg_cell *cell, int Nfield, double *vertVal, double *nodeVal){
    register int i,fld,sk=0;
    double *r = dg_cell_r(cell);
    double *s = dg_cell_s(cell);
    const int Np = dg_cell_Np(cell);
    for (i=0;i<Np;++i) {
        double ri=r[i];
        double si=s[i];
        for(fld=0;fld<Nfield;fld++){
            nodeVal[sk++] = 0.5*( -vertVal[0*Nfield+fld]*(ri + si)
                                  + vertVal[1*Nfield+fld]*(1.+ri)
                                  + vertVal[2*Nfield+fld]*(1.+si));
        }
    }
}