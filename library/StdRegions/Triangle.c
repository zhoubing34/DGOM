#include <stdlib.h>
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

#include "StdRegions.h"

/* private functions of triangle element */
void xytors(int Np, double *x, double *y, double *r, double *s);
void rstoad(int Np, double *r, double *s, double *a, double *b);
void GradSimplex2DP(int Np, double *a, double *b, int id, int jd, double *dmodedr, double *dmodeds);
void GetTriDeriV2d(int N, int Np, double *r, double *s, double **Vr, double **Vs);
void GetTriDeriM(int N, int Np, double *r, double *s, double **V, double **Dr, double **Ds);
void BuildTriFmask(int N, int **Fmask);
void Warpfactor(int, double *, int, double *);
void GetTriCoord(int, double*, double*);
void GetTriV(int N, int Nr, double *r, double *s, double **V);
void GetTriM(int Np, double **V, double **M);

/**
 * @brief
 * Generation of standard triangle element
 *
 * @details
 *
 * @param[in] N polynomial order
 *
 * @return
 * name     | type     | description of value
 * -------- |----------|----------------------
 * tri | StdRegions2d* | structure of standard triangle element
 *
 */
StdRegions2d* StdTriEle_create(int N){
    StdRegions2d *tri = (StdRegions2d *) calloc(1, sizeof(StdRegions2d));

    /* basic info */
    tri->N      = N;
    tri->Np     = ((N+1)*(N+2)/2);
    tri->Nv     = 3;
    tri->Nfaces = 3;
    tri->Nfp    = N+1;

    /* nodes at faces, Fmask */
    tri->Fmask = IntMatrix_create(tri->Nfaces, tri->Nfp);
    BuildTriFmask(tri->N, tri->Fmask);

    /* coordinate, r and s */
    tri->r = Vector_create(tri->Np);
    tri->s = Vector_create(tri->Np);

    GetTriCoord(N, tri->r, tri->s);

    /* vandermonde matrix, V */
    tri->V = Matrix_create(tri->Np, tri->Np);
    GetTriV(N, tri->Np, tri->r, tri->s, tri->V);

    /* mass matrix, M */
    tri->M = Matrix_create(tri->Np, tri->Np);
    GetTriM(tri->Np, tri->V, tri->M);

    /* Derivative Matrix, Dr and Ds */
    tri->Dr = Matrix_create(tri->Np, tri->Np);
    tri->Ds = Matrix_create(tri->Np, tri->Np);
    GetTriDeriM(tri->N, tri->Np, tri->r, tri->s, tri->V, tri->Dr, tri->Ds);

    /* suface LIFT matrix, LIFT */
    double **Mes = Matrix_create(tri->Np, tri->Nfaces * tri->Nfp);
    tri->LIFT = Matrix_create(tri->Np, tri->Nfaces * tri->Nfp);
    GetSurfLinM(tri->N, tri->Nfaces, tri->Fmask, Mes);
    GetLIFT2d(tri, Mes, tri->LIFT);

    Matrix_free(Mes);

    /* integration coeff, ws and wv */
    tri->ws = Vector_create(tri->Nfp);
    double *r = Vector_create(tri->Nfp);
    zwglj(r, tri->ws, tri->Nfp, 0, 0);
    Vector_free(r);

    tri->wv = Vector_create(tri->Np);
    int i,j;
    for(i=0;i<tri->Np;i++){
        for(j=0;j<tri->Np;j++){
            tri->wv[j] += tri->M[i][j];
        }
    }

    /* float version */
    size_t sz = tri->Np*(tri->Nfp)*(tri->Nfaces)*sizeof(real);
    tri->f_LIFT = (real *) malloc(sz);
    sz = tri->Np*tri->Np*sizeof(real);
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
            tri->f_Dr[sk] = (real) tri->Dr[n][m];
            tri->f_Ds[sk] = (real) tri->Ds[n][m];
            ++sk;
        }
    }

    return tri;
}

/**
 * @brief
 * Get LIFT matrix
 *
 * @param [StdRegions2d] shape               standard element
 * @param [double]       Mes[Np][Nfp*Nfaces] mass matrix of surface integral
 *
 * @return
 * return values:
 * name     | type     | description of value
 * -------- |----------|----------------------
 * LIFT   | double[Np][Nfp*Nfaces] | LIFT = M^{-1}*Mes
 *
 */
void GetLIFT2d(StdRegions2d *shape, double **Mes, double **LIFT){
    unsigned int Np = (unsigned)shape->Np;
    int Nfp = shape->Nfp, Nfaces=shape->Nfaces;
    double *vc = Vector_create(Np * Np); /* vandermonde matrix */
    double *vct = Vector_create(Np * Np); /* transform of vandermonde matrix */
    double *me = Vector_create(Np * Nfp * Nfaces);
    double *temp1 = Vector_create(Np * Np); /* inverse of mass matrix */
    double *temp2 = Vector_create(Np * Nfp * Nfaces);

    int i,j,sk=0;
    /* assignment */
    for(i=0;i<Np;i++){
        for(j=0;j<Np;j++){
            /* row counts first */
            vc[sk] = shape->V[i][j];
            vct[sk++] = shape->V[j][i];
        }
    }

    sk = 0;
    for(i=0;i<Np;i++){
        for(j=0;j<Nfp*Nfaces;j++) {
            me[sk++] = Mes[i][j];
        }
    }

    /* get the inverse mass matrix M^{-1} = V*V' */
    dgemm_(Np, Np, Np,  vc, vct, temp1);
    /* LIFT = M^{-1}*Mes */
    dgemm_(Np, Np, (unsigned)Nfp*Nfaces, temp1, me, temp2);

    sk = 0;
    for(i=0;i<Np;i++){
        for(j=0;j<Nfp*Nfaces;j++)
            LIFT[i][j] = temp2[sk++];
    }

    Vector_free(vc);
    Vector_free(vct);
    Vector_free(me);
    Vector_free(temp1);
    Vector_free(temp2);
}

/**
 * @brief
 * Get surface mass matrix
 * @details
 * @param [int] N order
 * @param [int] Fmask[Nfaces][Nfp] nodes list at faces
 * @return
 * name  | type     | description of value
 * ----- |----------|----------------------
 * Mes   | double[Np][Nfaces*Nfp] | surface mass matrix
 *
 */
void GetSurfLinM(int N, int Nfaces, int **Fmask, double **Mes){

    unsigned Nfp=(unsigned)N+1;
    double *r= Vector_create(Nfp);
    double *w= Vector_create(Nfp);

    double *invt = Vector_create(Nfp * Nfp);
    double *inv = Vector_create(Nfp * Nfp);
    double *m = Vector_create(Nfp * Nfp);
    /* get surface mass matrix */
    zwglj(r, w, Nfp, 0, 0); /* get coordinate */
    int i,j;
    for(i=0;i<Nfp;i++){
        /* get vandermonde matrix */
        jacobiP(Nfp, r, w, i, 0, 0);
        for(j=0;j<Nfp;j++){
            inv[j*Nfp+i] = w[j];
        }
    }
    invM(inv, Nfp);
    /* transform of vandermonde matrix */
    for(i=0;i<Nfp;i++){
        for(j=0;j<Nfp;j++)
            invt[j+Nfp*i] = inv[j*Nfp+i];
    }
    /* get M = inv(V)'*inv(V) */
    dgemm_(Nfp, Nfp, Nfp, invt, inv, m);

    int k, sr, sk;
    for(i=0;i<Nfaces;i++){
        for(j=0;j<Nfp;j++){         /* row index of M */
            for(k=0;k<Nfp;k++){     /* column index of M */
                sr = Fmask[i][j];   /* row index of Mes */
                sk = i*Nfp + k;     /* columns index of Mes */
                Mes[sr][sk] = m[j*Nfp + k];
            }
        }
    }

    Vector_free(invt);
    Vector_free(inv);
    Vector_free(m);
    Vector_free(r);
    Vector_free(w);
}

/**
 * @brief
 * Get the Gradient matrix of Lagrange basis at (r,s) at order N
 * @details
 * The Gradient matrix \f$ \mathbf{Dr} \f$ and \f$ \mathbf{Ds} \f$ is obtained through
 * \f[ \mathbf{Dr} \cdot \mathbf{V} = \mathbf{Vr},
 * \quad \mathbf{Ds} \cdot \mathbf{V} = \mathbf{Vs} \f]
 * where
 * \f[ Dr_{(ij)} = \left. \frac{\partial l_j}{\partial r} \right|_{ \mathbf{r}_i },
 * \quad Ds_{(ij)} = \left. \frac{\partial l_j}{\partial s} \right|_{ \mathbf{r}_i } \f]
 *
 * @param [unsigned] N
 * @param [unsigned] Np number of points
 * @param [double]   r[Np] coordinate
 * @param [double]   s[Np] coordinate
 * @param [double]   V[Np][Np] vandermonde matrix
 *
 * @return
 * name     | type     | description of value
 * -------- |----------|----------------------
 * Dr | double[Np][Np]   |
 * Ds | object[Np][Np]   |
 */
void GetTriDeriM(int N, int Np, double *r, double *s, double **V, double **Dr, double **Ds){
    double **Vr = Matrix_create(Np, Np);
    double **Vs = Matrix_create(Np, Np);
    double *temp = Vector_create(Np * Np);
    double *vtemp = Vector_create(Np * Np);
    double *dtemp = Vector_create(Np * Np);

    int i,j,sk=0;

    /* inverse of vandermonde matrix */
    for(i=0;i<Np;i++){
        for(j=0;j<Np;j++){
            /* row counts first */
            temp[sk++] = V[i][j];
        }
    }
    invM(temp, Np);

    /* get the gradient of the modal basis */
    GetTriDeriV2d(N, Np, r, s, Vr, Vs);

    /* matrix multiply */
    unsigned n = (unsigned)Np;
    sk = 0;
    for(i=0;i<Np;i++){
        for(j=0;j<Np;j++){
            /* row counts first */
            vtemp[sk++] = Vr[i][j];
        }
    }
    /* \f$ \mathbf{Dr} = \mathbf{Vr}*\mathbf{V}^{-1} \f$ */
    dgemm_(n, n, n, vtemp, temp, dtemp);
    sk = 0;
    for(i=0;i<Np;i++){
        for(j=0;j<Np;j++){
            /* row counts first */
            Dr[i][j] = dtemp[sk++];
        }
    }

    sk = 0;
    for(i=0;i<Np;i++){
        for(j=0;j<Np;j++){
            /* row counts first */
            vtemp[sk++] = Vs[i][j];
        }
    }
    /* \f$ \mathbf{Ds} = \mathbf{Vs}*\mathbf{V}^{-1} \f$ */
    dgemm_(n, n, n, vtemp, temp, dtemp);
    sk = 0;
    for(i=0;i<Np;i++){
        for(j=0;j<Np;j++){
            /* row counts first */
            Ds[i][j] = dtemp[sk++];
        }
    }

    /* dealloc mem */
    Matrix_free(Vr);
    Matrix_free(Vs);
    Vector_free(temp);
    Vector_free(vtemp);
    Vector_free(dtemp);
}

/**
 * @brief
 * get the gradient matrix of the modal basis (i,j) at (r,s) at order N
 * @details
 *
 * @param [unsigned] N order
 * @param [int]      Np number of points
 * @param [double]   r[Np] coordinate
 * @param [double]   s[Np] coordinate
 *
 * @return
 * name     | type     | description of value
 * -------- |----------|----------------------
 * Vr | double[Np][Np] | gradient matrix of r
 * Vs | double[Np][Np] | gradient matrix of s
 *
 */
void GetTriDeriV2d(int N, int Np, double *r, double *s, double **Vr, double **Vs){
    double *a = Vector_create(Np);
    double *b = Vector_create(Np);
    double *dpdr, *dpds;
    int i,j,k,sk=0;

    dpdr = Vector_create(Np); dpds = Vector_create(Np);

    rstoad(Np, r, s, a, b);
    for(i=0;i<N+1;i++){
        for(j=0;j<N-i+1;j++){
            GradSimplex2DP(Np, a, b, i, j, dpdr, dpds);
            /* assignment */
            for(k=0;k<Np;k++){
                Vr[k][sk] = dpdr[k];
                Vs[k][sk] = dpds[k];
            }
            sk++;
        }
    }
    Vector_free(dpdr);
    Vector_free(dpds);
    Vector_free(a);
    Vector_free(b);

}

/**
 * @brief
 * Get the gradient of the modal basis (id,jd) on the 2D simplex at (a,b).
 * @details
 * Coordinate (a,b) is calculated from (r,s) by
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
 * @param [int]      Np number of coordinate
 * @param [double]   a[Np]
 * @param [double]   b[Np]
 * @param [int]      id
 * @param [int]      jd
 *
 * @return
 * name     | type     | description of value
 * -------- |----------|----------------------
 * car_id   | int      | 车源编号
 * car_info | object   | json对象格式的车源信息
 *
 */
void GradSimplex2DP(int Np, double *a, double *b, int id, int jd, double *dmodedr, double *dmodeds){
    double *fa  = Vector_create(Np);
    double *dfa = Vector_create(Np);
    double *gb  = Vector_create(Np);
    double *dgb = Vector_create(Np);
    double *temp = Vector_create(Np);
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

    Vector_free(fa);
    Vector_free(dfa);
    Vector_free(gb);
    Vector_free(dgb);
    Vector_free(temp);
}


/**
 * @brief
 * Generate the mass matrix
 * @details
 * The mass matrix is calculated with
 * \f[ \mathbf{M} = (\mathbf{V}^T)^{-1} \cdot \mathbf{V}^{-1} \f]
 *
 * @param [unsigned] Np
 * @param [double]   V[Np][Np] Vandermonde matrix
 *
 * @return
 * return values:
 * name     | type     | description of value
 * -------- |----------|----------------------
 * M   | double[Np][Np] | Mass Matrix
 */
void GetTriM(int Np, double **V, double **M){
    double *temp = Vector_create(Np * Np);
    double *invt = Vector_create(Np * Np);
    double *Mv = Vector_create(Np * Np);
    const unsigned n = (unsigned)Np;

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

    dgemm_(n, n, n, invt, temp, Mv);

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
 * The orthogonal basis function \f$ \varphi(\mathbf{r}) \f$ is obtained by
 * \f[ \varphi_n(\mathbf{r}) = \sqrt{2}P_i(a)P_j^{(2i+1, 0)}(b)(1-b)^i \f]
 *
 * @param [int]      Np number of points
 * @param [double]   a[Np] coordinate
 * @param [double]   b[Np] coordinate
 * @param [int]      i order
 * @param [int]      j order
 *
 * @return
 * name     | type     | description of value
 * -------- |----------|----------------------
 * poly   | double[Np] | Orthogonal Basis Function
 *
 */
void Simplex2dP(int Np, double *a, double *b, int i, int j, double *poly){
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
 * @brief
 * Generate the Vandermonde matrix of triangle element
 *
 * @details
 * Vandermonde matrix \f$ \mathbf{V} \f$ is derivated by
 * \f[ \mathbf{V}_{mn} = \varphi_n \left( \mathbf{r}_m \right) \f]
 * where \f$ \varphi_n \left( \mathbf{r} \right) \f$ is the modal basis obtained by
 * \f[ \varphi_n \left( \mathbf{r} \right) = \sqrt{2}P_i(a)P_j^{(2i+1, 0)}(b)(1-b)^i \f]
 *
 * where the modal basis is sorted as
 * \f[ n(i,j) = \left\{ \begin{array}{lllll}
 * (0,0) & (0,1) & \cdots & (0,N-1) & (0, N) \cr
 * (1,0) & (1,1) & \cdots & (1,N-1) \cr
 * \cdots \cr
 * (0,N) \end{array} \right\} \f]
 *
 *
 * @param [int]      N order
 * @param [int]      Nr number of coordinate
 * @param [doule]    r[Nr]
 * @param [doule]    s[Nr]
 *
 * @return
 * name     | type     | description of value
 * -------- |----------|----------------------
 * V  | double[Np][Np] | vandermonde matrix
 *
 */
void GetTriV(int N, int Nr, double *r, double *s, double **V){
    int i,j,k,sk=0;

    double *temp = Vector_create(Nr);
    double *a= Vector_create(Nr);
    double *b= Vector_create(Nr);

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

    Vector_free(a);
    Vector_free(b);
    Vector_free(temp);
}


/**
 * @brief
 * Get the nature coordinate of triangle element
 *
 * @details
 *
 * @param [int] N polynomial order
 *
 * @return
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

    L1 = Vector_create(Np);
    L2 = Vector_create(Np);
    L3 = Vector_create(Np);
    dL = Vector_create(Np);
    x  = Vector_create(Np);
    y  = Vector_create(Np);

    warpf1 = Vector_create(Np);

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
    Vector_free(L1);
    Vector_free(L2);
    Vector_free(L3);
    Vector_free(x);
    Vector_free(y);
    Vector_free(warpf1);
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
 * @param [int]     N order of degree
 * @param [double]  r[Nr] input points
 * @param [int]     Nr number of input points
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
 * @brief
 * Build the index matrix of nodes on faces
 * @details
 * Three faces of the standard triangle element is
 * \f[ s=-1, \quad r+s=0, \quad r=-1 \f]
 *
 * @param [int] N order
 *
 * @return
 * name     | type     | description of value
 * -------- |----------|----------------------
 * Fmask   | int[Nfaces][Nfp] |
 *
 */
void BuildTriFmask(int N, int **Fmask){
    int Nfp = N+1;
    int *temp = IntVector_create(Nfp);
    int *nrp = IntVector_create(Nfp); /* # of nodes from s=-1 to s=1 */
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

    IntVector_free(temp);
    IntVector_free(nrp);
}


/**
 * @brief
 *
 * @details
 * transfer coordinate (x,y) on equilateral triangle to natural coordinate (r,s)
 * on right triangle
 *
 * @param [int]     Np number of points
 * @param [double]  x[Np] coordinate of equilateral triangle
 * @param [double]  y[Np] coordinate of equilateral triangle
 *
 * @return
 *
 * |name     | type     | description of value |
 * | --- | --- | --- |
 * |r | double[Np] | coordinate |
 * |s | double[Np] | coordinate |
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
 * @details
 *
 * @param [in] Np number of points
 * @param [in] r coordinate
 * @param [in] s coordinate
 *
 * @return
 *
 * | name | type | description of value
 * | --- | --- | ---
 * |a | double[Np] | collapse coordinate
 * |b | double[Np] | collapse coordinate
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