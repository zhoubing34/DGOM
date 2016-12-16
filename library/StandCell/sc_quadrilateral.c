/**
 * @file
 * Quadrilateral.c
 *
 * @brief
 * Functions related with standard quadrilateral elements
 *
 * @author
 * li12242, Tianjin University, li12242@tju.edu.cn
 */

#include "sc_stdcell.h"

/* Private functions */
void BuildQuadFmask(int N, int **Fmask);
void GetQuadCoord(int N, double *r, double *s);
void GetQuadV(int N, int Np, double *r, double *s, double **V);
void GetQuadM(int Np, double **V, double **M);
void GetQuadDeriV2d(int N, int Np, double *r, double *s, double **Vr, double **Vs);
void GetQuadDeriM(int N, int Np, double *r, double *s, double **V, double **Dr, double **Ds);

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
 * tri | stdCell* |
 *
 */
stdCell* sc_createQuad(int N){
    stdCell *quad = (stdCell *) calloc(1, sizeof(stdCell));

    int Np = (N+1)*(N+1);
    /* basic info */
    quad->N      = N;
    quad->Np     = Np;
    quad->Nfaces = 4;
    quad->Nfp    = N+1;
    quad->Nv     = 4;

    /* nodes at faces, Fmask */
    quad->Fmask = IntMatrix_create(quad->Nfaces, quad->Nfp);
    BuildQuadFmask(quad->N, quad->Fmask);

    /* coordinate */
    quad->r = (double *) calloc(Np, sizeof(double));
    quad->s = (double *) calloc(Np, sizeof(double));
    GetQuadCoord(N, quad->r, quad->s);

    /* vandermonde matrix */
    quad->V = Matrix_create(Np, Np);
    GetQuadV(N, Np, quad->r, quad->s, quad->V);

    /* mass matrix */
    quad->M = Matrix_create(Np, Np);
    GetQuadM(Np, quad->V, quad->M);

    /* Derivative Matrix, Dr and Ds */
    quad->Dr = Matrix_create(quad->Np, quad->Np);
    quad->Ds = Matrix_create(quad->Np, quad->Np);
    GetQuadDeriM(N, Np, quad->r, quad->s, quad->V, quad->Dr, quad->Ds);

    /* suface LIFT matrix, LIFT */
    double **Mes = Matrix_create(quad->Np, quad->Nfaces * quad->Nfp);
    quad->LIFT = Matrix_create(quad->Np, quad->Nfaces * quad->Nfp);
    GetSurfLinM(quad->N, quad->Nfaces, quad->Fmask, Mes);
    GetLIFT2d(quad, Mes, quad->LIFT);

//    PrintMatrix_test("Mes", Mes, quad->Np, quad->Nfaces*quad->Nfp);
    Matrix_free(Mes);

    /* integration coeff, ws and wv */
    /* integration coeff, ws and wv */
    quad->ws = Vector_create(quad->Nfp);
    double *r = Vector_create(quad->Nfp);
    zwglj(r, quad->ws, quad->Nfp, 0, 0);
    Vector_free(r);

    quad->wv = Vector_create(quad->Np);
    int i,j;
    for(i=0;i<quad->Np;i++){
        for(j=0;j<quad->Np;j++){
            quad->wv[j] += quad->M[i][j];
        }
    }

    /* float version */
    size_t sz = quad->Np*(quad->Nfp)*(quad->Nfaces)*sizeof(real);
    quad->f_LIFT = (real *) malloc(sz);
    sz = quad->Np*quad->Np*sizeof(real);
    quad->f_Dr = (real*) malloc(sz);
    quad->f_Ds = (real*) malloc(sz);

    int sk = 0, n, m;
    for(n=0;n<quad->Np;++n){
        for(m=0;m<quad->Nfp*quad->Nfaces;++m){
            quad->f_LIFT[sk++] = (real) quad->LIFT[n][m];
        }
    }

    sk = 0;
    for(n=0;n<quad->Np;++n){
        for(m=0;m<quad->Np;++m){
            quad->f_Dr[sk] = (real) quad->Dr[n][m];
            quad->f_Ds[sk] = (real) quad->Ds[n][m];
            ++sk;
        }
    }


    return quad;
}


void GetQuadDeriM(int N, int Np, double *r, double *s, double **V, double **Dr, double **Ds){
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
    Matrix_Inverse(temp, Np);

    /* get the gradient of the modal basis */
    GetQuadDeriV2d(N, Np, r, s, Vr, Vs);

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
    Matrix_Multiply(n, n, n, vtemp, temp, dtemp);
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
    Matrix_Multiply(n, n, n, vtemp, temp, dtemp);
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

void GetQuadDeriV2d(int N, int Np, double *r, double *s, double **Vr, double **Vs){
    double *h1 = Vector_create(Np);
    double *h2 = Vector_create(Np);

    int i,j,k,sk;
    for(i=0;i<(N+1);i++){
        GradjacobiP(Np, r, h1, i, 0, 0);
        for(j=0;j<(N+1);j++){
            jacobiP(Np, s, h2, j, 0, 0);
            sk = i*(N+1)+j;  /* column index */
            for(k=0;k<Np;k++){
                Vr[k][sk] = h1[k]*h2[k];
            }
        }
    }

    for(i=0;i<(N+1);i++){
        jacobiP(Np, r, h1, i, 0.0, 0.0);
        for(j=0;j<(N+1);j++){
            GradjacobiP(Np, s, h2, j, 0.0, 0.0);
            sk = i*(N+1)+j;  /* column index */
            for(k=0;k<Np;k++){
                Vs[k][sk] = h1[k]*h2[k];
            }
        }
    }

    Vector_free(h1);
    Vector_free(h2);
}


void GetQuadM(int Np, double **V, double **M){
    double *temp = Vector_create(Np * Np);
    double *invt = Vector_create(Np * Np);
    double *Mv   = Vector_create(Np * Np);
    const unsigned n = (unsigned)Np;

    int i,j,sk=0;

    for(i=0;i<Np;i++){
        for(j=0;j<Np;j++){
            /* row counts first */
            temp[sk++] = V[i][j];
        }
    }

    Matrix_Inverse(temp, Np);

    for(i=0;i<Np;i++){
        for(j=0;j<Np;j++){
            invt[j*Np + i] = temp[i*Np + j];
        }
    }

    Matrix_Multiply(n, n, n, invt, temp, Mv);

    for(i=0;i<Np;i++){
        for(j=0;j<Np;j++){
            M[i][j] = Mv[i*Np + j];
        }
    }
}

/**
 * @brief
 * Generate the Vandermonde matrix of quadrilateral element
 *
 * @details
 * Vandermonde matrix \f$ \mathbf{V} \f$ is derivated by
 * \f[ \mathbf{V}_{mn} = \varphi_n \left( \mathbf{r}_m \right) \f]
 * where \f$ \varphi_n \left( \mathbf{r} \right) \f$ is the modal basis obtained by
 * \f[ \varphi_n \left( \mathbf{r} \right) = P_i(r)P_j(s) \f]
 *
 * where the modal basis is arranged to count \f[ j \f] first.
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
void GetQuadV(int N, int Np, double *r, double *s, double **V){
    double *h1 = Vector_create(Np);
    double *h2 = Vector_create(Np);

    int i,j,k,sk=0;
    for(i=0;i<(N+1);i++){
        jacobiP(Np, r, h1, i, 0.0, 0.0);
        for(j=0;j<(N+1);j++){
            jacobiP(Np, s, h2, j, 0.0, 0.0);
            sk = i*(N+1)+j;  /* column index */
            for(k=0;k<Np;k++){
                V[k][sk] = h1[k]*h2[k];
            }
        }
    }
}

/**
 * @brief
 * Get the nature coordinate of quadrilateral element
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
 * r and s shuold be allocated before calling GetQuadCoord
 *
 */
void GetQuadCoord(int N, double *r, double *s){
    int Np = (N+1);
    double *t = Vector_create(Np);
    double *w = Vector_create(Np);

    /* get Gauss-Lobatto-Jacobi zeros and weights */
    zwglj(t, w, Np, 0, 0);

    int i,j,sk=0;
    for(i=0;i<Np;i++){
        for(j=0;j<Np;j++){
            r[sk]   = t[j];
            s[sk++] = t[i];
        }
    }

    Vector_free(t);
    Vector_free(w);
}


/**
 * @brief
 * Build the index matrix of nodes on faces
 * @details
 * Four faces of the standard quadrilateral element is
 * \f[ s=-1, \quad r=1, \quad s=1, \quad r=-1 \f]
 *
 * @param [int] N order
 *
 * @return
 * name     | type     | description of value
 * -------- |----------|----------------------
 * Fmask   | int[Nfaces][Nfp] |
 *
 */
void BuildQuadFmask(int N, int **Fmask){
    int Nfp = N+1;
    int i, std, td;

    /* face 1, s=-1 */
    for(i=0;i<Nfp;i++)
        Fmask[0][i] = i;

    /* face 2, r=+1 */
    std = Nfp-1; /* start index */
    td  = Nfp;
    for(i=0;i<Nfp;i++) {
        Fmask[1][i] = std;
        std += td;
    }

    /* face 3, s=+1 */
    std = Nfp*Nfp - 1;
    td  = -1;
    for(i=0;i<Nfp;i++) {
        Fmask[2][i] = std;
        std += td;
    }

    /* face 4, r=-1 */
    std = (Nfp - 1)*Nfp;
    td  = -Nfp;
    for(i=0;i<Nfp;i++) {
        Fmask[3][i] = std;
        std += td;
    }
}