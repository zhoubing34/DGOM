/**
 * @file
 * Standard quadrilateral element
 *
 * @brief
 * Functions related with standard quadrilateral elements
 *
 * @author
 * li12242, Tianjin University, li12242@tju.edu.cn
 */

#include "sc_stdcell.h"
#include "Polylib/polylib.h"

/* build the nodes index matrix on each faces */
int** sc_fmask_quad(stdCell *quad);

/* get the nature coordinate of interpolation nodes in standard quadrilateral element */
void sc_coord_quad(stdCell *quad);

/* transform the index of orthogonal function to [ti,tj] for quadrilateral elements */
void sc_transInd_quad(int N, int ind, int *ti, int *tj);

/* get orthogonal function value at interpolation nodes */
void sc_orthogFunc_quad(stdCell *quad, int ind, double *func);

/* calculate the value of derivative function at interpolation points */
void sc_deriOrthogFunc_quad(stdCell *quad, int ind, double *dr, double *ds);

/**
 * @brief create standard quadrilateral element
 * @param[in] N polynomial order
 * @return quad standard quadrilateral element
 */
stdCell* sc_createQuad(int N){
    stdCell *quad = (stdCell *) calloc(1, sizeof(stdCell));

    /* cell type */
    quad->type = QUADRIL;

    const int Np = (N+1)*(N+1);
    const int Nfaces = 4;
    const int Nfp = N+1;
    const int Nv = 4;
    /* basic info */
    quad->N      = N;
    quad->Np     = Np;
    quad->Nfaces = Nfaces;
    quad->Nfp    = Nfp;
    quad->Nv     = Nv;

    /* nodes at faces, Fmask */
    quad->Fmask = sc_fmask_quad(quad);

    /* coordinate, r and s */
    sc_coord_quad(quad);

    /* Vandermonde matrix, V */
    quad->V = sc_VandMatrix2d(quad, sc_orthogFunc_quad);

    /* mass matrix, M */
    quad->M = sc_massMatrix(quad);

    /* Derivative Matrix, Dr and Ds */
    sc_deriMatrix2d(quad, sc_deriOrthogFunc_quad);

    /* suface LIFT matrix, LIFT */
    quad->LIFT = sc_liftMatrix2d(quad);

    /* integration coefficients, ws and wv */
    sc_GaussQuadrature2d(quad);

    /* float version */
    size_t sz = Np*Nfp*Nfaces*sizeof(real);
    quad->f_LIFT = (real *) malloc(sz);
    sz = Np*Np*sizeof(real);
    quad->f_Dr = (real*) malloc(sz);
    quad->f_Ds = (real*) malloc(sz);

    int sk = 0, n, m;
    for(n=0;n<Np;++n){
        for(m=0;m<Nfp*Nfaces;++m){
            quad->f_LIFT[sk++] = (real) quad->LIFT[n][m];
        }
    }

    sk = 0;
    for(n=0;n<Np;++n){
        for(m=0;m<Np;++m){
            quad->f_Dr[sk] = (real) quad->Dr[n][m];
            quad->f_Ds[sk] = (real) quad->Ds[n][m];
            ++sk;
        }
    }

    return quad;
}

/**
 * @brief transform the index of orthogonal function to [ti,tj] for quadrilateral elements
 * @details the index is arranged as
 * \f[ [i,j] = \left\{ \begin{array}{lll}
 * (0,0) & \cdots & (0, N) \cr
 * (1,0) & \cdots & (0, N) \cr
 * \cdots \cr
 * (1,0) & \cdots & (0, N) \end{array} \right\} \f]
 *
 * @param [in] N polynomial order
 * @param [in] ind orthogonal polynomial index
 * @param [out] ti
 * @param [out] tj
 */
void sc_transInd_quad(int N, int ind, int *ti, int *tj){
    int i,j,sk=0;
    for(i=0;i<N+1;i++){
        for(j=0;j<N+1;j++){
            if(sk == ind){
                *ti = i;
                *tj = j;
            }
            sk++;
        }
    }
}

/**
 * @brief calculate the value of derivative function at interpolation points.
 * @details
 * @param [in] quad standard quadrilateral element
 * @param [in] ind index of orthogonal function
 * @param [out] dr derivative function with r
 * @param [out] ds derivative function with s
 * @note
 * precondition: dr and ds should be allocated before calling @ref sc_deriOrthogFunc_quad
 */
void sc_deriOrthogFunc_quad(stdCell *quad, int ind, double *dr, double *ds){
    const int N = quad->N;
    const int Np = quad->Np;

    double *r = quad->r;
    double *s = quad->s;
    double temp[Np];
    /* transform to the index [ti,tj] */
    int ti,tj;
    sc_transInd_quad(N, ind, &ti, &tj);

    /* Vr */
    GradjacobiP(Np, r, temp, ti, 0, 0);
    jacobiP(Np, s, dr, tj, 0, 0);
    int i;
    for(i=0;i<Np;i++)
        dr[i] *= temp[i];

    /* Vs */
    jacobiP(Np, r, temp, ti, 0.0, 0.0);
    GradjacobiP(Np, s, ds, tj, 0.0, 0.0);
    for(i=0;i<Np;i++)
        ds[i] *= temp[i];
}

/**
 * @brief get orthogonal function value at interpolation nodes
 * @details
 * The orthogonal function on quadrilateral is derived by the
 * 1d orthogonal function
 *
 * @param [in] cell quadrilateral cell
 * @param [in] ind index of orthogonal function
 * @param [out] func value of orthogonal function
 */
void sc_orthogFunc_quad(stdCell *quad, int ind, double *func){
    const int N = quad->N;
    const int Np = quad->Np;

    double temp[Np];
    /* transform to the index [ti,tj] */
    int ti,tj;
    sc_transInd_quad(N, ind, &ti, &tj);

    double *r = quad->r;
    double *s = quad->s;

    jacobiP(Np, r, temp, ti, 0.0, 0.0);
    jacobiP(Np, s, func, tj, 0.0, 0.0);

    int i;
    for(i=0;i<Np;i++)
        func[i] *= temp[i];
}

/**
 * @brief get the nature coordinate of interpolation nodes in standard quadrilateral element.
 * @param [int,out] quad standard quadrilateral element
 * @note
 * the nodes is arranged along r coordinate first
 */
void sc_coord_quad(stdCell *quad){
    const int Nfp = quad->Nfp;
    const int Np = quad->Np;
    double t[Nfp],w[Nfp];

    quad->r = Vector_create(Np);
    quad->s = Vector_create(Np);

    /* get Gauss-Lobatto-Jacobi zeros and weights */
    zwglj(t, w, Nfp, 0, 0);

    int i,j,sk=0;
    for(i=0;i<Nfp;i++){
        for(j=0;j<Nfp;j++){
            quad->r[sk]   = t[j];
            quad->s[sk++] = t[i];
        }
    }
}


/**
 * @brief build the nodes index matrix on each faces
 * @details
 * Four faces of the standard quadrilateral element is
 * \f[ s=-1, \quad r=1, \quad s=1, \quad r=-1 \f]
 *
 * @param [int] quad standard quadralteral element
 * @return Fmask[Nfaces][Nfp]
 */
int** sc_fmask_quad(stdCell *quad){
    const int Nfp = quad->Nfp;
    const int Nfaces = quad->Nfaces;

    int **Fmask = IntMatrix_create(Nfaces, Nfp);
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

    return Fmask;
}


/**
 * @brief
 * Project the quadrilateral vertex value to interpolation nodes.
 *
 * @param[in] vertVal value of vertex
 * @param[in] nodeVal value of nodes
 *
 */
void sc_vertProj_quad(stdCell *cell, double *vertVal, double *nodeVal){
    int i;
    double ri, si;
    for (i=0;i<cell->Np;++i) {
        ri = cell->r[i];
        si = cell->s[i];
        nodeVal[i] = 0.25 * (vertVal[0] * (1. - ri) * (1. - si) + vertVal[1] * (1. + ri) * (1. - si)
                             + vertVal[2] * (1. + ri) * (1. + si) + vertVal[3] * (1. - ri)*(1. + si));
    }
}