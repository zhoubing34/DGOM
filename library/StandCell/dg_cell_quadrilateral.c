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

#include "dg_cell.h"
#include "Polylib/polylib.h"

#define DEBUG 0
#if DEBUG
#include "Utility/unit_test.h"
#endif
/* transform the index of orthogonal function to [ti,tj] for quadrilateral elements */
void sc_transInd_quad(int N, int ind, int *ti, int *tj);

void dg_quad_info(dg_cell *cell, int N){
    /* cell type */
    cell->type = QUADRIL;

    const int Np = (N+1)*(N+1);
    const int Nfaces = 4;
    const int Nfp = N+1;
    const int Nv = 4;
    /* basic info */
    cell->N      = N;
    cell->Np     = Np;
    cell->Nfaces = Nfaces;
    cell->Nv     = Nv;
    cell->Nfptotal = Nfaces*Nfp;
    cell->Nfp    = (int *)calloc(Nfaces, sizeof(int));
    int f;
    for(f=0;f<Nfaces;f++){
        cell->Nfp[f] = Nfp;
    }

    return;
}

/**
 * @brief calculate the Gauss quadrature weights
 * @param [in] cell 2d standard elements
 */
void dg_quad_gauss_weight(dg_cell *cell){
    const int Nfp = dg_cell_N(cell)+1;
    const int Np = cell->Np;
    const int Nfaces = dg_cell_Nfaces(cell);
    double r[Nfp], ws[Nfp];
    // Gauss quadrature weights for face2d
    cell->ws = matrix_double_create(Nfaces, Nfp);
    zwglj(r, ws, Nfp, 0, 0);
    int i,j;
    for(i=0;i<Nfaces;i++){
        for(j=0;j<dg_cell_Nfp(cell, i);j++){
            cell->ws[i][j] = ws[j];
        }
    }

    // Gauss quadrature weights for volume
    cell->wv = vector_double_create(Np);
    for(i=0;i<Np;i++){
        for(j=0;j<Np;j++){
            cell->wv[j] += cell->M[i][j];
        }
    }
    return;
}
/**
 * @brief
 * generation of surface mass matrix
 * @param cell
 * @return surface mass matrix
 */
double ** dg_quad_surf_mass_matrix(dg_cell *cell){
    const int Np = dg_cell_Np(cell);
    const int Nfaces = dg_cell_Nfaces(cell);
    const int Nfp = dg_cell_N(cell)+1;

    double **Mes = matrix_double_create(Np, Nfaces * Nfp);
    int **Fmask = cell->Fmask;
    /* coefficients for faces */
    double r[Nfp], w[Nfp];
    double invt[Nfp*Nfp], inv[Nfp*Nfp], m[Nfp*Nfp];
    /* get mass matrix of line */
    zwglj(r, w, Nfp, 0, 0); /* get coordinate */
    int i,j;
    for(i=0;i<Nfp;i++){
        /* get vandermonde matrix of line */
        jacobiP(Nfp, r, w, i, 0, 0);
        for(j=0;j<Nfp;j++){
            inv[j*Nfp+i] = w[j];
        }
    }
    matrix_inverse(inv, Nfp);
    /* transform of vandermonde matrix */
    for(i=0;i<Nfp;i++){
        for(j=0;j<Nfp;j++)
            invt[j+Nfp*i] = inv[j*Nfp+i];
    }
    /* get M = inv(V)'*inv(V) */
    matrix_multiply(Nfp, Nfp, Nfp, invt, inv, m);

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
    return Mes;
}

void dg_quad_free(dg_cell *cell){

    matrix_int_free(cell->Fmask);
    vector_double_free(cell->r);
    vector_double_free(cell->s);
    /* vandermonde matrix */
    matrix_double_free(cell->V);
    /* mass matrix */
    matrix_double_free(cell->M);
    /* Derivative Matrix */
    matrix_double_free(cell->Dr);
    matrix_double_free(cell->Ds);
    matrix_double_free(cell->Dt);

    /* LIFT */
    matrix_double_free(cell->LIFT);

    /* Gauss quadrature */
    matrix_double_free(cell->ws);
    vector_double_free(cell->wv);

    /* float version */
    free(cell->f_LIFT);
    free(cell->f_Dr);
    free(cell->f_Ds);

    free(cell);
    return;
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
void dg_quad_deri_orthog_func(dg_cell *quad, int ind, double *dr, double *ds, double *dt){
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
    for(i=0;i<Np;i++) {dr[i] *= temp[i];}

    /* Vs */
    jacobiP(Np, r, temp, ti, 0.0, 0.0);
    GradjacobiP(Np, s, ds, tj, 0.0, 0.0);
    for(i=0;i<Np;i++) {ds[i] *= temp[i];}
    return;
}

/**
 * @brief
 * get orthogonal function value at interpolation nodes
 * @details
 * The orthogonal function on quadrilateral is derived by the
 * 1d orthogonal function
 *
 * @param [in] cell quadrilateral cell
 * @param [in] ind index of orthogonal function
 * @param [out] func value of orthogonal function
 */
void dg_quad_orthog_func(dg_cell *quad, int ind, double *func){
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
    for(i=0;i<Np;i++) {func[i] *= temp[i];}
    return;
}

/**
 * @brief
 * get the nature coordinate of interpolation nodes in standard quadrilateral element.
 * @param [int,out] quad standard quadrilateral element
 * @note
 * the nodes is arranged along r coordinate first
 */
void dg_quad_coord(dg_cell *quad){
    const int Np = dg_cell_Np(quad);
    const int Nfp = dg_cell_N(quad)+1;
    double t[Nfp],w[Nfp];

    quad->r = vector_double_create(Np);
    quad->s = vector_double_create(Np);

    /* get Gauss-Lobatto-Jacobi zeros and weights */
    zwglj(t, w, Nfp, 0, 0);
    int i,j,sk=0;
    for(i=0;i<Nfp;i++){
        for(j=0;j<Nfp;j++){
            quad->r[sk]   = t[j];
            quad->s[sk++] = t[i];
        }
    }
    return;
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
int** dg_quad_fmask(dg_cell *quad){
    const int Nfp = dg_cell_N(quad)+1;
    const int Nfaces = dg_cell_Nfaces(quad);
    const int Nfptotal = dg_cell_Nfptotal(quad);

    // allocation
    int **Fmask = (int **)calloc(Nfaces, sizeof(int *));
    Fmask[0] = calloc((size_t) Nfptotal, sizeof(int));
    int n;
    for(n=1;n<Nfaces;++n){
        Fmask[n] = Fmask[n-1]+ dg_cell_Nfp(quad, n-1);
    }
    // point index on each face
    int i, std, td;
    /* face2d 1, s=-1 */
    for(i=0;i<Nfp;i++)
        Fmask[0][i] = i;
    /* face2d 2, r=+1 */
    std = Nfp-1; /* start index */
    td  = Nfp;
    for(i=0;i<Nfp;i++) {
        Fmask[1][i] = std;
        std += td;
    }
    /* face2d 3, s=+1 */
    std = Nfp*Nfp - 1;
    td  = -1;
    for(i=0;i<Nfp;i++) {
        Fmask[2][i] = std;
        std += td;
    }
    /* face2d 4, r=-1 */
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
 * @param[in] vertVal value of vertex
 * @param[in] nodeVal value of nodes
 *
 */
void dg_quad_proj_vert2node(dg_cell *cell, double *vertVal, double *nodeVal){
    int i;
    for (i=0;i<cell->Np;++i) {
        double ri = cell->r[i];
        double si = cell->s[i];
        nodeVal[i] = 0.25 * (vertVal[0] * (1. - ri) * (1. - si)
                             + vertVal[1] * (1. + ri) * (1. - si)
                             + vertVal[2] * (1. + ri) * (1. + si)
                             + vertVal[3] * (1. - ri)*(1. + si));
    }
    return;
}