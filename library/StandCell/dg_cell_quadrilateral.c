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

#define DEBUG 1
#if DEBUG
#include "Utility/unit_test.h"
#endif
/* transform the index of orthogonal function to [ti,tj] for quadrilateral elements */
static void dg_quad_transInd(int N, int ind, int *ti, int *tj);


dg_cell_info* dg_cell_quad_info(int N){
    dg_cell_info *info = (dg_cell_info *)calloc(1, sizeof(dg_cell_info));
    const int Nv = 4;
    const int Nfaces = 4;
    const int Nfv = 2;
    info->N = N;
    info->Nv = Nv;
    info->Nfaces = Nfaces;
    info->type = QUADRIL;
    info->face_type = calloc(Nfaces, sizeof(dg_cell_type));
    info->FToV = matrix_int_create(Nfaces, Nfv);
    info->vr = (double *)calloc(Nv, sizeof(double));
    info->vs = (double *)calloc(Nv, sizeof(double));
    info->vt = (double *)calloc(Nv, sizeof(double));
    /* faces */
    info->face_type[0] = LINE; info->FToV[0][0] = 0; info->FToV[0][1] = 1;
    info->face_type[1] = LINE; info->FToV[1][0] = 1; info->FToV[1][1] = 2;
    info->face_type[2] = LINE; info->FToV[2][0] = 2; info->FToV[2][1] = 3;
    info->face_type[3] = LINE; info->FToV[3][0] = 3; info->FToV[3][1] = 0;
    /* vertex */
    info->vr[0] = -1; info->vr[1] =  1; info->vr[2] = 1; info->vr[3] = -1;
    info->vs[0] = -1; info->vs[1] = -1; info->vs[2] = 1; info->vs[3] =  1;
    return info;
}

/**
 * @brief
 * get the nature coordinate of interpolation nodes in standard quadrilateral element.
 * @param [int,out] quad standard quadrilateral element
 * @note
 * the nodes is arranged along r coordinate first
 */
void dg_cell_quad_set_nood(dg_cell *cell, int *Np, double **r, double **s, double **t){

    const int N = dg_cell_N(cell);
    const int Npt = (N+1)*(N+1);
    const int Nfp = N+1;
    double x[Nfp],w[Nfp];
    double *rt = vector_double_create(Npt);
    double *st = vector_double_create(Npt);

    /* get Gauss-Lobatto-Jacobi zeros and weights */
    zwglj(x, w, Nfp, 0, 0);
    int i,j,sk=0;
    for(i=0;i<Nfp;i++){
        for(j=0;j<Nfp;j++){
            rt[sk] = x[j];
            st[sk++] = x[i];
        }
    }
    // assignment
    *Np = Npt;
    *r = rt;
    *s = st;
    *t = vector_double_create(Npt);
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
static void dg_quad_transInd(int N, int ind, int *ti, int *tj){
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
void dg_cell_quad_orthog_func(int N, int ind, int Np, double *r, double *s, double *t, double *fun){

    double temp[Np];
    /* transform to the index [ti,tj] */
    int ti,tj;
    dg_quad_transInd(N, ind, &ti, &tj);

    jacobiP(Np, r, temp, ti, 0.0, 0.0);
    jacobiP(Np, s, fun, tj, 0.0, 0.0);

    int i;
    for(i=0;i<Np;i++) {fun[i] *= temp[i];}
    return;
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
void dg_cell_quad_deri_orthog_func(int N, int ind, int Np, double *r, double *s, double *t,
                                   double *dr, double *ds, double *dt){

    double temp[Np];
    /* transform to the index [ti,tj] */
    int ti,tj;
    dg_quad_transInd(N, ind, &ti, &tj);

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
 * Project the quadrilateral vertex value to interpolation nodes.
 * @param[in] vertVal value of vertex
 * @param[in] nodeVal value of nodes
 *
 */
void dg_cell_quad_proj(dg_cell *cell, double *vertVal, double *nodeVal){
    int i;
    double *r = dg_cell_r(cell);
    double *s = dg_cell_s(cell);
    for (i=0;i<dg_cell_Np(cell);++i) {
        double ri = r[i];
        double si = s[i];
        nodeVal[i] = 0.25 * (vertVal[0] * (1. - ri) * (1. - si)
                             + vertVal[1] * (1. + ri) * (1. - si)
                             + vertVal[2] * (1. + ri) * (1. + si)
                             + vertVal[3] * (1. - ri)*(1. + si));
    }
    return;
}