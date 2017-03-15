//
// Created by li12242 on 17/3/14.
//

#include "dg_cell_volume.h"

static double** dg_cell_volume_vand_matrix(int N, int Np, double *r, double *s, double *t, Orthogonal_Func);
static double** dg_cell_volume_mass_matrix(int Np, double **V);
static void dg_cell_deri_matrix(int N, int Np, double *r, double *s, double *t, double **V,
                                Derivative_Orthogonal_Func, double ***Dr, double ***Ds, double ***Dt);

/**
 * @brief
 * create dg_cell_volume structure.
 * @param cell
 * @param creator
 * @return
 */
dg_cell_volume *dg_cell_volume_create(dg_cell *cell, const dg_cell_creator *creator){
    const int N = dg_cell_N(cell);
    dg_cell_volume *volume = (dg_cell_volume *)calloc(1, sizeof(dg_cell_volume));
    /* get coordinate */
    int Np;
    double *r, *s, *t;
    creator->set_node(cell, &Np, &r, &s, &t);
    volume->Np = Np;
    volume->r = r;
    volume->s = s;
    volume->t = t;
    /* Vandermonde matrix */
    volume->V = dg_cell_volume_vand_matrix(N, Np, r, s, t, creator->orthogonal_func);
    /* mass matrix */
    volume->M = dg_cell_volume_mass_matrix(Np, volume->V);
    /* derivative matrix */
    double **Dr, **Ds, **Dt;
    dg_cell_deri_matrix(N, Np, r, s, t, volume->V, creator->deri_orthgonal_func, &Dr, &Ds, &Dt);
    volume->Dr = Dr;
    volume->Ds = Ds;
    volume->Dt = Dt;
    /* volume quadrature weights */
    volume->w = vector_double_create(Np);
    int i,j;
    for(i=0;i<Np;i++){
        for(j=0;j<Np;j++){
            volume->w[j] += volume->M[i][j];
        }
    }
    return volume;
}

void dg_cell_volume_free(dg_cell_volume *volume){
    vector_double_free(volume->r);
    vector_double_free(volume->s);
    vector_double_free(volume->t);

    matrix_double_free(volume->V);
    matrix_double_free(volume->Dr);
    matrix_double_free(volume->Ds);
    matrix_double_free(volume->Dt);
    free(volume);
    return;
}

/**
 * @brief calculate the Vandermonde matrix
 * @details
 * @param [in] cell standard cell
 * @param [in] orthfunc orthogonal function handle
 * @return Vandermonde matrix
 */
static double** dg_cell_volume_vand_matrix(int N, int Np,
                                           double *r, double *s, double *t,
                                           Orthogonal_Func orthogonal_func) {
    // allocation
    double **V = matrix_double_create(Np, Np);
    // assignment
    int dim1,dim2;
    double temp[Np];
    for(dim2=0;dim2<Np;dim2++){
        orthogonal_func(N, dim2, Np, r, s, t, temp);
        for(dim1=0;dim1<Np;dim1++){
            V[dim1][dim2] = temp[dim1];
        }
    }
    return V;
}

/**
 * @brief calculate the mass matrix
 * @details
 * The mass matrix is calculated with
 * \f[ \mathbf{M} = (\mathbf{V}^T)^{-1} \cdot \mathbf{V}^{-1} \f]
 *
 * @param [in] cell standard cell
 * @param [out] M Vandermonde matrix
 * @return mass matrix
 */
static double** dg_cell_volume_mass_matrix(int Np, double **V){

    double inv[Np*Np], invt[Np*Np];
    int i,j,sk=0;
    // allocation
    double **M = matrix_double_create(Np, Np);
    // copy matrix to a vector
    for(i=0;i<Np*Np;i++) {inv[sk++] = V[0][i];}

    matrix_inverse(inv, Np);
    // transpose of invV
    for(i=0;i<Np;i++){
        for(j=0;j<Np;j++){
            invt[j*Np + i] = inv[i*Np + j];
        }
    }
    matrix_multiply(Np, Np, Np, invt, inv, M[0]);
    return M;
}

/**
 * @brief get the derivative Vandermonde matrix
 * @param [in] cell standard element
 * @param [in] derorthfunc derivative of orthogonal function
 * @param [out] Vr the derivative of Vandermonde matrix on r coordinate
 * @param [out] Vs the derivative of Vandermonde matrix on s coordinate
 * @param [out] Vt the derivative of Vandermonde matrix on t coordinate
 * @note
 * Vr, Vs and Vt should be allocated before calling.
 */
static void dg_deri_vand_matrix(
        int N, int Np, double *r, double *s, double *t,
        Derivative_Orthogonal_Func deri_orthgonal_func,
        double **Vr, double **Vs, double **Vt) {
    double dr[Np],ds[Np],dt[Np];

    int dim1, dim2;
    for(dim2=0;dim2<Np;dim2++){
        deri_orthgonal_func(N, dim2, Np, r, s, t, dr, ds, dt);
        for(dim1=0;dim1<Np;dim1++){
            Vr[dim1][dim2] = dr[dim1];
            Vs[dim1][dim2] = ds[dim1];
            Vt[dim1][dim2] = dt[dim1];
        }
    }
    return;
}


/**
 * @brief
 * get the gradient matrix Dr and Ds of Lagrange basis at (r,s) at order N
 * @detail
 * The Gradient matrix \f$ \mathbf{Dr} \f$ and \f$ \mathbf{Ds} \f$ is obtained through
 * \f[ \mathbf{Dr} \cdot \mathbf{V} = \mathbf{Vr}, \quad \mathbf{Ds} \cdot \mathbf{V} = \mathbf{Vs} \f]
 * where
 * \f[ Dr_{(ij)} = \left. \frac{\partial l_j}{\partial r} \right|_{ \mathbf{r}_i },
 * \quad Ds_{(ij)} = \left. \frac{\partial l_j}{\partial s} \right|_{ \mathbf{r}_i } \f]
 *
 * @param [in,out] cell standard element
 * @param [in] derorthfunc derivative orthogonal function handle
 */
static void dg_cell_deri_matrix(
        int N, int Np, double *r, double *s, double *t, double **V,
        Derivative_Orthogonal_Func deri_orthgonal_func,
        double ***Dr, double ***Ds, double ***Dt) {
    // allocation
    double **dr = matrix_double_create(Np, Np);
    double **ds = matrix_double_create(Np, Np);
    double **dt = matrix_double_create(Np, Np);

    double **Vr = matrix_double_create(Np, Np);
    double **Vs = matrix_double_create(Np, Np);
    double **Vt = matrix_double_create(Np, Np);
    double inv[Np*Np];

    // inverse of Vandermonde matrix
    int i, sk=0; // copy matrix to a vector
    for(i=0;i<Np*Np;i++) {inv[sk++] = V[0][i];}

    matrix_inverse(inv, Np);
    // get derivative Vandermonde matrix
    dg_deri_vand_matrix(N, Np, r, s, t, deri_orthgonal_func, Vr, Vs, Vt);

    /* \f$ \mathbf{Dr} = \mathbf{Vr}*\mathbf{V}^{-1} \f$ */
    matrix_multiply(Np, Np, Np, Vr[0], inv, dr[0]);
    matrix_multiply(Np, Np, Np, Vs[0], inv, ds[0]);
    matrix_multiply(Np, Np, Np, Vt[0], inv, dt[0]);

    matrix_double_free(Vr);
    matrix_double_free(Vs);
    matrix_double_free(Vt);
    /* assignment */
    *Dr = dr;
    *Ds = ds;
    *Dt = dt;

    return;
}
