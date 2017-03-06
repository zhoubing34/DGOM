#include "dg_cell.h"
#include "dg_cell_quadrilateral.h"
#include "dg_cell_triangle.h"

typedef struct dg_cell_creator{
    void (*set_info)(dg_cell *cell, int N); ///< infomation property
    void (*set_coor)(dg_cell *cell); ///< coordinate function
    int   **(*fmask)(dg_cell *cell); ///< return fmask matrix
    void (*orthogonal_func)(dg_cell *cell, int ind, double *fun); ///< orthogonal function
    void (*deri_orthgonal_func)(dg_cell *cell, int ind,
                                double *dr, double *ds, double *dt); ///< derivative of orthogonal func
    double **(*surf_mass_matrix)(dg_cell *cell); ///< return surf mass matrix
    void (*gauss_weight)(dg_cell *cell); ///< gauss quadrature weights
    void (*free_func)(dg_cell *cell); ///< free strucutre
    void (*proj_vert2node)(dg_cell *cell, double *vert_val, double *node_val); ///< project vertex value to node value
} dg_cell_creator;


static const dg_cell_creator tri_creator = {
        dg_tri_info,
        dg_tri_coord,
        dg_tri_fmask,
        dg_tri_orthog_func,
        dg_tri_deriorthog_func,
        dg_tri_surf_mass_matrix,
        dg_tri_gauss_weight,
        dg_tri_free,
        dg_tri_proj_vert2node,
};

static const dg_cell_creator quad_creator = {
        dg_quad_info,
        dg_quad_coord,
        dg_quad_fmask,
        dg_quad_orthog_func,
        dg_quad_deri_orthog_func,
        dg_quad_surf_mass_matrix,
        dg_quad_gauss_weight,
        dg_quad_free,
        dg_quad_proj_vert2node,
};

/// declaration of local functions
static double** dg_cell_lift_matrix(dg_cell *cell, double **(*surf_mass_matrix)(dg_cell *));
static double** dg_cell_vand_matrix(dg_cell *cell, void (*orthfunc)(dg_cell *, int ind, double *func));
static double** dg_cell_mass_matrix(dg_cell *cell);
static void dg_cell_deri_matrix(dg_cell *cell,
                                void (*derorthfunc)(dg_cell *, int ind, double *dr, double *ds, double *dt));
static void dg_cell_d2f(dg_cell *cell);

/***
 * @brief
 * @param N
 * @param type
 * @return
 */
dg_cell *dg_cell_creat(int N, dg_cell_type type){
    dg_cell *cell = (dg_cell *) calloc(1, sizeof(dg_cell));
    const dg_cell_creator *creator;
    switch (type){
        case TRIANGLE:
            creator = &tri_creator; break;
        case QUADRIL:
            creator = &quad_creator; break;
        default:
            fprintf(stderr, "%s (%d): Unknown cell type %d\n", __FUNCTION__, __LINE__, type);
            exit(-1);
    }
    creator->set_info(cell, N); // basic information
    creator->set_coor(cell); // coordinate
    cell->Fmask = creator->fmask(cell); // fmask
    cell->V = dg_cell_vand_matrix(cell, creator->orthogonal_func);
    cell->M = dg_cell_mass_matrix(cell);
    dg_cell_deri_matrix(cell, creator->deri_orthgonal_func);
    cell->LIFT = dg_cell_lift_matrix(cell, creator->surf_mass_matrix);
    dg_cell_d2f(cell);
    creator->gauss_weight(cell);

    cell->free_func = creator->free_func;
    cell->proj_vert2node = creator->proj_vert2node;
    return cell;
}

/**
 * @brief
 * copy double type to user specific precision.
 * @param [in,out] cell dg_cell structure
 */
static void dg_cell_d2f(dg_cell *cell){
    const int Np = cell->Np;
    const int Nfaces = cell->Nfaces;
    const int Nfp = cell->Nfp;
    /* float version */
    size_t sz = (size_t) Np*Nfp*Nfaces;
    cell->f_LIFT = (dg_real *) calloc(sz, sizeof(dg_real));

    int sk = 0, n, m;
    for(n=0;n<Np;++n){
        for(m=0;m<Nfp*Nfaces;++m){
            cell->f_LIFT[sk++] = (dg_real) cell->LIFT[n][m];
        }
    }

    sz = (size_t) Np*Np;
    cell->f_Dr = (dg_real*) calloc(sz, sizeof(dg_real));
    cell->f_Ds = (dg_real*) calloc(sz, sizeof(dg_real));
    sk = 0;
    for(n=0;n<Np;++n){
        for(m=0;m<Np;++m){
            cell->f_Dr[sk  ] = (dg_real) cell->Dr[n][m];
            cell->f_Ds[sk++] = (dg_real) cell->Ds[n][m];
        }
    }
    return;
}

/**
 * @brief
 * free stdCell object
 * @param[in] cell dg_cell structure
 */
void dg_cell_free(dg_cell *cell){
    cell->free_func(cell);
}

/**
 * @brief Project the vertex value to interpolation nodes.
 * @param[in] vertVal value on vertex
 * @param[in] nodeVal value on nodes
 */
void dg_cell_proj_vert2node(dg_cell *cell, double *vertVal, double *nodeVal){
    cell->proj_vert2node(cell, vertVal, nodeVal);
    return;
}

/**
 * @brief get the derivative Vandermonde matrix
 * @param [in] cell standard element
 * @param [in] derorthfunc derivative of orthogonal function
 * @param [out] Vr the derivative of Vandermonde matrix on r coordinate
 * @param [out] Vs the derivative of Vandermonde matrix on s coordinate
 * @param [out] Vt the derivative of Vandermonde matrix on t coordinate
 * @note
 * Vr, Vs and Vt should be allocated before calling sc_deriVandMatrix2d
 */
static void dg_deri_vand_matrix(dg_cell *cell,
                                void (*derorthfunc)(dg_cell *, int ind, double *dr, double *ds, double *dt),
                                double **Vr, double **Vs, double **Vt)
{
    const int Np = cell->Np;
    double dr[Np],ds[Np],dt[Np];

    int dim1, dim2;
    for(dim2=0;dim2<Np;dim2++){
        derorthfunc(cell, dim2, dr, ds, dt);
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
static void dg_cell_deri_matrix(dg_cell *cell,
                                void (*derorthfunc)(dg_cell *, int ind, double *dr, double *ds, double *dt))
{
    const int Np = cell->Np;
    // allocation
    cell->Dr = matrix_double_create(Np, Np);
    cell->Ds = matrix_double_create(Np, Np);
    cell->Dt = matrix_double_create(Np, Np);

    double **V = cell->V;
    double **Vr = matrix_double_create(Np, Np);
    double **Vs = matrix_double_create(Np, Np);
    double **Vt = matrix_double_create(Np, Np);
    double inv[Np*Np];

    // inverse of Vandermonde matrix
    int i, sk=0; // copy matrix to a vector
    for(i=0;i<Np*Np;i++) {inv[sk++] = V[0][i];}

    matrix_inverse(inv, Np);
    // get derivative Vandermonde matrix
    dg_deri_vand_matrix(cell, derorthfunc, Vr, Vs, Vt);

    /* \f$ \mathbf{Dr} = \mathbf{Vr}*\mathbf{V}^{-1} \f$ */
    matrix_multiply(Np, Np, Np, Vr[0], inv, cell->Dr[0]);
    matrix_multiply(Np, Np, Np, Vs[0], inv, cell->Ds[0]);
    matrix_multiply(Np, Np, Np, Vt[0], inv, cell->Dt[0]);

    matrix_double_free(Vr);
    matrix_double_free(Vs);
    matrix_double_free(Vt);
    return;
}

/**
 * @brief calculate the Vandermonde matrix
 * @details
 * @param [in] cell standard cell
 * @param [in] orthfunc orthogonal function handle
 * @return Vandermonde matrix
 */
static double** dg_cell_vand_matrix(dg_cell *cell,
                                    void (*orthfunc)(dg_cell *, int ind, double *func))
{
    const int Np = cell->Np;
    // allocation
    double **V = matrix_double_create(Np, Np);

    // assignment
    int dim1,dim2;
    double temp[Np];
    for(dim2=0;dim2<Np;dim2++){
        orthfunc(cell, dim2, temp);
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
static double** dg_cell_mass_matrix(dg_cell *cell){
    const int Np = cell->Np;
    double **V = cell->V;

    double inv[Np*Np];
    double invt[Np*Np];

    int i,j,sk=0;
    // allocation
    double **M = matrix_double_create(Np, Np);

    for(i=0;i<Np*Np;i++) // copy matrix to a vector
        inv[sk++] = V[0][i];

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
 * @brief
 * create LIFT matrix
 * @param cell dg_cell structure
 * @param surf_mass_matrix function that return surface mass matrix
 * @return
 */
static double** dg_cell_lift_matrix(dg_cell *cell,
                                    double **(*surf_mass_matrix)(dg_cell *))
{
    const int Np = cell->Np;
    const int Nfaces = cell->Nfaces;
    const int Nfp = cell->Nfp;
    // allocation
    double **LIFT = matrix_double_create(Np, Nfaces * Nfp);

    double **V = cell->V;
    double vt[Np*Np]; //transpose Vandermonde matrix
    int i, j, sk=0;
    for(i=0;i<Np;i++){
        for(j=0;j<Np;j++){
            vt[sk++] = V[j][i];
        }
    }
    double invM[Np*Np];
    // get the inverse mass matrix M^{-1} = V*V'
    matrix_multiply(Np, Np, Np, V[0], vt, invM);
    // get surface mass matrix Mes
    double **Mes = surf_mass_matrix(cell);
    /* LIFT = M^{-1}*Mes */
    matrix_multiply(Np, Np, Nfp * Nfaces, invM, Mes[0], LIFT[0]);

    // free surface mass matrix
    matrix_double_free(Mes);
    return LIFT;
}