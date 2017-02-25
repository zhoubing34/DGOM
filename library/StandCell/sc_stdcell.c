#include "sc_stdcell.h"
#include "sc_quadrilateral.h"
#include "sc_triangle.h"

/**
 * @brief return a pointer to a new stand cell object
 * @param[in] N order
 * @param[in] cellType cell type
 */
stdCell* sc_create(int N, sc_cellType type){
    stdCell *std;
    switch (type){
        case TRIANGLE:
            std = sc_create_tri(N);
            std->dim = 2;
            break;
        case QUADRIL:
            std = sc_create_quad(N);
            std->dim = 2;
            break;
        default:
            fprintf(stderr, "%s (%d): Unknown cell type %d\n",
                    __FUNCTION__, __LINE__, type);
            exit(-1);
    }
    return std;
}

/**
 * @brief free stdCell object
 * @param[in] stdcell stdCell object
 */
void sc_free(stdCell *stdcell){
    /* fmaks */
    matrix_int_free(stdcell->Fmask);
    /* coordinate */
    switch (stdcell->dim){
        case 2:
            vector_double_free(stdcell->r);
            vector_double_free(stdcell->s);
            break;
        case 3:
            vector_double_free(stdcell->r);
            vector_double_free(stdcell->s);
            vector_double_free(stdcell->t);
            break;
        default: // TRIANGLE and QUADRIL
            printf("StandCell (sc_free): Wrong dimensions %d\n", stdcell->dim);
            exit(-1);
    }
    /* vandermonde matrix */
    matrix_double_free(stdcell->V);
    /* mass matrix */
    matrix_double_free(stdcell->M);
    /* Derivative Matrix */
    matrix_double_free(stdcell->Dr);
    matrix_double_free(stdcell->Ds);

    /* LIFT */
    matrix_double_free(stdcell->LIFT);

    /* Gauss quadrature */
    vector_double_free(stdcell->ws);
    vector_double_free(stdcell->wv);

    /* float version */
    free(stdcell->f_LIFT);
    free(stdcell->f_Dr);
    free(stdcell->f_Ds);

    free(stdcell);
}

/**
 * @brief Project the vertex value to interpolation nodes.
 * @param[in] vertVal value of vertex
 * @param[in] nodeVal value of nodes
 */
void sc_proj_vert2node(stdCell *cell, double *vertVal, double *nodeVal){
    switch (cell->type){
        case TRIANGLE:
            sc_proj_vert2node_tri(cell, vertVal, nodeVal); break;
        case QUADRIL:
            sc_proj_vert2node_quad(cell, vertVal, nodeVal); break;
        default:
            printf("%s (%d): Unknown cell type %d\n", __FUNCTION__, __LINE__, cell->type);
            exit(-1);
    }
    return;
}

/**
 * @brief calculate the Vandermonde matrix
 * @details
 * @param [in] cell standard cell
 * @param [in] orthfunc orthogonal function handle
 * @return Vandermonde matrix
 */
double** sc_VandMatrix(stdCell *cell, void (*orthfunc)(stdCell *, int ind, double *func))
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
double** sc_massMatrix(stdCell *cell){
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
 * @brief create LIFT matrix for 2d elements (triangle or quadrilateral)
 * @param [in,out] cell stand element
 */
double** sc_liftMatrix(stdCell *cell, void (*surf_mass_matrix)(stdCell *, double **))
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
    double **Mes = matrix_double_create(Np, Nfaces * Nfp);
    surf_mass_matrix(cell, Mes);
    /* LIFT = M^{-1}*Mes */
    matrix_multiply(Np, Np, Nfp * Nfaces, invM, Mes[0], LIFT[0]);

    // free surface mass matrix
    matrix_double_free(Mes);
    return LIFT;
}