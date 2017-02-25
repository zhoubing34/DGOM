//
// Created by li12242 on 17/2/25.
//

#include "sc_stdcell2d.h"
#include "Polylib/polylib.h"

static void sc_deriVandMatrix2d(stdCell *cell,
                         void (*derorthfunc)(stdCell *, int ind, double *dr, double *ds),
                         double **Vr, double **Vs);
/**
 * @brief get the gradient matrix Dr and Ds of Lagrange basis at (r,s) at order N
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
void sc_deriMatrix2d(stdCell *cell, void (*derorthfunc)(stdCell *, int ind, double *dr, double *ds))
{
    const int Np = cell->Np;
    // allocation
    cell->Dr = matrix_double_create(Np, Np);
    cell->Ds = matrix_double_create(Np, Np);

    double **V = cell->V;
    double **Vr = matrix_double_create(Np, Np);
    double **Vs = matrix_double_create(Np, Np);
    double inv[Np*Np];

    // inverse of Vandermonde matrix
    int i, sk=0;
    for(i=0;i<Np*Np;i++) // copy matrix to a vector
        inv[sk++] = V[0][i];

    matrix_inverse(inv, Np);

    // get derivative Vandermonde matrix
    sc_deriVandMatrix2d(cell, derorthfunc, Vr, Vs);

    /* \f$ \mathbf{Dr} = \mathbf{Vr}*\mathbf{V}^{-1} \f$ */
    matrix_multiply(Np, Np, Np, Vr[0], inv, cell->Dr[0]);
    /* \f$ \mathbf{Ds} = \mathbf{Vs}*\mathbf{V}^{-1} \f$ */
    matrix_multiply(Np, Np, Np, Vs[0], inv, cell->Ds[0]);

    matrix_double_free(Vr);
    matrix_double_free(Vs);
}

/**
 * @brief get the derivative Vandermonde matrix
 * @param [in] cell standard element
 * @param [in] derorthfunc derivative orthogonal function handle
 * @param [out] Vr the derivative of Vandermonde matrix on r coordinate
 * @param [out] Vs the derivative of Vandermonde matrix on s coordinate
 * @note
 * Vr and Vs should be allocated before calling sc_deriVandMatrix2d
 */
 static void sc_deriVandMatrix2d(stdCell *cell,
                                 void (*derorthfunc)(stdCell *, int ind, double *dr, double *ds),
                                 double **Vr, double **Vs)
{
    const int Np = cell->Np;
    double dr[Np],ds[Np];

    int dim1, dim2;
    for(dim2=0;dim2<Np;dim2++){
        derorthfunc(cell, dim2, dr, ds);
        for(dim1=0;dim1<Np;dim1++){
            Vr[dim1][dim2] = dr[dim1];
            Vs[dim1][dim2] = ds[dim1];
        }
    }
}

/**
 * @brief calculate the Gauss quadrature weights for faces (ws) and volume (wv) integral
 * @param [in] cell 2d standard elements
 */
void sc_GaussQuadrature2d(stdCell *cell){
    const int Nfp = cell->Nfp;
    const int Np = cell->Np;
    double r[Nfp];
    // Gauss quadrature weights for face
    cell->ws = vector_double_create(Nfp);
    zwglj(r, cell->ws, Nfp, 0, 0);

    // Gauss quadrature weights for volume
    cell->wv = vector_double_create(Np);
    int i,j;
    for(i=0;i<Np;i++){
        for(j=0;j<Np;j++){
            cell->wv[j] += cell->M[i][j];
        }
    }
}

/**
 * @brief create mass matrix of edges
 * @param [in] cell standard element
 * @param [out] Mes mass matrix of edges
 */
void sc_surfMassMatrix2d(stdCell *cell, double **Mes){
    const int Nfp = cell->Nfp;
    const int Nfaces = cell->Nfaces;

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
}