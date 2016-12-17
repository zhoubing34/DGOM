#include "sc_stdcell.h"

/* Quadrilateral.c */
stdCell* sc_createQuad(int N);
/* Triangle.c */
stdCell* sc_createTri(int N);

/**
 * @brief create stand cell
 * @param[in] N order
 * @param[in] cellType cell enum type
 */
stdCell* sc_create(int N, cellType type){
    stdCell *std;
    switch (type){
        case TRIANGLE:
            std = sc_createTri(N); break;
        case QUADRIL:
            std = sc_createQuad(N); break;
        default:
            printf("Unknown cell type %d\n", type);
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
    IntMatrix_free(stdcell->Fmask);
    /* coordinate */
    switch (stdcell->type){
        case TETRA:
        case TRIPRISM:
        case HEXA:
            Vector_free(stdcell->t);
        default: // TRIANGLE and QUADRIL
            Vector_free(stdcell->r);
            Vector_free(stdcell->s);
    }
    /* vandermonde matrix */
    Matrix_free(stdcell->V);
    /* mass matrix */
    Matrix_free(stdcell->M);
    /* Derivative Matrix */
    Matrix_free(stdcell->Dr);
    Matrix_free(stdcell->Ds);

    /* LIFT */
    Matrix_free(stdcell->LIFT);

    /* float version */
    free(stdcell->f_LIFT);
    free(stdcell->f_Dr);
    free(stdcell->f_Ds);
}

/**
 * @brief calculate the Vandermonde matrix
 * @details
 * @param [in] cell standard cell
 * @param [in] orthfunc orthogonal function handle
 * @return Vandermonde matrix
 */
double** sc_VandMatrix2d(stdCell *cell, void (*orthfunc)(stdCell*, int ind, double *func))
{
    int dim1,dim2;
    const int Np = cell->Np;

    // allocation
    double **V = Matrix_create(Np, Np);

    // assignment
    double temp[cell->Np];
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
    double **M = Matrix_create(Np, Np);

    for(i=0;i<Np*Np;i++) // copy matrix to a vector
        inv[sk++] = V[0][i];

    Matrix_Inverse(inv, Np);

    // transpose of invV
    for(i=0;i<Np;i++){
        for(j=0;j<Np;j++){
            invt[j*Np + i] = inv[i*Np + j];
        }
    }

    Matrix_Multiply(Np, Np, Np, invt, inv, M[0]);
    return M;
}

/**
 * @brief get the derivative Vandermonde matrix
 * @param [in] cell standard element
 * @param [in] derorthfunc derivative orthogonal function handle
 * @param [out] Vr the derivative of Vandermonde matrix on r coordinate
 * @param [out] Vs the derivative of Vandermonde matrix on s coordinate
 * @note
 * precondition: Vr and Vs should be allocated before calling @ref sc_deriVandMatrix2d
 */
void sc_deriVandMatrix2d(stdCell *cell, void (*derorthfunc)
        (stdCell *, int ind, double *dr, double *ds), double **Vr, double **Vs){
    const int Np = cell->Np;
    double dr[Np],ds[Np];

    int dim1, dim2;
    for(dim2=0;dim2<Np;dim2++){
        derorthfunc(cell, dim2, dr, ds);
        for(dim1=0;dim1<Np;dim1++){
            Vr[dim1][dim2] = dr[dim1];
            Vs[dim1][dim2] = dr[dim1];
        }
    }
}

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
void sc_deriMatrix2d(stdCell *cell, void (*derorthfunc)
        (stdCell *, int ind, double *dr, double *ds)){
    const int Np = cell->Np;
    // allocation
    cell->Dr = Matrix_create(Np, Np);
    cell->Ds = Matrix_create(Np, Np);

    double **V = cell->V;
    double **Vr = Matrix_create(Np, Np);
    double **Vs = Matrix_create(Np, Np);
    double inv[Np*Np];

    // inverse of Vandermonde matrix
    int i, sk=0;
    for(i=0;i<Np*Np;i++) // copy matrix to a vector
        inv[sk++] = V[0][i];

    Matrix_Inverse(inv, Np);

    // get derivative Vandermonde matrix
    sc_deriVandMatrix2d(cell, derorthfunc, Vr, Vs);

    /* \f$ \mathbf{Dr} = \mathbf{Vr}*\mathbf{V}^{-1} \f$ */
    Matrix_Multiply(Np, Np, Np, Vr[0], inv, cell->Dr[0]);
    /* \f$ \mathbf{Ds} = \mathbf{Vs}*\mathbf{V}^{-1} \f$ */
    Matrix_Multiply(Np, Np, Np, Vs[0], inv, cell->Ds[0]);

    Matrix_free(Vr);
    Matrix_free(Vs);
}


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
    Matrix_Inverse(inv, Nfp);
    /* transform of vandermonde matrix */
    for(i=0;i<Nfp;i++){
        for(j=0;j<Nfp;j++)
            invt[j+Nfp*i] = inv[j*Nfp+i];
    }
    /* get M = inv(V)'*inv(V) */
    Matrix_Multiply(Nfp, Nfp, Nfp, invt, inv, m);

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

/**
 * @brief create LIFT matrix for 2d elements (triangle or quadrilateral)
 * @param [in,out] cell stand element
 *
 */
double** sc_liftMatrix2d(stdCell *cell){
    const int Np = cell->Np;
    const int Nfaces = cell->Nfaces;
    const int Nfp = cell->Nfp;
    // allocation
    double **LIFT = Matrix_create(Np, Nfaces*Nfp);

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
    Matrix_Multiply(Np, Np, Np, V[0], vt, invM);
    // get surface mass matrix Mes
    double **Mes = Matrix_create(Np, Nfaces*Nfp);
    sc_surfMassMatrix2d(cell, Mes);
    /* LIFT = M^{-1}*Mes */
    Matrix_Multiply(Np, Np, Nfp*Nfaces, invM, Mes[0], LIFT[0]);

    // free surface mass matrix
    Matrix_free(Mes);
    return LIFT;
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
    cell->ws = Vector_create(Nfp);
    zwglj(r, cell->ws, Nfp, 0, 0);

    // Gauss quadrature weights for volume
    cell->wv = Vector_create(Np);
    int i,j;
    for(i=0;i<Np;i++){
        for(j=0;j<Np;j++){
            cell->wv[j] += cell->M[i][j];
        }
    }
}