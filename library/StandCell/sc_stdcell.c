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
 * @param [in] orthfunc orthogonal function
 * @return Vandermonde matrix
 */
double** sc_VandMatrix2d(stdCell *cell, void (*orthfunc)
        (stdCell *cell, int ind, double *func))
{
    int dim1,dim2;
    const int Np = cell->Np;
    const int N = cell->N;
    double *r = cell->r;
    double *s = cell->s;

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

    double temp[Np*Np];
    double invt[Np*Np];

    // allocation
    double **M = Matrix_create(Np, Np);

    int i,j,sk=0;

    for(i=0;i<Np;i++){
        for(j=0;j<Np;j++){
            /* row counts first */
            temp[sk++] = V[i][j];
        }
    }

    Matrix_Inverse(temp, Np);

    // transpose of invV
    for(i=0;i<Np;i++){
        for(j=0;j<Np;j++){
            invt[j*Np + i] = temp[i*Np + j];
        }
    }

    Matrix_Multiply(Np, Np, Np, invt, temp, M[0]);

    return M;
}