/**
 * @file
 * Matrix and vector operation functions
 *
 * @brief
 *
 * @author
 * li12242, Tianjin University, li12242@tju.edu.cn
 */

#include "LibUtilities.h"

/* matrix operations */

/**
 * @brief
 * Inverse of square matrix A
 *
 * @param[int]          N row of square matrix
 * @param[doublereal]   A[N*N] reshape matrix M
 *
 * @note
 * Row counts first to generalize the vector A, which means that
 * A[i][j] = A[i*N+j]
 */
void invM(double* A, int N){
    int lda = N;
    int ipiv[N];
    int info;

    info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR,N,N,A,lda,ipiv);
    if( info > 0 ) {
        printf( "The algorithm failed to compute LU decomposition.\n" );
        exit( 1 );
    }
    info = LAPACKE_dgetri(LAPACK_ROW_MAJOR,N,A,lda,ipiv);
    if( info > 0 ) {
        printf( "The algorithm failed to compute matrix inverse.\n" );
        exit( 1 );
    }
}


/**
 * @brief
 * Matrix Multiply
 *
 * @details
 * \f[ \mathbf{C} = \mathbf{A}* \mathbf{B} \f]
 *
 * @param[in] lda the leading dimension of the matrix
 * @param[in] A is M-by-K matrix
 * @param[in] B is K-by-N matrix
 * @param[inout] C is M-by-N matrix
 *
 * @note
 * Row counts first to generalize the vector A, B and C, which means that
 * A[i][j] = A[i*N+j]
 */

void dgemm_(const unsigned lda,
             const unsigned M, const unsigned N, const unsigned K,
             const double *A, const double *B, double *C) {
    unsigned i, j, k;

    for (i = 0; i < M; ++i) {
        const double *Ai_ = A + i*K;
        for (j = 0; j < N; ++j) {
            const double *B_j = B + j;

            double cij = 0.0;

            for (k = 0; k < K; ++k) {
                cij += *(Ai_ + k) * *(B_j + k*N);
            }

            *(C + j + i*N) = cij;
        }
    }
}


/* some very basic memory allocation routines */
/* row major storage for a 2D matrix array */
double **BuildMatrix(int Nrows, int Ncols){
    int n;
    double **A = (double**) calloc(Nrows, sizeof(double*));

    A[0] = (double*) calloc(Nrows*Ncols, sizeof(double));

    for(n=1;n<Nrows;++n){
        A[n] = A[n-1]+ Ncols;
    }

    return A;
}
    
double *BuildVector(int Nrows){

    double *A = (double*) calloc(Nrows, sizeof(double));

    return A;
}

/* row major storage for a 2D matrix array */
int **BuildIntMatrix(int Nrows, int Ncols){
  int n;
  int **A = (int**) calloc(Nrows, sizeof(int*));

  A[0] = (int*) calloc(Nrows*Ncols, sizeof(int));

  for(n=1;n<Nrows;++n){
    A[n] = A[n-1]+ Ncols;
  }

  return A;
}

int *BuildIntVector(int Nrows){

  int *A = (int*) calloc(Nrows, sizeof(int));

  return A;
}

double *DestroyVector(double *v){
  free(v);
  return NULL;
}

double **DestroyMatrix(double **A){
  free(A[0]);
  free(A);

  return NULL;
}

int *DestroyIntVector(int *v){
  free(v);
  return NULL;
}

int **DestroyIntMatrix(int **A){
  free(A[0]);
  free(A);

  return NULL;
}
