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
 * @param[in] N row of square matrix
 * @param[in] A[N*N] reshape matrix M
 *
 * @note
 * Row counts first to generalize the vector A, which means that
 * A[i][j] = A[i*N+j]
 */
void Matrix_Inverse(double *A, int N){
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
 * @param[in] A is M-by-K matrix
 * @param[in] B is K-by-N matrix
 * @param[in,out] C is M-by-N matrix
 *
 * @note
 * Row counts first to generalize the vector A, B and C, which means that
 * A[i][j] = A[i*N+j]
 */

void Matrix_Multiply(const unsigned M, const unsigned K, const unsigned N,
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

/**
 * @brief
 * Transfer string to integer.
 *
 * @param[in] str           the input string
 * @param[in] errmessage    error message
 *
 * @return
 * return values:
 * name     | type     | description of value
 * -------- |----------|----------------------
 * N        | int      | the integer result
 *
 * Usages:
 *
 *      str2int(argv[1], &N , "Wrong degree input");
 */
void str2int(char *str, int *N, char* errmessage){
    int info = sscanf(str,"%d",N);
    if (info!=1) {
        fprintf(stderr, "%s:%s \n", errmessage, str);
        exit(-1);
    }
}

/* some very basic memory allocation routines */
/* row major storage for a 2D matrix array */
double **Matrix_create(int Nrows, int Ncols){
    int n;
    double **A = (double**) calloc(Nrows, sizeof(double*));

    A[0] = (double*) calloc(Nrows*Ncols, sizeof(double));

    for(n=1;n<Nrows;++n){
        A[n] = A[n-1]+ Ncols;
    }

    return A;
}
    
double *Vector_create(int Nrows){

    double *A = (double*) calloc(Nrows, sizeof(double));

    return A;
}

/* row major storage for a 2D matrix array */
int **IntMatrix_create(int Nrows, int Ncols){
  int n;
  int **A = (int**) calloc(Nrows, sizeof(int*));

  A[0] = (int*) calloc(Nrows*Ncols, sizeof(int));

  for(n=1;n<Nrows;++n){
    A[n] = A[n-1]+ Ncols;
  }

  return A;
}

int *IntVector_create(int Nrows){

  int *A = (int*) calloc(Nrows, sizeof(int));

  return A;
}

double *Vector_free(double *v){
  free(v);
  return NULL;
}

double **Matrix_free(double **A){
  free(A[0]);
  free(A);

  return NULL;
}

int *IntVector_free(int *v){
  free(v);
  return NULL;
}

int **IntMatrix_free(int **A){
  free(A[0]);
  free(A);

  return NULL;
}
