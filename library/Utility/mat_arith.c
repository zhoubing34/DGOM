//
// Created by li12242 on 17/1/29.
//

#include "mat_arith.h"

/**
 * @brief Inverse of square matrix A
 *
 * @param[in] N row of square matrix
 * @param[in] A matrix
 *
 * @note
 * Row counts first to generalize the vector A, which means that A[i][j] = A[i*N+j]
 */
void Matrix_inverse(double *A, int N){
#ifdef MKL_LIB
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
#else
    integer W = (integer) N;
    integer  LDA = W;
    integer  IPIV[W];
    integer  ERR_INFO;
    integer  LWORK = W * W;
    double Workspace[LWORK];
    // - Compute the LU factorization of a M by N matrix W
    dgetrf_(&W, &W, A, &LDA, IPIV, &ERR_INFO);
    // - Generate inverse of the matrix given its LU decompsotion
    dgetri_(&W, A, &LDA, IPIV, Workspace, &LWORK, &ERR_INFO);
#endif
    return;
}


/**
 * @brief matrix multiply
 *
 * @param[in] A is M-by-K matrix
 * @param[in] B is K-by-N matrix
 * @param[out] C is M-by-N matrix
 *
 * @note
 * Row counts first to generalize the vector A, B and C, which means that A[i*N+j] = A[i][j]
 */
void Matrix_multiply(const int M, const int K, const int N,
                     const double *A, const double *B, double *C) {
    int i, j, k;

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