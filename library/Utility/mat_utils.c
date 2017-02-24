#include "utility.h"

/**
 * @brief Inverse of square matrix A
 * @param[in] N row of square matrix
 * @param[in] A matrix
 * @note
 * Row counts first of A, which is A[i][j] = A[i*N+j]
 */
void matrix_inverse(double *A, int N){

//    int lda = N;
//    int ipiv[N];
//    int info;
//
//    info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR,N,N,A,lda,ipiv);
//    if( info > 0 ) {
//        printf( "The algorithm failed to compute LU decomposition.\n" );
//        exit( 1 );
//    }
//    info = LAPACKE_dgetri(LAPACK_ROW_MAJOR,N,A,lda,ipiv);
//    if( info > 0 ) {
//        printf( "The algorithm failed to compute matrix inverse.\n" );
//        exit( 1 );
//    }

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
    return;
}

/**
 * @brief matrix multiply
 * @param[in] A M-by-K matrix
 * @param[in] B K-by-N matrix
 * @param[out] C M-by-N matrix
 *
 * @note
 * Row counts first of A, B and C, which is A[i*N+j] = A[i][j]
 */
void matrix_multiply(const int M, const int K, const int N,
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

/* some very basic memory allocation routines */
/* row major storage for a 2D matrix array */

#define __T__ double
#define __MATRIX_CREATE_FUNC matrix_double_create
#include "mat_utils.hh"

#define __VECTOR_CREATE_FUNC vector_double_create
#include "mat_utils.hh"

#define __MATRIX_FREE_FUNC matrix_double_free
#include "mat_utils.hh"

#define __VECTOR_FREE_FUNC vector_double_free
#include "mat_utils.hh"
#undef __T__

#define __T__ float
#define __MATRIX_CREATE_FUNC matrix_float_create
#include "mat_utils.hh"

#define __VECTOR_CREATE_FUNC vector_float_create
#include "mat_utils.hh"

#define __MATRIX_FREE_FUNC matrix_float_free
#include "mat_utils.hh"

#define __VECTOR_FREE_FUNC vector_float_free
#include "mat_utils.hh"
#undef __T__

#define __T__ dg_real
#define __MATRIX_CREATE_FUNC matrix_real_create
#include "mat_utils.hh"

#define __VECTOR_CREATE_FUNC vector_real_create
#include "mat_utils.hh"

#define __MATRIX_FREE_FUNC matrix_real_free
#include "mat_utils.hh"

#define __VECTOR_FREE_FUNC vector_real_free
#include "mat_utils.hh"
#undef __T__

#define __T__ int
#define __MATRIX_CREATE_FUNC matrix_int_create
#include "mat_utils.hh"

#define __VECTOR_CREATE_FUNC vector_int_create
#include "mat_utils.hh"

#define __MATRIX_FREE_FUNC matrix_int_free
#include "mat_utils.hh"

#define __VECTOR_FREE_FUNC vector_int_free
#include "mat_utils.hh"
#undef __T__
