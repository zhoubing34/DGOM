/**
 * @file Matrix and vector operation functions
 *
 * @author
 * li12242, Tianjin University, li12242@tju.edu.cn
 */

#include "LibUtilities.h"
#include "mkl_lapacke.h"

/* matrix operations */

/**
 * @brief vector division
 *
 * @param[in] N number of elements
 * @param[in] v1 dividend vector
 * @param[in] v2 divisor vector
 */
void Vector_division(int N, double *v1, double *v2){
    int i;
    for(i=0;i<N;i++)
        v1[i] /= v2[i];
    return;
}

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
void Matrix_inverse(double *A, int N){
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

/**
 * @brief
 * Transfer string to integer.
 *
 * @param[in] str the input string
 * @param[out] N  output integer
 * @param[in] errmessage    error message
 *
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

void str2double(char *str, double *scal, char* errmessage){
    int info = sscanf(str,"%lf", scal);
    if (info!=1) {
        fprintf(stderr, "%s:%s \n", errmessage, str);
        exit(-1);
    }
}

/* some very basic memory allocation routines */
/* row major storage for a 2D matrix array */

#define __T__ double
#define __MATRIX_CREATE_FUNC matrix_double_create
#include "Utils.h"

#define __VECTOR_CREATE_FUNC vector_double_create
#include "Utils.h"

#define __MATRIX_FREE_FUNC matrix_double_free
#include "Utils.h"

#define __VECTOR_FREE_FUNC vector_double_free
#include "Utils.h"
#undef __T__

#define __T__ float
#define __MATRIX_CREATE_FUNC matrix_float_create
#include "Utils.h"

#define __VECTOR_CREATE_FUNC vector_float_create
#include "Utils.h"

#define __MATRIX_FREE_FUNC matrix_float_free
#include "Utils.h"

#define __VECTOR_FREE_FUNC vector_float_free
#include "Utils.h"
#undef __T__

#define __T__ real
#define __MATRIX_CREATE_FUNC matrix_real_create
#include "Utils.h"

#define __VECTOR_CREATE_FUNC vector_real_create
#include "Utils.h"

#define __MATRIX_FREE_FUNC matrix_real_free
#include "Utils.h"

#define __VECTOR_FREE_FUNC vector_real_free
#include "Utils.h"
#undef __T__

#define __T__ int
#define __MATRIX_CREATE_FUNC matrix_int_create
#include "Utils.h"

#define __VECTOR_CREATE_FUNC vector_int_create
#include "Utils.h"

#define __MATRIX_FREE_FUNC matrix_int_free
#include "Utils.h"

#define __VECTOR_FREE_FUNC vector_int_free
#include "Utils.h"
#undef __T__
