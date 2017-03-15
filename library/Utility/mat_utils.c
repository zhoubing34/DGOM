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
#define ARRAY_CREATE_FUNC_TEMPLATE(__T__, __FUNC_NAME)           \
__T__** __FUNC_NAME(int Nrows, int *Ncols){             \
    int n, Ntol = 0;                                    \
    for(n=0;n<Nrows;n++){                               \
        Ntol += Ncols[n];                               \
    }                                                   \
    __T__ **A = (__T__ **) calloc((size_t) Nrows, sizeof(__T__ *));     \
    A[0] = (__T__ *) calloc((size_t) Ntol, sizeof(__T__));              \
    for(n=1;n<Nrows;++n){                               \
        A[n] = A[n-1]+ Ncols[n-1];                      \
    }                                                   \
    return A;                                           \
}

ARRAY_CREATE_FUNC_TEMPLATE(int, array_int_create)
ARRAY_CREATE_FUNC_TEMPLATE(double, array_double_create)
ARRAY_CREATE_FUNC_TEMPLATE(float, array_float_create)
ARRAY_CREATE_FUNC_TEMPLATE(dg_real, array_real_create)

/**
 * @brief Template for creating matrix variable
 * @details
 * use marco __T__ to determine the type of matrix to create
 *
 *    MATRIX_CREATE_FUNC_TEMPLATE(int, matrix_int_create)
 *
 * @param[in] Nrows rows number
 * @param[in] Ncols cols number
 * @return double pointer __T__ **A
 */
#define MATRIX_CREATE_FUNC_TEMPLATE(__T__, __FUNC_NAME)     \
__T__** __FUNC_NAME(int Nrows, int Ncols){                  \
    int n;                                                  \
    __T__ **A = (__T__ **) calloc((size_t) Nrows, sizeof(__T__ *));     \
    A[0] = (__T__ *) calloc((size_t) Nrows*Ncols, sizeof(__T__));       \
    for(n=1;n<Nrows;++n){                                   \
        A[n] = A[n-1]+ Ncols;                               \
    }                                                       \
    return A;                                               \
}

MATRIX_CREATE_FUNC_TEMPLATE(int, matrix_int_create)
MATRIX_CREATE_FUNC_TEMPLATE(double, matrix_double_create)
MATRIX_CREATE_FUNC_TEMPLATE(float, matrix_float_create)
MATRIX_CREATE_FUNC_TEMPLATE(dg_real, matrix_real_create)


/**
 * @brief Template for creating vector variable
 * @details
 * use marco __T__ to determine the type of matrix to create
 *
 *    #include __T__ int
 *    #define __VECTOR_CREATE_FUNC IntVector_create
 *    #include "Utils.h"
 *
 * @param[in] Nrows rows number
 * @return pointer __T__ *A
 */
#define VECTOR_CREATE_FUNC_TEMPLATE(__T__, __FUNC_NAME)     \
__T__ *__FUNC_NAME(int Nrows){                     \
    __T__ *A = (__T__*) calloc((size_t) Nrows, sizeof(__T__));  \
    return A;                                               \
}

VECTOR_CREATE_FUNC_TEMPLATE(int, vector_int_create)
VECTOR_CREATE_FUNC_TEMPLATE(double, vector_double_create)
VECTOR_CREATE_FUNC_TEMPLATE(float, vector_float_create)
VECTOR_CREATE_FUNC_TEMPLATE(dg_real, vector_real_create)


/**
 * @brief Template for free matrix memory
 * @details
 * use marco __T__ to determine the type of matrix to create
 *
 *    #include __T__ int
 *    #define __MATRIX_FREE_FUNC IntMatrix_free
 *    #include "Utils.h"
 *
 * @param[in] A matrix pointer
 * @return NULL
 */
#define MATRIX_FREE_FUNC_TEMPLATE(__T__, __FUNC_NAME)   \
__T__ **__FUNC_NAME(__T__ **A){         \
    free(A[0]);                         \
    free(A);                            \
    return NULL;                        \
}
MATRIX_FREE_FUNC_TEMPLATE(int, array_int_free)
MATRIX_FREE_FUNC_TEMPLATE(double, array_double_free)
MATRIX_FREE_FUNC_TEMPLATE(float, array_float_free)
MATRIX_FREE_FUNC_TEMPLATE(dg_real, array_real_free)

MATRIX_FREE_FUNC_TEMPLATE(int, matrix_int_free)
MATRIX_FREE_FUNC_TEMPLATE(double, matrix_double_free)
MATRIX_FREE_FUNC_TEMPLATE(float, matrix_float_free)
MATRIX_FREE_FUNC_TEMPLATE(dg_real, matrix_real_free)

/**
 * @brief Template for free vector memory
 * @details
 * use marco __T__ to determine the type of vector to create
 *
 *    #include __T__ int
 *    #define __MATRIX_FREE_FUNC IntVector_free
 *    #include "Utils.h"
 *
 * @param[in] A vector pointer
 * @return NULL
 */
#define VECTOR_FREE_FUNC_TEMPLATE(__T__, __FUNC_NAME)   \
__T__ *__FUNC_NAME(__T__ *v){   \
    free(v);                    \
    return NULL;                \
}                               \

VECTOR_FREE_FUNC_TEMPLATE(int, vector_int_free)
VECTOR_FREE_FUNC_TEMPLATE(double, vector_double_free)
VECTOR_FREE_FUNC_TEMPLATE(float, vector_float_free)
VECTOR_FREE_FUNC_TEMPLATE(dg_real, vector_dg_real_free)
