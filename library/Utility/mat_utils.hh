/**
 * @brief Template for creating matrix variable
 * @details
 * use marco __T__ to determine the type of matrix to create
 *
 *    #include __T__ int
 *    #define __MATRIX_CREATE_FUNC IntMatrix_create
 *    #include "Utils.h"
 *
 * @param[in] Nrows rows number
 * @param[in] Ncols cols number
 * @return double pointer __T__ **A
 */
#ifdef __MATRIX_CREATE_FUNC
__T__** __MATRIX_CREATE_FUNC(int Nrows, int Ncols){
    int n;
    __T__ **A = (__T__ **) calloc((size_t) Nrows, sizeof(__T__ *));

    A[0] = (__T__ *) calloc((size_t) Nrows*Ncols, sizeof(__T__));

    for(n=1;n<Nrows;++n){
        A[n] = A[n-1]+ Ncols;
    }

    return A;
}

#undef __MATRIX_CREATE_FUNC
#endif

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
 * @param[in] Ncols cols number
 * @return pointer __T__ *A
 */
#ifdef __VECTOR_CREATE_FUNC
__T__ *__VECTOR_CREATE_FUNC(int Nrows){

    __T__ *A = (__T__*) calloc((size_t) Nrows, sizeof(__T__));

    return A;
}
#undef __VECTOR_CREATE_FUNC
#endif

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
#ifdef __MATRIX_FREE_FUNC

__T__ **__MATRIX_FREE_FUNC(__T__ **A){
  free(A[0]);
  free(A);

  return NULL;
}
#undef __MATRIX_FREE_FUNC
#endif

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
#ifdef __VECTOR_FREE_FUNC
__T__ *__VECTOR_FREE_FUNC(__T__ *v){
  free(v);
  return NULL;
}
#undef __VECTOR_FREE_FUNC
#endif