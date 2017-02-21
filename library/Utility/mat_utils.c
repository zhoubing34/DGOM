/**
 * @file Matrix and vector operation functions
 *
 * @author
 * li12242, Tianjin University, li12242@tju.edu.cn
 */

#include "utility.h"

/* matrix operations */

/**
 * @brief vector division
 *
 * @param[in] N number of elements
 * @param[in] v1 dividend vector
 * @param[in] v2 divisor vector
 */
void Vector_division(int N, double *v1, double *v2){
    register int i;
    for(i=0;i<N;i++)
        v1[i] /= v2[i];
    return;
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
#include "mat_utils.h"

#define __VECTOR_CREATE_FUNC vector_double_create
#include "mat_utils.h"

#define __MATRIX_FREE_FUNC matrix_double_free
#include "mat_utils.h"

#define __VECTOR_FREE_FUNC vector_double_free
#include "mat_utils.h"
#undef __T__

#define __T__ float
#define __MATRIX_CREATE_FUNC matrix_float_create
#include "mat_utils.h"

#define __VECTOR_CREATE_FUNC vector_float_create
#include "mat_utils.h"

#define __MATRIX_FREE_FUNC matrix_float_free
#include "mat_utils.h"

#define __VECTOR_FREE_FUNC vector_float_free
#include "mat_utils.h"
#undef __T__

#define __T__ real
#define __MATRIX_CREATE_FUNC matrix_real_create
#include "mat_utils.h"

#define __VECTOR_CREATE_FUNC vector_real_create
#include "mat_utils.h"

#define __MATRIX_FREE_FUNC matrix_real_free
#include "mat_utils.h"

#define __VECTOR_FREE_FUNC vector_real_free
#include "mat_utils.h"
#undef __T__

#define __T__ int
#define __MATRIX_CREATE_FUNC matrix_int_create
#include "mat_utils.h"

#define __VECTOR_CREATE_FUNC vector_int_create
#include "mat_utils.h"

#define __MATRIX_FREE_FUNC matrix_int_free
#include "mat_utils.h"

#define __VECTOR_FREE_FUNC vector_int_free
#include "mat_utils.h"
#undef __T__
