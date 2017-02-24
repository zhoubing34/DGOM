/**
 * @file
 * 矩阵与向量基本操作
 * @details
 * 1. `matrix_inverse` calculate the inverse of matrix;
 * 2. `matrix_multiply` multiply the matrix;
 * 3. `matrix_T_create` create T type matrix;
 * 4. `matrix_T_free` free matrix of T type;
 * 5. `vector_T_create` create T type vector;
 * 6. `vector_T_free` free vector of T type;
 *
 * @author
 * li12242, Tianjin University, li12242@tju.edu.cn
 */

#ifndef DGOM_MAT_UTILS_H
#define DGOM_MAT_UTILS_H

#include "utility.h"
#include "f2c.h"
#include "blaswrap.h"
#include "clapack.h"

/* matrix operation */
void matrix_inverse(double *A, int N);
void matrix_multiply(const int M, const int K, const int N,
                     const double *A, const double *B, double *C);

/* allocation */
double **matrix_double_create(int Nrows, int Ncols);
double * vector_double_create(int Nrows);
int    **matrix_int_create(int Nrows, int Ncols);
int    * vector_int_create(int Nrows);
float  **matrix_float_create(int Nrows, int Ncols);
float  * vector_float_create(int Nrows);
dg_real   **matrix_real_create(int Nrows, int Ncols);
dg_real   * vector_real_create(int Nrows);

double **matrix_double_free(double **);
double * vector_double_free(double *);
int    **matrix_int_free(int **);
int    * vector_int_free(int *);
float  **matrix_float_free(float **);
float  * vector_float_free(float *);
dg_real   **matrix_real_free(dg_real **);
dg_real   * vector_real_free(dg_real *);

#endif