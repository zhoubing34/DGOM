#ifndef LIBUTILITIES_H
#define LIBUTILITIES_H

#include <strings.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "pnetcdf.h"

/* float type */
#ifdef DOUBLE_PRECISION
typedef double real; ///> user defined float type
#define MPI_TYPE MPI_DOUBLE ///> variable type for MPI subroutines
#define NC_TYPE NC_DOUBLE ///> variable type for netCDF subroutines
#else
typedef float real; ///> user defined float type
#define MPI_SIZE MPI_FLOAT ///> variable type for MPI subroutines
#define NC_TYPE NC_FLOAT ///> variable type for netCDF subroutines
#endif

/* min and max function */
#define max(a,b)  ( (a>b)?a:b )
#define min(a,b)  ( (a<b)?a:b )

/* string to int */
void str2int(char *str, int *N, char* errmessage);
void str2double(char *str, double *scal, char* errmessage);

/* matrix operation */
void Matrix_inverse(double *A, int N);
void Matrix_multiply(const int M, const int K, const int N,
                     const double *A, const double *B, double *C);

/* allocate mem */
double **matrix_double_create(int Nrows, int Ncols);
double * vector_double_create(int Nrows);
int    **matrix_int_create(int Nrows, int Ncols);
int    * vector_int_create(int Nrows);
float  **matrix_float_create(int Nrows, int Ncols);
float  * vector_float_create(int Nrows);
real   **matrix_real_create(int Nrows, int Ncols);
real   * vector_real_create(int Nrows);

double **matrix_double_free(double **);
double * vector_double_free(double *);
int    **matrix_int_free(int **);
int    * vector_int_free(int *);
float  **matrix_float_free(float **);
float  * vector_float_free(float *);
real   **matrix_real_free(real **);
real   * vector_real_free(real *);

#endif //LIBUTILITIES_H