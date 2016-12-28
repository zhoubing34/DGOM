#ifndef LIBUTILITIES_H
#define LIBUTILITIES_H

#include <strings.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "pnetcdf.h"

/* float type */
#ifdef DOUBLE_PRECISION
typedef double real;
#define MPI_TYPE MPI_DOUBLE
#define NC_TYPE NC_DOUBLE
#else
typedef float real;
#define MPI_SIZE MPI_FLOAT
#define NC_TYPE NC_FLOAT
#endif

/* min and max function */
#define max(a,b)  ( (a>b)?a:b )
#define min(a,b)  ( (a<b)?a:b )

/* string to int */
void str2int(char *str, int *N, char* errmessage);
void str2double(char *str, double *scal, char* errmessage);

/* matrix operation */
void Matrix_inverse(double *A, int N);
void Matrix_multiply(const unsigned M, const unsigned K, const unsigned N,
                     const double *A, const double *B, double *C);

/* allocate mem */
double **Matrix_create(int Nrows, int Ncols);
double *Vector_create(int Nrows);
int    **IntMatrix_create(int Nrows, int Ncols);
int    *IntVector_create(int Nrows);

double **Matrix_free(double **);
double *Vector_free(double *);
int    **IntMatrix_free(int **);
int    *IntVector_free(int *);

#endif //LIBUTILITIES_H