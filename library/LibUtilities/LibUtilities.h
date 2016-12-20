#ifndef LIBUTILITIES_H
#define LIBUTILITIES_H

#include <strings.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* float type */
#ifdef DOUBLE_PRECISION
typedef double real;
#define MPI_SIZE MPI_DOUBLE
#else
typedef float real;
#define MPI_SIZE MPI_FLOAT
#endif

/* min and max function */
#define max(a,b)  ( (a>b)?a:b )
#define min(a,b)  ( (a<b)?a:b )

/* string to int */
void str2int(char *str, int *N, char* errmessage);

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

/* ParallelPairs.c */
void ParallelPairs(void *objs, int Nmyobjs, int sizeobj,
                   int  (*numget)(const void *),
                   void (*numset)(const void *, int ),
                   int  (*procget)(const void *),
                   void (*marry)(const void *, const void *),
                   int (*compare_objs)(const void *, const void *));

#endif //LIBUTILITIES_H