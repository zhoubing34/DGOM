#ifndef LIBUTILITIES_H
#define LIBUTILITIES_H

#include <strings.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mkl_lapacke.h"
//#include "f2c.h"
//#include "clapack.h"

#include "mpi.h"

#define max(a,b)  ( (a>b)?a:b )
#define min(a,b)  ( (a<b)?a:b )

#define TOTALERR 1.0e-8
#define RELATIVEERROR 1.0e-10

/* matrix operation */
void invM(double* A, int N);
void dgemm_(const unsigned lda,
            const unsigned M, const unsigned N, const unsigned K,
            const double *A, const double *B, double *C);

/* allocate mem */
double **BuildMatrix(int Nrows, int Ncols);
double  *BuildVector(int Nrows);
int    **BuildIntMatrix(int Nrows, int Ncols);
int     *BuildIntVector(int Nrows);

double **DestroyMatrix(double **);
double  *DestroyVector(double *);
int    **DestroyIntMatrix(int **);
int     *DestroyIntVector(int *);

/* ParallelPairs.c */
void ParallelPairs(void *objs, int Nmyobjs, int sizeobj,
                   int  (*numget)(const void *),
                   void (*numset)(const void *, int ),
                   int  (*procget)(const void *),
                   void (*marry)(const void *, const void *),
                   int (*compare_objs)(const void *, const void *));

/* UTest.c */
void PrintVector(char *message, double *A, int Ncols);
void PrintMatrix(char *message, double **A, int Nrows, int Ncols);
void SaveMatrix(char *filename, double **A, int Nrows, int Ncols);

int CreateVectorTest(char *message, double *A, double *ExactA, int Ncols);
int CreateMatrixTest(char *message, double **A, double **ExactA, int Nrows, int Ncols);
FILE* CreateLog(char *funname, int nprocs, int rank);

#endif