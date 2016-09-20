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

/* variable type */
#ifdef DOUBLE_PRECISION
typedef double real;
#else
typedef float real;
#endif

#define max(a,b)  ( (a>b)?a:b )
#define min(a,b)  ( (a<b)?a:b )

#define TOTALERR 1.0e-6
#define RELATIVEERROR 1.0e-8

/* string to int */
void str2int(char *str, int *N, char* errmessage);

/* GenUniformMesh */
void GenUniformTriMesh(int Mx, int My, double xmin, double xmax, double ymin, double ymax, int type,
                       int *Ne, int *Nv, int ***newEToV, double **newVX, double **newVY);

void GenUniformQuadMesh(int Mx, int My, double xmin, double xmax, double ymin, double ymax,
                        int *Ne, int *Nv, int ***newEToV, double **newVX, double **newVY);

void GenParallelUniformTriMesh(int Mx, int My, double xmin, double xmax, double ymin, double ymax, int type,
                               int procid, int nprocs, int *parK, int *newNv,
                               int ***parEToV, double **VX, double **VY);

void GenParallelUniformQuadMesh(int Mx, int My, double xmin, double xmax, double ymin, double ymax,
                                int procid, int nprocs, int *parK, int *newNv,
                                int ***parEToV, double **VX, double **VY);

/* matrix operation */
void invM(double* A, int N);
void dgemm_(const unsigned M, const unsigned K, const unsigned N,
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
void PrintIntMatrix(char *message, int **A, int Nrows, int Ncols);
void SaveMatrix(char *filename, double **A, int Nrows, int Ncols);
void PrintIntMatrix2File(FILE *fp, char *message, int **Mat, int row, int col);
void PrintMatrix2File(FILE *fp, char *message, double **Mat, int row, int col);
void PrintIntVector2File(FILE *fp, char *message, int *Mat, int len);
void PrintVector2File(FILE *fp, char *message, double *Mat, int len);

int CreateVectorTest(char *message, double *A, double *ExactA, int Ncols);
int CreateMatrixTest(char *message, double **A, double **ExactA, int Nrows, int Ncols);
FILE* CreateLog(char *funname, int nprocs, int rank);

int CreateIntMatrixTest(char *message, int **A, int **ExactA, int Nrows, int Ncols);

#endif