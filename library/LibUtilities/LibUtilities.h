#ifndef LIBUTILITIES_H
#define LIBUTILITIES_H

#include <strings.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "mpi.h"

double **BuildMatrix(int Nrows, int Ncols);
double  *BuildVector(int Nrows);
int    **BuildIntMatrix(int Nrows, int Ncols);
int     *BuildIntVector(int Nrows);

double **DestroyMatrix(double **);
double  *DestroyVector(double *);
int    **DestroyIntMatrix(int **);
int     *DestroyIntVector(int *);

void PrintMatrix(char *message, double **A, int Nrows, int Ncols);
void SaveMatrix(char *filename, double **A, int Nrows, int Ncols);

FILE* CreateLog(char *funname, int nprocs, int rank);

#define max(a,b)  ( (a>b)?a:b )
#define min(a,b)  ( (a<b)?a:b )

void ParallelPairs(void *objs, int Nmyobjs, int sizeobj,
                   int  (*numget)(const void *),
                   void (*numset)(const void *, int ),
                   int  (*procget)(const void *),
                   void (*marry)(const void *, const void *),
                   int (*compare_objs)(const void *, const void *));


#endif