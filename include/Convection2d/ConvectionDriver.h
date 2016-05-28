#include <strings.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#if 1 // triangle shape
#define Tri
#else // quad shape
#endif

#include "Convection2d/Mesh2d.h"

Mesh* ReadTriMesh();
Mesh* ReadQuadMesh();

/* prototypes for storage functions (Utils.c) */
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

void PrintMeshTri ( Mesh * );
void FacePairTri( Mesh * );
void PrintMeshConnectionTri( Mesh * );


void ParallelPairs(void *objs, int Nmyobjs, int sizeobj,
                   int  (*numget)(const void *),
                   void (*numset)(const void *, int ),
                   int  (*procget)(const void *),
                   void (*marry)(const void *, const void *),
                   int (*compare_objs)(const void *, const void *));

void SetupTri(Mesh *);


void BuildTriMaps(Mesh *);

void NormalsTri(Mesh *, int, double *, double *, double *);

void GeometricFactorsTri(Mesh *, int ,
                         double *, double *, double *, double *,
                         double *);

double InitTriMeshInfo(Mesh *, int);

// Initial Condition
void InitData(Mesh *);

void ConvectionRun2d(Mesh *, double, double);

void ConvectionRHS2d(Mesh *, float, float, float);