/** \file   Convection2d.h
 *  \brief
 *  header file for 2d convection problem
 */

#include <strings.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* apply slope limit */
//#define LIMIT

/* print log file */
//#define DEBUG

/* element shape */
//#define TRI
#define QUAD

/** number of unknown variables */
#define p_Nfields 1

#include "Convection2d/Mesh2d.h"
#include "NcOutput.h"

/* Mesh2d.c */
Mesh* ReadTriMesh();
Mesh* ReadQuadMesh();
void GeometricFactors(Mesh *, int ,
                      double *, double *, double *, double *,
                      double *);
void Normals(Mesh *, int, double *, double *, double *);
void PrintMesh ( Mesh * );

/* PairFace.c */
void FacePair( Mesh * );
void PrintMeshConnection( Mesh * );

/* LoadBalance.c */
void LoadBalance(Mesh *);

/* ParallelPairs.c */
void ParallelPairs(void *, int, int ,
                   int  (*)(const void *),
                   void (*)(const void *, int ),
                   int  (*)(const void *),
                   void (*)(const void *, const void *),
                   int (*)(const void *, const void *));

/* Setup.c */
void SetMeshCoeff(Mesh *mesh);
void SetupTriCoeff(Mesh *);
void SetupQuadCoeff(Mesh *);

/* BuildMaps.c */
void BuildMaps(Mesh*);

/* InitMeshInfo.c */
double InitMeshInfo(Mesh *, int);

/* InitialCondition.c */
void InitData(Mesh *);

/* ConvectionRun2d.c */
void ConvectionRun2d(Mesh *, Ncfile * , double, double);

/* ConvectionRHS2d.c */
void ConvectionRHS2d(Mesh *, float, float, float);

/* ConvectionDriver2d.c */
void ConvectionFinish(Mesh *, Ncfile *);

/* DsiDetector.c */
void DisDetector(Mesh *);

/* LimiterBJ2d.c */
void LimiterBJ2d(Mesh *);

/* Utils.c */
FILE* CreateLog(char *, int, int);

/* Postprocess.c */
void Postprocess(Mesh *);

/* prototypes for storage functions */
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