/**
 * @file
 * StdRegions.h
 *
 * @brief
 * Basic element structure definition
 *
 */

#ifndef STDREGION_H
#define STDREGION_H

#include "LibUtilities/LibUtilities.h"
#include "Polylib/polylib.h"
#include <math.h>

struct StdReg2d{
    /** Polynomial order */
    int N;
    /** number of points */
    int Np;
    /** number of vertex */
    int Nv;
    /** number of faces */
    int Nfaces;
    /** number of points per face */
    int Nfp;

    /** index of node at faces */
    int **Fmask;

    /** face integration coeff */
    double *wj;
    /** volume integration coeff */
    double *wv;

    /* natural coordinate of nodes */
    double *r, *s;
    /** Vandermonde matrix */
    double **V;
    /** mass matrix */
    double **M;

    /* local nodal derivative matrices */
    double **Dr, **Ds;
    /** local lift matrix */
    double **LIFT;
};

typedef struct StdReg2d StdRegions2d;


/* Triangle.c */

void Warpfactor(int, double *, int, double *);
void GetTriCoord(int, double*, double*);

#endif