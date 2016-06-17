/**
 * @file
 * Public function from standard element library
 *
 * @brief
 * Basic element structure definition
 *
 * @author
 * li12242, Tianjin University, li12242@tju.edu.cn
 */

#ifndef STDREGION_H
#define STDREGION_H

#include "LibUtilities/LibUtilities.h"
#include "Polylib/polylib.h"

struct StdReg2d{
    /** Polynomial order */
    unsigned N;
    /** number of points */
    unsigned Np;
    /** number of vertex */
    unsigned Nv;
    /** number of faces */
    unsigned Nfaces;
    /** number of points per face */
    unsigned Nfp;

    /** index of node at faces */
    int **Fmask;

    /** face integration coeff */
    double *ws;
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

/* Quadrilateral.c */

/* Triangle.c - public functions */
StdRegions2d* GenStdTriEle(const unsigned N);
void FreeStdRegions2d(StdRegions2d * );

#endif