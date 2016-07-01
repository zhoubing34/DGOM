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

    /* float version coefficient */
    real *f_Dr, *f_Ds;
    real *f_LIFT;

};
typedef struct StdReg2d StdRegions2d;

/* Quadrilateral.c */
StdRegions2d* GenStdQuadEle(int N);

/* Triangle.c - public functions */
StdRegions2d* GenStdTriEle(int N);
void FreeStdRegions2d(StdRegions2d * );
void GetLIFT2d(StdRegions2d *shape, double **Mes, double **LIFT);
void GetSurfLinM(int N, int Nfaces, int **Fmask, double **Mes);

#endif