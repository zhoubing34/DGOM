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

#ifndef DGOM_STDCELL_H
#define DGOM_STDCELL_H

#include "LibUtilities/LibUtilities.h"

typedef enum {
    TRIANGLE, ///< triangle
    QUADRIL,  ///< quadrilateral
    TETRA,    ///< tetrahedron
    TRIPRISM, ///< triangle prism
    HEXA      ///< hexahedron
} cellType;

typedef struct {
    cellType type; ///< cell enum type

    int N; ///< polynomial order
    int Np; ///< number of points
    int Nv; ///< number of vertex
    int Nfaces; ///< number of faces
    int Nfp; ///< number of points on each face

    int **Fmask; ///< index of node at faces
    double *ws; ///< surface Gauss quadrature
    double *wv; ///< volume Gauss quadrature

    double *r, *s, *t; ///< coordinate
    double **V; ///< Vandermonde matrix
    double **M; ///< mass matrix

    double **Dr, **Ds; ///< nodal basis derivative matrix
    double **LIFT; ///< lift matrix

    /* float version coefficient */
    real *f_Dr, *f_Ds, *f_LIFT; ///< user specific version

} stdCell;

/* functions for standard elements */
stdCell* sc_create(int N, cellType type);
void sc_free(stdCell *);
double** sc_massMatrix(stdCell *cell);

/* functions for 2d elements (tri and quad) */
void sc_deriMatrix2d(stdCell *cell, void (*derorthfunc)(stdCell *, int ind, double *dr, double *ds));
double** sc_VandMatrix2d(stdCell *cell, void (*orthfunc)(stdCell *, int ind, double *func));
double** sc_liftMatrix2d(stdCell *cell);
void sc_GaussQuadrature2d(stdCell *cell);

#endif //DGOM_STDCELL_H