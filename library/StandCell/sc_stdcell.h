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
    TRIANGLE=0, ///< triangle
    QUADRIL=1,  ///< quadrilateral
    TETRA=2,    ///< tetrahedron
    TRIPRISM=3, ///< triangle prism
    HEXA=4,     ///< hexahedron
} sc_cellType;

typedef struct {
    sc_cellType type; ///< cell enum type

    int dim; ///< dimension

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

/* ======================== functions for standard elements ======================== */
/* create stand cell object */
stdCell* sc_create(int N, sc_cellType type);
/* free stdCell object */
void sc_free(stdCell *);
/* calculate the mass matrix */
double** sc_massMatrix(stdCell *cell);
/* project the vertex value to interpolation nodes */
void sc_vertProj(stdCell *cell, double *vertVal, double *nodeVal);


/* ======================== functions for 2d elements (tri and quad) ======================== */
/* get the gradient matrix Dr and Ds of Lagrange basis at (r,s) at order N */
void sc_deriMatrix2d(stdCell *cell, void (*derorthfunc)(stdCell *, int ind, double *dr, double *ds));
/* calculate the Vandermonde matrix */
double** sc_VandMatrix2d(stdCell *cell, void (*orthfunc)(stdCell *, int ind, double *func));
/* create LIFT matrix for 2d elements (triangle or quadrilateral) */
double** sc_liftMatrix2d(stdCell *cell);
/* calculate the Gauss quadrature weights for faces (ws) and volume (wv) integral */
void sc_GaussQuadrature2d(stdCell *cell);


#endif //DGOM_STDCELL_H