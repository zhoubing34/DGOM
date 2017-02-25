/**
 * @file
 * Public function for standard element
 *
 * @brief
 * Basic element structure definition
 *
 * @author
 * li12242, Tianjin University, li12242@tju.edu.cn
 */

#ifndef DGOM_STDCELL_H
#define DGOM_STDCELL_H

#include "Utility/utility.h"

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

    double **Dr; ///< nodal basis derivative matrix
    double **Ds; ///< nodal basis derivative matrix
    double **LIFT; ///< lift matrix

    /* float version coefficient */
    dg_real *f_Dr, *f_Ds, *f_LIFT; ///< user specific version

} stdCell;

/* ======================== functions for standard elements ======================== */
/* create stand cell object */
stdCell* sc_create(int N, sc_cellType type);
/* free stdCell object */
void sc_free(stdCell *);
/* project the vertex value to interpolation nodes */
void sc_proj_vert2node(stdCell *cell, double *vertVal, double *nodeVal);

/* calculate the Vandermonde matrix */
double** sc_VandMatrix(stdCell *cell, void (*orthfunc)(stdCell *, int ind, double *func));
/* calculate the mass matrix */
double** sc_massMatrix(stdCell *cell);
/* create LIFT matrix */
double** sc_liftMatrix(stdCell *cell, void (*surf_mass_matrix)(stdCell *, double **));

#endif //DGOM_STDCELL_H