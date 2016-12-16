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
#include "Polylib/polylib.h"

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

stdCell* sc_create(int N, cellType type);
void sc_free(stdCell *);

double** sc_VandMatrix2d(stdCell *cell, void (*orthfunc)(stdCell *cell, int ind, double *func));
double** sc_massMatrix(stdCell *cell);


void GetLIFT2d(stdCell *shape, double **Mes, double **LIFT);
void GetSurfLinM(int N, int Nfaces, int **Fmask, double **Mes);

#endif //DGOM_STDCELL_H