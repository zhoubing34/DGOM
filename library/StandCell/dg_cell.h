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
} dg_cell_type;

typedef struct dg_cell{
    dg_cell_type type; ///< cell enum type
    int N; ///< polynomial order
    int Np; ///< number of points
    int Nv; ///< number of vertex
    int Nfaces; ///< number of faces
    int *Nfp; ///< number of points on each face
    int Nfptotal; ///< total number of points on faces

    int **Fmask; ///< index of node at faces
    double **ws; ///< surface Gauss quadrature at each face
    double *wv; ///< volume Gauss quadrature

    double *r, *s, *t; ///< coordinate
    double **V; ///< Vandermonde matrix
    double **M; ///< mass matrix

    double **Dr, **Ds, **Dt; ///< nodal basis derivative matrix
    double **LIFT; ///< lift matrix

    /* float version coefficient */
    dg_real *f_Dr, *f_Ds, *f_LIFT; ///< user specific version

    void (*free_func)(struct dg_cell *cell);
    void (*proj_vert2node)(struct dg_cell *cell, double *vertVal, double *nodeVal);
} dg_cell;

/* create stand cell object */
dg_cell* dg_cell_creat(int N, dg_cell_type type);
/* free dg_cell object */
void dg_cell_free(dg_cell *);
/* project the vertex value to interpolation nodes */
void dg_cell_proj_vert2node(dg_cell *cell, double *vertVal, double *nodeVal);

#define dg_cell_celltype(cell) cell->type
#define dg_cell_N(cell) cell->N
#define dg_cell_Np(cell) cell->Np
#define dg_cell_Nv(cell) cell->Nv
#define dg_cell_Nfaces(cell) cell->Nfaces
#define dg_cell_Nfptotal(cell) cell->Nfptotal
#define dg_cell_Nfp(cell, f) cell->Nfp[f]  ///< number of points on fth face

#endif //DGOM_STDCELL_H