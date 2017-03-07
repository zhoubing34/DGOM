/**
 * @file
 * Mesh object for multi-process
 *
 * @brief
 * The basic mesh information
 *
 * @author
 * li12242, Tianjin University, li12242@tju.edu.cn
 */

#ifndef DGOM_MULTIREGIONS_H
#define DGOM_MULTIREGIONS_H

#include "MultiRegions/Grid/dg_grid.h"

typedef struct dg_region{

    int procid; ///< index of local process
    int nprocs; ///< number of processes

    dg_cell *cell; ///< standard element
    dg_grid *grid; ///< geometry grid
    double **x; ///< node coordinate
    double **y; ///< node coordinate
    double **z; ///< node coordinate

    /* info of each element */
    double **drdx, **drdy, **drdz; ///< transformation of partial derivative on each nodes
    double **dsdx, **dsdy, **dsdz; ///< transformation of partial derivative on each nodes
    double **dtdx, **dtdy, **dtdz; ///< transformation of partial derivative on each nodes
    double **J; ///< jacobi-coefficients on each nodes on each nodes

    double **nx, **ny, **nz; ///< outward normal vector of each faces
    double **sJ; ///< jacobi-coefficients of each faces

    double *len; ///< radius of each elements
    double *size; ///< area of 2d elements or volume of 3d elements

    void (*free_func)(struct dg_region *region);
}dg_region;

/* create of multi-region object */
dg_region* mr_reg_create(dg_grid *grid);
/* free the memory of dg_region */
void mr_reg_free(dg_region *region);
/* integral in the specific element */
double mr_reg_integral(dg_region *region, int ind, double *nodalVal);

#endif //DGOM_MULTIREGIONS_H
