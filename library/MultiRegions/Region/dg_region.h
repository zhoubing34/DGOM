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
dg_region* dg_region_create(dg_grid *grid);
/* free the memory of dg_region */
void dg_region_free(dg_region *region);
/* integral in the specific element */
double dg_region_integral(dg_region *region, int ind, double *nodalVal);

#define dg_region_procid(region) region->procid
#define dg_region_nprocs(region) region->nprocs
#define dg_region_x(region) region->x
#define dg_region_y(region) region->y
#define dg_region_z(region) region->z
#define dg_region_drdx(region) region->drdx
#define dg_region_drdy(region) region->drdy
#define dg_region_drdz(region) region->drdz
#define dg_region_dsdx(region) region->dsdx
#define dg_region_dsdy(region) region->dsdy
#define dg_region_dsdz(region) region->dsdz
#define dg_region_dtdx(region) region->dtdx
#define dg_region_dtdy(region) region->dtdy
#define dg_region_dtdz(region) region->dtdz
#define dg_region_nx(region) region->nx
#define dg_region_ny(region) region->ny
#define dg_region_nz(region) region->nz
#define dg_region_J(region) region->J
#define dg_region_sJ(region) region->sJ

#endif //DGOM_MULTIREGIONS_H
