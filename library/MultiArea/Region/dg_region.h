/**
 * @file
 * dg_region structure.
 *
 * @author
 * li12242, Tianjin University, li12242@tju.edu.cn
 */

#ifndef DGOM_MULTIREGIONS_H
#define DGOM_MULTIREGIONS_H

#include "../Grid/dg_grid.h"

typedef struct dg_region{

    dg_grid *grid; ///< geometry grid;
    double **x; ///< node coordinate;
    double **y; ///< node coordinate;
    double **z; ///< node coordinate;

    /* info of each element */
    double **drdx, **drdy, **drdz; ///< transformation of partial derivative on each nodes;
    double **dsdx, **dsdy, **dsdz; ///< transformation of partial derivative on each nodes;
    double **dtdx, **dtdy, **dtdz; ///< transformation of partial derivative on each nodes;
    double **J; ///< jacobi-coefficients on each nodes on each nodes;

    double **nx, **ny, **nz; ///< outward normal vector of each faces;
    double **sJ; ///< jacobi-coefficients of each faces;

    double *len; ///< radius of each elements;
    double *size; ///< area of 2d elements or volume of 3d elements;
    double **face_size; ///< length or area of each faces;

    void (*vol_integral)(struct dg_region *region, int Nfield, int k, dg_real *f_Q, dg_real *c_Q);
    void (*face_integral)(struct dg_region *region, int Nfield, int k, dg_real *f_Q, dg_real *face_Q);
}dg_region;

/* create of multi-region object */
dg_region* dg_region_create(dg_grid *grid);
/* free the memory of dg_region */
void dg_region_free(dg_region *region);

#define dg_region_grid(region) region->grid
#define dg_region_cell(region) dg_grid_cell(region->grid)
#define dg_region_procid(region) region->grid->procid
#define dg_region_nprocs(region) region->grid->nprocs
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
#define dg_region_size(region) (region->size)
#define dg_region_len(region) (region->len)
#define dg_region_face_size(region) (region->face_size)

#endif //DGOM_MULTIREGIONS_H
