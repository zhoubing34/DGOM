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

#include <StandCell/sc_stdcell.h>
#include "mr_grid.h"

typedef struct {

    int procid; ///< index of local process
    int nprocs; ///< number of processes
    int dim; ///< dimension

    stdCell *cell; ///< standard element
    geoGrid *grid; ///< geometry grid
    double **x; ///< node coordinate
    double **y; ///< node coordinate
    double **z; ///< node coordinate

    /* info of each element */
    int Nvgeo; ///< number of volume geometric factors
    real *vgeo; ///< volume geometric factors
    double **J; ///< jacobi-coefficients on each nodes

    double *len; ///< radius of each elements
    double *size; ///< area of 2d elements or volume of 3d elements

}multiReg;

/* create of multi-region object */
multiReg* mr_reg_create(geoGrid *grid);

/* free the memory of multiReg */
void mr_reg_free(multiReg *region);

/* integral in the specific element */
double mr_reg_integral(multiReg *region, int ind, double *nodalVal);

#endif //DGOM_MULTIREGIONS_H
