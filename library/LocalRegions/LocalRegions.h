#ifndef LOCAL_REGIONS_H
#define LOCAL_REGIONS_H

#include "LibUtilities/LibUtilities.h"
#include "StdRegions/StdRegions.h"


/* public functions */
void Normals(int Nv, double *GX, double *GY, double *nx, double *ny, double *sJ);
void GeometricFactors(int Np, double *x, double *y,
                      double **Dr, double **Ds,
                      double *drdx, double *dsdx,
                      double *drdy, double *dsdy, double *J);

void MapTriCoor(StdRegions2d *shape,
                double *GX, double *GY,
                double *x, double *y);

void MapQuadCoor(StdRegions2d *shape,
                 double *GX, double *GY,
                 double *x, double *y);
#endif