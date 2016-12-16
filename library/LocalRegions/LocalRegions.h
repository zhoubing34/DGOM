#ifndef LOCAL_REGIONS_H
#define LOCAL_REGIONS_H

#include "LibUtilities/LibUtilities.h"
#include "StdRegions/StdRegions.h"

/* public functions */
void Normals2d(int Nv, double *GX, double *GY, double *nx, double *ny, double *sJ);
void GeoFactor2d(int Np, double *x, double *y,
                 double **Dr, double **Ds,
                 double *drdx, double *dsdx,
                 double *drdy, double *dsdy, double *J);

void MapTriCoor(stdCell *shape,
                double *GX, double *GY,
                double *x, double *y);

void MapQuadCoor(stdCell *shape,
                 double *GX, double *GY,
                 double *x, double *y);

#endif // LOCAL_REGIONS_H