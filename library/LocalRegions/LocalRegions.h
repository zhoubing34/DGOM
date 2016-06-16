#ifndef LOCAL_REGIONS_H
#define LOCAL_REGIONS_H

#include "LibUtilities/LibUtilities.h"


/* public functions */
void Normals(int Nv, double *GX, double *GY, double *nx, double *ny, double *sJ);
void GeometricFactors(int Np, const double *x, const double *y,
                      double **Dr, double **Ds,
                      double *drdx, double *dsdx,
                      double *drdy, double *dsdy, double *J);

#endif