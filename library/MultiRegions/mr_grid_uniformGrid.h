#ifndef DGOM_UNIFORMMESH_H
#define DGOM_UNIFORMMESH_H

#include "mr_grid.h"

/* generation of uniform triangle mesh */
geoGrid* mr_grid_createUniformGrid_tri
        (stdCell *shape, int Mx, int My,
         double xmin, double xmax, double ymin, double ymax, int type);

/* generation of uniform quadrilateral mesh */
geoGrid* mr_grid_createUniformGrid_quad
        (stdCell *shape, int Mx, int My,
         double xmin, double xmax, double ymin, double ymax);

#endif // DGOM_UNIFORMMESH_H
