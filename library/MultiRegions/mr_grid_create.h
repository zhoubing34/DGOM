#ifndef DGOM_UNIFORMMESH_H
#define DGOM_UNIFORMMESH_H

#include "mr_grid.h"

/* generation of uniform triangle mesh */
geoGrid* mr_grid_create_uniform_tri(
        stdCell *shape, int Mx, int My, double xmin, double xmax, double ymin, double ymax, int type);

/* generation of uniform quadrilateral mesh */
geoGrid* mr_grid_create_uniform_quad(
        stdCell *shape, int Mx, int My, double xmin, double xmax, double ymin, double ymax);

geoGrid* mr_grid_read_file2d(stdCell *shape, char *casename);
#endif // DGOM_UNIFORMMESH_H
