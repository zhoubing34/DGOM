#ifndef DGOM_UNIFORMMESH_H
#define DGOM_UNIFORMMESH_H

#include "mr_grid.h"

/* generation of uniform triangle mesh */
dg_grid* mr_grid_create_uniform_tri(
        dg_cell *shape, int Mx, int My, double xmin, double xmax, double ymin, double ymax, int type);

/* generation of uniform quadrilateral mesh */
dg_grid* mr_grid_create_uniform_quad(
        dg_cell *shape, int Mx, int My, double xmin, double xmax, double ymin, double ymax);

dg_grid* mr_grid_read_file2d(dg_cell *shape, char *casename);
#endif // DGOM_UNIFORMMESH_H
