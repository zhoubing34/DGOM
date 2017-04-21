#ifndef DGOM_UNIFORMMESH_H
#define DGOM_UNIFORMMESH_H

#include "dg_grid.h"

/* generation of uniform triangle mesh */
dg_grid* dg_grid_create_uniform_tri(dg_cell *cell, int Mx, int My, double xmin, double xmax,
                                    double ymin, double ymax, int type);

/* generation of uniform quadrilateral mesh */
dg_grid* dg_grid_create_uniform_quad(dg_cell *cell, int Mx, int My, double xmin, double xmax,
                                     double ymin, double ymax);

dg_grid* dg_grid_create_from_file2d(dg_cell *cell, char *casename);

#endif // DGOM_UNIFORMMESH_H
