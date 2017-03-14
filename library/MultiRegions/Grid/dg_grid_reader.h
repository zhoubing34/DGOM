#ifndef DGOM_UNIFORMMESH_H
#define DGOM_UNIFORMMESH_H

#include "dg_grid.h"

/* generation of uniform triangle mesh */
dg_grid* dg_grid_uniform_tri(
        dg_cell *cell, int Mx, int My, double xmin, double xmax, double ymin, double ymax, int type);

/* generation of uniform quadrilateral mesh */
dg_grid* dg_grid_uniform_quad(
        dg_cell *cell, int Mx, int My, double xmin, double xmax, double ymin, double ymax);

dg_grid* dg_grid_read_file2d(dg_cell *cell, char *casename);
void dg_grid_read_BSfile2d(dg_grid *grid, char *casename);
void dg_grid_add_EToR(dg_grid *grid, int cell_num, int **RgionID);

#endif // DGOM_UNIFORMMESH_H
