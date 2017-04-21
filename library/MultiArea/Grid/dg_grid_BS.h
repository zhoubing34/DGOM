//
// Created by li12242 on 17/3/9.
//

#ifndef DGOM_DG_GRID_BC_H
#define DGOM_DG_GRID_BC_H

#include "dg_grid.h"

void dg_grid_init_BS(dg_grid *grid);
void dg_grid_set_EToBS2d(dg_grid *grid, int Nsurface, int **SFToV);
void dg_grid_set_EToBS3d(dg_grid *grid, int Nsurface, int **SFToV);

#endif //DGOM_DG_GRID_BC_H
