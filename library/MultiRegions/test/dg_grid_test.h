//
// Created by li12242 on 12/18/16.
//

#ifndef DGOM_MR_GRID_TEST_H
#define DGOM_MR_GRID_TEST_H

#include "multiregion_test.h"

int dg_grid_Ktol_test(dg_grid *grid, int verbose);
int dg_grid_EToV_test(dg_grid *grid, int verbose);
int dg_grid_vertex_test(dg_grid *grid, int verbose);
int dg_grid_connect_test(dg_grid *grid, int verbose);
int dg_grid_EToBS_test(dg_grid *grid, int verbose);

#endif //DGOM_MR_GRID_TEST_H
