//
// Created by li12242 on 17/3/14.
//

#ifndef DGOM_DG_CELL_VOLUME_H
#define DGOM_DG_CELL_VOLUME_H

#include "dg_cell.h"

dg_cell_volume *dg_cell_volume_create(dg_cell *cell, const dg_cell_creator *creator);
void dg_cell_volume_free(dg_cell_volume *volume);
#endif //DGOM_DG_CELL_VOLUME_H
