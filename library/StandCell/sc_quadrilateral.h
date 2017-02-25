#ifndef DGOM_SC_QUADRILATERAL_H
#define DGOM_SC_QUADRILATERAL_H

#include "sc_stdcell.h"
stdCell* sc_create_quad(int N);
void sc_proj_vert2node_quad(stdCell *cell, double *vertVal, double *nodeVal);

#endif