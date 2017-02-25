#ifndef DGOM_SC_TRIANGLE_H
#define DGOM_SC_TRIANGLE_H

stdCell* sc_create_tri(int N);
void sc_proj_vert2node_tri(stdCell *cell, double *vertVal, double *nodeVal);

#endif