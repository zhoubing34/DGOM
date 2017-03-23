#ifndef DGOM_SC_TRIANGLE_H
#define DGOM_SC_TRIANGLE_H

#include "dg_cell.h"

dg_cell_info* dg_cell_tri_info(int N);
void dg_cell_tri_set_node(dg_cell *cell, int *Np, double **r, double **s, double **t);
void dg_cell_tri_orthog_func(int N, int ind, int Np, double *r, double *s, double *t, double *func);
void dg_cell_tri_deriorthog_func(int N, int ind, int Np,
                                 double *r, double *s, double *t,
                                 double *dr, double *ds, double *dt);

void dg_cell_tri_proj(dg_cell *cell, int Nfield, double *vertVal, double *nodeVal);

#endif