//
// Created by li12242 on 17/3/14.
//

#ifndef DGOM_DG_CELL_LINE_H
#define DGOM_DG_CELL_LINE_H

#include "dg_cell.h"
dg_cell_info* dg_cell_line_info(int N);
void dg_cell_line_set_nood(dg_cell *cell, int *Np, double **r, double **s, double **t);
void dg_cell_line_orthog_func(int N, int ind, int Np, double *r, double *s, double *t, double *fun);
void dg_cell_line_deri_orthog_func(int N, int ind, int Np, double *r, double *s, double *t,
                                   double *dr, double *ds, double *dt);
void dg_cell_line_proj(dg_cell *cell, double *vertVal, double *nodeVal);
#endif //DGOM_DG_CELL_LINE_H
