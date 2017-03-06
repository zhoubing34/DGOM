#ifndef DGOM_SC_QUADRILATERAL_H
#define DGOM_SC_QUADRILATERAL_H

void dg_quad_proj_vert2node(dg_cell *cell, double *vertVal, double *nodeVal);

void dg_quad_info(dg_cell *cell, int N);
void dg_quad_coord(dg_cell *quad);
int** dg_quad_fmask(dg_cell *quad);
void dg_quad_orthog_func(dg_cell *quad, int ind, double *func);
void dg_quad_deri_orthog_func(dg_cell *quad, int ind, double *dr, double *ds, double *dt);
double ** dg_quad_surf_mass_matrix(dg_cell *cell);
void dg_quad_gauss_weight(dg_cell *cell);
void dg_quad_free(dg_cell *cell);
void dg_quad_proj_vert2node(dg_cell *cell, double *vertVal, double *nodeVal);

#endif