#ifndef DGOM_SC_TRIANGLE_H
#define DGOM_SC_TRIANGLE_H

void dg_tri_info(dg_cell *cell, int N);
void dg_tri_coord(dg_cell *tri);
int** dg_tri_fmask(dg_cell *tri);
void dg_tri_orthog_func(dg_cell *cell, int ind, double *func);
void dg_tri_deriorthog_func(dg_cell *tri, int ind, double *dr, double *ds, double *dt);
double ** dg_tri_surf_mass_matrix(dg_cell *cell);
void dg_tri_gauss_weight(dg_cell *cell);
void dg_tri_free(dg_cell *cell);
void dg_tri_proj_vert2node(dg_cell *cell, double *vertVal, double *nodeVal);

#endif