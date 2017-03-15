#ifndef DGOM_SC_QUADRILATERAL_H
#define DGOM_SC_QUADRILATERAL_H


dg_cell_info* dg_cell_quad_info(int N);
void dg_cell_quad_set_nood(dg_cell *cell, int *Np, double **r, double **s, double **t);
void dg_cell_quad_orthog_func(int N, int ind, int Np, double *r, double *s, double *t, double *fun);
void dg_cell_quad_deri_orthog_func(int N, int ind, int Np, double *r, double *s, double *t,
                                   double *dr, double *ds, double *dt);
void dg_cell_quad_Fmask(dg_cell *cell, int **Fmask);


int** dg_quad_fmask(dg_cell *quad);
double ** dg_quad_surf_mass_matrix(dg_cell *cell);
void dg_quad_gauss_weight(dg_cell *cell);
void dg_quad_free(dg_cell *cell);
void dg_cell_quad_proj(dg_cell *cell, double *vertVal, double *nodeVal);

#endif