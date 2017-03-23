//
// Created by li12242 on 17/3/14.
//

#include "dg_cell_point.h"

dg_cell_info* dg_cell_point_info(int N){
    dg_cell_info *info = (dg_cell_info *)calloc(1, sizeof(dg_cell_info));
    const int Nv = 1;
    const int Nfaces = 0;
    const int Nfv = 0;
    info->N = N;
    info->Nv = Nv;
    info->Nfaces = Nfaces;
    info->type = POINT;
    info->face_type = calloc(Nfaces, sizeof(dg_cell_type));
    info->FToV = matrix_int_create(Nfaces, Nfv);
    info->vr = calloc(Nv, sizeof(double));
    info->vs = calloc(Nv, sizeof(double));
    info->vt = calloc(Nv, sizeof(double));
    return info;
}

void dg_cell_point_set_node(dg_cell *cell, int *Np, double **r, double **s, double **t){
    const int Npt = 1;
    *Np = Npt;
    *r = vector_double_create(1);
    *s = vector_double_create(1);
    *t = vector_double_create(1);

    *r[0] = 0;
    return;
}

void dg_cell_point_orthog_func(int N, int ind, int Np,
                               double *r, double *s, double *t, double *fun){
    fun[0] = 1;
    return;
}

void dg_cell_point_deri_orthog_func(int N, int ind, int Np, double *r, double *s, double *t,
                                    double *dr, double *ds, double *dt){
    /* Vr */
    dr[0] = 0;
    return;
}

void dg_cell_point_proj(dg_cell *cell, int Nfield, double *vertVal, double *nodeVal){
    nodeVal[0] = vertVal[0];
    return;
}