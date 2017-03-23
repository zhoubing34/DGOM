//
// Created by li12242 on 17/3/14.
//

#include "dg_cell_line.h"
#include "Polylib/polylib.h"

#define DEBUG 0

dg_cell_info* dg_cell_line_info(int N){

    dg_cell_info *info = (dg_cell_info *)calloc(1, sizeof(dg_cell_info));
    const int Nv = 2;
    const int Nfaces = 2;
    const int Nfv = 1;
    info->N = N;
    info->Nv = Nv;
    info->Nfaces = Nfaces;
    info->type = LINE;
    info->face_type = calloc(Nfaces, sizeof(dg_cell_type));
    info->FToV = matrix_int_create(Nfaces, Nfv);
    info->vr = calloc(Nv, sizeof(double));
    info->vs = calloc(Nv, sizeof(double));
    info->vt = calloc(Nv, sizeof(double));

    /* faces */
    info->face_type[0] = POINT; info->FToV[0][0] = 0;
    info->face_type[1] = POINT; info->FToV[1][0] = 1;
    /* vertex */
    info->vr[0] = -1; info->vr[1] = 1;

    return info;
}

void dg_cell_line_set_nood(dg_cell *cell, int *Np, double **r, double **s, double **t){
    const int Npt = dg_cell_N(cell)+1;

    *Np = Npt;
    *r = vector_double_create(Npt);
    *s = vector_double_create(Npt);
    *t = vector_double_create(Npt);

    double w[Npt];
    zwglj(*r, w, Npt, 0, 0);
    return;
}

void dg_cell_line_orthog_func(int N, int ind, int Np,
                              double *r, double *s, double *t, double *fun){
    jacobiP(Np, r, fun, ind, 0.0, 0.0);
    return;
}


void dg_cell_line_deri_orthog_func(int N, int ind, int Np, double *r, double *s, double *t,
                                   double *dr, double *ds, double *dt){
    /* Vr */
    GradjacobiP(Np, r, dr, ind, 0, 0);
    return;
}

/**
 * @brief Project from node value from vertex value.
 * @param cell
 * @param Nfield
 * @param vertVal
 * @param nodeVal
 */
void dg_cell_line_proj(dg_cell *cell, int Nfield, double *vertVal, double *nodeVal){
    register int i,fld,sk=0;
    const int Np = dg_cell_Np(cell);
    double *r = dg_cell_r(cell);
    for (i=0;i<Np;++i) {
        double ri=r[i];
        for(fld=0;fld<Nfield;fld++){
            nodeVal[sk++] = 0.5*( vertVal[0*Nfield+fld]*(1.-ri) + vertVal[1*Nfield+fld]*(1.+ri) );
        }
    }
    return;
}