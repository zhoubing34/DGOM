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
    info->vs = NULL;
    info->vt = NULL;

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
    *s = NULL;
    *t = NULL;

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

void dg_cell_line_Fmask(dg_cell *cell, int **Fmask){
    const int N = dg_cell_N(cell);
    const int Np = dg_cell_Np(cell);
    const int Nfaces = dg_cell_Nfaces(cell);
    const dg_cell_type *face_type = dg_cell_facetype(cell);

    int m,n,f;
    for(f=0;f<Nfaces;f++){
        dg_cell *face_cell = dg_cell_creat(N, face_type[f]);
        /* assignment of Fmask */
        int Nface_node = dg_cell_Np(face_cell);
        int Nface_vert = dg_cell_Nv(face_cell);
        double f_vr[Nface_vert];
        for(n=0;n<Nface_vert;n++){
            f_vr[n] = dg_cell_vr(cell)[ dg_cell_FToV(cell)[f][n] ];
        }
        double f_r[Nface_node]; // face node coordinate
        face_cell->proj_vert2node(face_cell, f_vr, f_r);
        for(n=0;n<Nface_node;n++){
            double x1 = f_r[n];
            for(m=0;m<Np;m++){
                double x2 = dg_cell_r(cell)[m];
                double d12 = (x1-x2)*(x1-x2);
                if(d12 < EPS) {Fmask[f][n] = m; break;}
            }
        }
    }
    return;
}


void dg_cell_line_proj(dg_cell *cell, double *vertVal, double *nodeVal){
    register int i;
    double *r = dg_cell_r(cell);
    for (i = 0; i < dg_cell_Np(cell); ++i) {
        double ri=r[i];
        nodeVal[i] = 0.5*( vertVal[0]*(1.-ri) + vertVal[1]*(1.+ri) );
    }
    return;
}