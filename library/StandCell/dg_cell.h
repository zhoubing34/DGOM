/**
 * @file
 * Public function for standard element
 *
 * @brief
 * Basic element structure definition
 *
 * @author
 * li12242, Tianjin University, li12242@tju.edu.cn
 */

#ifndef DGOM_STDCELL_H
#define DGOM_STDCELL_H

#include "Utility/utility.h"

typedef enum {
    POINT = 0,
    LINE = 1,
    TRIANGLE=2, ///< triangle
    QUADRIL=3,  ///< quadrilateral
    TETRA=4,    ///< tetrahedron
    TRIPRISM=5, ///< triangle prism
    HEXA=6,     ///< hexahedron
} dg_cell_type;

typedef struct dg_cell_info{
    int N; ///< order
    int Nv; ///< number of vertex
    int Nfaces; ///< number of faces
    dg_cell_type type; ///< cell type
    dg_cell_type *face_type; ///< face type
    double *vr,*vs,*vt; ///< vertex coordinate
    int **FToV; ///< Face to vertex list
}dg_cell_info;

typedef struct dg_cell_volume{
    int Np;
    double *r,*s,*t;
    double *w;
    double **M, **V;
    double **Dr,**Ds,**Dt;
} dg_cell_volume;

typedef struct dg_cell_face{
    int *Nfp;
    int Nfptotal;
    int **Fmask;
    double **ws;
    double **LIFT;
} dg_cell_face;

typedef struct dg_cell{
    dg_cell_info *info;
    dg_cell_volume *volume;
    dg_cell_face *face;
    /* user defined precision */
    dg_real *f_Dr,*f_Ds,*f_Dt,*f_LIFT; ///< user specific version
    void (*proj_vert2node)(struct dg_cell *cell, double *vertVal, double *nodeVal);
} dg_cell;

/** function of create cell_info structure */
typedef dg_cell_info* (*Cell_Info_Create_Func)(int N);
/** function of return node number and the address of pointer */
typedef void (*Cell_Node_Func)(dg_cell *cell, int *Np, double **r, double **s, double **t);
/** function of orthogonal function value at nodes (r,s,t) */
typedef void (*Orthogonal_Func)(int N, int ind, int Np, double *r, double *s, double *t, double *fun);
/** function of derivatives orthogonal function value at nodes (r,s,t) */
typedef void (*Derivative_Orthogonal_Func)(int N, int ind, int Np, double *r, double *s, double *t,
                                           double *dr, double *ds, double *dt);
/** function of creating Fmask matrix */
typedef void (*Fmask_Func)(dg_cell *cell, int **Fmask);
/** function of projecting vertex to nodes */
typedef void (*Proj_Func)(dg_cell *cell, double *vertVal, double *nodeVal);

typedef struct dg_cell_creator{
    Cell_Info_Create_Func cell_info_create; ///< generate dg_cell_info struct
    Cell_Node_Func set_node;
    Orthogonal_Func orthogonal_func; ///< orthogonal function
    Derivative_Orthogonal_Func deri_orthgonal_func; ///<derivative of orthogonal func
    Fmask_Func face_Fmask;
    Proj_Func proj_func;
} dg_cell_creator;

/** create stand cell object */
dg_cell* dg_cell_creat(int N, dg_cell_type type);
/** free dg_cell object */
void dg_cell_free(dg_cell *);
/** project the vertex value to interpolation nodes */
void dg_cell_proj_vert2node(dg_cell *cell, double *vertVal, double *nodeVal);

/* basic info */
#define dg_cell_celltype(cell) cell->info->type
#define dg_cell_Nv(cell) cell->info->Nv
#define dg_cell_N(cell) cell->info->N
#define dg_cell_Nfaces(cell) cell->info->Nfaces
#define dg_cell_vr(cell) cell->info->vr
#define dg_cell_vs(cell) cell->info->vs
#define dg_cell_vt(cell) cell->info->vt
#define dg_cell_FToV(cell) cell->info->FToV
#define dg_cell_facetype(cell) cell->info->face_type

/* volume properties */
#define dg_cell_Np(cell) cell->volume->Np
#define dg_cell_r(cell) cell->volume->r
#define dg_cell_s(cell) cell->volume->s
#define dg_cell_t(cell) cell->volume->t
#define dg_cell_w(cell) cell->volume->w
#define dg_cell_V(cell) cell->volume->V
#define dg_cell_M(cell) cell->volume->M
#define dg_cell_Dr(cell) cell->volume->Dr
#define dg_cell_Ds(cell) cell->volume->Ds
#define dg_cell_Dt(cell) cell->volume->Dt

/* face properties */
#define dg_cell_Nfptotal(cell) cell->face->Nfptotal
#define dg_cell_Nfp(cell) cell->face->Nfp ///< number of points on fth face
#define dg_cell_ws(cell) cell->face->ws
#define dg_cell_Fmask(cell) cell->face->Fmask
#define dg_cell_LIFT(cell) cell->face->LIFT

#endif //DGOM_STDCELL_H