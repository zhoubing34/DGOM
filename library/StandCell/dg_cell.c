#include "dg_cell.h"
#include "dg_cell_point.h"
#include "dg_cell_line.h"
#include "dg_cell_quadrilateral.h"
#include "dg_cell_triangle.h"
#include "dg_cell_face.h"
#include "dg_cell_volume.h"

#define DEBUG 0
#if DEBUG
#include "Utility/unit_test.h"
#endif

/** creator for creating standard point cell */
static const dg_cell_creator point_creator = {
        dg_cell_point_info,
        dg_cell_point_set_node,
        dg_cell_point_orthog_func,
        dg_cell_point_deri_orthog_func,
        dg_cell_point_proj,
};
/** creator for creating standard line cell */
static const dg_cell_creator line_creator = {
        dg_cell_line_info,
        dg_cell_line_set_nood,
        dg_cell_line_orthog_func,
        dg_cell_line_deri_orthog_func,
        dg_cell_line_proj,
};
/** creator for creating standard triangle cell */
static const dg_cell_creator tri_creator = {
        dg_cell_tri_info,
        dg_cell_tri_set_node,
        dg_cell_tri_orthog_func,
        dg_cell_tri_deriorthog_func,
        dg_cell_tri_proj,
};
/** creator for creating standard quadrilateral cell */
static const dg_cell_creator quad_creator = {
        dg_cell_quad_info,
        dg_cell_quad_set_nood,
        dg_cell_quad_orthog_func,
        dg_cell_quad_deri_orthog_func,
        dg_cell_quad_proj,
};

/** copy the Dr,Dr,Dt and LIFT from double to user specific precision */
static void dg_cell_d2f(dg_cell *cell);

/***
 * @brief creating standard dg_cell structure.
 * @param N order;
 * @param type cell type;
 * @return
 * pointer of dg_cell type.
 */
dg_cell *dg_cell_creat(int N, dg_cell_type type){
    dg_cell *cell = (dg_cell *) calloc(1, sizeof(dg_cell));
    const dg_cell_creator *creator;
    switch (type){
        case POINT:
            creator = &point_creator; break;
        case LINE:
            creator = &line_creator; break;
        case TRIANGLE:
            creator = &tri_creator; break;
        case QUADRIL:
            creator = &quad_creator; break;
        default:
            fprintf(stderr, "%s (%d): Unknown cell type %d\n", __FUNCTION__, __LINE__, type);
            exit(-1);
    }
    cell->info = creator->cell_info_create(N);
    cell->volume = dg_cell_volume_create(cell, creator);
    cell->face = dg_cell_face_create(cell);
    cell->proj_vert2node = creator->proj_func;
    dg_cell_d2f(cell);
    return cell;
}

/**
 * @brief deallocate the memory of dg_cell_info structure.
 * @param info pointer of dg_cell_info type;
 */
static void dg_cell_info_free(dg_cell_info *info){
    free(info->face_type);
    matrix_int_free(info->FToV);
    free(info->vr);
    free(info->vs);
    free(info->vt);
    free(info);
    return;
}

/**
 * @brief deallocate the memory of dg_cell structure.
 * @param cell pointer of dg_cell type;
 */
void dg_cell_free(dg_cell *cell){
    dg_cell_info_free(cell->info);
    dg_cell_volume_free(cell->volume);
    dg_cell_face_free(cell->face);
    free(cell);
    return;
}

/**
 * @brief
 * copy double type of Dr,Ds,Dt and LIFT to user specific precision.
 * @param[in,out] cell pointer of dg_cell type;
 */
static void dg_cell_d2f(dg_cell *cell){
    const int Np = dg_cell_Np(cell);
    int NfpTotal = dg_cell_Nfptotal(cell);
    /* float version */
    size_t sz = (size_t) NfpTotal*Np;
    cell->f_LIFT = (dg_real *) calloc(sz, sizeof(dg_real));

    int sk=0,n,m;
    for(n=0;n<Np;++n){
        for(m=0;m<NfpTotal;++m){
            cell->f_LIFT[sk++] = (dg_real) dg_cell_LIFT(cell)[n][m];
        }
    }

    sz = (size_t) Np*Np;
    cell->f_Dr = (dg_real*) calloc(sz, sizeof(dg_real));
    cell->f_Ds = (dg_real*) calloc(sz, sizeof(dg_real));
    cell->f_Dt = (dg_real*) calloc(sz, sizeof(dg_real));
    sk = 0;
    for(n=0;n<Np;++n){
        for(m=0;m<Np;++m){
            cell->f_Dr[sk  ] = (dg_real) dg_cell_Dr(cell)[n][m];
            cell->f_Ds[sk  ] = (dg_real) dg_cell_Ds(cell)[n][m];
            cell->f_Dt[sk++] = (dg_real) dg_cell_Ds(cell)[n][m];
        }
    }
    return;
}
