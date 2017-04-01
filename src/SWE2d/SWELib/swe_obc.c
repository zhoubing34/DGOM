//
// Created by li12242 on 17/3/29.
//

#include "swe_lib.h"

static void obc_zero_gradient(dg_real *f_M, dg_real *f_P);
static void obc_clamped_all(dg_real *f_M, dg_real *f_ext, dg_real *f_P);
static void obc_clamped_h(dg_real *f_M, dg_real *f_ext, dg_real *f_P);
static void obc_clamped_flux(dg_real *f_M, dg_real *f_ext, dg_real *f_P);
static void obc_flather_flow(dg_real nx, dg_real ny, dg_real *f_M, dg_real *f_ext, dg_real *f_P);
static void obc_flather_h(dg_real nx, dg_real ny, dg_real *f_M, dg_real *f_ext, dg_real *f_P);

int swe_obc(dg_real nx, dg_real ny, dg_real *f_M, dg_real *f_ext, int obc_ind, dg_real *f_P){
    switch (obc_ind){
        case ZERO_GRADIENT:
            obc_zero_gradient(f_M, f_P); break;
        case CLAMPED_ALL:
            obc_clamped_all(f_M, f_ext, f_P); break;
        case CLAMPED_H:
            obc_clamped_h(f_M, f_ext, f_P); break;
        case CLAMPED_FLOW:
            obc_clamped_flux(f_M, f_ext, f_P); break;
        case FLATHER_H:
            obc_flather_h(nx, ny, f_M, f_ext, f_P); break;
        case FLATHER_FLOW:
            obc_flather_flow(nx, ny, f_M, f_ext, f_P); break;
        default:
            fprintf(stderr, "Unknown open boundary type %d\n", obc_ind);
            return -1;
    }
    return 0;
}

static void obc_flather_h(dg_real nx, dg_real ny, dg_real *f_M, dg_real *f_ext, dg_real *f_P){
    extern SWE_Solver solver;

    dg_real hM = f_M[0];
    dg_real hM_ext =  f_ext[0];
    dg_real qx_ext =  f_ext[1];
    dg_real qy_ext =  f_ext[2];
    dg_real qn_ext =  qx_ext*nx + qy_ext*ny;
    //dg_real qv_ext = -qx_ext*ny + qy_ext*nx;
    dg_real qn = qn_ext - sqrt(solver.gra*hM)*(hM_ext - hM);
    //dg_real qv = qv_ext;

    dg_real hP = hM_ext - (qn_ext - qn)/sqrt(solver.gra*hM);
    f_P[0] = 2*hP - hM;
    f_P[1] = f_M[1];
    f_P[2] = f_M[2];
    f_P[3] = f_M[3];
    return;
}

static void obc_flather_flow(dg_real nx, dg_real ny, dg_real *f_M, dg_real *f_ext, dg_real *f_P){
    extern SWE_Solver solver;

    dg_real hM = f_M[0];
    dg_real hM_ext =  f_ext[0];
    dg_real qx_ext =  f_ext[1];
    dg_real qy_ext =  f_ext[2];
    dg_real qn_ext =  qx_ext*nx + qy_ext*ny;
    dg_real qv_ext = -qx_ext*ny + qy_ext*nx;
    dg_real qn = qn_ext - sqrt(solver.gra*hM)*(hM_ext - hM);
    dg_real qv = qv_ext;

    f_P[0] = f_M[0];
    f_P[1] = qn*nx - qv*ny;
    f_P[2] = qn*ny + qv*nx;
    f_P[3] = f_M[3];
    return;
}

static void obc_clamped_flux(dg_real *f_M, dg_real *f_ext, dg_real *f_P){
    f_P[0] = f_M[0];
    f_P[1] = 2*f_ext[1] - f_M[1];
    f_P[2] = 2*f_ext[2] - f_M[2];
    f_P[3] = f_M[3];
    return;
}

static void obc_clamped_h(dg_real *f_M, dg_real *f_ext, dg_real *f_P){
    f_P[0] = 2*f_ext[0] - f_M[0];
    f_P[1] = f_M[1];
    f_P[2] = f_M[2];
    f_P[3] = f_M[3];
    return;
}

/**
 * @brief
 * @param f_M
 * @param f_ext
 * @param f_P
 */
static void obc_clamped_all(dg_real *f_M, dg_real *f_ext, dg_real *f_P){
    f_P[0] = 2*f_ext[0] - f_M[0];
    f_P[1] = 2*f_ext[1] - f_M[1];
    f_P[2] = 2*f_ext[2] - f_M[2];
    f_P[3] = f_M[3];
    return;
}

/**
 * @brief
 * @param f_M
 * @param f_P
 */
static void obc_zero_gradient(dg_real *f_M, dg_real *f_P){
    f_P[0] = f_M[0];
    f_P[1] = f_M[1];
    f_P[2] = f_M[2];
    f_P[3] = f_M[3];
    return;
}