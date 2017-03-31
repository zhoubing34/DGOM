//
// Created by li12242 on 17/3/29.
//

#include "swe_lib.h"

static void obc_zero_gradient(dg_real *f_M, dg_real *f_P);
static void obc_clamped_all(dg_real *f_M, dg_real *f_ext, dg_real *f_P);
static void obc_clamped_h(dg_real *f_M, dg_real *f_ext, dg_real *f_P);
static void obc_clamped_flux(dg_real *f_M, dg_real *f_ext, dg_real *f_P);

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
        case FLATHER:
        default:
            fprintf(stderr, "Unknown open boundary type %d\n", obc_ind);
            return -1;
    }
    return 0;
}

static void obc_clamped_flux(dg_real *f_M, dg_real *f_ext, dg_real *f_P){
    f_P[0] = f_M[0];
    f_P[1] = f_ext[1];
    f_P[2] = f_ext[2];
    f_P[3] = f_M[3];
    return;
}

static void obc_clamped_h(dg_real *f_M, dg_real *f_ext, dg_real *f_P){
    f_P[0] = f_ext[0];
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
    f_P[0] = f_ext[0];
    f_P[1] = f_ext[1];
    f_P[2] = f_ext[2];
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