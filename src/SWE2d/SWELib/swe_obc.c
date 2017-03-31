//
// Created by li12242 on 17/3/29.
//

#include "swe_lib.h"

int swe_obc(dg_real nx, dg_real ny, dg_real *f_M, dg_real *f_ext, int obc_ind, dg_real *f_P){
    switch (obc_ind){
        case ZERO_GRADIENT:
        case CLAMPED_ALL:
            f_P[0] = f_ext[0]; f_P[1] = f_ext[1]; f_P[2] = f_ext[2]; f_P[3] = f_M[3];
            break;
        case CLAMPED_H:
        case CLAMPED_FLOW:
        case FLATHER:
        default:
            fprintf(stderr, "Unknown open boundary type %d\n", obc_ind);
            return -1;
    }
    return 0;
}