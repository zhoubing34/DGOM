//
// Created by li12242 on 17/3/29.
//

#include "swe_lib.h"

void swe_obc(dg_real nx, dg_real ny, dg_real *f_M, int obc_ind, dg_real *f_P){
    switch (obc_ind){
        case ZERO_GRADIENT:
        case CLAMPED_ALL:
        case CLAMPED_H:
        case CLAMPED_FLOW:
        case FLATHER:
        default:
            fprintf(stderr, "Unknown open boundary type %d\n", obc_ind);
            exit(-1);
    }
    return;
}