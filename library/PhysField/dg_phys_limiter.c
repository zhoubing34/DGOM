//
// Created by li12242 on 17/3/23.
//

#include "dg_phys_limiter.h"
#include "PhysField/Limiter/dg_phys_limiter_BJ2d.h"

static void dg_phys_set_limiter(dg_phys_limiter *dg_limit, Limiter_Type type);

static dg_phys_limiter BJ2d_limiter = {BJ_LIMITER, dg_phys_set_limiter, dg_phys_limiter_BJ2d};

dg_phys_limiter* dg_phys_limiter_create(){
    dg_phys_limiter *dg_limit = &BJ2d_limiter; // default limiter - BJ
    return dg_limit;
}

static void dg_phys_set_limiter(dg_phys_limiter *dg_limit, Limiter_Type type){
    switch (type){
        case BJ_LIMITER:
            dg_limit = &BJ2d_limiter; break;
        default:
            fprintf(stderr, "%s (%d): Unknown limiter\n", __FUNCTION__, __LINE__);
            exit(-1);
    }
    return;
}
