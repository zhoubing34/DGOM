//
// Created by li12242 on 17/3/23.
//

#include "dg_phys_limiter.h"

static void dg_phys_limit_func(dg_phys_info *phys_info, double parameter);

dg_phys_limiter* dg_phys_limiter_create(){
    dg_phys_limiter *dg_limit = calloc(1, sizeof(dg_phys_limiter));
    dg_limit->type = BJ_LIMITER;
    dg_limit->limit = dg_phys_limit_func;
    return dg_limit;
}

static void dg_phys_limit_func(dg_phys_info *phys_info, double parameter){

}