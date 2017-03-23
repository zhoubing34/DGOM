//
// Created by li12242 on 17/3/23.
//

#ifndef DGOM_DG_PHYS_LIMITER_H
#define DGOM_DG_PHYS_LIMITER_H

#include "dg_phys_info.h"

/** type of limiter */
typedef enum{
    BJ_LIMITER = 0,
} Limiter_Type;

/** structure of limiter */
typedef struct dg_phys_limiter{
    Limiter_Type type;
    void (*limit)(dg_phys_info *info, double parameter);
}dg_phys_limiter;

#endif //DGOM_DG_PHYS_LIMITER_H
