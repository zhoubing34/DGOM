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

typedef enum{
    EDGE_INDICATOR = 0,
} Indicator_Type;

/** structure of limiter */
typedef struct dg_phys_limiter{
    int *indicator;
    void (*set_indicator)(struct dg_phys_limiter *limiter, Indicator_Type type);
    void (*set_limiter)(struct dg_phys_limiter *limiter, Limiter_Type type);
    void (*indicator_func)(dg_phys_info *info, int *indicator);
    void (*limiter_func)(dg_phys_info *info, int *indicator, double parameter);
    void (*limit_trouble_cell)(struct dg_phys_limiter *limiter, dg_phys_info *info, double parameter);
}dg_phys_limiter;

dg_phys_limiter* dg_phys_limiter_create(dg_phys_info *info);
void dg_phys_limiter_free(dg_phys_limiter *limiter);

#endif //DGOM_DG_PHYS_LIMITER_H
