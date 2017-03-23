//
// Created by li12242 on 17/3/22.
//

#ifndef DGOM_DG_PHYS_TEMPORAL_DISCRETIZE_H
#define DGOM_DG_PHYS_TEMPORAL_DISCRETIZE_H

#include "dg_phys_info.h"

typedef enum{
    TIME_RK45 = 0,
    TIME_AB = 1,
} dg_phys_time_discretize_type;

typedef struct dg_phys_time{
    dg_phys_info *info;

    dg_real *f_rhsQ; ///< RHS data;
    dg_real *f_resQ; ///< residual data;

    dg_phys_time_discretize_type time_discretize_type;
} dg_phys_time;

#endif //DGOM_DG_PHYS_TEMPORAL_DISCRETIZE_H
