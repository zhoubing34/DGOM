//
// Created by li12242 on 17/1/1.
//

#ifndef DGOM_PHYS_STRONG_VISCOSITY_LDG_FLUX2D_H
#define DGOM_PHYS_STRONG_VISCOSITY_LDG_FLUX2D_H

#include "dg_phys.h"

void dg_phys_LDG_solve(dg_phys *phys, Wall_Condition slipwall_func, Wall_Condition non_slipwall_func, OBC_Fun obc_fun);

#endif //DGOM_PHYS_STRONG_VISCOSITY_LDG_FLUX2D_H
