//
// Created by li12242 on 16/12/30.
//

#ifndef DGOM_PHYS_ADD_VISCOSITY_SOLVER_H
#define DGOM_PHYS_ADD_VISCOSITY_SOLVER_H

#include "dg_phys_info.h"

typedef struct dg_phys_LDG{
    dg_phys_info *info;
    dg_real *px;
    dg_real *py;
    dg_real *pz;
    dg_real *vissqrt;

    void (*set_vis)(struct dg_phys_LDG *solver, dg_real *vis);
} dg_phys_LDG;

dg_phys_LDG* dg_phys_LDG_create(dg_phys_info *info);
void dg_phys_LDG_free(dg_phys_LDG *ldg);

#endif //DGOM_PHYS_ADD_VISCOSITY_SOLVER_H
