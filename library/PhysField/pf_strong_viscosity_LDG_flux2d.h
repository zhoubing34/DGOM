//
// Created by li12242 on 17/1/1.
//

#ifndef DGOM_PHYS_STRONG_VISCOSITY_LDG_FLUX2D_H
#define DGOM_PHYS_STRONG_VISCOSITY_LDG_FLUX2D_H

#include "pf_phys.h"

typedef int (* viscosity_wall_condition_func)(dg_real nx, dg_real ny, dg_real *varM, dg_real *varP,
                                    dg_real *pxM, dg_real *pxP, dg_real *pyM, dg_real *pyP);

void pf_strong_viscosity_LDG_flux2d(physField *phys,
                                    viscosity_wall_condition_func slipwall_condition,
                                    viscosity_wall_condition_func non_slipwall_condition,
                                    dg_real c11, dg_real c12, dg_real c22);

#endif //DGOM_PHYS_STRONG_VISCOSITY_LDG_FLUX2D_H
