//
// Created by li12242 on 17/1/1.
//

#ifndef DGOM_PHYS_STRONG_VISCOSITY_LDG_FLUX2D_H
#define DGOM_PHYS_STRONG_VISCOSITY_LDG_FLUX2D_H

#include "pf_phys.h"

typedef int (* viscosity_wall_condition_func)(real nx, real ny, real *varM, real *varP,
                                    real *pxM, real *pxP, real *pyM, real *pyP);

void pf_strong_viscosity_LDG_flux2d(physField *phys,
                                    viscosity_wall_condition_func slipwall_condition,
                                    viscosity_wall_condition_func non_slipwall_condition,
                                    real c11, real c12, real c22);

#endif //DGOM_PHYS_STRONG_VISCOSITY_LDG_FLUX2D_H
