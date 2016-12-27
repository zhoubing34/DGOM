//
// Created by li12242 on 12/27/16.
//

#ifndef DGOM_PHYS_STRONG_SURFACE_INTEGRAL2D_H
#define DGOM_PHYS_STRONG_SURFACE_INTEGRAL2D_H

#include "phs_physField.h"

typedef int (* wall_condition_func)(real nx, real ny, real *varM, real *varP);
typedef int (*nodal_flux_func)(real *var, real *Eflux, real *Gflux);
typedef int (*numerical_flux_func)(real nx, real ny, real *varM, real *varP, real *Fhs);

void phys_strong_surface_integral2d
        (physField *phys,
         wall_condition_func slipwall_condition,  // slip wall condition
         wall_condition_func non_slipwall_condition, // non-slip wall condition
         nodal_flux_func nodal_flux, // flux term
         numerical_flux_func numerical_flux);

#endif //DGOM_PHYS_STRONG_SURFACE_INTEGRAL2D_H
