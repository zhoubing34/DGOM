//
// Created by li12242 on 12/27/16.
//

#ifndef DGOM_PHYS_STRONG_SURFACE_INTEGRAL2D_H
#define DGOM_PHYS_STRONG_SURFACE_INTEGRAL2D_H

#include "pf_phys.h"

typedef int (* wall_condition_func)(dg_real nx, dg_real ny, dg_real *varM, dg_real *varP);
typedef int (*nodal_flux_func)(dg_real *var, dg_real *Eflux, dg_real *Gflux);
typedef int (*numerical_flux_func)(dg_real nx, dg_real ny, dg_real *varM, dg_real *varP, dg_real *Fhs);

void pf_strong_surface_flux2d
        (physField *phys,
         wall_condition_func slipwall_condition,  // slip wall condition
         wall_condition_func non_slipwall_condition, // non-slip wall condition
         nodal_flux_func nodal_flux, // flux term
         numerical_flux_func numerical_flux);

#endif //DGOM_PHYS_STRONG_SURFACE_INTEGRAL2D_H
