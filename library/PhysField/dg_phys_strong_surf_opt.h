//
// Created by li12242 on 12/27/16.
//

#ifndef DGOM_PHYS_STRONG_SURFACE_INTEGRAL2D_H
#define DGOM_PHYS_STRONG_SURFACE_INTEGRAL2D_H

#include "dg_phys.h"

typedef int (*Wall_Condition)(dg_real nx, dg_real ny, dg_real *varM, dg_real *varP);
typedef int (*Nodal_Flux_Fun)(dg_real *var, dg_real *Eflux, dg_real *Gflux);
typedef int (*Numerical_Flux)(dg_real nx, dg_real ny, dg_real *varM, dg_real *varP, dg_real *Fhs);

void dg_phys_strong_surf_opt2d
        (dg_phys *phys,
         Wall_Condition slipwall_func,  // slip wall condition
         Wall_Condition non_slipwall_func, // non-slip wall condition
         Nodal_Flux_Fun nodal_flux, // flux term
         Numerical_Flux numerical_flux);

#endif //DGOM_PHYS_STRONG_SURFACE_INTEGRAL2D_H
