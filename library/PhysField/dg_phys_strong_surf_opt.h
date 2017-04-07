//
// Created by li12242 on 12/27/16.
//

#ifndef DGOM_PHYS_STRONG_SURFACE_INTEGRAL2D_H
#define DGOM_PHYS_STRONG_SURFACE_INTEGRAL2D_H

#include "dg_phys.h"


void dg_phys_strong_surf_opt2d
        (dg_phys *phys,
         Wall_Condition slipwall_func,  // slip wall condition
         Wall_Condition non_slipwall_func, // non-slip wall condition
         OBC_Fun obc_func,
         Nodal_Flux_Fun nodal_flux, // flux term
         Numerical_Flux numerical_flux);

#endif //DGOM_PHYS_STRONG_SURFACE_INTEGRAL2D_H
